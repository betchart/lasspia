import baofast
import math
import numpy as np
from scipy.sparse import csr_matrix
from astropy.io import fits
import sys

DEGTORAD = math.pi / 180

class ThetaChunker(object):

    @staticmethod
    def trigProducts(centersDec, centersRA):
        sinDec = np.sin(centersDec)
        cosDec = np.cos(centersDec)
        return {'sinsinDec': np.multiply.outer(sinDec,sinDec),
                'coscosDec': np.multiply.outer(cosDec,cosDec),
                'cosDeltaRA': np.cos(np.subtract.outer(centersRA,centersRA))}

    def __init__(self, trigPrd, binsDecRA1, binsDecRA2):
        for k,v in trigPrd.iteritems(): setattr(self, k, v)
        self.binsDec1, self.binsRA1 = binsDecRA1
        self.binsDec2, self.binsRA2 = binsDecRA2

    def __call__(self, slice1, slice2):
        cosTheta = (
            self.cosDeltaRA[self.binsRA1[slice1]][ :,self.binsRA2[slice2]] *
            self.coscosDec[self.binsDec1[slice1]][:,self.binsDec2[slice2]] +
            self.sinsinDec[self.binsDec1[slice1]][:,self.binsDec2[slice2]])
        np.clip(cosTheta,-1,1,out=cosTheta)
        return np.arccos(cosTheta,out=cosTheta)


class combinatorial(baofast.routine):

    def __call__(self):
        self.trig = ThetaChunker.trigProducts(
            DEGTORAD * self.getPre("centerDec").data["binCenter"],
            DEGTORAD * self.getPre("centerRA").data["binCenter"] )

        self.hdus.append( self.binCentersTheta() )
        self.hdus.extend( self.fguOfTheta() )
        self.writeToFile()

    def binCentersTheta(self):
        centers = np.array( baofast.utils.centers(self.config.edgesTheta()),
                            dtype = [("binCenter", np.float64)])
        return fits.BinTableHDU(centers, name="centerTheta")

    def fguOfTheta(self):
        ang = self.getPre("ANG").data
        angzd = csr_matrix(self.getPre("ANGZD").data)
        cType = np.int64 if np.issubdtype(type(ang["countR"][0]), np.integer) else np.float64
        frq = np.zeros(len(self.config.edgesTheta())-1, dtype=cType)
        binsDecRA = (ang['binDec'], ang['binRA'])
        thetaChunk = ThetaChunker(self.trig, binsDecRA, binsDecRA)

        splits = range(0, len(ang), self.config.chunkSize())
        slices = [slice(i,j) for i,j in zip(splits,splits[1:]+[None])] # full array of indices possible in place of slices: regioning
        chunks = [(slices[i],jSlice)
                  for i in range(len(slices))
                  for jSlice in slices[i:]][self.iJob::self.nJobs]

        for slice1,slice2 in chunks:
            chunkT = thetaChunk(slice1, slice2)
            countcount = np.multiply.outer(ang["countR"][slice1],
                                           ang["countR"][slice2]).astype(cType)
            if slice1 != slice2: countcount *= 2 # fill histogram with twice-weights
            frq += np.histogram( chunkT, weights = countcount,
                                 **self.config.binningTheta())[0]

            # thetas = np.array(len(slicez) * [chunkT]).T -- 3D
            # just histogram bin index for z: int16 sufficient
            # zs = np.array(len(chunkT) * [range(len(zcenters))])
            # weights = outer( ang['countR'][slice1], angzD['countD'][slice2])
            # mask zero weights and unravel
            # gthetaz = np.histogram2D(thetas[mask].flat, zs[mask].flat, weights[mask].flat, binning)

        if self.iJob is None:
            frq /= 2
        fTheta = np.array(frq, dtype = [('count',cType)])
        hdu = fits.BinTableHDU(fTheta, name="fTheta")
        return [hdu]

    @property
    def inputFileName(self):
        return self.config.stageFileName('preprocessing')

    def getPre(self, name):
        hdulist = fits.open(self.inputFileName)
        return hdulist[name]

    def combineOutput(self):
        jobFiles = [self.outputFileName + self.jobString(iJob)
                    for iJob in range(self.nJobs)]

        jobHDUs = [fits.open(f) for f in jobFiles]

        assert all([baofast.utils.identicalHDUs("centerTheta", jobHDUs[0], h)
                    for h in jobHDUs])
        self.hdus.append(jobHDUs[0]["centerTheta"])

        fTheta = jobHDUs[0]['fTheta']
        fTheta.data['count'] = np.zeros(fTheta.data['count'].shape,
                                   dtype=fTheta.data['count'].dtype)
        for hdu in jobHDUs:
            fTheta.data['count'] += hdu['fTheta'].data['count']
        fTheta.data['count'] /= 2
        self.hdus.append(fTheta)

        self.writeToFile()
