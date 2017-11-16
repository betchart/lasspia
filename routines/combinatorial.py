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
            print '.'
            chunkT = thetaChunk(slice1, slice2)
            countR1 = ang["countR"][slice1]
            countR2 = ang["countR"][slice2]
            zD1 = angzd[slice2]
            zD2 = angzd[slice2]

            ccR = np.multiply.outer(countR1, countR2).astype(cType)
            if slice1 != slice2: ccR *= 2 # fill histogram with twice-weights
            frq += np.histogram( chunkT, weights = ccR,
                                 **self.config.binningTheta())[0]

            iAng, iZ = zD2.nonzero()
            thetas = chunkT[:,iAng]
            iZs = np.multiply.outer(np.ones(len(countR1), dtype=np.int16), iZ).astype(np.int16)
            weight = np.multiply.outer(countR1, zD2.data) # relies on zD2.data in order with zD2.nonzero()

            zBins = zD1.shape[1]
            binningThetaZ = self.config.binningDD([self.config.binningTheta(),
                                                   {"bins": zBins, "range":(0,zBins)}])

            grq, x, y = np.histogram2d( thetas.ravel(),
                                        iZs.ravel(),
                                        weights=weight.ravel(),
                                        **binningThetaZ)

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
