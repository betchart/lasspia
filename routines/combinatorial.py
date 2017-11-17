import baofast
import math
import numpy as np
from scipy.sparse import csr_matrix
from astropy.io import fits

DEGTORAD = math.pi / 180

class ThetaChunker(object):

    def __init__(self, centersDec, centersRA, binsDec, binsRA):
        sinDec = np.sin(centersDec)
        cosDec = np.cos(centersDec)
        self.sinsinDec = np.multiply.outer(sinDec,sinDec)
        self.coscosDec = np.multiply.outer(cosDec,cosDec)
        self.cosDeltaRA = np.cos(np.subtract.outer(centersRA,centersRA))
        self.binsDec = binsDec
        self.binsRA = binsRA

    def __call__(self, slice1, slice2):
        cosTheta = (
            self.cosDeltaRA[self.binsRA[slice1]][ :,self.binsRA[slice2]] *
            self.coscosDec[self.binsDec[slice1]][:,self.binsDec[slice2]] +
            self.sinsinDec[self.binsDec[slice1]][:,self.binsDec[slice2]])
        np.clip(cosTheta,-1,1,out=cosTheta)
        return np.arccos(cosTheta,out=cosTheta)


class combinatorial(baofast.routine):

    def __call__(self):
        self.hdus.append( self.binCentersTheta() )
        self.hdus.extend( self.fguOfTheta() )
        self.writeToFile()

    def binCentersTheta(self):
        centers = np.array( baofast.utils.centers(self.config.edgesTheta()),
                            dtype = [("binCenter", np.float64)])
        return fits.BinTableHDU(centers, name="centerTheta")

    def chunks(self, size):
        # Full array of indices possible in place of slices: regioning.
        # However, fancy indices makes a copy, while slice makes a view
        splits = range(0, size, self.config.chunkSize())
        slices = [slice(i,j) for i,j in zip(splits,splits[1:]+[None])]
        return [(slices[i],jSlice)
                for i in range(len(slices))
                for jSlice in slices[i:]][self.iJob::self.nJobs]

    def fguOfTheta(self):
        ang = self.getPre("ANG").data
        angzd = csr_matrix(self.getPre("ANGZD").data)
        cType = np.int64 if np.issubdtype(ang["countR"].dtype, np.integer) else np.float64
        frq = np.zeros(len(self.config.edgesTheta())-1, dtype=cType)

        thetaChunk = ThetaChunker(
            DEGTORAD * self.getPre("centerDec").data["binCenter"],
            DEGTORAD * self.getPre("centerRA").data["binCenter"],
            ang['binDec'], ang['binRA'])

        for slice1,slice2 in self.chunks(len(ang)):
            equalSlices = slice1 == slice2
            print '.'
            chunkT = thetaChunk(slice1, slice2)
            countR1 = ang["countR"][slice1]
            countR2 = ang["countR"][slice2]
            zD1 = angzd[slice2]
            zD2 = angzd[slice2]

            ccR = np.multiply.outer(countR1, countR2).astype(frq.dtype)
            if not equalSlices: ccR *= 2 # fill histogram with twice-weights
            frq += np.histogram( chunkT, weights = ccR,
                                 **self.config.binningTheta())[0]

            iAng, iZ = zD2.nonzero()
            thetas = chunkT[:,iAng]
            iZs = np.multiply.outer(np.ones(len(countR1), dtype=np.int16), iZ.astype(np.int16))
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
        fTheta = np.array(frq, dtype = [('count',frq.dtype)])
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
