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

    def fguInitialized(self, typeR, typeD, zBins):
        tR = np.int64 if np.issubdtype(typeR, np.integer) else np.float64
        tD = np.int64 if np.issubdtype(typeD, np.integer) else np.float64
        tDR = np.int64 if (np.issubdtype(typeR, np.integer) and
                           np.issubdtype(typeD, np.integer)) else np.float64

        tBins = self.config.binningTheta()['bins']
        return (np.zeros(tBins, dtype=tR),
                csr_matrix((tBins, zBins), dtype=tDR),
                csr_matrix((tBins, zBins*zBins), dtype=tD))

    def fguOfTheta(self):
        ang = self.getPre("ANG").data
        angzd = csr_matrix(self.getPre("ANGZD").data)

        thetaChunk = ThetaChunker(
            DEGTORAD * self.getPre("centerDec").data["binCenter"],
            DEGTORAD * self.getPre("centerRA").data["binCenter"],
            ang['binDec'], ang['binRA'])

        zBins = angzd.shape[1]
        invThetaBinWidth = ( self.config.binningTheta()['bins']
                             / (self.config.binningTheta()['range'][1] -
                                self.config.binningTheta()['range'][0]))

        fTheta, gThetaZ, uThetaZZ = self.fguInitialized(ang['countR'].dtype,
                                                        angzd.dtype,
                                                        zBins)

        for slice1,slice2 in self.chunks(len(ang)):
            print '.'
            thetas = thetaChunk(slice1, slice2)
            iThetas = (thetas * invThetaBinWidth).astype(np.int16)
            countR1 = ang["countR"][slice1].astype(fTheta.dtype)
            countR2 = ang["countR"][slice2].astype(fTheta.dtype)
            zD1 = angzd[slice2].astype(gThetaZ.dtype)
            zD2 = angzd[slice2].astype(gThetaZ.dtype)

            def ft(dbCnt):
                ccR = np.multiply.outer(dbCnt * countR1, countR2)
                return np.histogram( thetas, weights=ccR, **self.config.binningTheta())[0]

            def gtz(iTh, countR, zD):
                iAng, iZ = zD.nonzero()
                iZs = np.multiply.outer(np.ones(len(countR), dtype=iZ.dtype), iZ)
                weight = np.multiply.outer(countR, zD.data)
                return csr_matrix((weight.flat, (iTh[:,iAng].flat, iZs.flat)),
                                  shape=gThetaZ.shape)

            def utzz(dbCnt):
                iAng1, iZ1 = zD1.nonzero()
                iAng2, iZ2 = zD2.nonzero()
                iTh = iThetas[iAng1][:,iAng2]
                iZ1s = np.multiply.outer(iZ1, np.ones(len(iZ2), dtype=iZ2.dtype))
                iZ2s = np.multiply.outer(np.ones(len(iZ1), dtype=iZ1.dtype), iZ2)
                weight = np.multiply.outer(dbCnt * zD1.data, zD2.data)

                iZZs = zBins*iZ1s + iZ2s
                return csr_matrix((weight.flat, (iTh.flat, iZZs.flat)),
                                  shape=uThetaZZ.shape)

            ii = slice1 == slice2
            dbl = 1 if ii else 2
            fTheta += ft(dbl)
            uThetaZZ += utzz(dbl)
            gThetaZ += gtz(iThetas, countR1, zD2)
            if not ii:
                gThetaZ += gtz(iThetas.T, countR2, zD1)
            pass

        if self.iJob is None:
            fTheta /= 2
            uThetaZZ /= 2

        fThetaRec = np.array(fTheta, dtype = [('count',fTheta.dtype)])
        hdu = fits.BinTableHDU(fThetaRec, name="fTheta")
        hdu2 = fits.ImageHDU(gThetaZ.toarray(), name="gThetaZ")
        return [hdu, hdu2]

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
