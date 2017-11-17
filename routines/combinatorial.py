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

        tBins = len(self.config.edgesTheta()) - 1
        zBins = angzd.shape[1]
        binningiZ = {"bins": zBins, "range":(0,zBins)}
        binningTheta = self.config.binningTheta()
        binningThetaZ = self.config.binningDD([binningTheta, binningiZ])
        binningThetaZdZ = self.config.binningDD([binningTheta, binningiZ, binningiZ])

        typeR = np.int64 if np.issubdtype(ang["countR"].dtype, np.integer) else np.float64
        typeD = np.int64 if np.issubdtype(angzd.dtype, np.integer) else np.float64
        typeDR = np.int64 if np.issubdtype(typeR, np.integer) and np.issubdtype(typeD, np.integer) else np.float64
        fTheta = np.zeros(tBins, dtype=typeR)
        gThetaZ = np.zeros((tBins, zBins), dtype=typeDR)
        uThetaZdZ = np.zeros((tBins, zBins, zBins), dtype=typeD)

        thetaChunk = ThetaChunker(
            DEGTORAD * self.getPre("centerDec").data["binCenter"],
            DEGTORAD * self.getPre("centerRA").data["binCenter"],
            ang['binDec'], ang['binRA'])

        for slice1,slice2 in self.chunks(len(ang)):
            equalSlices = slice1 == slice2
            print '.'
            chunkT = thetaChunk(slice1, slice2)
            countR1 = ang["countR"][slice1].astype(fTheta.dtype)
            countR2 = ang["countR"][slice2].astype(fTheta.dtype)
            zD1 = angzd[slice2].astype(gThetaZ.dtype)
            zD2 = angzd[slice2].astype(gThetaZ.dtype)

            def ft(dbCnt):
                ccR = np.multiply.outer(dbCnt * countR1, countR2)
                return np.histogram( chunkT, weights=ccR, **binningTheta)[0]

            def gtz(cT, countR, zD):
                iAng, iZ = zD.nonzero()
                thetas = cT[:,iAng]
                iZs = np.multiply.outer(np.ones(len(countR), dtype=np.int16),
                                        iZ.astype(np.int16))
                weight = np.multiply.outer(countR, zD.data) # zD.data aligned w zD.nonzero()?
                return np.histogram2d( thetas.ravel(), iZs.ravel(),
                                       weights=weight.ravel(), **binningThetaZ)[0]

            def utzdz(dbCnt):
                iAng1, iZ1 = zD1.nonzero()
                iAng2, iZ2 = zD2.nonzero()
                thetas = chunkT[iAng1][:,iAng2]
                iZ1s = np.multiply.outer(iZ1.astype(np.int16),
                                         np.ones(len(iZ2), dtype=np.int16))
                iZ2s = np.multiply.outer(np.ones(len(iZ1), dtype=np.int16),
                                         iZ2.astype(np.int16))
                weight = np.multiply.outer(dbCnt * zD1.data, zD2.data)
                triplets = np.vstack([thetas.ravel(),
                                      np.minimum(iZ1s,iZ2s).ravel(),
                                      np.abs(iZ1s-iZ2s).ravel()]).T
                return np.histogramdd(triplets, weights=weight.ravel(),
                                      **binningThetaZdZ)[0]

            dbl = 1 if equalSlices else 2
            fTheta += ft(dbl)
            gThetaZ += gtz(chunkT, countR1, zD2)
            if not equalSlices:
                gThetaZ += gtz(chunkT.T, countR2, zD1)
            uThetaZdZ += utzdz(dbl)

        if self.iJob is None:
            fTheta /= 2
            uThetaZdZ /= 2
        fThetaRec = np.array(fTheta, dtype = [('count',fTheta.dtype)])
        hdu = fits.BinTableHDU(fThetaRec, name="fTheta")
        hdu2 = fits.ImageHDU(gThetaZ, name="gThetaZ")
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
