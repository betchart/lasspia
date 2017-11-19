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

class Chunker(object):
    def __init__(self, comb):
        self.ang = comb.getPre("ANG").data
        self.angzd = csr_matrix(comb.getPre("ANGZD").data)

        self.thetaChunk = ThetaChunker(
            DEGTORAD * comb.getPre("centerDec").data["binCenter"],
            DEGTORAD * comb.getPre("centerRA").data["binCenter"],
            self.ang['binDec'], self.ang['binRA'])

        self.typeR = self.largeType(self.ang['countR'])
        self.typeD = self.largeType(self.angzd)
        self.typeDR = type(self.typeR()*self.typeD())

        self.zBins = self.angzd.shape[1]
        self.binningTheta = comb.config.binningTheta()
        self.invThetaBinWidth = (
            lambda d: d['bins'] / (d['range'][1] - d['range'][0])
        )(self.binningTheta)

    @staticmethod
    def largeType(a):
        return np.int64 if np.issubdtype(a.dtype, np.integer) else np.float64

    def set(self, slice1, slice2):
        self.countR1 = self.ang["countR"][slice1].astype(self.typeR)
        self.countR2 = self.ang["countR"][slice2].astype(self.typeR)

        self.zD1 = self.angzd[slice1].astype(self.typeDR)
        self.zD2 = self.angzd[slice2].astype(self.typeDR)

        self.thetas = self.thetaChunk(slice1, slice2)
        self.iThetas = (self.invThetaBinWidth * self.thetas).astype(np.int16)
        self.ii = slice1 == slice2
        self.dbl = 1 if self.ii else 2


class combinatorial(baofast.routine):

    def __call__(self):
        self.hdus.append( self.binCentersTheta() )
        self.hdus.append( self.getPre('centerZ'))
        self.hdus.append( self.getPre('pdfZ'))
        self.hdus.extend( self.fguHDU(*self.fguLoop()) )
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

    def fguInit(self, typeR, typeD, zBins):
        tBins = self.config.binningTheta()['bins']
        return (np.zeros(tBins, dtype=typeR),
                csr_matrix((tBins, zBins), dtype=type(typeR()*typeD())),
                csr_matrix((tBins, zBins*zBins), dtype=typeD))

    @staticmethod
    def ft(ch):
        ccR = np.multiply.outer(ch.dbl * ch.countR1, ch.countR2)
        return np.histogram( ch.thetas, weights=ccR, **ch.binningTheta)[0]

    @staticmethod
    def gtz(iTh, countR, zD, shp):
        iAng, iZ = zD.nonzero()
        iZs = np.multiply.outer(np.ones(len(countR), dtype=iZ.dtype), iZ)
        weight = np.multiply.outer(countR, zD.data)
        return csr_matrix((weight.flat, (iTh[:,iAng].flat, iZs.flat)),
                          shape=shp)

    @staticmethod
    def utzz(ch, shp):
        iAng1, iZ1 = ch.zD1.nonzero()
        iAng2, iZ2 = ch.zD2.nonzero()
        iTh = ch.iThetas[iAng1][:,iAng2]
        iZ1s = np.multiply.outer(iZ1, np.ones(len(iZ2), dtype=iZ2.dtype))
        iZ2s = np.multiply.outer(np.ones(len(iZ1), dtype=iZ1.dtype), iZ2)
        weight = np.multiply.outer(ch.dbl * ch.zD1.data, ch.zD2.data)

        iZ = np.minimum(iZ1s, iZ2s)
        diZ = np.abs(iZ1s - iZ2s)

        iZdZs = ch.zBins*iZ + diZ
        return csr_matrix((weight.flat, (iTh.flat, iZdZs.flat)),
                          shape=shp)

    def fguLoop(self):
        ch = Chunker(self)

        (fTheta,
         gThetaZ,
         uThetaZZ) = self.fguInit(ch.typeR, ch.typeD, ch.zBins)

        for chunk in self.chunks(len(ch.ang)):
            ch.set(*chunk)
            fTheta += self.ft(ch)
            uThetaZZ += self.utzz(ch, uThetaZZ.shape)
            gThetaZ += self.gtz(ch.iThetas, ch.countR1, ch.zD2, gThetaZ.shape)
            if not ch.ii:
                gThetaZ += self.gtz(ch.iThetas.T, ch.countR2, ch.zD1, gThetaZ.shape)
            pass

        if self.iJob is None:
            fTheta /= 2
            uThetaZZ /= 2

        return (fTheta, gThetaZ, uThetaZZ)

    def fguHDU(self, fTheta, gThetaZ, uThetaZZ):

        fThetaRec = np.array(fTheta, dtype = [('count',fTheta.dtype)])
        iTheta, iZdZ = uThetaZZ.nonzero()
        c1 = fits.Column(name='binTheta', array = iTheta.astype(np.int16), format='I')
        c2 = fits.Column(name='binZdZ', array = iZdZ, format='J')
        c3 = fits.Column(name='count', array = uThetaZZ.data.astype(np.float32), format='E')

        return [fits.BinTableHDU(fThetaRec, name="fTheta"),
                fits.ImageHDU(gThetaZ.toarray(), name="gThetaZ"),
                fits.BinTableHDU.from_columns([c1,c2,c3], name="uThetaZZ")
        ]

    @property
    def inputFileName(self):
        return self.config.stageFileName('preprocessing')

    def getPre(self, name):
        hdulist = fits.open(self.inputFileName)
        return hdulist[name]

    def combineOutput(self):
        jobFiles = [self.outputFileName + self.jobString(iJob)
                    for iJob in range(self.nJobs)]

        with fits.open(jobFiles[0]) as h0:
            (fTheta,
             gThetaZ,
             uThetaZZ) = self.fguInit(h0['fTheta'].data['count'].dtype.type,
                                      h0['uThetaZZ'].data['count'].dtype.type,
                                      len(h0['centerZ'].data))

            for jF in jobFiles:
                with fits.open(jF) as h:
                    for name in ['centerTheta','centerZ','pdfZ']:
                        assert baofast.utils.identicalHDUs(name, h0, h)
                    fTheta += h['fTheta'].data['count']
                    gThetaZ += csr_matrix(h['gThetaZ'].data, shape=gThetaZ.shape)
                    u = h['uThetaZZ'].data
                    uThetaZZ += csr_matrix((u['count'], (u['binTheta'],u['binZdZ'])),
                                           shape=uThetaZZ.shape)

            self.hdus.append(h0["centerTheta"])
            self.hdus.append(h0["centerZ"])
            self.hdus.append(h0["pdfZ"])
            fTheta /= 2
            uThetaZZ /= 2
            self.hdus.extend(self.fguHDU(fTheta, gThetaZ, uThetaZZ))
            self.writeToFile()
