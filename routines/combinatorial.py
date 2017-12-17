import lasspia as La
import math
import numpy as np
from scipy.sparse import csr_matrix
from astropy.io import fits
from lasspia.timing import timedHDU

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
        self.zBins = len(comb.getPre('centerz').data['binCenter'])
        self.binningTheta = comb.config.binningTheta()
        self.diZmax = comb.config.maxDeltaZ()
        if self.diZmax: self.diZmax *= La.utils.invBinWidth(comb.config.binningZ())

        self.ang = comb.getPre("ANG").data
        angzd = comb.getPre('ANGZD').data
        self.angzd = csr_matrix((angzd["count"], (angzd['iAlign'],angzd['iZ'])), shape=(len(self.ang), self.zBins))
        self.angzd_e2 = csr_matrix((angzd["err2"], (angzd['iAlign'],angzd['iZ'])), shape=(len(self.ang), self.zBins))

        self.thetaChunk = ThetaChunker(
            DEGTORAD * comb.getPre("centerDec").data["binCenter"],
            DEGTORAD * comb.getPre("centerRA").data["binCenter"],
            self.ang['binDec'], self.ang['binRA'])

        self.typeR = self.largeType(self.ang['countR'])
        self.typeD = self.largeType(self.angzd)
        self.typeDR = type(self.typeR()*self.typeD())

    @staticmethod
    def largeType(a):
        return np.int64 if np.issubdtype(a.dtype, np.integer) else np.float64

    def set(self, slice1, slice2):
        self.countR1 = self.ang["countR"][slice1].astype(self.typeR)
        self.countR2 = self.ang["countR"][slice2].astype(self.typeR)

        self.zD1 = self.angzd[slice1].astype(self.typeDR)
        self.zD2 = self.angzd[slice2].astype(self.typeDR)

        self.zD1_e2 = self.angzd_e2[slice1].astype(self.typeDR)
        self.zD2_e2 = self.angzd_e2[slice2].astype(self.typeDR)

        self.thetas = self.thetaChunk(slice1, slice2)
        self.iThetas = La.utils.toBins(self.thetas, self.binningTheta, dtype=np.int16)
        self.ii = np.all(slice1 == slice2)
        self.dbl = 1 if self.ii else 2


class combinatorial(La.routine):

    def __call__(self):
        self.hdus.append( self.binCentersTheta() )
        self.hdus.append( self.getPre('centerZ'))
        self.hdus.append( self.getPre('pdfZ'))
        self.hdus.extend( self.fguHDU() )
        self.writeToFile()

    @timedHDU
    def binCentersTheta(self):
        centers = np.array( La.utils.centers(self.config.edgesTheta()),
                            dtype = [("binCenter", np.float64)])
        return fits.BinTableHDU(centers, name="centerTheta")

    def chunks(self, size):
        if any([self.config.maxDeltaRA(),
                self.config.maxDeltaDec()]):
            return self.__orderedChunks__()
        return self.__sliceChunks__(size)

    def __sliceChunks__(self, size):
        slices = La.utils.slices(size, self.config.chunkSize())
        return [(slices[i],jSlice)
                for i in range(len(slices))
                for jSlice in slices[i:]][self.iJob::self.nJobs]

    def __indexChunks__(self, size):
        ang = self.getPre('ang').data
        def regionIndices(ra,dc):
            mask = reduce(np.logical_and, [ra.start <= ang['binRA'],
                                           ang['binRA'] < ra.stop,
                                           dc.start <= ang['binDec'],
                                           ang['binDec'] < dc.stop],
                          0 <= ang['binDec'])
            return mask.nonzero()[0]

        regions = [[(i,j) for i in self.config.binRegionsRA()]
                   for j in self.config.binRegionsDec()]

        rpairs = sum([
            [(r,r) for row in regions for r in row],                                             # self
            [(r1,r2) for row in regions for r1,r2 in zip(row,row[1:])],                          # right
            [(r1,r2) for row1,row2 in zip(regions,regions[1:]) for r1,r2 in zip(row1,row2)],     # below
            [(r1,r2) for row1,row2 in zip(regions,regions[1:]) for r1,r2 in zip(row1,row2[1:])], # below right
            [(r1,r2) for row1,row2 in zip(regions,regions[1:]) for r1,r2 in zip(row1[1:],row2)], # below left
            ], [])
        cz = self.config.chunkSize()
        for r1,r2 in rpairs:
            i1 = regionIndices(*r1)
            i2 = regionIndices(*r2)
            schnks = ( self.__sliceChunks__(len(i1)) if r1==r2 else
                       [(s1,s2)
                        for s1 in La.utils.slices(len(i1),cz)
                        for s2 in La.utils.slices(len(i2),cz)][self.iJob::self.nJobs])
            for s1,s2 in schnks:
                yield i1[s1], i2[s2]
        return

    def __orderedChunks__(self):
        ang = self.getPre('ang').data
        sp =  self.getPre('slicePoints').data['bin']

        def sliceMnMx(z):
            slc = slice(*z)
            mnRA,mxRA = min(ang['binRA'][slc]), max(ang['binRA'][slc])
            mnD, mxD = min(ang['binDec'][slc]), max(ang['binDec'][slc])
            return slc, mnRA, mxRA, mnD, mxD

        def inRange((iS,imnRA,imxRA,imnD,imxD),
                    (jS,jmnRA,jmxRA,jmnD,jmxD)):
            iDeltaRA = La.utils.invBinWidth(self.config.binningRA()) * self.config.maxDeltaRA()
            iDeltaDec = La.utils.invBinWidth(self.config.binningDec()) * self.config.maxDeltaDec()
            return (max(imnRA,jmnRA) - min(imxRA,jmxRA) < iDeltaRA and
                    max(imnD,jmnD) - min(imxD,jmxD) < iDeltaDec)

        slices = [(sliceMnMx(z)) for z in zip(sp,sp[1:])]
        return [(iSlcMnMx[0], jSlcMnMx[0])
                for i,iSlcMnMx in enumerate(slices)
                for jSlcMnMx in slices[i:]
                if inRange(iSlcMnMx,jSlcMnMx)][self.iJob::self.nJobs]

    def fgueInit(self, typeR, typeD, zBins):
        tBins = self.config.binningTheta()['bins']
        return (np.zeros(tBins, dtype=typeR),
                csr_matrix((tBins, zBins), dtype=type(typeR()*typeD())),
                csr_matrix((tBins, zBins*zBins), dtype=typeD),
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
        weight_e2 = np.multiply.outer(ch.dbl * ch.zD1_e2.data, ch.zD2_e2.data)

        iZ = np.minimum(iZ1s, iZ2s)
        diZ = np.abs(iZ1s - iZ2s)
        iZdZs = ch.zBins*iZ + diZ
        mask = slice(None) if not ch.diZmax else (diZ < ch.diZmax)

        def calcU(w):
            unbinned = 0 if not ch.diZmax else w[mask==False].sum()
            return ( csr_matrix((w[mask].flat, (iTh[mask].flat, iZdZs[mask].flat)),
                                shape=shp)
                     + csr_matrix(([unbinned], ([shp[0]-1], [shp[1]-1])), shape=shp))

        return calcU(weight), calcU(weight_e2)

    def fgueLoop(self):
        ch = Chunker(self)

        (fTheta,
         gThetaZ,
         uThetaZZ,
         uThetaZZe2) = self.fgueInit(ch.typeR, ch.typeD, ch.zBins)

        for chunk in self.chunks(len(ch.ang)):
            ch.set(*chunk)
            fTheta += self.ft(ch)
            dU, dUe2 = self.utzz(ch, uThetaZZ.shape)
            uThetaZZ += dU
            uThetaZZe2 += dUe2
            gThetaZ += self.gtz(ch.iThetas, ch.countR1, ch.zD2, gThetaZ.shape)
            if not ch.ii:
                gThetaZ += self.gtz(ch.iThetas.T, ch.countR2, ch.zD1, gThetaZ.shape)
            pass

        if self.iJob is None:
            fTheta /= 2
            uThetaZZ /= 2
            uThetaZZe2 /= 2

        return (fTheta, gThetaZ, uThetaZZ, uThetaZZe2)

    @timedHDU
    def fguHDU(self, fgue=None):
        fTheta, gThetaZ, uThetaZZ, uThetaZZe2 = fgue if fgue else self.fgueLoop()

        fThetaRec = np.array(fTheta, dtype = [('count',fTheta.dtype)])
        iTheta, iZdZ = uThetaZZ.nonzero()
        c1 = fits.Column(name='binTheta', array = iTheta.astype(np.int16), format='I')
        c2 = fits.Column(name='binZdZ', array = iZdZ, format='J')
        c3 = fits.Column(name='count', array = uThetaZZ.data.astype(np.float32), format='E')
        c4 = fits.Column(name='err2', array = uThetaZZe2.data.astype(np.float32), format='E')

        return [fits.BinTableHDU(fThetaRec, name="fTheta"),
                fits.ImageHDU(gThetaZ.toarray(), name="gThetaZ"),
                fits.BinTableHDU.from_columns([c1,c2,c3,c4], name="uThetaZZ")
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
             uThetaZZ,
             uThetaZZe2) = self.fgueInit(h0['fTheta'].data['count'].dtype.type,
                                         h0['uThetaZZ'].data['count'].dtype.type,
                                         len(h0['centerZ'].data))
            cputime = 0.

            for jF in jobFiles:
                with fits.open(jF) as h:
                    for name in ['centerZ','pdfZ']:
                        assert La.utils.hduDiff(name, h0, h).identical
                    assert La.utils.hduDiff('centerTheta', h0, h).diff_data.identical
                    fTheta += h['fTheta'].data['count']
                    gThetaZ += csr_matrix(h['gThetaZ'].data, shape=gThetaZ.shape)
                    u = h['uThetaZZ'].data
                    uThetaZZ += csr_matrix((u['count'], (u['binTheta'],u['binZdZ'])),
                                           shape=uThetaZZ.shape)
                    uThetaZZe2 += csr_matrix((u['err2'], (u['binTheta'],u['binZdZ'])),
                                             shape=uThetaZZ.shape)
                    cputime += h['uThetaZZ'].header['cputime']

            self.hdus.append(h0["centerTheta"])
            self.hdus.append(h0["centerZ"])
            self.hdus.append(h0["pdfZ"])
            fTheta /= 2
            uThetaZZ /= 2
            uThetaZZe2 /= 2
            self.hdus.extend(self.fguHDU((fTheta, gThetaZ, uThetaZZ, uThetaZZe2)))
            self.hdus[-1].header['cputime'] = cputime
            self.writeToFile()
