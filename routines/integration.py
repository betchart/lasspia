import math
import numpy as np
import lasspia as La
from astropy.io import fits
from scipy.integrate import quad
from lasspia.timing import timedHDU

class integration(La.routine):

    def __call__(self):
        self.hdus.append(self.tpcf())
        self.writeToFile()
        return

    @timedHDU
    def tpcf(self):
        self.pdfz = self.getInput('pdfZ').data['probability']

        binsS = self.config.binningS()['bins']

        binsTh = len(self.getInput('centertheta').data)
        slcT = ( slice(None) if self.iJob is None else
                 La.utils.slices(binsTh,
                                 binsTh/self.nJobs)[self.iJob] )

        s = np.sqrt(sum(np.power(a,2) for a in self.sigmaPiGrids(slcT)))
        hdu = fits.BinTableHDU.from_columns([
            fits.Column(name='s', array=self.centersS(), format='E'),
            fits.Column(name='RR', array=self.calcRR(s, slcT), format='D'),
            fits.Column(name='DR', array=self.calcDR(s, slcT), format='D'),
            fits.Column(name='DD', array=self.calcDD(s, slcT), format='D')],
                                            name="TPCF")
        hdu.header.add_comment("Two-point correlation function for pairs of galaxies,"+
                               " by distance s.")
        return hdu

    def sigmaPiGrids(self, slcT):
        '''A cubic grid of sigma (pi) values
        for pairs of galaxies with coordinates (iTheta, iZ1, iZ2).'''
        Iz = self.zIntegral()
        rOfZ = Iz * (self.config.lightspeed()/self.config.H0())
        tOfZ = rOfZ * (1 + self.config.omegasMKL()[1]/6 * Iz**2)

        thetas = self.getInput('centertheta').data['binCenter'][slcT]
        sinT2 = np.sin(thetas/2)
        cosT2 = np.cos(thetas/2)

        sigmas = sinT2[:,None,None] * (tOfZ[None,:,None] + tOfZ[None,None,:])
        pis = cosT2[:,None,None] * (rOfZ[None,:,None] - rOfZ[None,None,:])
        return sigmas, pis

    def calcRR(self,s, slcT):
        ft = self.getInput('fTheta').data['count'][slcT]
        counts = ft[:,None,None] * self.pdfz[None,:,None] * self.pdfz[None,None,:]
        rr = np.histogram(s, weights=counts, **self.config.binningS())[0]
        if self.iJob is None:
            rr /= sum(rr)
        del counts
        return rr

    def calcDR(self,s, slcT):
        gtz = self.getInput('gThetaZ').data
        counts = gtz[slcT,:,None] * self.pdfz[None,None,:]
        dr = np.histogram(s, weights=counts, **self.config.binningS())[0]
        if self.iJob is None:
            dr /= sum(dr)
        del counts
        return dr

    def calcDD(self,s, slcT):
        utzz = self.getInput('uThetaZZ').data
        overflow = utzz['binZdZ'][-1]+1 == s.shape[1]**2
        slc = slice(-1 if overflow else None)

        iThetas = utzz['binTheta'][slc]
        mask = (slice(None) if slcT==slice(None) else
                np.logical_and(slcT.start <= iThetas, iThetas < slcT.stop))

        counts = utzz['count'][slc][mask]
        iTh = iThetas[mask] - (slcT.start or 0)
        iZdZ = utzz['binZdZ'][slc][mask]
        iZ = iZdZ / s.shape[1]
        diZ = iZdZ % s.shape[1]
        iZ2 = iZ + diZ

        dd = np.histogram(s[iTh,iZ,iZ2], weights=counts, **self.config.binningS())[0]
        if overflow and self.iJob in [0,None]:
            dd[-1] = dd[-1] + utzz['count'][-1]
        if self.iJob is None:
            dd /= sum(dd)
        return dd

    def centersS(self):
        return La.utils.centers(self.config.edgesFromBinning(self.config.binningS()))

    def zIntegral(self):
        zCenters = self.getInput('centerz').data['binCenter']
        zz = zip(np.hstack([[0.],zCenters]), zCenters)
        dIz = [quad(self.integrand, z1, z2, args=self.config.omegasMKL())[0]
               for z1,z2 in zz]
        return np.cumsum(dIz)

    @staticmethod
    def integrand(z, omegaM, omegaK, omegaLambda):
        return 1./math.sqrt(omegaM * (1+z)**3 +
                            omegaK * (1+z)**2 +
                            omegaLambda)

    @property
    def inputFileName(self):
        return self.config.stageFileName('combinatorial')

    def getInput(self, name):
        hdulist = fits.open(self.inputFileName)
        return hdulist[name]


    def combineOutput(self):
        jobFiles = [self.outputFileName + self.jobString(iJob)
                    for iJob in range(self.nJobs)]

        with fits.open(jobFiles[0]) as h0:
            hdu = h0['TPCF']
            cputime = hdu.header['cputime']
            for jF in jobFiles[1:]:
                with fits.open(jF) as jfh:
                    assert np.all( hdu.data['s'] == jfh['TPCF'].data['s'])
                    cputime += jfh['TPCF'].header['cputime']
                    for col in ['RR','DR','DD']:
                        hdu.data[col] += jfh['TPCF'].data[col]
            for col in ['RR','DR','DD']:
                hdu.data[col] /= sum(hdu.data[col])
            hdu.header['cputime'] = cputime
            self.hdus.append(hdu)
            self.writeToFile()
        return
