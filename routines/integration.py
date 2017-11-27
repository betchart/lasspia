import math
import numpy as np
import baofast as bf
from astropy.io import fits
from scipy.integrate import quad

class integration(bf.routine):

    def __call__(self):
        self.pdfz = self.getInput('pdfZ').data['probability']
        s = np.sqrt(sum(np.power(a,2) for a in self.sigmaPiGrids()))

        hdu = fits.BinTableHDU.from_columns([
            fits.Column(name='s', array=self.centersS(), format='E'),
            fits.Column(name='RR', array=self.calcRR(s), format='E'),
            fits.Column(name='DR', array=self.calcDR(s), format='E'),
            fits.Column(name='DD', array=self.calcDD(s), format='E')],
                                            name="TPCF")
        hdu.header.add_comment("Two-point correlation function for pairs of galaxies,"+
                               " by distance s.")
        self.hdus.append(hdu)
        self.writeToFile()
        return

    def sigmaPiGrids(self):
        '''A cubic grid of sigma (pi) values
        for pairs of galaxies with coordinates (iTheta, iZ1, iZ2).'''
        Iz = self.zIntegral()
        rOfZ = Iz * (self.config.lightspeed()/self.config.H0())
        tOfZ = rOfZ * (1 + self.config.omegasMKL()[1]/6 * Iz**2)

        thetas = self.getInput('centertheta').data['binCenter']
        sinT2 = np.sin(thetas/2)
        cosT2 = np.cos(thetas/2)

        sigmas = sinT2[:,None,None] * (tOfZ[None,:,None] + tOfZ[None,None,:])
        pis = cosT2[:,None,None] * (rOfZ[None,:,None] - rOfZ[None,None,:])
        return sigmas, pis

    def calcRR(self,s):
        ft = self.getInput('fTheta').data['count']
        counts = ft[:,None,None] * self.pdfz[None,:,None] * self.pdfz[None,None,:]
        rr = np.histogram(s, weights=counts, **self.config.binningS())[0]
        rr /= np.sum(rr)
        del counts
        return rr

    def calcDR(self,s):
        gtz = self.getInput('gThetaZ').data
        counts = gtz[:,:,None] * self.pdfz[None,None,:]
        dr = np.histogram(s, weights=counts, **self.config.binningS())[0]
        dr /= np.sum(dr)
        del counts
        return dr

    def calcDD(self,s):
        utzz = self.getInput('uThetaZZ').data
        counts = utzz['count']
        iThetas = utzz['binTheta']
        iZdZ = utzz['binZdZ']
        iZ = iZdZ / s.shape[1]
        diZ = iZdZ % s.shape[1]
        iZ2 = iZ + diZ
        dd = np.histogram(s[iThetas,iZ,iZ2], weights=counts, **self.config.binningS())[0]
        dd /= np.sum(dd)
        return dd

    def centersS(self):
        return bf.utils.centers(self.config.edgesFromBinning(self.config.binningS()))

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
