import baofast
import math
from astropy.io import fits
import numpy as np
from scipy.integrate import quad

class integration(baofast.routine):

    def __call__(self):
        Iz = self.zIntegral()
        rOfZ = Iz * (self.config.lightspeed()/self.config.H0())
        tOfZ = rOfZ * (1 + self.config.omegasMKL()[1]/6 * Iz**2)
        zCenters = self.getInput('centerZ').data['binCenter']
        thetas = self.getInput('centertheta').data['binCenter']
        sinT2 = np.sin(thetas/2)
        cosT2 = np.cos(thetas/2)

        ft = self.getInput('fTheta').data['count']
        gtz = self.getInput('gThetaZ').data
        utzz = self.getInput('uThetaZZ').data
        pdfz = self.getInput('pdfZ').data['probability']

        sigmas = sinT2[:,None,None] * (tOfZ[None,:,None] + tOfZ[None,None,:])
        pis = cosT2[:,None,None] * (rOfZ[None,:,None] - rOfZ[None,None,:])
        s = np.sqrt(sigmas**2 + pis**2)

        def calcRR():
            counts = ft[:,None,None] * pdfz[None,:,None] * pdfz[None,None,:]
            rr = np.histogram(s, weights=counts, **self.config.binningS())[0]
            rr /= np.sum(rr)
            del counts
            return rr

        def calcDR():
            counts = gtz[:,:,None] * pdfz[None,None,:]
            dr = np.histogram(s, weights=counts, **self.config.binningS())[0]
            dr /= np.sum(dr)
            del counts
            return dr

        def calcDD():
            counts = utzz['count']
            iThetas = utzz['binTheta']
            iZdZ = utzz['binZdZ']
            iZ = iZdZ / len(zCenters)
            diZ = iZdZ % len(zCenters)
            iZ2 = iZ + diZ
            dd = np.histogram(s[iThetas,iZ,iZ2], weights=counts, **self.config.binningS())[0]
            dd /= np.sum(dd)
            return dd

        centers = baofast.utils.centers(self.config.edgesFromBinning(self.config.binningS()))
        RR = calcRR()
        DR = calcDR()
        DD = calcDD()
        for x,y,z,a in zip(centers,RR,DR,DD): print x,y,z,a
        return

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
