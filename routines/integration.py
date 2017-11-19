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

    def RR(self):
        ft = None
        pdfz = None
        binningTheta = None
        binningZ = None

    def DR(self):
        gtz = None
        pdfz = None
        binningTheta = None
        binningZ = None

    def DD(self):
        utzz = None
        pdfz = None
        binningTheta = None
        binningZ = None


    def zIntegral(self):
        zCenters = self.getPre('centerz').data['binCenter']
        zz = zip(np.hstack([[0.],zCenters]), zCenters)
        dIz = [quad(self.integrand, z1, z2, args=self.config.omegasMKL())[0]
               for z1,z2 in zz]
        return np.cumsum(dIz)

    @staticmethod
    def integrand(z, omegaM, omegaK, omegaLambda):
        return math.sqrt(omegaM * (1+z)**3 +
                         omegaK * (1+z)**2 +
                         omegaLambda)

    @property
    def inputFileName(self):
        return self.config.stageFileName('combinatorial')

    def getPre(self, name):
        hdulist = fits.open(self.inputFileName)
        return hdulist[name]
