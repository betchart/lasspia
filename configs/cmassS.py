import baofast
import numpy as np

class cmassS(baofast.configuration):

    def dataDir(self):
        """Directory of catalog files."""
        return '../data/'

    def outputLocation(self):
        return self.dataDir()

    def filesRandom(self):
        return [self.dataDir() + "randoms_DR9_CMASS_South.fits"]

    def filesObserved(self):
        return [self.dataDir() + "galaxies_DR9_CMASS_South.fits"]

    def catalogRandom(self): return baofast.wrapRandomSDSS(self.filesRandom())
    def catalogObserved(self): return baofast.wrapObservedSDSS(self.filesObserved())

    def binsZ(self): return np.arange(0.4, 0.72, 0.01)
    def binsRA(self): return np.arange(100, 270.001, 1./30) # FIXME periodic boundary conditions
    def binsDec(self): return np.arange(-10, 20, 0.1)
