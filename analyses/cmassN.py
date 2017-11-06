import baofast
import numpy as np

class cmassN(baofast.analysis):

    def dataDir(self):
        """Directory of catalog files."""
        return '../data/'

    def filesRandom(self):
        return [self.dataDir() + "randoms_DR9_CMASS_North.fits"]

    def filesObserved(self):
        return [self.dataDir() + "galaxies_DR9_CMASS_North.fits"]

    def catalogRandom(self): return baofast.wrapRandomSDSS(self.filesRandom())
    def catalogObserved(self): return baofast.wrapObservedSDSS(self.filesObserved())

    def binsZ(self): return np.arange(0.4, 0.72, 0.01)
    def binsRA(self): return np.arange(100, 270.001, 1./30)
    def binsDec(self): return np.arange(-10, 60, 0.1)
