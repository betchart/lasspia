import baofast
import numpy as np

class cmassS(baofast.configuration):

    def dataDir(self):
        """Directory of catalog files."""
        return 'data/'

    def outputLocation(self):
        return self.dataDir()

    def inputFilesRandom(self):
        return [self.dataDir() + "randoms_DR9_CMASS_South.fits"]

    def inputFilesObserved(self):
        return [self.dataDir() + "galaxies_DR9_CMASS_South.fits"]

    def catalogRandom(self):
        return baofast.wrapRandomSDSS(self.inputFilesRandom(), shiftRA=True)

    def catalogObserved(self):
        return baofast.wrapObservedSDSS(self.inputFilesObserved(), shiftRA=True)

    def binsZ(self): return np.arange(0.4, 0.72, 0.01)
    def binsRA(self): return np.arange(-50, 50, 1./30)
    def binsDec(self): return np.arange(-10, 20, 0.1)

    def chunkSize(self) : return 7000
