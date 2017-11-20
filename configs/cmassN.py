import baofast
import math

class cmassN(baofast.configuration):

    def dataDir(self):
        """Directory of catalog files."""
        return 'data/'

    def outputLocation(self):
        return self.dataDir()

    def inputFilesRandom(self):
        return [self.dataDir() + "randoms_DR9_CMASS_North.fits"]

    def inputFilesObserved(self):
        return [self.dataDir() + "galaxies_DR9_CMASS_North.fits"]

    def catalogRandom(self): return baofast.wrapRandomSDSS(self.inputFilesRandom())
    def catalogObserved(self): return baofast.wrapObservedSDSS(self.inputFilesObserved())

    def binningZ(self): return {"bins":900, "range":(0.43,0.7)}
    def binningRA(self): return {"bins":3500 , "range":(105,265)}
    def binningDec(self): return {"bins":1600, "range":(-4,57)}
    def binningTheta(self): return {"bins":3142, "range":(0,math.pi)}

    def chunkSize(self) : return 2000
