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

    def binningZ(self): return {"bins":300, "range":(0.4,0.7)}
    def binningRA(self): return {"bins":5100 , "range":(100,270)}
    def binningDec(self): return {"bins":700, "range":(-10,60)}
    def binningTheta(self): return {"bins":1000, "range":(0,math.pi)}

    def chunkSize(self) : return 3000
