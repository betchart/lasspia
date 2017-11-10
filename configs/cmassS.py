import baofast
import math

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

    def binningZ(self): return {"bins":300, "range":(0.4,0.7)}
    def binningRA(self): return {"bins": 3000, "range":(-50,50)}
    def binningDec(self): return {"bins":300, "range":(-10,20)}
    def binningTheta(self): return {"bins":1000, "range":(0,math.pi)}
    
    def chunkSize(self) : return 7500
