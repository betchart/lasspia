import lasspia as La
import math
from cmassS import cmassS

class cmassN(cmassS):

    def inputFilesRandom(self):
        return [self.dataDir() + "randoms_DR9_CMASS_North.fits"]

    def inputFilesObserved(self):
        return [self.dataDir() + "galaxies_DR9_CMASS_North.fits"]

    def catalogRandom(self):
        return La.wrapRandomSDSS(self.inputFilesRandom())

    def catalogObserved(self):
        return La.wrapObservedSDSS(self.inputFilesObserved())

    def binningRA(self): return {"bins":3500 , "range":(105,265)}
    def binningDec(self): return {"bins":1600, "range":(-4,57)}
