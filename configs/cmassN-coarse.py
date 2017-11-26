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

    def binningZ(self): return {"bins":300, "range":(0.43,0.7)}
    def binningRA(self): return {"bins":1200 , "range":(105,265)}
    def binningDec(self): return {"bins":600, "range":(-4,57)}
    def binningTheta(self): return {"bins":1000, "range":(0,math.pi)}

    def chunkSize(self) : return 2000


    '''Configuration affecting only the "integration" routine.'''

    def omegasMKL(self):
        '''Cosmology parameters (\Omega_M, \Omega_K, \Omega_\Lambda).'''
        return (0.274, 0, 0.726)

    def H0(self):
        '''Hubble constant in (h km/s / Mpc)'''
        return 100.

    def lightspeed(self):
        '''Speed of light in km/s'''
        return 299792

    def binningS(self): return {"bins":1200, "range":(0,6000)}
