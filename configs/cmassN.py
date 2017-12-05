import lasspia as La
import math

class cmassN(La.configuration):

    def dataDir(self):
        """Directory of catalog files."""
        return '../data/'

    def outputLocation(self):
        return self.dataDir()

    def inputFilesRandom(self):
        return [self.dataDir() + "randoms_DR9_CMASS_North.fits"]

    def inputFilesObserved(self):
        return [self.dataDir() + "galaxies_DR9_CMASS_North.fits"]

    def catalogRandom(self): return La.wrapRandomSDSS(self.inputFilesRandom())
    def catalogObserved(self): return La.wrapObservedSDSS(self.inputFilesObserved())

    def binningZ(self): return {"bins":900, "range":(0.43,0.7)}
    def binningRA(self): return {"bins":3500 , "range":(105,265)}
    def binningDec(self): return {"bins":1600, "range":(-4,57)}
    #def binningTheta(self): return {"bins":3142, "range":(0,math.pi)}
    def binningTheta(self):
        ra = self.maxDeltaRA() * math.pi/180
        dc = self.maxDeltaDec() * math.pi/180
        amx = math.acos(math.cos(2*ra) + math.sin(dc)*math.sin(-dc))
        res = 0.001
        return {"bins":int(amx/res), "range":(0,amx)}

    def chunkSize(self) : return 2000


    '''Parameters for avoiding unnecessary combinatorial calculations at large s.
    Galaxies farther apart than these parameters may not be included in result.'''

    def maxDeltaRA(self): return 25
    def maxDeltaDec(self): return 25
    def maxDeltaZ(self): return 0.1

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
