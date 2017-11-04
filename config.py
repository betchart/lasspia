import configSSDS
import numpy as np

"""Configuration: Stage 1 (Preprocessing)"""

def dataDir():
    """Data Directory"""
    return '../data/'

def filesRandom():
    dd = dataDir()
    return [dd + s for s in
            ['randoms_DR9_CMASS_North.fits',
             #'randoms_DR9_CMASS_South.fits'
            ]]

def filesObserved():
    dd = dataDir()
    return [dd + s for s in
            ['galaxies_DR9_CMASS_North.fits',
             #'galaxies_DR9_CMASS_South.fits'
            ]]

def catalogRandom(): return configSSDS.wrapRandomSDSS(filesRandom())
def catalogObserved(): return configSSDS.wrapObservedSDSS(filesObserved())

def binsZ():
    """Left bin edges and rightmost edge of P_z"""
    return np.arange(0.4, 0.701, 0.01) # DR9_CMASS

def binsRA():
    """Left bin edges and rightmost edge of R_ang(ra,dec): ra"""
    return np.arange(100, 270.001, 1./30) # DR9_CMASS_North

def binsDec():
    """Lower bin edges and uppermost edge of R_ang(ra,dec): dec"""
    return np.arange(-10, 60, 0.1) # DR9_CMASS_North
