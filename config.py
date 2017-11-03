"""Configuration: Stage 1 (Preprocessing)"""

import configSSDS

def dataDir():
    """Data Directory"""
    return '../data/'

def filesRandom():
    dd = dataDir()
    return [dd + s for s in
            ['randoms_DR9_CMASS_North.fits',
             'randoms_DR9_CMASS_South.fits'
            ]]

def filesObserved():
    dd = dataDir()
    return [dd + s for s in
            ['galaxies_DR9_CMASS_North.fits',
             'galaxies_DR9_CMASS_South.fits'
            ]]


def catalogRandom(): return configSSDS.wrapRandomSDSS(filesRandom())
def catalogObserved(): return configSSDS.wrapObservedSDSS(filesObserved())

