from astropy.io import fits
import numpy as np

def openSDSS(files):
    """Return galaxy catalog from Sloan Digital Sky Survey files."""
    return np.concatenate([
        fits.open(f)[1].data
        for f in files
    ])

class wrapSDSS(object):
    """A catalog wrapper interface for SDSS.

    Wrappers must provide the following properties:
    * z
    * weight
    """
    def __init__(self, files):
        self.ctlg = openSDSS(files)

    @property
    def z(self): return self.ctlg['z']


class wrapRandomSDSS(wrapSDSS):
    def __init__(self, files):
        super(wrapRandomSDSS, self).__init__(files)

    @property
    def weight(self): return self.ctlg['weight']
        
class wrapObservedSDSS(wrapSDSS):
    def __init__(self, files):
        super(wrapObservedSDSS, self).__init__(files)

    @property
    def weight(self):
        ct = self.ctlg
        return (ct['weight_sdc'] * ct['weight_fkp'] *
                (ct['weight_noz'] + ct['weight_cp'] - 1))
