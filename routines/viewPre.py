import baofast
from astropy.io import fits

class viewPre(baofast.routine):

    def __call__(self):
        self.openPre()
    
    @property
    def inputFileName(self):
        return self.config.stageFileName('preprocessing')
    
    def openPre(self):
        hdulist = fits.open(self.inputFileName)
        hdulist.info()
        for hdu in hdulist:
            print
            print repr(hdu.header)
        print
        return
