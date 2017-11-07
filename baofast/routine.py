from astropy.io import fits

class routine(object):

    def __init__(self, config):
        self.config = config  # subclass of baofast.config
        self.hdus = [fits.PrimaryHDU()]

    @property
    def outputFileName(self):
        return self.config.stageFileName( self.__class__.__name__)

    def writeToFile(self):
        hdulist = fits.HDUList(self.hdus)
        hdulist.writeto(self.outputFileName, clobber=True) # clobber is overwrite in astropy v2
        print "Wrote %s" % self.outputFileName

    def __call__(self):
        '''Defined in subclasses.'''
        pass
