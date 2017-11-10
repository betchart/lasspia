from astropy.io import fits

class routine(object):

    def __init__(self, config, nJobs=None, iJob=None):
        self.config = config  # subclass of baofast.config
        self.nJobs = nJobs
        self.iJob = iJob
        self.hdus = [fits.PrimaryHDU()]

    def jobString(self):
        if self.iJob is not None:
            return "-%03d_%d" % (self.iJob, self.nJobs)
        return ""

    @property
    def outputFileName(self):
        return self.config.stageFileName( self.__class__.__name__) + self.jobString()

    def writeToFile(self):
        hdulist = fits.HDUList(self.hdus)
        hdulist.writeto(self.outputFileName, clobber=True) # clobber is overwrite in astropy v2
        print "Wrote %s" % self.outputFileName

    def __call__(self):
        '''Defined in subclasses.'''
        pass

    def combineJobOutput(self):
        pass
