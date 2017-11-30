from astropy.io import fits
import os

class routine(object):

    def __init__(self, config, nJobs=None, iJob=None):
        self.config = config  # subclass of lasspia.config
        self.nJobs = nJobs
        self.iJob = iJob
        self.hdus = [fits.PrimaryHDU()]

    def jobString(self, iJob=None):
        if iJob is None:
            iJob = self.iJob
        if iJob is not None:
            return "-%03d_%d" % (iJob, self.nJobs)
        return ""

    @property
    def outputFileName(self):
        return self.config.stageFileName( self.__class__.__name__) + self.jobString()

    def writeToFile(self):
        hdulist = fits.HDUList(self.hdus)
        hdulist.writeto(self.outputFileName, clobber=True) # clobber is overwrite in astropy v2
        print "Wrote %s" % self.outputFileName

    def showFitsHeaders(self):
        if not os.path.exists(self.outputFileName):
            print 'Not found:', self.outputFileName
            print 'Perhaps you need to first create',
            print 'it by running the routine.'
            return

        with fits.open(self.outputFileName) as hdus:
            hdus.info()
            for h in hdus[1:]:
                print
                print repr(h.header)
            print
        return

    def __call__(self):
        '''Defined in subclasses.'''
        pass

    def combineOutput(self):
        '''Defined in subclasses for paralellizable routines.'''
        raise Exception("The class '%s(lasspia.routine)' does not define 'combineOutput'." % self.__class__.__name__)