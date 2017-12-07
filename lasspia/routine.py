from astropy.io import fits
import sys
import os

class routine(object):

    def __init__(self, config, nJobs=None, iJob=None):
        self.config = config  # subclass of lasspia.config
        self.nJobs = nJobs
        self.iJob = iJob
        self.hdus = [fits.PrimaryHDU()]
        self.outstream = (sys.stdout if not config.txtToFile else
                          open(self.outputFileName.replace('fits','txt'),'w'))

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
        print >> self.outstream, "Writing %s" % self.outputFileName
        hdulist = fits.HDUList(self.hdus)
        if os.path.exists(self.outputFileName):
            os.remove(self.outputFileName)
        for iTry in range(self.config.nWriteAttempts):
            try:
                hdulist.writeto(self.outputFileName)
            except IOError as e:
                print >> self.outstream, e
            if os.path.exists(self.outputFileName):
                print >> self.outstream, "Wrote %s" % self.outputFileName
                break
            else:
                print >> self.outstream, "Attempt %d failed: hdulist.writeto(%s)" % (iTry, self.outputFileName)
                from time import sleep
                sleep(10)

    def showFitsHeaders(self):
        if not os.path.exists(self.outputFileName):
            print >> self.outstream, 'Not found:', self.outputFileName
            print >> self.outstream, 'Perhaps you need to first create',
            print >> self.outstream, 'it by running the routine.'
            return

        with fits.open(self.outputFileName) as hdus:
            hdus.info(self.outstream)
            for h in hdus[1:]:
                print >> self.outstream
                print >> self.outstream, repr(h.header)
            print >> self.outstream
        return

    def __call__(self):
        '''Defined in subclasses.'''
        pass

    def combineOutput(self):
        '''Defined in subclasses for paralellizable routines.'''
        raise Exception("The class '%s(lasspia.routine)' does not define 'combineOutput'." % self.__class__.__name__)
