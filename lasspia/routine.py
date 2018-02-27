from __future__ import print_function
from astropy.io import fits
import sys
import os

class routine(object):

    def __init__(self, config, nJobs=None, iJob=None):
        self.config = config  # subclass of lasspia.config
        self.nJobs = nJobs
        self.iJob = iJob
        self.hdus = [fits.PrimaryHDU()]
        self.hdus[0].header['cmd'] = ' '.join(sys.argv)
        self.streamFile = (None if not config.txtToFile else
                           self.outputFileName.replace('fits','txt'))
        if self.streamFile and os.path.exists(self.streamFile):
            os.remove(self.streamFile)
        config.checkConsistency(self.out)
        config.info(self.out)
        self.closeStreams()

    def destruct(self):
        self.closeStreams()
        del self.hdus[:]
        return

    def jobString(self, iJob=None):
        if iJob is None:
            iJob = self.iJob
        if iJob is not None:
            return "-%03d_%d" % (iJob, self.nJobs)
        return ""

    @property
    def outputFileName(self):
        return self.config.stageFileName( self.__class__.__name__) + self.jobString()

    @property
    def out(self):
        if not hasattr(self, '_ostream') or self._ostream.closed:
            self._ostream = sys.stdout if not self.streamFile else open(self.streamFile, 'a')
        return self._ostream

    def closeStreams(self):
        if hasattr(self, '_ostream'):
            if self._ostream == sys.stdout:
                self._ostream = None
            elif self._ostream:
                self._ostream.close()
            del self._ostream
        return

    def writeToFile(self):
        hdulist = fits.HDUList(self.hdus)
        if os.path.exists(self.outputFileName):
            os.remove(self.outputFileName)
        hdulist.writeto(self.outputFileName)
        print("Wrote %s" % self.outputFileName, file=self.out)
        self.destruct()

    def showFitsHeaders(self):
        if not os.path.exists(self.outputFileName):
            print('Not found:', self.outputFileName, file=self.out)
            print('Perhaps you need to first create',
                  'it by running the routine.',
                  file = self.out)
            return

        with fits.open(self.outputFileName) as hdus:
            hdus.info(self.out)
            for h in hdus:
                print(file=self.out)
                print(repr(h.header), file=self.out)
            print(file=self.out)
        return

    def __call__(self):
        '''Defined in subclasses.'''
        pass

    def combineOutput(self):
        '''Defined in subclasses for paralellizable routines.'''
        raise Exception("The class '%s(lasspia.routine)' does not define 'combineOutput'." % self.__class__.__name__)

    def combineOutputZ(self):
        '''Defined in subclasses for routines which can combine output from multiple slices in z: see zSlicing.'''
        raise Exception('\n'.join(["The class '%s(lasspia.routine)' does not define 'combineOutputZ'."
                                   % self.__class__.__name__,
                                   "Perhaps you want to specify --iSliceZ"]))

    def plot(self):
        '''Defined in subclasses, to be run when lasspia.py receives the --plot flag.'''
        print( 'No plot() method defined for %s.' % self.__class__.__name__)
