from __future__ import print_function
import sys
import math
import numpy as np
from lasspia import utils
from lasspia import slicing

class configuration(object):

    def outputLocation(self):
        return "."

    def inputFilesRandom(self):
        '''List of random catalog file names.'''
        pass
    
    def inputFilesObserved(self):
        '''List of observed catalog file names.'''
        pass

    def binningZ(self): pass
    def binningRA(self): pass
    def binningDec(self): pass
    def binningTheta(self): pass

    def chunkSize(self): pass

    '''Parameters for avoiding unnecessary combinatorial calculations at large s.
    Galaxies farther apart than these parameters may not be included in result.'''
    def maxDeltaRA(self): return None
    def maxDeltaDec(self): return None
    def maxDeltaZ(self): return None
    def regionBasedCombinations(self): return False

    def edgesZ(self): return self.edgesFromBinning(self.binningZ())
    def edgesRA(self): return self.edgesFromBinning(self.binningRA())
    def edgesDec(self): return self.edgesFromBinning(self.binningDec())
    def edgesTheta(self): return self.edgesFromBinning(self.binningTheta())

    def binningS(self): pass

    def __init__(self, txtToFile=False):
        self.txtToFile = txtToFile

    @property
    def name(self) : return self.__class__.__name__
        
    def stageFileName(self, stage):
        return '/'.join([self.outputLocation().rstrip('/'),
                         "%s_%s.fits" % (self.name, str(stage))])

    @staticmethod
    def binningDD(binnings):
        binning = {"bins": tuple([b['bins'] for b in binnings])}
        if "range" in binnings[0]:
            binning["range"] = [b["range"] for b in binnings]
        return binning

    @staticmethod
    def edgesFromBinning(binning):
        _,edges = np.histogram([], **binning)
        return edges

    def checkConsistency(self, output=None):
        self.checkThetaRange(output=output)


    def checkThetaRange(self, output=None):
        thetaLo, thetaHi = self.binningTheta()['range']

        if 0 < thetaLo:
            text = ["Warning: %s.binningTheta() defines" % self.name,
                    "a range minimum greater than zero,",
                    "which will cause routines/combinatorial.py to crash."]
            print("\n         ".join(text), "\n", file=output)

        def deltaRA():
            raLo, raHi = [ra/180. for ra in self.binningRA()['range']]
            if self.maxDeltaRA():
                aMax = self.maxDeltaRA()
                return min(2*aMax, raHi-raLo)
            return raHi-raLo

        def decLOHI():
            dLo, dHi = [dec/180. for dec in self.binningDec()['range']]
            if self.maxDeltaDec():
                dMax = self.maxDeltaDec()
                return (dLo,dHi) if (dHi-dLo) < 2*dMax else (-dMax,dMax)

        def maxTheta():
            dLo,dHi = decLOHI()
            ss = math.sin(dLo)*math.sin(dHi)
            cc = math.cos(dLo)*math.cos(dHi)
            return math.acos(math.cos(deltaRA()) * cc + ss)

        maxT = maxTheta()
        if thetaHi < maxT:
            text = ["Warning: %s.maxDeltaRA() and %s.maxDeltaDec()" % (self.name, self.name),
                    "imply possible values of theta up to %f" % maxT,
                    "but %s.binningTheta() defines a range maximum of only %f," % (self.name, thetaHi),
                    "which may cause routines/combinatorial.py to crash."]
            print("\n         ".join(text), "\n", file=output)
