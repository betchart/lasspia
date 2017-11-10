import baofast
import math
import numpy as np
from astropy.io import fits
import sys

DEGTORAD = math.pi / 180

class combinatorial(baofast.routine):

    def __call__(self):
        binsDec = DEGTORAD * self.getPre("centerDec").data["binCenter"]
        binsRA = DEGTORAD * self.getPre("centerRA").data["binCenter"]
        self.rAng = self.getPre("RANG").data

        self.sinDec = np.sin(binsDec)
        self.cosDec = np.cos(binsDec)
        self.sinsinDec = np.multiply.outer(self.sinDec,self.sinDec)
        self.coscosDec = np.multiply.outer(self.cosDec,self.cosDec)
        self.cosDeltaRA = np.cos(np.subtract.outer(binsRA,binsRA))

        splits = range(0, len(self.rAng), self.config.chunkSize())
        slices = [slice(i,j) for i,j in zip(splits,splits[1:]+[None])]

        # Full arrays of indices can work in place of slices.
        # A slight penalty is incurred, but makes regioning possible.
        # Regioning would be important in the case of a full-sky dataset
        # slices = [range(len(self.rAng))[slice(i,j)]
        #           for i,j in zip(splits,splits[1:]+[None])]

        typeRR = np.int64
        RR = np.zeros(len(self.config.edgesTheta())-1, dtype=typeRR)

        chunks = [(slices[i],jSlice) for i in range(len(slices)) for jSlice in slices[i:]]
        trunc = 47
        print "There are %d chunks" % len(chunks)
        if trunc: print "Truncating after %d chunks" % trunc

        for slice1,slice2 in chunks[:trunc]:
            chunkT = self.thetaChunk(slice1, slice2)
            countcount = np.multiply.outer(self.rAng["count"][slice1],
                                           self.rAng["count"][slice2]).astype(typeRR)
            if slice1 != slice2: countcount *= 2 # fill histogram with twice-weights
            frq,_ = np.histogram( chunkT, weights = countcount,
                                  **self.config.binningTheta())
            RR += frq
            print '.',
            sys.stdout.flush()
        print
        RR /= 2
        outFile = open("points.txt","w")
        for i,j in zip(self.config.edgesTheta(), RR): print>>outFile, i,j
                
    def thetaChunk(self, slice1, slice2):
        cosTheta = (
            self.cosDeltaRA[self.rAng["binRA"][slice1]][ :,self.rAng["binRA"][slice2]] *
            self.coscosDec[self.rAng["binDec"][slice1]][:,self.rAng["binDec"][slice2]] +
            self.sinsinDec[self.rAng["binDec"][slice1]][:,self.rAng["binDec"][slice2]])
        np.clip(cosTheta,-1,1,out=cosTheta)
        return np.arccos(cosTheta,out=cosTheta) # clip and arccos add 30% time cost compared to returning cosTheta

    @property
    def inputFileName(self):
        return self.config.stageFileName('preprocessing')

    def getPre(self, name):
        hdulist = fits.open(self.inputFileName)
        return hdulist[name]
        
