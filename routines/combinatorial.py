import baofast
import math
import numpy as np
from astropy.io import fits

DEGTORAD = math.pi / 180

class combinatorial(baofast.routine):

    def __call__(self):
        binsDec = DEGTORAD * self.getPre("centerDec").data["binCenter"]
        binsRA = DEGTORAD * self.getPre("centerRA").data["binCenter"]
        self.rAng = self.getPre("RANG").data

        # calculate and histogram first row dTheta12, weighted with R1*R2

        self.sinDec = np.sin(binsDec)
        self.cosDec = np.cos(binsDec)
        self.prodSinDec = np.multiply.outer(self.sinDec,self.sinDec)
        self.prodCosDec = np.multiply.outer(self.cosDec,self.cosDec)
        self.cosDeltaRA = np.cos(np.subtract.outer(binsRA,binsRA))

        print len(self.rAng)
        print len(self.rAng)**2
        chunk = self.cosThetaChunk(slice(self.config.chunkSize()),
                                   slice(self.config.chunkSize()))
        print len(chunk.flat)
        print len(chunk.flat) / float(len(self.rAng)**2)

    def cosThetaChunk(self, slice1, slice2):
        return (self.cosDeltaRA[self.rAng["binRA"][slice1]][ :,self.rAng["binRA"][slice2]] *
                self.prodCosDec[self.rAng["binDec"][slice1]][:,self.rAng["binDec"][slice2]] +
                self.prodSinDec[self.rAng["binDec"][slice1]][:,self.rAng["binDec"][slice2]] )

    @property
    def inputFileName(self):
        return self.config.stageFileName('preprocessing')

    def getPre(self, name):
        hdulist = fits.open(self.inputFileName)
        return hdulist[name]
        
