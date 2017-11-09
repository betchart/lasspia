import baofast
import math
import numpy as np
from astropy.io import fits

DEGTORAD = math.pi / 180

class combinatorial(baofast.routine):

    def __call__(self):
        binsDec = DEGTORAD * self.getPre("centerDec").data["binCenter"]
        binsRA = DEGTORAD * self.getPre("centerRA").data["binCenter"]
        rAng = self.getPre("RANG").data

        # calculate and histogram first row dTheta12, weighted with R1*R2

        self.sinDec = np.sin(binsDec)
        self.cosDec = np.cos(binsDec)
        self.prodSinDec = np.multiply.outer(self.sinDec,self.sinDec)
        self.prodCosDec = np.multiply.outer(self.cosDec,self.cosDec)
        self.cosDeltaRA = np.cos(np.subtract.outer(binsRA,binsRA))

        for i in range(len(rAng)):
            cT = self.cosThetaRow(i, rAng)
            if not (len(cT)%100): print len(cT), cT[:6], cT[-6:]

    def cosThetaRow(self, i, rAng):
        return (self.cosDeltaRA[rAng["binRA"][i], rAng["binRA"][i:]] *
                self.prodCosDec[rAng["binDec"][i], rAng["binDec"][i:]] +
                self.prodSinDec[rAng["binDec"][i], rAng["binDec"][i:]] )

    @property
    def inputFileName(self):
        return self.config.stageFileName('preprocessing')

    def getPre(self, name):
        hdulist = fits.open(self.inputFileName)
        return hdulist[name]
        
