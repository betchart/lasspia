import baofast as bf
import numpy as np
from astropy.io import fits
from scipy.sparse import csr_matrix

class xcheckPreprocessing(bf.routine):

    def __call__(self):

        def getAngZD(shape):
            a = self.getInput("ANGZD").data
            return csr_matrix((a['count'], (a['iAlign'], a['iZ'])),
                              shape=shape)

        zCenter = self.getInput("centerz").data['binCenter']
        ang = self.getInput("ANG").data
        angzd = getAngZD((len(ang),len(zCenter)))
        assert len(ang['countD'].nonzero()[0]) == len(angzd.sum(axis=1).nonzero()[0])
        assert 1e-6 > abs(angzd.sum() - ang['countD'].sum()) / ang['countD'].sum()
        print 'Cross-checks pass!'

    @property
    def inputFileName(self):
        return self.config.stageFileName('preprocessing')

    def getInput(self, name):
        hdulist = fits.open(self.inputFileName)
        return hdulist[name]
