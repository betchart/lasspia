import baofast
import numpy as np
from astropy.io import fits

class preprocessing(baofast.routine):
    """Preprocessing for fast 2-point correlations.

    Open input catalogs, create histograms, save to file.
    """
    def __call__(self):
        ctlg = self.config.catalogRandom()
        self.hdus.append( self.pdfZ(ctlg) )
        self.writeToFile()

    def pdfZ(self, ctlg):
        frq, edges = np.histogram(ctlg.z, self.config.binsZ(), # use scipy.stats.binned_statistic if want to keep mean value
                                  weights = ctlg.weight / sum(ctlg.weight))
        pdfz = np.array(zip(edges, frq),
                        dtype = [("lowEdge", np.float64),
                                 ("probability", np.float32)])
        hdu = fits.BinTableHDU(pdfz)
        hdu.header['name'] = ("pdfZ",
                              "Redshift probability density.")
        for i,f in enumerate(self.config.inputFilesRandom()):
            hdu.header['prov%d'%i] = (f.split('/')[-1],
                                      "File from which pdfZ was constructed.")
        return hdu

