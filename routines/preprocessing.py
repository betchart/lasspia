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
        self.hdus.append( self.Rang(ctlg) )
        self.writeToFile()

    @staticmethod
    def addProvenance(hdu, inputFiles):
        for i,f in enumerate(inputFiles):
            hdu.header['prov%d'%i] = (f.split('/')[-1],
                                      "Source data file.")

    def pdfZ(self, ctlg):
        frq, edges = np.histogram(ctlg.z, self.config.binsZ(),
                                  weights = ctlg.weight / sum(ctlg.weight))

        pdfz = np.array(zip(edges, frq),
                        dtype = [("lowEdge", np.float64),
                                 ("probability", np.float32)])

        hdu = fits.BinTableHDU(pdfz, name="pdfZ")
        #hdu.header['name'] = ("pdfZ", "Redshift probability density.")
        self.addProvenance(hdu, self.config.inputFilesRandom())

        return hdu

    def Rang(self, ctlg):

        frq, xedges, yedges = np.histogram2d(ctlg.ra, ctlg.dec,
                                             [self.config.binsRA(),
                                              self.config.binsDec()])
        xx, yy = np.meshgrid(xedges, yedges)

        rang = np.array(zip(xx.flat, yy.flat, frq.flat),
                        dtype = [("lowEdgeRA", np.float64),
                                 ("lowEdgeDec", np.float64),
                                 ("count", np.float32)])

        hdu = fits.BinTableHDU(rang[frq.flat>0], name="Rang")
        #hdu.header['name'] = ("Rang", "2D angular distribution, random catalog.")
        self.addProvenance(hdu, self.config.inputFilesRandom())

        return hdu

# use scipy.stats.binned_statistic if want to keep mean value

# consider storing bin centers rather than edges
