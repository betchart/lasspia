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
        self.hdus.append( self.binCentersRA() )
        self.hdus.append( self.binCentersDec() )
        self.hdus.append( self.Rang(ctlg) )
        self.writeToFile()

    @staticmethod
    def addProvenance(hdu, inputFiles):
        for i,f in enumerate(inputFiles):
            hdu.header['prov%d'%i] = (f.split('/')[-1],
                                      "Source data file.")

    @staticmethod
    def centers(leftEdges):
        return leftEdges[:-1] + 0.5 * (leftEdges[1:] - leftEdges[:-1])

    @staticmethod
    def iType(iMax):
        return (np.int16 if iMax < np.iinfo(np.int16).max else
                np.int32 if iMax < np.iinfo(np.int32).max else
                np.int64)

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

    def binCentersRA(self):
        centers = np.array( self.centers(self.config.binsRA()),
                            dtype = [("binCenter", np.float64)])
        hdu = fits.BinTableHDU(centers, name="centerRA")
        return hdu

    def binCentersDec(self):
        centers = np.array( self.centers(self.config.binsDec()),
                            dtype = [("binCenter", np.float64)])
        hdu = fits.BinTableHDU(centers, name="centerDec")
        return hdu

    def Rang(self, ctlg):
        frq, xedges, yedges = np.histogram2d(ctlg.ra, ctlg.dec,
                                             [self.config.binsRA(),
                                              self.config.binsDec()])
        xx, yy = np.meshgrid(range(len(xedges)-1), range(len(yedges)-1))

        rang = np.array(zip(xx.flat, yy.flat, frq.T.flat),
                        dtype = [("binRA", self.iType(len(xedges))),
                                 ("binDec", self.iType(len(yedges))),
                                 ("count", self.iType(max(frq.flat)))])

        hdu = fits.BinTableHDU(rang[rang["count"]>0], name="Rang")
        #hdu.header['name'] = ("Rang", "Sparse 2D angular distribution, random catalog.")
        self.addProvenance(hdu, self.config.inputFilesRandom())

        return hdu

