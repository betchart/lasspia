import baofast as bf
import numpy as np
from astropy.io import fits
from scipy.sparse import csr_matrix

class preprocessing(bf.routine):
    """Preprocessing for fast 2-point correlations.

    Open input catalogs, create histograms, save to file.
    """
    def __call__(self):
        ctlgR = self.config.catalogRandom()
        ctlgD = self.config.catalogObserved()
        self.hdus.append( self.binCenters(self.config.edgesZ(), "centerZ") )
        self.hdus.append( self.binCenters(self.config.edgesRA(), "centerRA") )
        self.hdus.append( self.binCenters(self.config.edgesDec(), "centerDec") )
        self.hdus.append( self.pdfZ(ctlgR) )
        self.hdus.extend( self.ang(ctlgR, ctlgD) )
        self.writeToFile()

    def binCenters(self, edges, name):
        centers = np.array( bf.utils.centers(edges),
                            dtype = [("binCenter", np.float64)])
        return fits.BinTableHDU(centers, name=name)

    @staticmethod
    def addProvenance(hdu, inputFiles):
        for i,f in enumerate(inputFiles):
            hdu.header['prov%d'%i] = (f.split('/')[-1],
                                      "Source data file.")

    @staticmethod
    def iType(iMax):
        return (np.int16 if iMax < np.iinfo(np.int16).max else
                np.int32 if iMax < np.iinfo(np.int32).max else
                np.int64)

    def pdfZ(self, ctlg):
        frq, edges = np.histogram(ctlg.z,
                                  weights = ctlg.weightZ / sum(ctlg.weightZ),
                                  **self.config.binningZ())

        pdfz = np.array(zip(edges, frq),
                        dtype = [("lowEdge", np.float64),
                                 ("probability", np.float32)])

        hdu = fits.BinTableHDU(pdfz, name="pdfZ")
        hdu.header.add_comment("Redshift probability histogram.")
        self.addProvenance(hdu, self.config.inputFilesRandom())

        return hdu

    def ang(self, ctlgR, ctlgD):
        '''HDU with angular binning for both random and observed catalogs.

        Another HDU with aligned angular binning records the z-bin and weight of the data.'''

        binning2D = self.config.binningDD([self.config.binningRA(),
                                           self.config.binningDec()])

        angR, xedges, yedges = np.histogram2d(ctlgR.ra, ctlgR.dec, weights=ctlgR.weightNoZ,
                                              **binning2D)

        angD, xedges, yedges = np.histogram2d(ctlgD.ra, ctlgD.dec, weights=ctlgD.weight,
                                              **binning2D)

        xx, yy = np.meshgrid(range(len(xedges)-1), range(len(yedges)-1), indexing='ij')


        mask = np.logical_or(angR>0, angD>0)
        hdu = fits.BinTableHDU.from_columns([
            fits.Column(name="binRA", array=xx[mask], format='I'),
            fits.Column(name="binDec",array=yy[mask], format='I'),
            fits.Column(name="countR",array=angR[mask], format='I'),
            fits.Column(name="countD",array=angD[mask], format='E')],
                                            name="ang")
        hdu.header.add_comment("Unraveled angular (ra,dec) 2D histogram.")
        hdu.header.add_comment("Histogram for random catalog filled with z independent weights.")
        self.addProvenance(hdu,
                           self.config.inputFilesRandom() +
                           self.config.inputFilesObserved())

        binsZ = self.config.binningZ()['bins']
        binsRA = self.config.binningRA()['bins']
        binsDec = self.config.binningDec()['bins']

        iZ = bf.utils.toBins(ctlgD.z, self.config.binningZ())
        iRA = bf.utils.toBins(ctlgD.ra, self.config.binningRA())
        iDec= bf.utils.toBins(ctlgD.dec, self.config.binningDec())
        iAng = binsDec*iRA  + iDec

        frq = csr_matrix((ctlgD.weight, (iAng, iZ)), shape=(binsRA*binsDec, binsZ))

        angzD = frq[mask.ravel()]
        iA, iZ = angzD.nonzero()
        hdu2 = fits.BinTableHDU.from_columns([
            fits.Column(name="iA", array=iA, format='J'),
            fits.Column(name="iZ", array=iZ, format='I'),
            fits.Column(name='count', array=angzD.data, format='E')],
                                             name="angzD")

        hdu2.header.add_comment("Sparse 3D histogram (ra, dec, z) of observed galaxies.")
        hdu2.header.add_comment("Unraveled in (ra,dec) and masked to align with 'ANG' rows.")
        self.addProvenance(hdu2, self.config.inputFilesObserved())

        return [hdu, hdu2]
