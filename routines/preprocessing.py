import baofast
import numpy as np
from astropy.io import fits

class preprocessing(baofast.routine):
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
        centers = np.array( baofast.utils.centers(edges),
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

        binning3D = self.config.binningDD([self.config.binningRA(),
                                           self.config.binningDec(),
                                           self.config.binningZ()])

        angR, xedges, yedges = np.histogram2d(ctlgR.ra, ctlgR.dec, weights=ctlgR.weightNoZ,
                                              **binning2D)

        angD, xedges, yedges = np.histogram2d(ctlgD.ra, ctlgD.dec, weights=ctlgD.weight,
                                              **binning2D)

        xx, yy = np.meshgrid(range(len(xedges)-1), range(len(yedges)-1), indexing='ij')

        ang = np.array(zip(xx.flat, yy.flat, angR.flat, angD.flat),
                       dtype = [("binRA", self.iType(len(xedges))),
                                ("binDec", self.iType(len(yedges))),
                                ("countR", self.iType(max(angR.flat))),
                                ("countD", np.float32)])

        mask = np.logical_or(ang["countR"]>0, ang["countD"]>0)
        hdu = fits.BinTableHDU(ang[mask], name="ang")
        hdu.header.add_comment("Unraveled angular (ra,dec) 2D histogram.")
        hdu.header.add_comment("Histogram for random catalog filled with z independent weights.")
        self.addProvenance(hdu,
                           self.config.inputFilesRandom() +
                           self.config.inputFilesObserved())

        triplets = np.array(zip(ctlgD.ra, ctlgD.dec, ctlgD.z))
        frq, edges = np.histogramdd(triplets,
                                    weights=ctlgD.weight,
                                    **binning3D)

        nx,ny,nz = frq.shape
        angzD = frq.reshape(nx*ny,nz)[mask]

        # Floating point FITS images (which have BITPIX = -32 or -64)
        # usually contain too much 'noise' in the least significant
        # bits of the mantissa of the pixel values to be effectively
        # compressed with any lossless algorithm. Consequently,
        # floating point images are first quantized into scaled
        # integer pixel values (and thus throwing away much of the
        # noise) before being compressed with the specified algorithm
        # (either GZIP, RICE, or HCOMPRESS)
        hdu2 = fits.CompImageHDU(angzD.astype(np.float32), name="angzD")
        hdu2.header.add_comment("3D histogram (ra, dec, z) of observed galaxies.")
        hdu2.header.add_comment("Unraveled in (ra,dec) to align with 'ANG', rows are z dimension.")
        self.addProvenance(hdu2, self.config.inputFilesObserved())

        return [hdu, hdu2]
