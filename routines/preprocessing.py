import lasspia as La
import numpy as np
from astropy.io import fits
from scipy.sparse import csr_matrix
from lasspia.timing import timedHDU

class preprocessing(La.routine):
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

    @timedHDU
    def binCenters(self, edges, name):
        centers = np.array( La.utils.centers(edges),
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

    @timedHDU
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

    @timedHDU
    def ang(self, ctlgR, ctlgD):
        '''HDU with angular binning for both random and observed catalogs.

        Another HDU with aligned angular binning records the z-bin and weight of the data.'''

        binning2D = self.config.binningDD([self.config.binningRA(),
                                           self.config.binningDec()])

        angR = np.histogram2d(ctlgR.ra, ctlgR.dec, weights=ctlgR.weightNoZ,
                              **binning2D)[0]

        angD = np.histogram2d(ctlgD.ra, ctlgD.dec, weights=ctlgD.weight,
                              **binning2D)[0]

        slicePoints, iXs, iYs = La.utils.chunksWhere(np.logical_or(angR>0, angD>0),
                                                     1.3*self.config.chunkSize())

        hduSlc = fits.BinTableHDU(np.array( slicePoints,
                                            dtype = [("bin", np.int64)])
                                  , name="slicePoints")

        hduAng = fits.BinTableHDU.from_columns([
            fits.Column(name="binRA", array=iXs, format='I'),
            fits.Column(name="binDec",array=iYs, format='I'),
            fits.Column(name="countR",array=angR[iXs,iYs], format='I'),
            fits.Column(name="countD",array=angD[iXs,iYs], format='E')],
                                            name="ang")
        hduAng.header.add_comment("Unraveled angular (ra,dec) 2D histogram.")
        hduAng.header.add_comment("Histogram for random catalog filled with z independent weights.")
        self.addProvenance(hduAng,
                           self.config.inputFilesRandom() +
                           self.config.inputFilesObserved())

        binsZ = self.config.binningZ()['bins']
        binsRA = self.config.binningRA()['bins']
        binsDec = self.config.binningDec()['bins']

        jZ = La.utils.toBins(ctlgD.z, self.config.binningZ())
        jRA = La.utils.toBins(ctlgD.ra, self.config.binningRA())
        jDec= La.utils.toBins(ctlgD.dec, self.config.binningDec())
        jAng = binsDec*jRA  + jDec
        iRavelXY = binsDec*iXs + iYs

        def calcAngZD(weights):
            frq = csr_matrix((weights, (jAng, jZ)), shape=(binsRA*binsDec, binsZ))
            return frq[iRavelXY]

        angzD = calcAngZD(ctlgD.weight)
        err2 = calcAngZD(np.square(ctlgD.weight))

        iA, iZ = angzD.nonzero()
        hduAngzD = fits.BinTableHDU.from_columns([
            fits.Column(name="iAlign", array=iA, format='J'),
            fits.Column(name="iZ", array=iZ, format='I'),
            fits.Column(name='count', array=angzD.data, format='E'),
            fits.Column(name='err2', array=err2.data, format='E')],
                                             name="angzD")

        hduAngzD.header.add_comment("Sparse 3D histogram (ra, dec, z) of observed galaxies.")
        hduAngzD.header.add_comment("Unraveled in (ra,dec) and masked to align with 'ANG' rows.")
        self.addProvenance(hduAngzD, self.config.inputFilesObserved())

        return [hduSlc, hduAng, hduAngzD]


    def plot(self):
        from matplotlib import pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        infile = self.outputFileName

        def angPlot(pdf):
            ang = fits.getdata(infile, 'ANG')
            dc = fits.getdata(infile, 'centerDec').binCenter
            ra = fits.getdata(infile, 'centerRA').binCenter

            shp = len(ra), len(dc)
            h2d = csr_matrix((ang.countR, (ang.binRA, ang.binDec)), shape=shp)

            ddc = 0.5 * abs(dc[-1]-dc[0])/(len(dc)-1)
            dra = 0.5 * abs(ra[-1]-ra[0])/(len(ra)-1)
            ext = (ra[-1]-dra, ra[0]+dra, dc[0]-ddc, dc[-1]+ddc)

            plt.figure()
            plt.title(self.config.__class__.__name__)
            plt.imshow(np.fliplr(h2d.T.toarray()), origin='lower', extent=ext, interpolation='none', cmap='gray')
            plt.colorbar().set_label('random count')
            plt.xlabel(r'$\alpha$ [$\degree$]')
            plt.ylabel(r'$\delta$ [$\degree$]')
            pdf.savefig()
            plt.close()

        def zPlot(pdf):
            P = fits.getdata(infile, 'pdfz').probability
            z = fits.getdata(infile, 'centerz').binCenter
            dz = z[1]-z[0]

            plt.figure()
            plt.title(self.config.__class__.__name__)
            plt.bar(z, P, dz, color='blue', edgecolor='blue')
            plt.xlabel('redshift')
            plt.ylabel('probability')
            pdf.savefig()
            plt.close()

        with PdfPages(infile.replace('fits','pdf')) as pdf:
            angPlot(pdf)
            zPlot(pdf)
            print 'Wrote %s'% pdf._file.fh.name
        return
