from __future__ import print_function
import math
import numpy as np
import lasspia as La
from astropy.io import fits
from scipy.integrate import quad
from scipy.sparse import csr_matrix
from lasspia.timing import timedHDU

class integration(La.routine):

    def __call__(self):
        try:
            self.hdus.append( self.binCenters(self.config.binningS(), "centerS") )
            self.hdus.append( self.binCenters(self.config.binningSigma(), "centerSigma") )
            self.hdus.append( self.binCenters(self.config.binningPi(), "centerPi") )
            self.hdus.extend(self.tpcf())
            self.writeToFile()
        except MemoryError as e:
            print(e.__class__.__name__, e, file=self.out)
            print('\n'.join(['Use less memory by integrating via multiple jobs.',
                             'For example, use options: --nJobs 8 --nCores 1',
                             'Then combine job outputs: --nJobs 8']),
                            file=self.out)
        return

    def omegasMKL(self): return self.config.omegasMKL()
    def H0(self): return self.config.H0()

    @timedHDU
    def binCenters(self, binning, name):
        centers = np.array( La.utils.centers(self.config.edgesFromBinning(binning)),
                            dtype = [("binCenter", np.float64)])
        return fits.BinTableHDU(centers, name=name)

    @timedHDU
    def tpcf(self):
        self.pdfz = self.getInput('pdfZ').data['probability']
        self.zMask = np.ones(len(self.pdfz)**2, dtype=np.int).reshape(len(self.pdfz),len(self.pdfz))
        self.zMask[:self.config.nBinsMaskZ(),:self.config.nBinsMaskZ()] = 0

        slcT =( slice(None) if self.iJob is None else
                La.slicing.slices(len(self.getInput('centertheta').data),
                                  N=self.nJobs)[self.iJob] )

        def bundleHDU(name, addresses, binning, axes, dropZeros=False):
            rr, dr, dd, dde2 = self.calc(addresses, binning, slcT)
            mask = np.logical_or.reduce([a!=0 for a in [rr, dr, dd]]) if dropZeros else np.full(rr.shape, True, dtype=bool)
            grid = [fits.Column(name="i"+k, array=iK, format='I')
                    for k,iK in zip(axes, np.where(mask))]
            hdu = fits.BinTableHDU.from_columns(grid + [
                fits.Column(name='RR', array=rr[mask], format='D'),
                fits.Column(name='DR', array=dr[mask], format='D'),
                fits.Column(name='DD', array=dd[mask], format='D'),
                fits.Column(name='DDe2', array=dde2[mask], format='D')],
                                                name=name)

            hdu.header['NORMRR'] = self.getInput('fTheta').header['NORM']
            hdu.header['NORMDR'] = self.getInput('gThetaZ').header['NORM']
            hdu.header['NORMDD'] = self.getInput('uThetaZZ').header['NORM']
            hdu.header.add_comment("Two-point correlation function for pairs of galaxies,"+
                                   " by distance" + ("s" if len(axes)>1 else "") + " " +
                                   " and ".join(axes))
            return hdu

        sigmaPis = self.sigmaPiGrid(slcT)
        s = np.sqrt(np.power(sigmaPis,2).sum(axis=-1))

        b = self.config.binningDD([self.config.binningS()])
        b2 = self.config.binningDD([self.config.binningSigma(),
                                    self.config.binningPi()])

        hdu = bundleHDU("TPCF", s, b, ["S"])
        hdu2 = bundleHDU("TPCF2D", sigmaPis, b2, ["Sigma", "Pi"], dropZeros=True)

        return [hdu2, hdu]

    def sigmaPiGrid(self, slcT):
        '''A cubic grid of (sigma, pi) values
        for pairs of galaxies with coordinates (iTheta, iZ1, iZ2).'''
        Iz = self.zIntegral()
        rOfZ = Iz * (self.config.lightspeed() / self.H0())
        tOfZ = rOfZ * (1 + self.omegasMKL()[1]/6 * Iz**2)

        thetas = self.getInput('centertheta').data['binCenter'][slcT]
        sinT2 = np.sin(thetas/2)
        cosT2 = np.cos(thetas/2)

        sigmas = sinT2[:,None,None] * (tOfZ[None,:,None] + tOfZ[None,None,:])
        pis = cosT2[:,None,None] * (rOfZ[None,:,None] - rOfZ[None,None,:])
        return np.stack([sigmas, pis], axis=-1)

    def calc(self, *args):
        return (self.calcRR(*args),
                self.calcDR(*args),
                self.calcDD(*args, wName='count'),
                self.calcDD(*args, wName='err2'))

    def calcRR(self, addresses, binning, slcT):
        ft = self.getInput('fTheta').data['count'][slcT]
        counts = ft[:,None,None] * self.pdfz[None,:,None] * self.pdfz[None,None,:] * self.zMask[None,:]
        N = counts.size
        D = addresses.size // N
        rr = np.histogramdd(addresses.reshape(N,D), weights=counts.reshape(N), **binning)[0]
        del counts
        return rr

    def calcDR(self, addresses, binning, slcT):
        gtz = self.getInput('gThetaZ').data
        counts = gtz[slcT,:,None] * self.pdfz[None,None,:] * self.zMask[None,:]
        N = counts.size
        D = addresses.size // N
        dr = np.histogramdd(addresses.reshape(N,D), weights=counts.reshape(N), **binning)[0]
        del counts
        return dr

    def calcDD(self, addresses, binning, slcT, wName='count'):
        nZ = addresses.shape[1]

        utzz = self.getInput('uThetaZZ').data
        overflow = utzz['binZdZ'][-1]+1 == nZ**2
        slc = slice(-1 if overflow else None)

        iThetas = utzz['binTheta'][slc]
        mask = (slice(None) if slcT==slice(None) else
                np.logical_and(slcT.start <= iThetas, iThetas < slcT.stop))

        iTh = iThetas[mask] - (slcT.start or 0)
        iZdZ = utzz['binZdZ'][slc][mask]
        iZ = iZdZ // nZ
        diZ = iZdZ % nZ
        iZ2 = iZ + diZ
        counts = utzz[wName][slc][mask] * self.zMask[iZ,iZ2]

        dd = np.histogramdd(addresses[iTh,iZ,iZ2], weights=counts, **binning)[0]
        if overflow and self.iJob in [0,None]:
            dd[-1] = dd[-1] + utzz[wName][-1]
        return dd

    def zIntegral(self):
        zCenters = self.getInput('centerz').data['binCenter']
        zz = zip(np.hstack([[0.],zCenters]), zCenters)
        dIz = [quad(self.integrand, z1, z2, args=self.omegasMKL())[0]
               for z1,z2 in zz]
        return np.cumsum(dIz)

    @staticmethod
    def integrand(z, omegaM, omegaK, omegaLambda):
        return 1./math.sqrt(omegaM * (1+z)**3 +
                            omegaK * (1+z)**2 +
                            omegaLambda)

    @property
    def inputFileName(self):
        return self.config.stageFileName('combinatorial')

    def getInput(self, name):
        hdulist = fits.open(self.inputFileName)
        return hdulist[name]

    def combineOutput(self, jobFiles = None):
        if not jobFiles:
            jobFiles = [self.outputFileName + self.jobString(iJob)
                        for iJob in range(self.nJobs)]

        shape2D = (self.config.binningSigma()['bins'], self.config.binningPi()['bins'])

        with fits.open(jobFiles[0]) as h0:
            for h in ['centerS','centerSigma','centerPi']:
                self.hdus.append(h0[h])
            hdu = h0['TPCF']
            tpcf2d = AdderTPCF2D(h0['TPCF2D'], shape2D)
            cputime = hdu.header['cputime']

            for jF in jobFiles[1:]:
                with fits.open(jF) as jfh:
                    assert np.all( hdu.data['iS'] == jfh['TPCF'].data['iS'])
                    cputime += jfh['TPCF'].header['cputime']
                    tpcf2d += AdderTPCF2D(jfh['TPCF2D'], shape2D)
                    for col in ['RR','DR','DD','DDe2']:
                        hdu.data[col] += jfh['TPCF'].data[col]
            hdu.header['cputime'] = cputime
            self.hdus.append(tpcf2d.fillHDU(h0['TPCF2D']))
            self.hdus.append(hdu)
            self.writeToFile()
        return

    def combineOutputZ(self):
        zFiles = [self.outputFileName.replace(self.config.name,
                                              '_'.join([self.config.name,
                                                        self.config.suffixZ(iZ)]))
                  for iZ in range(len(self.config.binningsZ()))]
        self.combineOutput(zFiles)

    def plot(self):
        from matplotlib import pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        infile = self.outputFileName

        tpcf = fits.getdata(infile, 'TPCF')
        centerS = fits.getdata(infile, 'centerS').binCenter
        nRR,nDR,nDD = (lambda h:
                       (h['normrr'],
                        h['normdr'],
                        h['normdd']))(fits.getheader(infile, 'TPCF'))

        def tpcfPlot(pdf, binFactor):
            s = centerS[tpcf.iS]
            iStop = len(s) // binFactor
            plt.figure()
            plt.title(self.config.__class__.__name__)
            plt.step(s[:iStop], tpcf.RR[:iStop]/nRR, where='mid', label='RR', linewidth=0.4)
            plt.step(s[:iStop], tpcf.DR[:iStop]/nDR, where='mid', label='DR', linewidth=0.4)
            plt.step(s[:iStop], tpcf.DD[:iStop]/nDD, where='mid', label='DD', linewidth=0.4)
            plt.legend()
            plt.xlabel('s')
            plt.ylabel('probability')
            pdf.savefig()
            plt.close()

        def xissPlot(pdf, sMax):
            S = centerS[tpcf.iS]
            iStop = None if S[-1]<sMax else next(iter(np.where(S >= sMax)[0]))
            s = S[:iStop]
            xi = ( (tpcf.RR[:iStop]/nRR + tpcf.DD[:iStop]/nDD - 2*tpcf.DR[:iStop]/nDR)
                   / (tpcf.RR[:iStop]/nRR) )
            xie = ( (np.sqrt(tpcf.DDe2[:iStop])/nDD)
                    / (tpcf.RR[:iStop]/nRR) )
            ds = s[1]-s[0]

            plt.figure()
            plt.title(self.config.__class__.__name__)
            plt.errorbar(s, (xi*s*s), yerr=(xie*s*s), xerr=ds/2, fmt='.')
            plt.xlabel(r"$\mathrm{s\ [h^{-1} Mpc]}$")
            plt.ylabel(r"$\mathrm{\xi(s)s^2}$")
            plt.grid()
            pdf.savefig()
            plt.close()

        with PdfPages(infile.replace('fits','pdf')) as pdf:
            for i in range(5):
                tpcfPlot(pdf, 2**i)
            xissPlot(pdf, 200)
            print('Wrote %s'% pdf._file.fh.name, file=self.out)
        return

class AdderTPCF2D(object):
    def __init__(self, tpcf2d=None, shape2D=None):
        self.items = ['RR','DR','DD','DDe2']
        if not tpcf2d: return
        indices = (tpcf2d.data['iSigma'], tpcf2d.data['iPi'])
        for item in self.items:
            setattr(self, item, csr_matrix((tpcf2d.data[item], indices), shape2D))
        return

    def __add__(self, other):
        thesum = AdderTPCF2D()
        for item in self.items:
            setattr(thesum, item, getattr(self,item) + getattr(other, item))
        return thesum

    def fillHDU(self, hdu):
        allnonzero = sum(getattr(self, item) for item in self.items)
        iSigma, iPi = allnonzero.nonzero()
        hdu.data['iSigma'] = iSigma
        hdu.data['iPi'] = iPi
        for item in self.items:
            hdu.data[item] = getattr(self, item)[iSigma, iPi].A1
        return hdu
