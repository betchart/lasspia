from __future__ import print_function
import math
import numpy as np
import lasspia as La
from astropy.io import fits
from scipy.integrate import quad
from lasspia.timing import timedHDU

class integration(La.routine):

    def __call__(self):
        try:
            self.hdus.append(self.tpcf())
            self.writeToFile()
        except MemoryError as e:
            print(e.__class__.__name__, e, file=self.out)
            print('\n'.join(['Use less memory by integrating via multiple jobs.',
                             'For example, use options: --nJobs 8 --nCores 1',
                             'Then combine job outputs: --nJobs 8']),
                            file=self.out)
        return

    @timedHDU
    def tpcf(self):
        self.pdfz = self.getInput('pdfZ').data['probability']
        self.zMask = np.ones(len(self.pdfz)**2, dtype=np.int).reshape(len(self.pdfz),len(self.pdfz))
        self.zMask[:self.config.nBinsMaskZ(),:self.config.nBinsMaskZ()] = 0

        binsS = self.config.binningS()['bins']

        slcT =( slice(None) if self.iJob is None else
                La.slicing.slices(len(self.getInput('centertheta').data),
                                  N=self.nJobs)[self.iJob] )

        s = np.sqrt(sum(np.power(a,2) for a in self.sigmaPiGrids(slcT)))
        hdu = fits.BinTableHDU.from_columns([
            fits.Column(name='s', array=self.centersS(), format='E'),
            fits.Column(name='RR', array=self.calcRR(s, slcT), format='D'),
            fits.Column(name='DR', array=self.calcDR(s, slcT), format='D'),
            fits.Column(name='DD', array=self.calcDD(s, slcT, 'count'), format='D'),
            fits.Column(name='DDe2', array=self.calcDD(s, slcT, 'err2'), format='D')],
                                            name="TPCF")

        hdu.header['NORMRR'] = self.getInput('fTheta').header['NORM']
        hdu.header['NORMDR'] = self.getInput('gThetaZ').header['NORM']
        hdu.header['NORMDD'] = self.getInput('uThetaZZ').header['NORM']

        hdu.header.add_comment("Two-point correlation function for pairs of galaxies,"+
                               " by distance s.")
        return hdu

    def sigmaPiGrids(self, slcT):
        '''A cubic grid of sigma (pi) values
        for pairs of galaxies with coordinates (iTheta, iZ1, iZ2).'''
        Iz = self.zIntegral()
        rOfZ = Iz * (self.config.lightspeed() / self.config.H0())
        tOfZ = rOfZ * (1 + self.config.omegasMKL()[1]/6 * Iz**2)

        thetas = self.getInput('centertheta').data['binCenter'][slcT]
        sinT2 = np.sin(thetas/2)
        cosT2 = np.cos(thetas/2)

        sigmas = sinT2[:,None,None] * (tOfZ[None,:,None] + tOfZ[None,None,:])
        pis = cosT2[:,None,None] * (rOfZ[None,:,None] - rOfZ[None,None,:])
        return sigmas, pis

    def calcRR(self,s, slcT):
        ft = self.getInput('fTheta').data['count'][slcT]
        counts = ft[:,None,None] * self.pdfz[None,:,None] * self.pdfz[None,None,:] * self.zMask[None,:]
        rr = np.histogram(s, weights=counts, **self.config.binningS())[0]
        del counts
        return rr

    def calcDR(self,s, slcT):
        gtz = self.getInput('gThetaZ').data
        counts = gtz[slcT,:,None] * self.pdfz[None,None,:] * self.zMask[None,:]
        dr = np.histogram(s, weights=counts, **self.config.binningS())[0]
        del counts
        return dr

    def calcDD(self,s, slcT, wName='count'):
        utzz = self.getInput('uThetaZZ').data
        overflow = utzz['binZdZ'][-1]+1 == s.shape[1]**2
        slc = slice(-1 if overflow else None)

        iThetas = utzz['binTheta'][slc]
        mask = (slice(None) if slcT==slice(None) else
                np.logical_and(slcT.start <= iThetas, iThetas < slcT.stop))

        iTh = iThetas[mask] - (slcT.start or 0)
        iZdZ = utzz['binZdZ'][slc][mask]
        iZ = iZdZ // s.shape[1]
        diZ = iZdZ % s.shape[1]
        iZ2 = iZ + diZ
        counts = utzz[wName][slc][mask] * self.zMask[iZ,iZ2]

        dd = np.histogram(s[iTh,iZ,iZ2], weights=counts, **self.config.binningS())[0]
        if overflow and self.iJob in [0,None]:
            dd[-1] = dd[-1] + utzz[wName][-1]
        return dd

    def centersS(self):
        return La.utils.centers(self.config.edgesFromBinning(self.config.binningS()))

    def zIntegral(self):
        zCenters = self.getInput('centerz').data['binCenter']
        zz = zip(np.hstack([[0.],zCenters]), zCenters)
        dIz = [quad(self.integrand, z1, z2, args=self.config.omegasMKL())[0]
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

        with fits.open(jobFiles[0]) as h0:
            hdu = h0['TPCF']
            cputime = hdu.header['cputime']
            for jF in jobFiles[1:]:
                with fits.open(jF) as jfh:
                    assert np.all( hdu.data['s'] == jfh['TPCF'].data['s'])
                    cputime += jfh['TPCF'].header['cputime']
                    for col in ['RR','DR','DD','DDe2']:
                        hdu.data[col] += jfh['TPCF'].data[col]
            hdu.header['cputime'] = cputime
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
        nRR,nDR,nDD = (lambda h:
                       (h['normrr'],
                        h['normdr'],
                        h['normdd']))(fits.getheader(infile, 'TPCF'))

        def tpcfPlot(pdf, binFactor):
            iStop = len(tpcf.s) // binFactor
            plt.figure()
            plt.title(self.config.__class__.__name__)
            plt.step(tpcf.s[:iStop], tpcf.RR[:iStop]/nRR, where='mid', label='RR', linewidth=0.4)
            plt.step(tpcf.s[:iStop], tpcf.DR[:iStop]/nDR, where='mid', label='DR', linewidth=0.4)
            plt.step(tpcf.s[:iStop], tpcf.DD[:iStop]/nDD, where='mid', label='DD', linewidth=0.4)
            plt.legend()
            plt.xlabel('s')
            plt.ylabel('probability')
            pdf.savefig()
            plt.close()

        def xissPlot(pdf, sMax):
            iStop = None if tpcf.s[-1]<sMax else next(iter(np.where(tpcf.s >= sMax)[0]))
            s = tpcf.s[:iStop]
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
