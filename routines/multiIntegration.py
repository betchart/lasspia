from integration import integration

class multiIntegration(integration):

    '''Perform integration multiple times, scanning the integration parameters.

    This class is meant to be an example of how one might achieve such
    a scan.  The pattern used is to inherit from the 'integration'
    routine, and repeatedly call its methods after changing
    parameters.  Notice that the output file name also changes with
    parameters.

    Similar results could be achieved by creating a new instance of
    the 'integration' routine for every iteration of the parameters.
    One might also choose to inherit from the 'integration' routine,
    but micromanage the hdu list by calling more specific methods,
    like integration.tpcf.

    '''

    @property
    def variationsOmegasMKL(self):
        try:
            return self.config.variationsOmegasMKL()
        except:
            default = [(0.273, 0, 0.726), (0.274, 0, 0.726), (0.275, 0, 0.726)]
            if not (hasattr(self, 'warned') and self.warned):
                print('Warning: variationsOmegasMKL() is not defined in the configuration!')
                print('multiIntegration will proceed using the default variations:')
                print(default)
                self.warned = True
            return default

    @property
    def nVariations(self): return len(self.variationsOmegasMKL)

    def omegasMKL(self): return self.variationsOmegasMKL[self.iOmegas]
    def H0(self): return self.config.H0()

    @property
    def outputFileName(self):
        index = "_%03d" % self.iOmegas
        base = super(multiIntegration, self).outputFileName
        pre,suf = base.split('.fits')
        return pre + index + '.fits' + suf

    def __call__(self):
        for i in range(self.nVariations):
            self.iOmegas = i
            super(multiIntegration, self).__call__()

        return

    def combineOutput(self):
        for i in range(self.nVariations):
            self.iOmegas = i
            super(multiIntegration, self).combineOutput()

    def plot(self):
        for i in range(self.nVariations):
            self.iOmegas = i
            super(multiIntegration, self).plot()

    def showFitsHeaders(self):
        for i in range(self.nVariations):
            self.iOmegas = i
            super(multiIntegration, self).showFitsHeaders()
