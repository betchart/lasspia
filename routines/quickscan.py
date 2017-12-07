import lasspia as La

class quickscan(La.routine):

    def scan(self, label, ctlg):
        print >>self.outstream
        print >>self.outstream, "%s: %d" % (label, len(ctlg.z))
        print >>self.outstream, "  Sum of weights: %f" % sum(ctlg.weight)
        print >>self.outstream, "  Max weight: %f" % max(ctlg.weight)
        print >>self.outstream, "  Redshift range: [%f, %f] " % (min(ctlg.z), max(ctlg.z))
        print >>self.outstream, "  Declination range: [%f, %f]" % (min(ctlg.dec), max(ctlg.dec))
        print >>self.outstream, "  Right ascension range: [%f, %f]" % (min(ctlg.ra), max(ctlg.ra))
    
    def __call__(self):
        self.scan( "Observed Galaxies",
                   self.config.catalogObserved() )

        self.scan( "Random Galaxies",
                   self.config.catalogRandom() )
        print >>self.outstream
        return
