import lasspia as La

class quickscan(La.routine):

    def scan(self, label, ctlg):
        print >>self.out
        print >>self.out, "%s: %d" % (label, len(ctlg.z))
        print >>self.out, "  Sum of weights: %f" % sum(ctlg.weight)
        print >>self.out, "  Max weight: %f" % max(ctlg.weight)
        print >>self.out, "  Redshift range: [%f, %f] " % (min(ctlg.z), max(ctlg.z))
        print >>self.out, "  Declination range: [%f, %f]" % (min(ctlg.dec), max(ctlg.dec))
        print >>self.out, "  Right ascension range: [%f, %f]" % (min(ctlg.ra), max(ctlg.ra))

    def __call__(self):
        self.scan( "Observed Galaxies",
                   self.config.catalogObserved() )

        self.scan( "Random Galaxies",
                   self.config.catalogRandom() )
        print>>self.out
        return
