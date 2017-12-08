import lasspia as La

class quickscan(La.routine):

    def scan(self, label, ctlg):
        with open self.outStream() as f:
            print >>f
            print >>f, "%s: %d" % (label, len(ctlg.z))
            print >>f, "  Sum of weights: %f" % sum(ctlg.weight)
            print >>f, "  Max weight: %f" % max(ctlg.weight)
            print >>f, "  Redshift range: [%f, %f] " % (min(ctlg.z), max(ctlg.z))
            print >>f, "  Declination range: [%f, %f]" % (min(ctlg.dec), max(ctlg.dec))
            print >>f, "  Right ascension range: [%f, %f]" % (min(ctlg.ra), max(ctlg.ra))
    
    def __call__(self):
        self.scan( "Observed Galaxies",
                   self.config.catalogObserved() )

        self.scan( "Random Galaxies",
                   self.config.catalogRandom() )

        with self.outStream() as f:
            print>>f
        return
