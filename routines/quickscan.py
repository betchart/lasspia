import lasspia as La

class quickscan(La.routine):

    def scan(self, label, ctlg):
        print
        print "%s: %d" % (label, len(ctlg.z))
        print "  Sum of weights: %f" % sum(ctlg.weight)
        print "  Max weight: %f" % max(ctlg.weight)
        print "  Redshift range: [%f, %f] " % (min(ctlg.z), max(ctlg.z))
        print "  Declination range: [%f, %f]" % (min(ctlg.dec), max(ctlg.dec))
        print "  Right ascension range: [%f, %f]" % (min(ctlg.ra), max(ctlg.ra))
    
    def __call__(self):
        self.scan( "Observed Galaxies",
                   self.config.catalogObserved() )

        self.scan( "Random Galaxies",
                   self.config.catalogRandom() )
        print
        return
