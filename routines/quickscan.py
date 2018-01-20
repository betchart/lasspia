from __future__ import print_function
import lasspia as La

class quickscan(La.routine):

    def scan(self, label, ctlg):
        output = [''
                  ,"%s: %d" % (label, len(ctlg.z))
                  ,"  Sum of weights: %f" % sum(ctlg.weight)
                  ,"  Max weight: %f" % max(ctlg.weight)
                  ,"  Redshift range: [%f, %f] " % (min(ctlg.z), max(ctlg.z))
                  ,"  Declination range: [%f, %f]" % (min(ctlg.dec), max(ctlg.dec))
                  ,"  Right ascension range: [%f, %f]" % (min(ctlg.ra), max(ctlg.ra))
        ]
        print('\n'.join(output), file=self.out)

    def __call__(self):
        self.scan( "Observed Galaxies",
                   self.config.catalogObserved() )

        self.scan( "Random Galaxies",
                   self.config.catalogRandom() )
        print(file=self.out)
        return
