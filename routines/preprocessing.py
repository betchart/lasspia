import baofast
import numpy as np

class preprocessing(baofast.routine):
    """Preprocessing for fast 2-point correlations.

    Open input catalogs, create histograms, save to file.
    """
    def __call__(self):
        ctlg = self.config.catalogObserved()
        pdfZ = self.pdfZ(ctlg)
        
        print pdfZ['lowEdge']
        print pdfZ['probability']
        
        import matplotlib.pyplot as plt
        edg = pdfZ['lowEdge']
        plt.bar(edg, pdfZ['probability'], edg[1]-edg[0], align='edge')
        plt.show()

    def pdfZ(self, ctlg):
        frq, edges = np.histogram(ctlg.z, self.config.binsZ(),
                                  weights = ctlg.weight / sum(ctlg.weight))
        return np.array(zip(edges, frq),
                        dtype = [("lowEdge", np.float64),
                                 ("probability", np.float32)])
        
