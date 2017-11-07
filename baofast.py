#!/usr/bin/env python

def parseArgs():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Fast calculation of two-point correlations.")
    parser.add_argument('configFile', metavar='configFile', type=str, nargs=1,
                        help='A python file containing a subclass of baofast.configuration')
    args = parser.parse_args()
    return args

def getConfig(args):
    import sys
    path = args.configFile[0].split('/')
    name = path[-1].split('.')[0]
    sys.path.append('/'.join(path[:-1]))
    exec("from %s import %s " % (name, name))
    return eval(name)()

if __name__ == "__main__":
    args = parseArgs()
    config = getConfig(args)

    # temporary below
    import numpy as np
    ctlg = config.catalogObserved()
    frq, edges = np.histogram(ctlg.z, config.binsZ(), weights=ctlg.weight)

    print frq
    print edges

    import matplotlib.pyplot as plt
    plt.bar(edges[:-1], frq, edges[1]-edges[0], align='edge')
    plt.show()

    
