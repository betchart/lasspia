#!/usr/bin/env python

import sys

def parseArgs():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Fast calculation of two-point correlations.")
    parser.add_argument('configFile', metavar='configFile', type=str, nargs=1,
                        help='A python file containing a subclass of baofast.configuration')
    parser.add_argument('routineFile', metavar='routineFile', type=str, nargs=1,
                        help='A python file containing a subclass of baofast.routine')
    args = parser.parse_args()
    return args

def getInstance(argFile, init = ()):
    path = argFile[0].split('/')
    name = path[-1].split('.')[0]
    sys.path.append('/'.join(path[:-1]))
    exec("from %s import %s " % (name, name))
    return eval(name)(*init)

if __name__ == "__main__":
    args = parseArgs()
    config = getInstance(args.configFile)
    routine = getInstance(args.routineFile, (config,))
    routine()

    
