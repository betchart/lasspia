#!/usr/bin/env python

import sys

def parseArgs():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Fast calculation of two-point correlations.")

    parser.add_argument('configFile', metavar='configFile', type=str, nargs=1,
                        help='A python file containing a subclass of baofast.configuration')

    parser.add_argument('routineFile', metavar='routineFile', type=str, nargs=1,
                        help='A python file containing a subclass of baofast.routine')

    parser.add_argument('--nJobs', metavar='nJobs', type=int, nargs=1,
                        help='Divide the processing into nJobs portions.')

    parser.add_argument('--iJob', metavar='iJob', type=int, nargs=1,
                        help='Index of the job to process.')

    args = parser.parse_args()
    return args

def getInstance(argFile, args = (), kwargs={}):
    path = argFile[0].split('/')
    name = path[-1].split('.')[0]
    sys.path.append('/'.join(path[:-1]))
    exec("from %s import %s " % (name, name))
    return eval(name)(*args, **kwargs)

if __name__ == "__main__":
    args = parseArgs()
    config = getInstance(args.configFile)
    routine = getInstance(args.routineFile, (config,),
                          {"nJobs":args.nJobs[0], "iJob":args.iJob[0]})
    routine()

    
