#!/usr/bin/env python

import sys
from baofast import utils

def parseArgs():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Fast calculation of two-point correlations.")

    parser.add_argument('configFile', metavar='configFile', type=str, nargs=1,
                        help='A python file containing a subclass of baofast.configuration')

    parser.add_argument('routineFile', metavar='routineFile', type=str, nargs=1,
                        help='A python file containing a subclass of baofast.routine')

    parser.add_argument('--nJobs', metavar='nJobs', type=int, nargs=1,
                        help='Divide the processing into nJobs portions: process all jobs in parallel (with --nCores), or process just one job (with --iJob), or combine job outputs.')

    parser.add_argument('--iJob', metavar='iJob', type=int, nargs=1,
                        help='Index of the job to process (requires --nJobs).')

    parser.add_argument('--nCores', metavar='nCores', type=int, nargs=1,
                        help='Use nCores in parallel (requires --nJobs).')

    args = parser.parse_args()
    return args


def getInstance(argFile, args = (), kwargs={}):
    path = argFile[0].split('/')
    name = path[-1].split('.')[0]
    sys.path.append('/'.join(path[:-1]))
    exec("from %s import %s " % (name, name))
    return eval(name)(*args, **kwargs)

def getKWs(args):
    if args.nCores: return [{"nJobs":args.nJobs[0], "iJob":i} for i in range(args.nJobs[0])]
    if args.iJob: return {"nJobs":args.nJobs[0], "iJob":args.iJob[0]}
    if args.nJobs: return {"nJobs":args.nJobs[0]}
    return {}

if __name__ == "__main__":
    args = parseArgs()
    config = getInstance(args.configFile)
    kwargs = getKWs(args)

    if type(kwargs) is dict:
        routine = getInstance(args.routineFile, (config,), kwargs)
        if args.nJobs and not args.iJob:
            routine.combineOutput()
        else: routine()

    elif args.nCores:
        routines = [getInstance(args.routineFile, (config,), kw) for kw in kwargs]
        utils.callInParallel( args.nCores[0], routines )

    else:
        pass
