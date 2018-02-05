#!/usr/bin/env python

import sys
from lasspia import utils
from lasspia.zSlicing import SlicesZ

def parseArgs():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Large Scale Structure Probability Integration Algorithm")

    parser.add_argument('configFile', metavar='configFile', type=str, nargs=1,
                        help='A python file containing a subclass of lasspia.configuration of the same name.')

    parser.add_argument('routineFile', metavar='routineFile', type=str, nargs=1,
                        help='A python file containing a subclass of lasspia.routine of the same name.')

    parser.add_argument('--nJobs', metavar='nJobs', type=int, nargs=1,
                        help='Divide the processing into nJobs portions: process all jobs in parallel (with --nCores), or process just one job (with --iJob or --iJobEnv), or combine job outputs.')

    parser.add_argument('--iJob', metavar='iJob', type=int, nargs='+',
                        help='Index of the job to process (requires --nJobs; see also --iJobEnv).')

    parser.add_argument('--iJobEnv', metavar='iJobEnv', type=str, nargs=1,
                        help='Environment variable containing index of the job to process (requires --nJobs; alternative to --iJob).')

    parser.add_argument('--nCores', metavar='nCores', type=int, nargs=1,
                        help='Use nCores in parallel on your local machine. (requires --nJobs; not for use on batch systems)')

    parser.add_argument('--txtToFile', action='store_true',
                        help='Redirect messages to file.')

    parser.add_argument('--show', action='store_true',
                        help='Show info and HDU headers of the output file.')

    parser.add_argument('--plot', action='store_true',
                        help='Run the plot() method of the routine.')

    parser.add_argument('--iSliceZ', metavar='iSliceZ', type=int, nargs=1,
                        help='Set the index of the z-slice (only for classes inheriting from lasspia/zSlicing.py:SlicesZ.')

    args = parser.parse_args()
    parseEnv(args)
    return args

def parseEnv(args):
    if not args.iJob and args.iJobEnv:
        import os
        try:
            iJob = int(os.environ[args.iJobEnv[0]])
            args.iJob = [iJob]
        except KeyError:
            print( "No such environment variable %s" % args.iJobEnv[0])
            exit()
    return

def getInstance(argFile, args = (), kwargs={}):
    path = argFile[0].split('/')
    name = path[-1].split('.')[0]
    sys.path.append('/'.join(path[:-1]))
    exec("from %s import %s " % (name, name))
    return eval(name)(*args, **kwargs)

def getKWs(args):
    if args.nCores or args.iJob:
        n = args.nJobs[0] if args.nJobs else 1
        jobs = args.iJob if args.iJob else range(args.nJobs[0])
        return [{"nJobs": n, "iJob":i} for i in jobs]
    if args.nJobs: return {"nJobs":args.nJobs[0]}
    return {}

def getCfgArgs(args):
    cfgArgs = {'txtToFile':args.txtToFile}
    if args.iSliceZ: cfgArgs['iSliceZ'] = args.iSliceZ[0]
    return cfgArgs

if __name__ == "__main__":
    args = parseArgs()
    config = getInstance(args.configFile, kwargs=getCfgArgs(args))
    kwargs = getKWs(args)

    if type(kwargs) is dict:
        routine = getInstance(args.routineFile, (config,), kwargs)
        if args.show: routine.showFitsHeaders()
        elif args.plot: routine.plot()
        elif args.nJobs: routine.combineOutput()
        elif (issubclass(config.__class__, SlicesZ)
              and config.iSliceZ is None):
            routine.combineOutputZ()
        else: routine()

    elif type(kwargs) is list:
        routines = [getInstance(args.routineFile, (config,), kw) for kw in kwargs]
        utils.callInParallel( args.nCores[0] if args.nCores else 1, routines )

    else:
        pass
