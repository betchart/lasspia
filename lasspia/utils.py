from multiprocessing import JoinableQueue,Process
import sys,traceback
from astropy.io import fits
import numpy as np
import math

def centers(leftEdges):
    return leftEdges[:-1] + 0.5 * (leftEdges[1:] - leftEdges[:-1])


def invBinWidth(binning):
    return (binning['bins'] /
            float(binning['range'][1] -
                  binning['range'][0]))

def toBins(ary, binning, dtype=np.int32):
    return ((ary-binning['range'][0]) * invBinWidth(binning)).astype(dtype)

def binRegions(delta, binning):
    '''List of slices breaking binning["range"] into regions of delta.'''
    iDelta = int(delta * invBinWidth(binning))
    return slices(binning['bins'], iDelta)

def slices(size, step=None, N=None):
    step = step or int(math.ceil(size/N))
    splits = range(0, size, step if step else size)
    return [slice(i,j) for i,j in zip(splits,splits[1:]+[size])]


def callInParallel(nCores, itemsToCall):
    q = JoinableQueue()
    processes=[Process(target=qWorker(), args=(q,)) for _ in range(nCores)]
    for p in processes:
        p.daemon = True
        p.start()

    map(q.put, itemsToCall)
    q.join()
    for p in processes: p.terminate()
    return

class qWorker(object):
    def __call__(self,q):
        while True:
            try: q.get()()
            except Exception as e:
                traceback.print_tb(sys.exc_info()[2], limit=20, file=sys.stdout)
                print e.__class__.__name__,":", e
            q.task_done()
        return


def hduDiff(nameHDU, hdus1, hdus2):
    return fits.HDUDiff(hdus1[nameHDU], hdus2[nameHDU])
