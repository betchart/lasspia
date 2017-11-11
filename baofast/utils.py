from multiprocessing import JoinableQueue,Process
import sys,traceback
from astropy.io import fits

def centers(leftEdges):
    return leftEdges[:-1] + 0.5 * (leftEdges[1:] - leftEdges[:-1])


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


def identicalHDUs(nameHDU, hdus1, hdus2):
    return fits.HDUDiff(hdus1[nameHDU], hdus2[nameHDU]).identical
