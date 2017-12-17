import numpy as np
import utils

def xyClustersWhere(mask, limit):
    cc = groupedPoints(*np.where(mask), limit=limit)
    return indicesCat(cc)

def groupedPoints(iX,iY, limit=10):
    E = lambda x: sum(x)/float(len(x))
    if len(iX) < limit:
        return [(iX, iY)]
    EiX = E(iX)
    EiY = E(iY)
    vX = E(iX**2) - EiX**2
    vY = E(iY**2) - EiY**2
    pred = (iX < EiX if vX>vY else
            iY < EiY)
    _pred = np.invert(pred)
    return (groupedPoints(iX[pred], iY[pred], limit) +
            groupedPoints(iX[_pred], iY[_pred], limit))

def indicesCat(iXYs):
    iXs,iYs = zip(*iXYs)
    lens = map(len,iXs)
    assert lens == map(len,iYs)
    return np.cumsum([0]+lens), np.hstack(iXs), np.hstack(iYs)


def binRegions(delta, binning):
    '''List of slices breaking binning["range"] into regions of delta.'''
    iDelta = int(delta * utils.invBinWidth(binning))
    return slices(binning['bins'], iDelta)

def slices(size, step=None, N=None):
    step = step or int(math.ceil(size/N))
    splits = range(0, size, step if step else size)
    return [slice(i,j) for i,j in zip(splits,splits[1:]+[size])]
