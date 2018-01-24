from lasspia import utils
import numpy as np
import math

def overlapBinning(binning1, minOverlap, maxBinWidth2, minHi2):
    nBins1 = binning1['bins']
    lo1,hi1 = binning1['range']
    bw1 = (hi1 - lo1) / float(nBins1)
    nOverlap1 = math.ceil( minOverlap / bw1)
    hi1 = binning1['range'][1]
    lo2 = hi1 - nOverlap1 * bw1

    nOverlap2 = int(math.ceil( nOverlap1 * bw1 / maxBinWidth2 ))
    bw2 = nOverlap1 * bw1 / nOverlap2
    nBins2 = int(math.ceil((minHi2 - lo2) / bw2))
    hi2 = lo2 + nBins2*bw2

    return {'range':(lo2,hi2), 'bins':nBins2}, nOverlap2


def overlapBinnings(minOverlap, lo, hiTargets, maxBinWidths):
    '''Return list of pairs (binning, nBinsOverlap) for each region.

    Overlaping bins are at the beginning of the given range, and
    overlap with the previous binning definition in the list.
    '''
    assert len(hiTargets) == len(maxBinWidths)
    binning1 = {'range':(lo, hiTargets[0]),
                'bins': int(math.ceil((hiTargets[0]-lo) / maxBinWidths[0]))}

    binnings = [binning1]
    overlaps = [0]
    for hiT,maxBW in zip(hiTargets[1:], maxBinWidths[1:]):
        b,o = overlapBinning(binnings[-1], minOverlap, maxBW, hiT)
        binnings.append(b)
        overlaps.append(o)

    return zip(binnings, overlaps)


def checkOverlapBinning(binningsAndOverlaps):
    for (b1,_), (b2,nOverlap) in zip(binningsAndOverlaps,
                                     binningsAndOverlaps[1:]):
        ibw1 = utils.invBinWidth(b1)
        ibw2 = utils.invBinWidth(b2)

        N = utils.toBins(np.array(b1['range'][1]), b2, dtype=np.float32)
        assert N == nOverlap, "%f != %f" % (N, nOverlap)
        n3 = nOverlap * ibw1 / ibw2
        assert abs(round(n3) - n3) < 1e-6, "overlap not integer multiple of both binnings!"
    return
