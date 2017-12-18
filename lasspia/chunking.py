import numpy as np
import slicing

def fromSlices(slices, predFunc):
    return [(iSlice,jSlice)
            for i,iSlice in enumerate(slices)
            for jSlice in slices[i:]
            if predFunc(iSlice,jSlice)]


def fromAllSlices(size, step):
    slices = slicing.slices(size, step)
    return fromSlices(slices, (lambda *_: True))


def fromProximateSlices(slicePoints, X, Y, dX, dY):
    slices = [slice(*z) for z in zip(slicePoints,slicePoints[1:])]

    hsh = lambda s: (s.start, s.stop, s.step)
    rect = dict([(hsh(s), (min(X[s]), min(Y[s]),
                           max(X[s]), max(Y[s]))) for s in slices])

    def proximate(iS, jS):
        imnX,imnY, imxX,imxY = rect[hsh(iS)]
        jmnX,jmnY, jmxX,jmxY = rect[hsh(jS)]
        return (max(imnX,jmnX) - min(imxX,jmxX) < dX and
                max(imnY,jmnY) - min(imxY,jmxY) < dY)

    return fromSlices(slices, proximate)


def byRegions(X, Y, axisSlicesX, axisSlicesY, step):

    def regionIndices(x,y):
        mask = reduce(np.logical_and, [x.start <= X, X < x.stop,
                                       y.start <= Y, Y < y.stop], 0 <= Y)
        return mask.nonzero()[0]

    regions = [[(i,j) for i in axisSlicesX] for j in axisSlicesY]

    rpairs = sum([
        [(r,r) for row in regions for r in row],                                             # self
        [(r1,r2) for row in regions for r1,r2 in zip(row,row[1:])],                          # right
        [(r1,r2) for row1,row2 in zip(regions,regions[1:]) for r1,r2 in zip(row1,row2)],     # below
        [(r1,r2) for row1,row2 in zip(regions,regions[1:]) for r1,r2 in zip(row1,row2[1:])], # below right
        [(r1,r2) for row1,row2 in zip(regions,regions[1:]) for r1,r2 in zip(row1[1:],row2)], # below left
        ], [])

    for r1,r2 in rpairs:
        i1 = regionIndices(*r1)
        i2 = regionIndices(*r2)
        schnks = ( fromAllSlices(len(i1),step) if r1==r2 else
                   [(s1,s2)
                    for s1 in slicing.slices(len(i1),step)
                    for s2 in slicing.slices(len(i2),step)] )
        for s1,s2 in schnks:
            yield i1[s1], i2[s2]
    return
