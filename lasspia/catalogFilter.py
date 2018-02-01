from __future__ import print_function
import numpy as np
from copy import deepcopy
from functools import reduce

class FilteredCatalog(object):
    def addProp(p):
        return property(lambda s: getattr(s.cat, p)[s.indices])

    z = addProp('z')
    ra = addProp('ra')
    dec = addProp('dec')
    weight = addProp("weight")
    weightZ = addProp("weightZ")
    weightNoZ = addProp("weightNoZ")

    def addSumProp(p):
        '''Provide access to the unfiltered sum of weights.'''
        sumProp = 'sum'+p
        return property(lambda s: sum(getattr(s.cat, p)))

    sumweight = addSumProp('weight')
    sumweightZ = addSumProp('weightZ')
    sumweightNoZ = addSumProp('weightNoZ')

    def __init__(self, cat, indices):
        if type(cat) is FilteredCatalog:
            self.cat = cat.cat
            self.indices = cat.indices[indices]
        else:
            self.cat = cat
            self.indices = indices

    def __len__(self): return len(self.indices)

class catalogFilter(object):
    def __init__(self, binningRA=None, binningDec=None, binningZ=None, output=None):
        self.out = output
        self.cuts = []
        for atr, binning in zip(['ra','dec','z'],
                                ['binningRA','binningDec','binningZ']):
            bng = eval(binning)
            if not bng: continue
            self.cuts.append((atr, np.greater_equal, bng['range'][0]))
            self.cuts.append((atr, np.less, bng['range'][1]))

    def pred(self, cat):
        indicesCut = [compare(getattr(cat, atr), val) for atr, compare, val in self.cuts]
        return reduce(np.logical_and,
                      indicesCut,
                      np.full(len(cat), True, dtype=bool))

    def __call__(self, cat):
        indices = np.where(self.pred(cat))[0]
        if len(indices) < len(cat):
            print("Warning: Some catalog values are outside the configured ranges for %s and will be filtered." %
                  (' and '.join(set(a for a,c,v in self.cuts))),
                  file=self.out)
        return FilteredCatalog(cat, indices)
