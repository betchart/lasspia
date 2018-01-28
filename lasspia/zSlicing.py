from __future__ import print_function
from lasspia.configuration import configuration
from lasspia.overlap import overlapBinnings
from lasspia.utils import invBinWidth

class SlicesZ(object): pass

def zSlicing(Conf):
    '''Return a z-slicing class inheriting from the class Conf(lasspia.configuration).

    The returned class is functionally identical to Conf (a single
    slice in z).  You can enable z-slicing by inheriting the returned
    class and overriding 'zBreaks' and 'zMaxBinWidths'.
    Alternatively, you can override 'overlapBinningZ'.

    You could additionally define angular ranges dependent on iSliceZ
    by overriding 'binningRA', 'binningDec', and even 'binningTheta'.
    '''

    if not issubclass(Conf, configuration):
        raise Exception('"%s" is not a subclass of lasspia.configuration.' % Conf.__name__)


    class ParticularSlicesZ(Conf, SlicesZ):

        def zBreaks(self):
            return Conf.binningZ(self)['range']

        def zMaxBinWidths(self):
            return [1./invBinWidth(Conf.binningZ(self))]

        def overlapBinningZ(self):
            maxDeltaZ = Conf.maxDeltaZ(self)
            if maxDeltaZ is None:
                print("z-slicing only works with maxDeltaZ defined.")
                exit()

            zBreaks = self.zBreaks()
            if any(hi-lo < maxDeltaZ for lo,hi in zip(zBreaks,zBreaks[1:])):
                print("z-slice ranges must be greater than maxDeltaZ.")
                exit()

            return overlapBinnings(maxDeltaZ, zBreaks[0], zBreaks[1:], self.zMaxBinWidths())

        # Proceed with caution before overriding the functions defined below.

        def binningsZ(self):  return list(zip(*self.overlapBinningZ()))[0]
        def nOverlapsZ(self): return list(zip(*self.overlapBinningZ()))[1]

        def binningZ(self): return self.binningsZ()[self.iSliceZ]
        def nBinsMaskZ(self): return self.nOverlapsZ()[self.iSliceZ]

        @property
        def iSliceZ(self):
            if self._iSliceZ is None: pass
            elif self._iSliceZ < 0 or self._iSliceZ >= len(self.binningsZ()):
                print("--iSliceZ out of range (%d slices configured)"
                      % len(self.binningsZ()))
                exit()
            return self._iSliceZ

        def suffixes(self):
            return [] if self.iSliceZ is None else [self.suffixZ(self.iSliceZ)]

        def suffixZ(self, iZ): return "z%d" % iZ

        def __init__(self, iSliceZ=None, **kwargs):
            self._iSliceZ = iSliceZ
            Conf.__init__(self, **kwargs)

        def info(self, output=None):
            Conf.info(self, output)
            if self.iSliceZ is not None:
                print("Index %d of %d slices in z, %d bins on [%f, %f]"
                      % ((self.iSliceZ, len(self.binningsZ()),
                          self.binningZ()['bins'])
                         + self.binningZ()['range']),
                      file=output)

    return ParticularSlicesZ
