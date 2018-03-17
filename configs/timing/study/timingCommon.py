from cmassN import cmassN

class timingCommon(cmassN):

    def binningZ(self): return {"bins":621, "range":(0.43,0.7)}
    def binningDec(self): return {"bins":381, "range":(-4,8)}
    def chunkSize(self): return 4000


    '''Parameters for avoiding unnecessary combinatorial calculations at large s.
    Galaxies farther apart than these parameters may not be included in result.'''

    def maxDeltaRA(self): return 10.3
    def maxDeltaDec(self): return 10.3
    def maxDeltaZ(self): return 0.3
