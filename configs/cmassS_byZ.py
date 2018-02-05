
from cmassS import cmassS
from lasspia.zSlicing import zSlicing

class cmassS_byZ( zSlicing(cmassS) ):

    def zBreaks(self): return [0.43, 0.58, 0.7]
    def zMaxBinWidths(self): return [0.0003, 0.0003]
