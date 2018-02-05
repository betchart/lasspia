from cmassS_subsample import cmassS_subsample
from lasspia.zSlicing import zSlicing

class cmassS_subsample_byZ(zSlicing(cmassS_subsample)):

    # zSlicing
    def zBreaks(self): return [0.43, 0.58, 0.7]
    def zMaxBinWidths(self): return [0.0003, 0.0003]
