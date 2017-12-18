import math
from cmassN import cmassN

class cmassN_coarse(cmassN):

    def binningZ(self): return {"bins":300, "range":(0.43,0.7)}
    def binningRA(self): return {"bins":1200 , "range":(105,265)}
    def binningDec(self): return {"bins":600, "range":(-4,57)}
    def binningTheta(self): return {"bins":1000, "range":(0,math.pi)}

    def maxDeltaRA(self): return None
    def maxDeltaDec(self): return None
    def maxDeltaZ(self): return None
