import lasspia as La
import math
from cmassS import cmassS

class cmassS_coarse(cmassS):

    def binningZ(self): return {"bins":300, "range":(0.43,0.7)}
    def binningRA(self): return {"bins": 750, "range":(-50,50)}
    def binningDec(self): return {"bins":275, "range":(-10,20)}
    def binningTheta(self): return {"bins":1000, "range":(0,math.pi)}
