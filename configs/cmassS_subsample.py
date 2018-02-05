from cmassS import cmassS

class cmassS_subsample(cmassS):
    def binningRA(self): return {"bins": 220, "range":(0,10)}
    def binningDec(self): return {"bins":156, "range":(0,6)}
