from cmassS_byZ import cmassS_byZ

class cmassS_subsample_byZ(cmassS_byZ):
    def binningRA(self): return {"bins": 220, "range":(0,10)}
    def binningDec(self): return {"bins":156, "range":(0,6)}
