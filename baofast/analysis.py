
"""Configuration: Stage 1 (Preprocessing)"""

class analysis(object):

    def filesRandom(self):
        '''List of random catalog file names.'''
        pass
    
    def filesObserved(self):
        '''List of observed catalog file names.'''
        pass

    def binsZ(self):
        """Lower bin edges and uppermost edge of P_z"""
        pass

    def binsRA(self):
        """Lower bin edges and uppermost edge of R_ang(ra,dec): ra"""
        pass

    def binsDec(self):
        """Lower bin edges and uppermost edge of R_ang(ra,dec): dec"""
        pass


    def __init__(self):
        pass
    
    @property
    def name(self) : return self.__class__.__name__
        
