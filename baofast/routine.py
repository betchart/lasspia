
class routine(object):
    
    def __init__(self, config):
        self.config = config  # subclass of baofast.config

    def __call__(self):
        '''Defined in subclasses.'''
        pass
