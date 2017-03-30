import numpy as np


class ThresholdFunction(object):

    """ Function Threshold to smooth discontinuities (used on y)
    """

    def __init__(self,x):
        self.y = 1 - np.exp(-1000*(1-x))
        if np.sum(x>1)>0:
            coord = np.where(x>1)
            self.y[coord] = 0
