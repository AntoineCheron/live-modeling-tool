import numpy as np


class InitialConditions(object):

    def __init__(self, percentage_loaded, w, soil_depth, f):

        self.percentage_loaded = percentage_loaded
        self.Smax = f *np.reshape(w,(len(w),1)) * np.reshape(soil_depth,(len(soil_depth),1))
        if isinstance(percentage_loaded, float) or isinstance(percentage_loaded, int):
            self.sin = percentage_loaded * np.reshape(self.Smax, (len(self.Smax), 1))
        else:
            self.sin = np.reshape(self.percentage_loaded,(len(self.percentage_loaded),1)) * np.reshape(self.Smax, (len(self.Smax), 1))
        sup = np.where(self.sin>self.Smax)
        self.sin[sup[0]] = self.Smax[sup[0]]
