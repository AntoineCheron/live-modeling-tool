
# coding: utf-8

# In[ ]:

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from datetime import timedelta 


# In[ ]:

class time_properties()
   

    def __init__(self, tmin, tmax, Nt, unit):
        
        self.tmin = tmin
        self.tmax = tmax
        self.Nt = Nt
        self.t = self.tmin:(self.tmax-self.tmin)/(Nt-1):self.tmax
        self.unit = unit
        
        tmin = time_to_seconds(tmin, unit)
        tmax = time_to_seconds(tmax, unit)

    def get_tmin(self):
        tmin = self.tmin
    
    def get_tmax(self):
        tmax = self.tmax

    def get_unit(self):
        unit = self.unit

    def get_Nt(self):
        Nt=self.Nt
    
    def get_size(self):
        size_ = size(self.t)

    def [t, unit] = get_properties(self):
        t = self.t
        unit=self.unit

