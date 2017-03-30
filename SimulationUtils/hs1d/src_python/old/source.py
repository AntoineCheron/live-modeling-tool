
# coding: utf-8

# In[ ]:

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from datetime import timedelta 


# In[ ]:

class Source(object):
#period = 10

    def __init__(self, param[None,10]):
        self.period = param[1]
        self.pluvio_type = param[0]
    
    def source(self):

        if param[0] is None and param[1] is None:
            self.pluvio_type = 'steady'

        if (param[0] is not None or param[1] is not None) and param[0] == 'periodical':
            self.period = 10

        if param[0] == 'steady':
            self.period = 'inf'

        elif param[0] == 'random':
            self.period = 0

        elif param[0] == 'data_based':
            self.period = -1

        elif param[0] == 'periodical':
            
            self.period = TimeUnit.time_to_second(self.period)

        else:
            print("Unknown chronicle type")

        print(self.pluvio_type, self.period)

    source([None, None])
    source(['periodical', None])


# In[ ]:

def set_recharge_chronicle(t,recharge_rate):

    t = time_to_seconds(t)
    recharge_rate = 10
    
    if period == 'inf': 
        recharge_chronicle = recharge_rate * ones(time.get_Nt)
        
    elif period == 0:
        recharge_chronicle = set_random_recharge
        
    elif period == -1:
        recharge_chronicle = None
        
        print('Data has not been uploaded because they are not supplied in arguments with this function. \n Try to call function upload_recharge_chronicle \n')

        '''    else:
        [t_chronicle,~] = time.get_properties
        recharge_chronicle = recharge_rate*(1+cos(2*pi*t_chronicle/obj.period))'''


# In[ ]:

def compute_recharge_rate(t)

    if period == 'inf':
        Recharge = recharge_rate
        
    elif period == 1:
        
        [t_chronicle,~] = time.get_properties
        Recharge = interpn(t_chronicle, recharge_chronicle, t, 'linear')
    
    elif period == -1:
        
        [t_chronicle,~] = time.get_properties
        Recharge = interpn(t_chronicle, recharge_chronicle, t, 'nearest')
        
    else:
        
        Recharge = recharge_rate*(1+cos(2*pi*t/obj.period))


# In[ ]:

def set_random_recharge(t)
    recharge_chronicle=0


# In[ ]:

def upload_recharge_chronicle([Date, Pluvio, ratio_P_R])

    if(nargin<4):
        ratio_P_R = 1
        
    recharge_chronicle = ratio_P_R*Pluvio
    
    time=Date

