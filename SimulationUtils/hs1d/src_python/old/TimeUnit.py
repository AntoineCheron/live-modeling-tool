
# coding: utf-8

# In[ ]:

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from datetime import timedelta 


# In[ ]:

class TimeUnit(object):
    
    def time_to_second(t):
    # transforms travel times in time_unit to seconds

        day = timedelta(t)
        seconds = day.total_seconds()

        return seconds

