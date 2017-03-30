# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 11:24:21 2017

@author: Quentin Courtois
"""

class HydroObj(object):
    """
        Class containing hydrologic inputs for ModflowModel2D
    """
    def __init__(self, hk, sy, ss, soil_depth, percentage_loaded, recharge_initial, \
                 recharge_true, time_select):
        self.hk = hk
        self.sy = sy
        self.ss = ss
        self.soil_depth = soil_depth
        self.percentage_loaded = percentage_loaded
        self.recharge_initial = recharge_initial
        self.recharge_true = recharge_true
        self.time_select = time_select