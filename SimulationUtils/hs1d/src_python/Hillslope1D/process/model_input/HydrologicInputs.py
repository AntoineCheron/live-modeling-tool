# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:43:17 2017

@author: Quentin Courtois
"""

class HydrologicInputs(object):
    """
        Class storing all hydrologic inputs for BoussinesqSimulation

        @param
          tmin : minimal value of the time serie
          tmax : maximal value of the time serie
          Nt : number of values to consider in the time serie
          unit : unit of the time serie
          recharge_rate : either constant recharge rate or recharge chronicle
          time_custom : custom time serie
          period : period used to compute recharge chronicle if no known custom
            recharge
          recharge_type : determine the type of recharge to consider
          perc_loaded : initial state for each cell (% filled)
    """

    def __init__(self, tmin=0, tmax=35, Nt=35*24*10, unit='days', recharge_rate=30, \
                 time_custom=-1, period=None, recharge_type='square', \
                 perc_loaded=0):
        self.tmin = tmin
        self.tmax = tmax
        self.Nt = Nt
        self.unit = unit
        self.recharge_rate = recharge_rate
        self.time_custom = time_custom
        self.period = period
        self.recharge_type = recharge_type
        self.perc_loaded = perc_loaded