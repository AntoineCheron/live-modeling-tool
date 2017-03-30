# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 11:42:52 2017

@author: Quentin Courtois
"""

import numpy as np

class InputsFile(object):
    """
        Class loading input paramters for Hs1D model based on test cases
        @parameters
            - inputs_list : a list containing the path to each file to load
    """

    def __init__(self, inputs_list, interp):
        self.f, self.k, self.soil_depth = self.load_geol(inputs_list[0])
        self.time, self.recharge_chronicle = self.load_hydro(inputs_list[1])
        self.x, self.w, self.i, self.z_true, self.z_mod = self.load_morpho(inputs_list[2])
        if interp != 0:
            self.x_interp, self.w_interp, self.Smax_interp = self.load_spatial(inputs_list[3])

        #Reformating matrices and changing units
        self.k = self.k/3600                                   #from m/h to m/s
        self.x = np.reshape(self.x,(len(self.x),1))
        self.w = np.squeeze(np.reshape(self.w,(len(self.w),1)))
        self.i = np.squeeze(np.reshape(self.i,(len(self.i),1)))
        self.time = np.squeeze(np.reshape(self.time,(len(self.time),1)))
        self.recharge_chronicle = np.squeeze(np.reshape(self.recharge_chronicle,(len(self.recharge_chronicle),1)))


    def load_geol(self, file):
        """
            load geologic input file and read parameters
        """
        with open(file, 'r'):
            temporary = np.loadtxt(file, skiprows=3)
            f = temporary[0]
            k = temporary[1]
            soil_depth = temporary[2]

        return f, k, soil_depth

    def load_hydro(self, file):
        """
            load hydrologic input file and read parameters
        """
        with open(file, 'r'):
            temporary = np.loadtxt(file, skiprows=3)
            time = temporary[:, 0]
            recharge_chronicle = temporary[:, 1]

        return time, recharge_chronicle

    def load_morpho(self, file):
        """
            load morphologic input file and read parameters
        """
        with open(file, 'r'):
            temporary = np.loadtxt(file, skiprows=3)
            x = temporary[:, 0]
            w = temporary[:, 1]
            i = temporary[:, 2]
            z_true = temporary[:, 3]
            z_mod = temporary[:, 4]

        return x, w, i, z_true, z_mod

    def load_spatial(self, file):
        """
            load interpolated spatial parameters
        """
        with open(file, 'r'):
            temporary = np.loadtxt(file, skiprows=1)
            x_interp = temporary[:,0]
            w_interp = temporary[:,3]
            Smax_interp = temporary[:,2]

        return x_interp, w_interp, Smax_interp