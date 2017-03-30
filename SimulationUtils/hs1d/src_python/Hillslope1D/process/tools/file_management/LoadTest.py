# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 11:18:31 2016

@author: Quentin Courtois
"""

import os
import csv
import numpy as np

class LoadTest(object):

    def __init__(self, folder):
        #Load time data
        file_t_name = folder + "/" + "t_res"
        with open(file_t_name, 'r') as file:
            temp = file.read()
        temp_split = temp.split()
        temp_float = [float(i) for i in temp_split]
        self.t_res = np.zeros((len(temp_float), 1))
        i = 0
        while i < len(self.t_res):
            self.t_res[i] = temp_float[i]
            i += 1

        self.S = self.load_ref_values("S", folder)
        self.Q = self.load_ref_values("Q", folder)
        self.QS = self.load_ref_values("QS", folder)
        self.recharge_chronicle = self.load_ref_values("recharge", folder)
        self.t_recharge = self.load_ref_values("time", folder)
        self.x_S = self.load_ref_values("x_S", folder)
        self.x_Q = self.load_ref_values("x_Q", folder)
        self.Smax = self.load_ref_values("Smax", folder)
        self.w_node = self.load_ref_values("width_S", folder)
        self.w_edges = self.load_ref_values("width_Q", folder)

    def load_ref_values(self, var, folder):
        file = folder + "/" + var
        data_list = list(csv.reader(open(file, 'r'), delimiter='\t'))
        values = np.zeros((len(data_list), len(data_list[0])))
        i = 0
        while i < len(data_list):
            temp = data_list[i]
            temp_float = [float(j) for j in temp]
            values[i, :] = np.reshape(temp_float, (1, len(temp_float)))
            i += 1
        return values



