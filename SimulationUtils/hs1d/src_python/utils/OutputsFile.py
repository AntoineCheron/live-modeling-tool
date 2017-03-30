# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 11:43:09 2017

@author: Quentin Courtois
"""

import numpy as np

class OutputsFile(object):
    """
        Class loading outputs results for Hs1D model for a specified test case
        @parameters
            - outputs_list : a list containing the path to each file to load
    """

    def __init__(self, outputs_list):
        self.seepage, self.x_edges = self.load_seepage(outputs_list[0])
        self.storage, self.x_nodes = self.load_storage(outputs_list[1])
        self.subsurface_flux = self.load_subsurface_flux(outputs_list[2])
        self.q_river_tot, self.q_river_channel, self.runoff = self.load_river_runoff(outputs_list[3])
        self.perc_storage = self.load_relative_storage(outputs_list[4])

    def load_seepage(self, file):
        """
            load seepage results (x,t) (m²/s)
        """
        with open(file, 'r'):
            temporary = np.loadtxt(file, skiprows=2)
            seepage = temporary[:, 1:]
            x_edges = temporary[:, 0]

        return seepage, x_edges

    def load_storage(self, file):
        """
            load storage results (x,t) (m²)
        """
        with open(file, 'r'):
            temporary = np.loadtxt(file, skiprows=2)
            storage = temporary[:, 1:]
            x_nodes = temporary[:,0]

        return storage, x_nodes

    def load_subsurface_flux(self, file):
        """
            load subsurface flux results (x,t) (m3/s)
        """
        with open(file, 'r'):
            temporary = np.loadtxt(file, skiprows=2)
            subsurface_flux = temporary[:, 1:]

        return subsurface_flux

    def load_river_runoff(self, file):
        """
            load river  flow and runoff (t) (m3/s)
        """
        with open(file, 'r'):
            temporary = np.loadtxt(file, skiprows=1)
            q_river_tot = temporary[:, 1]
            q_river_channel = temporary[:, 2]
            runoff = temporary[:, 3]

        return q_river_tot, q_river_channel, runoff

    def load_relative_storage(self, file):
        """
            load relative storage results (x,t) (in %)
        """
        with open(file, 'r'):
            temporary = np.loadtxt(file, skiprows=2)
            perc_storage = temporary[:, 1:]

        return perc_storage
