# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 11:39:07 2016

@author: Quentin Courtois
"""

import numpy as np
import os


class WatershedReader(object):
    """
        Class designed to read files written using matlab in which are listed
        points composing a watershed identified by its name and its outlet.
        Files are extracted using numpy in order to create numpy arrays as
        attributes of self.
    """

    def __init__(self,name):
        self.name = name
        cwd = os.getcwd() + "\\" + name

        #Load X coordinates
        filename = cwd + "\\" + "coord_X_" + name
        with open(filename,'r'):
            coord_x = np.loadtxt(filename)
            self.coord_x = coord_x.astype(np.float)

        #Load Y coordinates
        filename = cwd + "\\" + "coord_Y_" + name
        with open(filename,'r'):
            coord_y = np.loadtxt(filename)
            self.coord_y = coord_y.astype(np.float)

        #Load Watershed's elevations
        filename = cwd + "\\" + name + "_elevation"
        with open(filename,'r'):
            elevation = np.loadtxt(filename)
            self.elevation = elevation.astype(np.float)

        #Load Watershed's outlet coordinates
        filename = cwd + "\\" + name + "_outlet_coord"
        with open(filename,'r'):
            outlet = np.loadtxt(filename)
            self.outlet = outlet.astype(np.float)

        #Load Stream's X coordinates
        filename = cwd + "\\" + "coord_X_" + name + "_stream_network"
        if os.path.isfile(filename):
            with open(filename,'r'):
                stream_x = np.loadtxt(filename)
                self.stream_x = stream_x.astype(np.float)

        #Load Stream's Y coordinates
        filename = cwd + "\\" + "coord_Y_" + name + "_stream_network"
        if os.path.isfile(filename):
            with open(filename,'r'):
                stream_y = np.loadtxt(filename)
                self.stream_y = stream_y.astype(np.float)

        #Load Stream's elevation
        filename = cwd + "\\" + name + "_stream_network_elevation"
        if os.path.isfile(filename):
            with open(filename,'r'):
                stream_elevation = np.loadtxt(filename)
                self.stream_elevation = stream_elevation.astype(np.float)


