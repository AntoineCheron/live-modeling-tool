# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:28:02 2017

@author: Quentin Courtois
"""

class MorphologicInputs(object):
    """
        Class storing morphologic inputs of BoussinesqSimulation:

            @param
              nx : number of cells along the hillslope
              xmin : minimal coordinate of the hillslope
              xmax : maximal coordinate of the hillslope
              discretization_type : type of wanted discretization
              x_custom : vector containing cells coordinates if discretization_type
                is custom
              angle : contains slope value along the hillslope
              z_custom : contains altitude values along the hillslope
              w : contains width values along the hillslope
    """

    def __init__(self, nx=100, xmin=0, xmax=100, discretization_type='linear', x_custom=-1, angle=0, \
                 z_custom=-1, w=0):
        self.nx = nx
        self.xmin = xmin
        self.xmax = xmax
        self.discretization_type = discretization_type
        self.x_custom = x_custom
        self.angle = angle
        self.z_custom = z_custom
        self.w = w