# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:38:46 2017

@author: Quentin Courtois
"""

class GeologicInputs(object):
    """
        Class storing geologic inputs for BoussinesqSimulation

        @param
          boundary_type : type of the boundary on bottom and top of the hillslope
          boundary_value : value for each boundary
          k : hydrauic conductivity
          f : porosity
          soil_depth : thickness of each cell

    """

    def __init__(self, boundary_type=['S','Q'], boundary_value=[0,0], k=1/3600, \
                 f=0.3, soil_depth=0):
        self.boundary_type = boundary_type
        self.boundary_value = boundary_value
        self.k = k
        self.f = f
        self.soil_depth = soil_depth