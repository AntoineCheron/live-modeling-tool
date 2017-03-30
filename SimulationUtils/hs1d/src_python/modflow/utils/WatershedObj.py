# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 11:19:40 2017

@author: Quentin Courtois
"""

class WatershedObj(object):
    """
        CLass containing geometric inputs for ModflowModel2D
    """

    def __init__(self, name, coord_x, coord_y, elevation, bottom, cell_size, outlet):
        self.model_name = name
        self.xcustom = coord_x
        self.ycustom = coord_y
        self.top_custom = elevation
        self.bot_custom = -1
        self.cell_size = cell_size
        self.outlet = outlet