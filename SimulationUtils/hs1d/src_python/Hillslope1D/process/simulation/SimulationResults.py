# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 13:38:42 2016

@author: Quentin Courtois
"""
import numpy as np

class SimulationResults(object):
    """
        Class storing the results of the integration (using IDA solver)

        @param:
          - t_res : times for which the system is solved (s)
          - y_res : the output of the solver for each time (m², m3/s, m²/s)
          - N_nodes : number of cells in the system (-)
          - N_edges : number of cell edges in the system (=N_nodes + 1) (-)
          - x_node : coordinates of cells centers (m)
          - x_edges : coordinates of cells edges (m)
          - recharge : recharge applied to the system (m/s)
          - dx : size of the cells along the hillslope

        @attributes:
          - t : times for which the system is solved (s)
          - S : stock in each cell as a function of time(lines) and x(columns) (m²)
          - Q : flowrate at each cell edge as a function of time(lines) and x(columns) (m3/s)
          - QS : seepage in each cell as a function of time(lines) and x(columns) (m3/s)
          - x_node : coordinates of cells centers (m)
          - x_edges : coordinates of cells edges (m)
    """

    def __init__(self, t_res, y_res, N_nodes, N_edges, x_node, x_edges, recharge, dx):
        self.t = t_res
        self.S = y_res[:, 0:N_nodes]
        self.Q = y_res[:, N_nodes : N_nodes + N_edges]
        self.QS = y_res[:, N_nodes + N_edges : 2*N_nodes + N_edges]
        #Calculating the seepage in m3/s
        for i in range(len(recharge)+1):
          self.QS[i,:] = self.QS[i,:]*np.reshape(dx,(1,len(dx)))
        self.x_node = x_node
        self.x_edges = x_edges
