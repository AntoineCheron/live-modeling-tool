import numpy as np

class BoundaryConditions(object):
    """
        Class defining boundary conditions in xmin and xmax as S, Q or free
        (imposed head,impose flux or free). This properties are stored in two 2
        valued arrays containing boundary types and boundary values. First
        Boundary condition is located downstream (on the river) and the other
        one is located upstream (on the ridge of the hillslope)
        #######################################################################
        INPUTS :
           boundary_type : a two-valued list containing the types corresponding
           to the two boundary_conditions : either S for a constant stock or
           Q for a forced flowrate
           boundary_value : a two-valued list containing the values
           corresponding to each limit. In m/s for Q
        #######################################################################
    """

    def __init__(self, boundary_type, boundary_value):
        self.boundary_type = boundary_type
        self.boundary_value = boundary_value
        self.edges = np.zeros(shape=(4, 1))
        self.edges_bool = np.zeros(shape=(4, 1))

    def fixed_edge_matrix_boolean(self):
        if self.boundary_type[0] == 'S':
            # Start with a fixed charge S for xmin and forced flux Q for xmax
            self.edges_bool[0] = 1
            self.edges_bool[2] = 1
        elif self.boundary_type[0] == 'Q':
            self.edges_bool[2] = 1

        if self.boundary_type[1] == 'S':
            self.edges_bool[1] = 1
            self.edges_bool[3] = 1
        elif self.boundary_type[1] == 'Q':
            self.edges_bool[3] = 1

        return self.edges_bool

    def fixed_edge_matrix_values(self):
        if self.boundary_type[0] == 'S':
            self.edges[0] = self.boundary_value[0]
        elif self.boundary_type[0] == 'Q':
            self.edges[2] = self.boundary_value[0]

        if self.boundary_type[1] == 'S':
            self.edges[1] = self.boundary_value[1]
        elif self.boundary_type[1] == 'Q':
            self.edges[3] = self.boundary_value[1]

        return self.edges
