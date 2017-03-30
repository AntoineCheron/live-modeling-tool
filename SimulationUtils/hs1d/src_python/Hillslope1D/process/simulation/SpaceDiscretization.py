import numpy as np
from scipy.interpolate import UnivariateSpline
import Hs1D as Hs
import copy

class SpaceDiscretization(object):
    """
        Class defining the discretization in space of the hillslope (linear, log, etc....)
        It's detirmined using X minimal and maximal values, soil depth, angle and width
        of the hillslope for each mesh

        @param :
          - xmin : minimal coordinate of the hillslope
          - xmax : maximal coordinate of the hillslope
          - nx : number of cells to define (if x_custom == -1)
          - discretization_type : type of discretization used : 'linear', 'logarithmic',
              'square' or 'custom'
          - x_custom : -1 or vector of coordinates if discretization_type == 'custom'
          - angle : float or vector of values corresponding to the slope
          - w : float or vector of values corresponding to the width
          - soil_depth : float or vector of values corresponding to the thickness of the layer
          - k : hydraulic conductivity
          - f: porosity
          - z_custom : -1 or a vector altitude of each cell (used to compute the slope
               if x_custom is a vector)

        @attributes:
          - xmin : minimal coordinate of the hillslope
          - xmax : maximal coordinate of the hillslope
          - discretization : type of discretization used
          - xcustom : -1 or vector of coordinates if discretization_type == 'custom'
          - b : computation matrix
          - a : computation matrix
          - omega : computation matrix
          - omega2 : computation matrix
          - x_node : coordinates of cells nodes
          - angle_node : slope of cells on nodes
          - w_node : width of cells on nodes
          - N_node : number of cells
          - soil_depth_node : layer thickness on nodes
          - dx_node : distance between two successive nodes
          - x_edges : coordinates of cells edges
          - N_edges : number of cells edges
          - dx_edges : distance between two successive edges
          - Hs : class Hs1D
    """

    def __init__(self, xmin=0, xmax=100, nx=100, discretization_type='linear', x_custom=-1, \
                 angle=0, w=0, soil_depth=0, k=1/3600, f=0.3, z_custom=-1):
        self.xmin = xmin
        self.xmax = xmax
        self.discretization = discretization_type # linear, log, or custom
        self.xcustom = x_custom
        self.x_edges = self.space_discretization()
        self.Hs = Hs.Hs1D(nx, angle, w, soil_depth, k, f, z_custom, self.x_edges)
        self.x_node = self.compute_x_node()
        self.set_matrix_properties()

    def space_discretization(self):

        if self.discretization == 'linear':
            self.x_edges = np.arange(self.xmin, self.xmax, (self.xmax - self.xmin) / self.N_edges)
        elif self.discretization == 'logarithmic':
            if self.xmin == 0:
                log_xmin = -2
                self.xmin = 10**(log_xmin)
            else:
                log_xmin = np.log10(self.xmin)
            log_xmax = np.log10(self.xmax)
            self.x_edges = np.logspace(log_xmin, log_xmax, self.N_edges)
        elif self.discretization == 'custom':
            self.x_edges = self.xcustom
            self.N_edges = len(self.x_edges)
            self.N_nodes = self.N_edges - 1
            self.xmin = min(self.x_edges)
            self.xmax = max(self.x_edges)
        elif self.discretization == 'square':
            self.x_edges = np.arange(self.xmin**0.5,self.xmax**0.5, \
                                     (self.xmax**0.5 - self.xmin**0.5)/self.N_edges)
            self.x_edges = self.x_edges**2
        else:
            print('no discretization type corresponding to ', self.discretization, '\n')

        return self.x_edges

    def resample_hs1D_spatial_variables(self):
        soil_depth = self.Hs.get_soil_depth_edges()
        w = self.Hs.get_w_edges()
        angle = self.Hs.get_angle_edges()
        self.x_node = self.compute_x_node()
         #UnivariateSpline does the same work as fit ('smoothingspline') in matlab
        smooth_width_function = UnivariateSpline(self.x_edges, w, k=4, s=0.9)
        self.w_node = smooth_width_function(self.x_node)
        smooth_slope_function = UnivariateSpline(self.x_edges, angle, s=0.9)
        self.angle_node = smooth_slope_function(self.x_node)
        self.soil_depth_node = np.interp(np.squeeze(self.x_node),np.squeeze(self.x_edges), soil_depth[0, :])
        return self.w_node, self.angle_node, self.soil_depth_node

    def get_angle_node(self):
        self.resample_hs1D_spatial_variables()
        return self.angle_node

    def get_w_node(self):
        self.resample_hs1D_spatial_variables()
        return self.w_node

    def get_soil_depth_node(self):
        self.resample_hs1D_spatial_variables()
        return self.soil_depth_node

    def set_matrix_properties(self):
        self.a = self.first_derivative_downstream()
        self.b = self.first_derivative_upstream()
        self.omega = self.weight_matrix()
        self.omega2 = self.weight_matrix_bis()
        return self.a
        return self.b
        return self.omega
        return self.omega2

    def compute_x_node(self):
        self.x_node = (self.x_edges[1:] + self.x_edges[0:-1])/2
        return self.x_node

    def compute_dx_edges(self):
        # length(dx) = Nx
        self.dx_edges = self.x_edges[1:]-self.x_edges[0:-1]
        return self.dx_edges

    def compute_dx_node(self):
        # length(dxS) = Nx - 1
        self.dx_node = self.x_node[1:]-self.x_node[0:-1]
        return self.dx_node

        ##########################################################################
        ###### stencil matrix(for basic derivation centered, downstream, upstream)
        ###### + mean matrix (Omega)######
        ##########################################################################

    def first_derivative_upstream(self):
        #######################################
        #         [ 0   ...   0   0   ...   0 ]
        #         [-1    1    0   0   ...   0 ]
        #    B =  [ 0  - 1    1   0   ...   0 ]
        #         [... ...   ... ...  ...  ...]
        #         [ 0  ...   ...  0   - 1   1 ]
        #         [ 0  ...    0   0   ...   0 ]
        #size B (Nx + 1) x(Nx)
        #######################################

        self.a = self.first_derivative_downstream()
        a = copy.copy(self.a[:-1, :-1])
        a[a > 0] = 1
        a[a < 0] = -1
        temp = np.zeros((1, (len(self.x_edges)-1)))
        self.b = np.append(temp, a, axis=0) #+ np.zeros(((len(self.x) - 1),1))
        self.b = np.append(self.b, temp, axis=0)
#        #B = [zeros(1, length(obj.x) - 1); A; zeros(1, length(obj.x) - 1)];
        dx_node = self.compute_dx_node(),
        dx_node = np.hstack((0, np.squeeze(dx_node), 0))
        self.b = np.dot(np.diag(1/dx_node), self.b)
        self.b[np.isnan(self.b)] = 0
        # self.B = np.eye(1 / self.dxS) * self.B
        #self.b = csr_matrix(self.b)
        return self.b

    def first_derivative_downstream(self):
        #########################################
        #       [-1   1   0    0    0   ...   0 ]
        # A =   [ 0 - 1   1    0    0   ...   0 ]
        #       [... ... ...  ...  ...  ...  ...]
        #       [ 0  ... ...  ...   0  - 1    1 ]
        #       size A = (Nx) x(Nx + 1)
        #########################################
        self.dx_edges = self.compute_dx_edges()
        self.a = -np.eye(len(self.x_edges)-1)
        temp = np.zeros(((len(self.x_edges)-1), 1))
        self.a = np.hstack((self.a, temp))
        a1 = np.eye(len(self.x_edges) - 1)
        a2 = np.zeros((len(self.x_edges) - 1, 1))
        abis = np.hstack((a2, a1))
        self.a = self.a + abis
        self.a = np.dot(np.diag(np.squeeze(1/self.dx_edges)), self.a)
#        print(np.size(self.a,0),'  ',np.size(self.a,1))
        # self.A = (np.eye(np.ones(((len(self.x)-1,1)))))*self.A
        #self.A = csr_matrix(self.A)

        return self.a

    def first_derivative_centered(self):
        ######################################################
        #           [ 0      1       0    0   ...           0]
        # C =       [ -1     0       1    0   ...           0]
        #           [ ...   ...     ...                     1]
        #           [ 0     ...     ...        0    - 1     0]
        ######################################################

        self.c = -np.diag((np.ones((len(self.x_edges) - 1, 1)), -1))
        self.c = self.c + np.diag(np.ones((len(self.x_edges) - 1, 1)), 1)
        self.c = np.diag(1 / self.dx_edges) * self.c

        return self.c

    def weight_matrix(self):
        self.omega = copy.copy(self.b)
        self.omega[(self.omega > 0)] = 0.5
        self.omega[(self.omega < 0)] = 0.5
        return self.omega

    def weight_matrix_bis(self):
        self.omega2 = copy.copy(self.b)
        self.omega2[self.omega2 < 0] = 0
        return self.omega2

