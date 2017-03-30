import numpy as np
import Source as So
import BoundaryConditions as BC
import ThresholdFunction as Tf
import SpaceDiscretization as SD
import InitialConditions as IC
import SimulationResults as SR
from assimulo.problem import Implicit_Problem
from assimulo.solvers.sundials import IDA

class BoussinesqSimulation(object):
#class storing the different object in order to make run a boussinesq
#simulation for different conditions
    """
    ###########################################################################
    Other classes called by BoussinesqSimulation

    - discretization : space_discretization object contains also the spatial
      properties of your hillslopes (x, width function w and soil_depth)
          named : SD
    - source_terms : source object containing the time properties and the
      recharge related to it
          named : So
    - boundary_conditions : boundary object containing the time of boundary on
      Morpho.xmin & Morpho.xmax (Q fixed, S fixed or free condtions)
          named : BC
    - initial_conditions : initial conditions for S, Q & QS
          named : IC
    - Integration results are contained in self.t_res, self.y_res, self.yd_res
          named : SR

    ###########################################################################
    INPUT format and units:

      @param:
        Inputs contained in 3 different classes
        - Geol : k, f and soil_depth (parameters)
        - Morpho : nx, xmin, xmax, discretization_type, x_custom, angle,
                   z_custom, w (geometry : spatial aspect)
        - Hydro : tmin, tmax, Nt, unit, recharge_rate, time_custom, period,
                  recharge_type, perc_loaded (temporal aspect)


            WARNING ! ALL CALCULATIONS ARE DONE IN M/S

    - k : hydraulic conductiviy (m/s) (ranges from e-8 to e-2 m/s)
    - f : porosity (in %) (ranges from e-5 to 0.3)
    - soil_depth : thickness of the layer (m) (ranges from 2 to 50)
    - period : used only if recharge_type is not data based : in the unit as t
    - recharge_type : either 'square','periodical','steady','random','databased'
    - recharge_rate : only if not data_based : in mm/day

    - tmin : minimal time of the serie
    - tmax : maximal time of the serie
    - Nt : Number of time steps in the time serie
    - unit : time unit of the serie

    - Geol.boundary_type : type of the two boundary conditions, first one corresponds
      to the boundary at the level of the river, second at the top of the
      hillslope
          Example : ['S','Q' ] with S : imposed level and Q : imposed flow rate
    - boundary_values : values corresponding to the two boundary conditions
          Example : [0,0] in m and in m/sec

    - Morpho.xmin : minimal coordinate of the hillslope (m)
    - Morpho.xmax : maximal coordinate of the hillslope (m)
    - discretization space along the hillslope (int) : default is 100
    - Morpho.discretization_type : type of the segmentation of the hillslope :
        'linear', 'logarithmic','square','custom'
    - Morpho.x_custom : used only for custom segmentation : must contain edges coordinates
    - Morpho.angle : Morpho.angle of the hillslope (rad)
    - w : width of the hillslope for each x
    - soil_depth : thickness of the media for each x
    - k : the kincematic porosity (in % : 30% = 0.3)
    - f : the hydraulic conductivity (m/second)

    - Id : the identification of the hillslope
    - percentage_loaded : initial charge in the hillslope (a percentage of the
      total capacity as 100% = 1)
    ###########################################################################
    """
#    def __init__(self, period=None, recharge_type='square', recharge_rate=30, \
#                 tmin=0, tmax=35, Nt=35*24*10, unit='days', time_custom=-1, \
#                 Geol.boundary_type=['S', 'Q'], boundary_value=[0, 0], \
#                 Morpho.xmin=0, Morpho.xmax=100, Morpho.nx=100, Morpho.discretization_type='linear', Morpho.x_custom=-1, \
#                 Morpho.z_custom=-1, Morpho.angle=0, w=0, soil_depth=0, k=1/3600, f=0.3, \
#                 Id=1, percentage_loaded=0):

    def __init__(self, Morpho, Geol, Hydro, Id):
        #List of INPUTS
        #1st row : INPUT for Source
        #2nd row : INPUT for Source.TimeProperties
        #3rd row : INPUT for BoundaryConditions
        #4th row : INPUT for SpaceDiscretization
        #5th row : INPUT for SpaceDiscretization.Hs1D
        #6th row : INPUT for InitialConditions
        Morpho.xmax = Morpho.xmax+Morpho.xmax/Morpho.nx
        if not isinstance(Morpho.x_custom, int):
            Morpho.nx = np.size(Morpho.x_custom) - 1
            Morpho.xmax = np.max(Morpho.x_custom)
            Morpho.xmin = np.min(Morpho.x_custom)

        #Build Source Class
        self.So = So.Source(Hydro.period, Hydro.recharge_type, Hydro.recharge_rate, \
                            Hydro.tmin, Hydro.tmax, Hydro.Nt, Hydro.unit, Hydro.time_custom)
        #Build SpaceDiscretization Class
        self.SD = SD.SpaceDiscretization(Morpho.xmin, Morpho.xmax, Morpho.nx, Morpho.discretization_type, Morpho.x_custom, \
                                         Morpho.angle, Morpho.w, Geol.soil_depth, Geol.k, Geol.f, Morpho.z_custom)

        #Get width and soild_depth over each node
        self.SD.resample_hs1D_spatial_variables()

        #Set Hydraulic Properties
        self.k = self.SD.Hs.get_k()
        self.f = self.SD.Hs.get_f()

        #Build Boundary Conditions
        h0 = self.f*self.SD.w_node[0]*self.SD.soil_depth_node[0]
        Geol.boundary_value[0] = h0
        self.BC = BC.BoundaryConditions(Geol.boundary_type, Geol.boundary_value)

        #Build Initial Conditions
        self.set_initial_conditions(percentage_loaded=Hydro.perc_loaded, w=self.SD.w_node, \
                                    soil_depth=self.SD.soil_depth_node, f=self.SD.Hs.f)

        #Identification of the hillslope
        self.Id = Id

        #Mass Matrix of the DAE
        self.m = self.compute_mass_matrix()

    def set_initial_conditions(self, percentage_loaded, w, soil_depth, f):
        """
            Set initial conditions of stock and flow over the hillslope based
            on percentage_loaded and boundary conditions
        """
        self.IC = IC.InitialConditions(percentage_loaded, w, soil_depth, f)

        boundary = self.BC
        edges = boundary.fixed_edge_matrix_values()
        edges_bool = boundary.fixed_edge_matrix_boolean()

        self.IC.sin[0] = edges_bool[0] * edges[0] + (1 - edges_bool[0]) * self.IC.sin[0]
        self.IC.sin[-1] = edges_bool[1] * edges[1] + (1 - edges_bool[1]) * self.IC.sin[-1]

        self.IC.qin = np.dot(-self.compute_q_from_s(self.IC.sin), self.IC.sin)
        self.IC.qin[0] = edges_bool[2] * edges[2] + [1 - edges_bool[2]] * self.IC.qin[0]
        self.IC.qin[-1] = edges_bool[3] * edges[3] + [1 - edges_bool[3]] * self.IC.qin[-1]

        self.IC.q_sin = np.dot(self.compute_qs_from_q(np.vstack((self.IC.sin, \
                        self.IC.qin, np.zeros((len(self.IC.qin), 1)))), self.So.TP.tmin), \
                        self.IC.qin)
        temp = np.where(self.IC.q_sin<0)
        self.IC.q_sin[temp[0],temp[1]] = 0
        return self

    def compute_c(self, y, t):
        """
            Compute a part of the DAE base on S, Q and QS values
        """
        ############################################
        #       [ 0      - alpha * beta * A      0 ]
        # C =   [P(y)       I                    0 ]
        #       [ 0      - (1 - alpha) * A     - I ]
        ############################################
        #Initialzing
        N_edges = self.SD.N_edges
        N_nodes = self.SD.N_nodes

        c = np.zeros((2 * N_nodes + N_edges, 2 * N_nodes + N_edges))

        dsdt_from_q = self.compute_dsdt_from_q(y, t)
        q_from_s = self.compute_q_from_s(y)
        qs_from_q = self.compute_qs_from_q(y, t)

            # first row block
        c[0:N_nodes, 0:N_nodes]                                 = 0
        c[0:N_nodes, N_nodes:(N_nodes + N_edges)]               = dsdt_from_q
        c[0:N_nodes, (N_nodes + N_edges):(2*N_nodes + N_edges)] = 0

            # second row block
        c[N_nodes :(N_nodes + N_edges), 0:N_nodes]                               = q_from_s
        c[N_nodes :(N_nodes + N_edges), N_nodes:N_nodes + N_edges]               = np.eye(N_edges)
        c[N_nodes :(N_nodes + N_edges), (N_nodes + N_edges):2*N_nodes + N_edges] = 0

            # third row block
        c[(N_nodes + N_edges):2*N_nodes + N_edges, 0:N_nodes]                               = 0
        c[(N_nodes + N_edges):2*N_nodes + N_edges, (N_nodes):N_nodes + N_edges]             = qs_from_q
        c[(N_nodes + N_edges):2*N_nodes + N_edges, (N_nodes + N_edges):2*N_nodes + N_edges] = -np.eye(N_nodes)

            # Introduce boundary conditions in the matrix
        edges = self.BC.fixed_edge_matrix_boolean()
        c[0, :] = (1-int(edges[0]))*c[0, :]
        c[N_nodes-1, :] = (1-edges[1])*c[N_nodes-1, :]

        return c

            #Darcy equation
    def compute_q_from_s(self, y):
        """
            Compute the flow rate in each cell (over the edges) using conversion
            matrices and stock value
        """
        N_nodes = self.SD.N_nodes
        tempo = y[:N_nodes]/(self.f*np.reshape(self.SD.w_node, (N_nodes, 1)))
        q_from_s = (self.k/self.f)*(np.cos(self.SD.Hs.angle_edges)*np.dot(self.SD.b, tempo) \
                    + np.sin(self.SD.Hs.angle_edges))
        q_from_s = np.diag(q_from_s[:, 0])
        q_from_s = np.dot(q_from_s, self.SD.omega)

        ## put boundary conditions on Q in the matrix
        edges = self.BC.fixed_edge_matrix_boolean()
        q_from_s[0, :] = (1-edges[2])*q_from_s[0, :]
        q_from_s[N_nodes, :] = (1-edges[3])*q_from_s[N_nodes, :]

        return q_from_s

    def compute_qs_from_q(self, y, t):
        """
            Compute Seepage in each cell (over nodes) based  on flow rate and
            conversion matrices
        """
        alpha = self.compute_alpha(y, t)
        alpha_complementar = np.diag(1 - alpha[:, 0])
        qs_from_q = np.dot(-alpha_complementar, self.SD.a)

        return qs_from_q

    def compute_dsdt_from_q(self, y, t):
        """
            Compute stock variation between two time steps based on flow rate and
            conversion matrices
        """
        alpha = np.diag(np.squeeze(self.compute_alpha(y, t)))
        beta = -np.diag(np.squeeze(self.compute_beta(y,t)))
        dsdt_from_q = np.dot(np.dot(beta, alpha), self.SD.a)

        return dsdt_from_q


    def test_derivative(self, y, t):
        """
            Test to determine if stock is still positive and seepage is occuring
            or not
        """
        N_edges = self.SD.N_edges
        N_nodes = self.SD.N_nodes

        recharge = self.So.compute_recharge_rate(t)
        recharge_rate_spatialized = recharge*np.reshape(self.SD.w_node, (len(self.SD.w_node), 1))

        test_deriv = np.dot(-self.SD.a, y[N_nodes: N_nodes + N_edges]) + recharge_rate_spatialized
        test_deriv = test_deriv >= 0

        return test_deriv


    def compute_alpha(self, y, t):
        """
            compute a matrix defining variations over time. Used to compute Q,
            S and QS
        """
        N_nodes = self.SD.N_nodes
        self.SD.w_node = np.reshape(self.SD.w_node, (len(self.SD.w_node), 1))
        self.SD.soil_depth_node = np.reshape(self.SD.soil_depth_node, \
                                             (len(self.SD.soil_depth_node), 1))
        xt = y[0:N_nodes]/(self.f* self.SD.soil_depth_node * self.SD.w_node)
        self.Tf = Tf.ThresholdFunction(xt)
        test_deriv = self.test_derivative(y, t)
        alpha = self.Tf.y * test_deriv + (1 - test_deriv)
        return alpha

    def compute_beta(self,y,t):
        """
            compute a matrix defining variations over time. Used to compute Q,
            S and QS
        """
        N_nodes = self.SD.N_nodes
        beta = np.ones((N_nodes, 1))
#        y_new = y[0:N_nodes] <= 0
#        beta=1-(1-self.test_derivative(y,t))*y_new
        return beta

    def compute_source_terms(self, y, t):
        """
            compute the recharge for a time step on each cell based on recharge
            defined by user
        """
        #Matrix describing source term (recharge)
        N_edges = self.SD.N_edges
        N_nodes = self.SD.N_nodes
        d = np.zeros((2*N_nodes + N_edges, 1))
        alpha = self.compute_alpha(y, t)
        beta = self.compute_beta(y,t)
        Recharge_rate = self.So.compute_recharge_rate(t)
        Recharge_rate_spatialized = Recharge_rate * self.SD.w_node

        d[0:N_nodes] = beta * alpha * Recharge_rate_spatialized
        d[(N_nodes + N_edges):2*N_nodes + N_edges] = (1 - alpha) * Recharge_rate_spatialized


        edges_bool = self.BC.fixed_edge_matrix_boolean()
        edges = self.BC.fixed_edge_matrix_values()

        d[0, :] = (1-edges_bool[0])*d[0, :]
        d[N_nodes-1, :] = (1-edges_bool[1])*d[N_nodes-1, :]
        d[N_nodes, :] = (1-edges_bool[2])*d[N_nodes, :]
        d[N_nodes + N_edges-1, :] = (1-edges_bool[3])*d[N_nodes + N_edges - 1, :]
#       Set D(edges) at the fixes value (for Q)
        d[N_nodes, :] = -edges_bool[2] * edges[2]
        d[N_nodes + N_edges-1, :] = -edges_bool[3] * edges[3]
        d[d == -0] = 0

        return d


    def rhs(t, y):
        """
            Compute differential equation dy/dt = C*dy + d
        """
        ##########################################
        #Differential equation : dy/dt = c*dy + d#
        ##########################################
        y = np.reshape(y, (len(y), 1))
        c = bouss_obj.compute_c(y, t)
        d = bouss_obj.compute_source_terms(y, t)
        dy = np.dot(c, y)+d
        dy = np.squeeze(np.reshape(dy, (np.size(dy), 1)))
        return dy

    def res(t, y, yd):
        """
            Compute algebraic differential equation m * dy/dt = C*dy + d
        """
        #####################################################
        #function describing the DAE as MM * dy/dt = c*y + d#
        #####################################################
        yd = np.reshape(yd, (len(yd), 1))
        y = np.reshape(y, (len(y), 1))
        dy_new = np.dot(bouss_obj.compute_c(y, t), y) + bouss_obj.compute_source_terms(y, t)
        res = np.dot(bouss_obj.m, yd) - np.reshape(dy_new, (len(dy_new), 1))
        res = np.squeeze(np.reshape(res, (len(res), 1)))
        return res

    def implicit_scheme_solver(self):
        """
            Resolution of the DAE using implicit problem solver DAE
        """
        #Building initial state
        y0 = np.vstack((self.IC.sin, self.IC.qin, self.IC.q_sin))
        y0 = np.squeeze(np.reshape(y0, (len(y0), 1)))

        f = BoussinesqSimulation.rhs
        # time range to solve
        t = self.So.TP.t

        global bouss_obj
        bouss_obj = self
        # Computation of the initial state
        yd0 = f(self.So.TP.t[0], y0)

        #Build the DAE problem to solve
        self.model = Implicit_Problem(BoussinesqSimulation.res, y0, yd0, t[0])
        self.model.name = 'Boussinesq Simulation'
        self.sim = IDA(self.model)
        self.sim.report_continuously = True
        self.sim.verbosity = 10
        self.sim.maxord = 5
        self.sim.maxsteps = 5000

        #Solving the DAE using IDA from Assimulo
        ncp = len(t)
        ncp_list = t
        self.success = 1
        try:
            t_res, y_res, yd_res = self.sim.simulate(t[-1], ncp, ncp_list)
        except:
            print("The solver didn't end normally")
            self.success = 0

        print("Success : ", self.success)
        bouss_obj = 0
        if self.success == 1:
            self.SR = SR.SimulationResults(t_res, y_res, self.SD.N_nodes, \
                                       self.SD.N_edges, self.SD.x_node, self.SD.x_edges,\
                                       self.So.recharge_chronicle, self.SD.dx_edges)

    def compute_mass_matrix(self):
        """
            Compute the matrix m of the DAE
        """
        #Mass matrix used as a multiplier for dy/dt in DAE
        N_edges = self.SD.N_edges
        N_nodes = self.SD.N_nodes
        n_unknowns = 2*N_nodes + N_edges
        m = np.zeros((n_unknowns, n_unknowns))
        m[0:N_nodes, 0:N_nodes] = np.eye(N_nodes)
        return m

    def output_simu(self, folder):
        """
            write Simulations Results (Q,S,QS and x_Q,x_S,t_res) in  .txt files
            (delimiter : tab) in the current working directory
        """
        x_S = self.SR.x_node
        t = self.SR.t
        x_Q = self.SD.x_edges

        S_sol = self.SR.S
        Q_sol = self.SR.Q
        QS_sol = self.SR.QS
        #Computation of the outgoing flow of the hillslope (m3/s)
        self.Q_hs = -self.SR.Q[1:,1] + np.sum(np.abs(self.SR.QS[1:,:]),1)* + self.So.recharge_chronicle*self.SD.w_node[0][0]*self.SD.dx_edges[0][0]

        #Save Stock integration results
        name_file = folder + "/S"
        with open(name_file, "wb") as f:
            np.savetxt(f, S_sol, fmt='%1.12e', delimiter="\t", newline='\n')
        #Save flux integration results
        name_file = folder + "/Q"
        with open(name_file, "wb") as f:
            np.savetxt(f, Q_sol, fmt='%1.12e', delimiter="\t", newline='\n')
        #Save Seepage integration results
        name_file = folder + "/QS"
        with open(name_file, "wb") as f:
            np.savetxt(f, QS_sol, fmt='%1.12e', delimiter="\t", newline='\n')
        #Save time requested points during integration
        name_file = folder +"/t_res"
        with open(name_file, "wb") as f:
            np.savetxt(f, t, fmt='%1.12e', delimiter="\t", newline='\n')
        #Save mesh's centers coordinates (Morpho.nx)
        name_file = folder +"/x_S"
        with open(name_file, "wb") as f:
            np.savetxt(f, x_S, fmt='%1.12e', delimiter="\t", newline='\n')
        #Save mesh's edges coordinates (Morpho.nx+1)
        name_file = folder +"/x_Q"
        with open(name_file, "wb") as f:
            np.savetxt(f, x_Q, fmt='%1.12e', delimiter="\t", newline='\n')
        #Save recharge chronicle
        name_file = folder +"/recharge"
        with open(name_file, "wb") as f:
            np.savetxt(f, self.So.recharge_chronicle, fmt='%1.12e', delimiter="\t", newline='\n')
        #Save recharge's time valeus
        name_file = folder +"/time"
        with open(name_file, "wb") as f:
            np.savetxt(f, self.So.TP.t, fmt='%1.12e', delimiter="\t", newline='\n')
        #Save maximal storga of each cell
        name_file = folder +"/Smax"
        with open(name_file, "wb") as f:
            np.savetxt(f, self.IC.Smax, fmt='%1.12e', delimiter="\t", newline='\n')
        #Save width of hillslope's cells
        name_file = folder +"/width_S"
        with open(name_file, "wb") as f:
            np.savetxt(f, self.SD.w_node, fmt='%1.12e', delimiter="\t", newline='\n')
        #Save wodth of hillslope's cells edges
        name_file = folder +"/width_Q"
        with open(name_file, "wb") as f:
            np.savetxt(f, self.SD.Hs.w_edges, fmt='%1.12e', delimiter="\t", newline='\n')
        #Save flow of the hillslope
        name_file = folder +"/Q_hillslope"
        with open(name_file, "wb") as f:
            np.savetxt(f, self.Q_hs, fmt='%1.12e', delimiter="\t", newline='\n')


