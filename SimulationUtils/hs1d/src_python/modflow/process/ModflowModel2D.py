# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 13:21:25 2016

@author: Quentin Courtois
"""
import numpy as np
import numpy.matlib as matlib
import flopy
import flopy.utils.binaryfile as bf
import os
import matplotlib.pyplot as plt
from scipy.integrate import trapz



class ModflowModel2D(object):
    """
        Class defining a 2D, one-layered model of a watershed defined in the
        input. High topographic is a no-flow boundary (crest line). Low
        topographic is a constant head boundary (river)
        #######################################################################
        INPUT :
            Warning ! All inputs must have consistent units
        - hk : hydraulic conductivity in m/d
        - sy : specific yield
        - soil_depth : thickness of the layer
        - LX : length of the model (if no known topography and coordinates)
        - LY : width of the model (same conditions as LX)
        - nper : number of stress period (if no known recharge)
        - model_name : name of the watershed modeled
        - percentage_loaded : % of cell's volume initially filled with water
        - xcustom : a list of W-E coordinates corresponding to the cells
            constituting the watershed
        - ycustom : same as xcustom for N-S coordinates
        - recharge_custom : list of recharge values for each stress period
        - top_custom : array of top elevation values for each cell
        - bot_custom : array of bottom elevation values for each cell
        - cell_size : length of a cell (assumed as being a square)
        - outlet : coordinates of the outlet of the watershed
        #######################################################################

    """

    def __init__(self, hydro, watershed, output=0, init=0):

        """Set the waterhsed's geometry and parameters"""
        self.nlay = 1
        if not isinstance(watershed.xcustom, int):
            self.build_grid(watershed, hydro)

        #Size of rows and columns
        self.delc = watershed.cell_size
        self.delr = watershed.cell_size

        #Hydrodynamic paramters
        self.hk = hydro.hk
        self.sy = hydro.sy
        self.ss = hydro.ss

        #Build the shape of the watershed by setting non belonging cells to
        #inactive cells
        self.ibound = np.ones((1, self.nrow, self.ncol), dtype=np.int32)
        print('Number of rows in the model : ', self.nrow, '\n')
        print('Number of columns in the model : ', self.ncol, '\n')

        coord_null = np.where(self.ztop > round(np.max(watershed.top_custom)) + 5)
        coord_null_x = coord_null[1]
        coord_null_y = coord_null[0]
        i = 0
        while i < len(coord_null_x):
            self.ibound[0, coord_null_y[i], coord_null_x[i]] = 0
            i += 1

#        if isinstance(watershed.outlet, int):
#            self.ibound[0, 0, 0] = -1
#        else:
#            self.ibound[0, self.outlet_y, self.outlet_x] = -1

        #Set drain on each cells
        self.lrcec = self.set_drain(hydro)
        #######################################################################
        """Build inital State"""
        if init != 0:
            #Set a constant recharge equal to 0.001 m/d
            self.rech_init = self.set_recharge(hydro.recharge_initial)

            #Set all states to transient
            nstp = np.ones(self.nper)
            steady = [True, False]
            for i in range(self.nper-2):
                steady.append(False)

            #Set initial state
            self.strt = np.reshape(self.zbot + hydro.percentage_loaded * \
                                   hydro.soil_depth, (1, np.size(self.zbot, 0), np.size(self.zbot, 1)))
            if isinstance(watershed.outlet, int):
                self.strt[0, 0, 0] = self.ztop[1, 0, 0]
            else:
                self.strt[0, self.outlet_y, self.outlet_x] = self.ztop[self.outlet_y, self.outlet_x]

            #Set boundary conditions
            outlet = [self.outlet_y, self.outlet_x]
#            self.bound_sp = self.set_boundaries(outlet)

            #Create the flopy's Model
            print('Building the Initial State model \n')
            mf = self.create_model(watershed, nstp, steady, self.strt, self.rech_init)

            #Run the modflow model
            print('Running Modflow \n')
            try:
                self.run_model(mf)
            except:
                print('Modflow did not terminate normally !')
            print('Run ended\n')

            #Getting the final state of the computation
            piezo_init = self.init_output(watershed.model_name)

       ########################################################################
        """Build the true model to be computed"""
       #Set a constant recharge equal to 0.001 m/d
        self.rech = self.set_recharge(hydro.recharge_true)

        #Set all states to transient
        nstp = np.ones(self.nper)
        steady = [False, False]
        for i in range(self.nper-2):
            steady.append(False)

        #Set initial state
        if init == 0:
            file = os.getcwd() + "/" + watershed.model_name + "/" + watershed.model_name + "_initial_state"
            with open(file, 'r'):
               piezo_init = np.loadtxt(file)
        self.strt = piezo_init
        if isinstance(watershed.outlet, int):
            self.strt[0, 0] = self.ztop[1, 0, 0]
        else:
            self.strt[self.outlet_y, self.outlet_x] = self.ztop[self.outlet_y, self.outlet_x]


        #Create the flopy's Model
        print('Building the model \n')
        mf = self.create_model(watershed, nstp, steady, self.strt, self.rech)

        #Run the modflow model
        print('Running Modflow \n')
        try:
            self.run_model(mf)
        except:
            print('Modflow did not terminate normally !')
        print('Run ended\n')

        #Getting the final state of the computation
        if output !=0:
            piezo = self.model_output(watershed.model_name, hydro.time_select)

    def build_grid(self, watershed, hydro):
        """
            Computes a grid containing elevation values of the watershed.
            Non watershed cells have a value of max_elevation + 10.
            Computes also the elevation of the bottom of the layer
        """
        #Defines Limits of the grid
        self.ncol = int(np.round((np.max(watershed.xcustom) - np.min(watershed.xcustom))/watershed.cell_size)+1)
        self.nrow = int(np.round((np.max(watershed.ycustom) - np.min(watershed.ycustom))/watershed.cell_size)+1)

        new_x = np.round((watershed.xcustom - np.min(watershed.xcustom))/watershed.cell_size)
        new_y = np.round((watershed.ycustom - np.min(watershed.ycustom))/watershed.cell_size)

        i = 0
        while i < len(new_x):
            new_x[i] = int(new_x[i])
            new_y[i] = int(new_y[i])
            i += 1
        #New matrix containing the elevation for each point in the watershed
        grid = np.zeros((self.nrow, self.ncol)) + round(np.max(watershed.top_custom)) + 10

        i = 0
        while i < len(new_x):
            grid[new_y[i], new_x[i]] = watershed.top_custom[i]
            i += 1
        self.ztop = grid
        self.cell_size = watershed.cell_size

        #Setting the real outlet
        test_x = np.abs(watershed.outlet[0] - watershed.xcustom)
        temp_x = new_x[np.where(test_x == np.min(test_x))]
        self.outlet_x = int(temp_x[0])
        test_y = np.abs(watershed.outlet[1] - watershed.ycustom)
        temp_y = new_y[np.where(test_y == np.min(test_y))]
        self.outlet_y = np.max(new_y) - int(temp_y[0])

        #Check if outlet coordinates belong  to the watershed
        test_min = np.where(self.ztop == np.min(self.ztop))
        self.test_min = test_min
        d_out = np.sqrt((self.outlet_x - test_min[1])**2 + (self.outlet_y - test_min[0])**2)
        coord_new_outlet = np.where(d_out == np.min(d_out))
        self.outlet_x = test_min[1][coord_new_outlet][0]
        self.outlet_y = test_min[0][coord_new_outlet][0]

        if isinstance(watershed.bot_custom, int):
            self.zbot = self.ztop - hydro.soil_depth
        else:
            i = 0
            while i < len(new_x):
                grid[new_y[i], new_x[i]] = watershed.top_custom[i]
                i += 1
            self.zbot = grid

        self.x = np.arange(len(np.unique(new_x)))
        self.y = np.arange(len(np.unique(new_y)))
        self.new_y = new_y

        return self

    def set_recharge(self, recharge_custom):
        """
            Computes recharge on each cell, mainly based on SURFEX values
        """
        if isinstance(recharge_custom, int):
            if recharge_custom == -1:
                self.nper = 35
                self.perlen = np.ones(self.nper)
                recharge_chronicle = np.ones((self.nper, 1))*0.1
                rech = {}
                for i in range(self.nper):
                    rech[i] = recharge_chronicle[i][0]
        elif isinstance(recharge_custom, float):
            self.nper = 5
            self.perlen = np.ones(self.nper)
            recharge_chronicle = np.ones((self.nper, 1))*recharge_custom
            rech = {}
            for i in range(self.nper):
                rech[i] = recharge_chronicle[i][0]
        else:
            self.nper = np.size(recharge_custom)
            self.perlen = np.ones(self.nper)
            self.unit = 'days'
            rech = {}
            for i in range(self.nper):
                    rech[i] = recharge_custom[i]
        return rech

    def set_boundaries(self, outlet):
        """
            Set boundaries for the model : constant level at the sea level and no flow
            up stream
        """
        bound_sp = {}
        if isinstance(outlet, int):
            for i in range(self.nper):
                conductance = self.hk * (self.ztop[0] - self.zbot[0]) * self.delr/self.delc
                bound = [0, 0, 0, self.strt[0, 0, 0], conductance]
                bound_sp[i] = bound
        else:
            for i in range(self.nper):
                conductance = self.hk * (self.ztop[self.outlet_y, self.outlet_x] - \
                                         self.zbot[self.outlet_y, self.outlet_x]) * self.delr/self.delc * 1000000
                bound = [0, self.outlet_y, self.outlet_x, self.ztop[self.outlet_y, self.outlet_x], conductance]
                bound_sp[i] = bound
        return bound_sp

    def set_drain(self,hydro):
        drain = []
        top_max = np.max(self.ztop)
        for ir in range(self.nrow):
            for ic in range(self.ncol):
#                if not(ic == self.outlet_x and ir == self.outlet_y) and self.ztop[ir,ic] < top_max:
                if self.ztop[ir,ic] < top_max:
                    drain.append([0, ir, ic, self.ztop[ir,ic], self.hk*self.delc*hydro.soil_depth/self.delr])
        self.drn = drain
        return drain


    def create_model(self, watershed, nstp, steady, strt, rech):
        """
            Creates the model based on the previous calculations and flopy
            properties
        """
        mf = flopy.modflow.Modflow(watershed.model_name, exe_name='MODFLOW-NWT', version='mfnwt', verbose=True)
        flopy.modflow.ModflowDis(mf, self.nlay, self.nrow, self.ncol, delr=self.delr, \
                                 delc=self.delc, top=self.ztop, botm=self.zbot, \
                                 nper=self.nper, perlen=self.perlen, nstp=nstp, \
                                 steady=steady, itmuni=3, lenuni=2)
        flopy.modflow.ModflowBas(mf, ibound=self.ibound, strt=strt, stoper=10)
        flopy.modflow.mfupw.ModflowUpw(mf, laytyp=1, hk=self.hk, sy=self.sy, ss=self.ss)
        flopy.modflow.ModflowNwt(mf, Continue=True, maxiterout=1000, linmeth=2, mxiterxmd=500)
        flopy.modflow.ModflowRch(mf, nrchop=3, rech=rech)
        flopy.modflow.ModflowDrn(mf, stress_period_data=self.lrcec)
#        flopy.modflow.ModflowGhb(mf, stress_period_data=self.bound_sp)

        spd = {}
        for i in range(self.nper):
            spd[(0, i)] = ['save head', 'save budget']
        flopy.modflow.ModflowOc(mf, stress_period_data=spd)
        return mf

    def run_model(self, mf):
        """
            Run the model using modflow 2005
        """
        mf.write_input()
        success, mfoutput = mf.run_model(silent=False, pause=False, report=False)
        if not success:
            raise Exception('MODFLOW did not terminate normally.')

    def model_output(self, model_name, time_select):
        """
            Save model ouptuts as a temporary variable
        """
        # Create the headfile object
        print('Writing Output\n')
        headobj = bf.HeadFile(model_name+'.hds')

        #Head stocked in temporary(n_time,n_lay,n_row,n_col)
        temporary = headobj.get_alldata()

        temp = np.zeros((self.nrow, self.ncol))
        temp = np.squeeze(temporary[time_select,:,:])
        piezo = temp

        mf_list = flopy.utils.MfListBudget(model_name + ".list", timeunit='hours')
        drain = []
        perc_disc = []
        for i in range(1,self.nper+1):
            data = mf_list.get_data(totim=i)
            drain.append(data[7])
            perc_disc.append(data[11])

        drain_tot = []
        for i in range(self.nper):
            drain_tot.append(drain[i][1])

        drain_tot = [0] + drain_tot
        drain_val = -np.diff(drain_tot)

        name_file = os.getcwd() + "/" + model_name + "/" + model_name + "_piezometry"
        with open(name_file,"wb") as f:
            np.savetxt(f,piezo, fmt='%1.8e', delimiter="\t", newline='\n')
        name_file = os.getcwd() + "/" + model_name + "/" + model_name + "_drain"
        with open(name_file,"wb") as f:
            np.savetxt(f,drain_val, fmt='%1.8e', delimiter="\t", newline='\n')
        return piezo

    def init_output(self, model_name):
        """
            Save Initial final state to be used as the initial state of the model
        """
                # Create the headfile object
        print('Writing Initial State Output\n')
        headobj = bf.HeadFile(model_name+'.hds')

        #Head stocked in temporary(n_time,n_lay,n_row,n_col)
        temporary = headobj.get_alldata()

        temp = np.zeros((self.nrow, self.ncol))
        self.temp = temporary
        temp = np.squeeze(temporary[-1,:,:])
#        temp[0, self.outlet_y, self.outlet_x,] = self.strt[0, self.outlet_y, self.outlet_x]
#        piezo = np.squeeze(temp)
        mf_list = flopy.utils.MfListBudget(model_name + ".list", timeunit='hours')
        self.list = mf_list
        drain = []
        perc_disc = []
        for i in range(1,self.nper+1):
            data = mf_list.get_data(totim=i)
            drain.append(data[7])
            perc_disc.append(data[11])

        drain_tot = []
        for i in range(self.nper):
            drain_tot.append(drain[i][1])

        drain_tot = [0] + drain_tot
        drain_val = -np.diff(drain_tot)

        name_file = os.getcwd() + "/" + model_name + "/" + model_name + "_drain_init"
        with open(name_file,"wb") as f:
            np.savetxt(f,drain_val, fmt='%1.8e', delimiter="\t", newline='\n')
        piezo = temp
        self.piezo_init_temp = piezo
        name_file = os.getcwd() + "/" + model_name + "/" + model_name + "_initial_state"
        with open(name_file,"wb") as f:
            np.savetxt(f,piezo, fmt='%1.8e', delimiter="\t", newline='\n')

        return piezo

    def plot_output(self, model_name, watershed, hydro):
        print("Ploting output")
        mf_list = flopy.utils.MfListBudget(model_name + ".list", timeunit='hours')
        self.list = mf_list
        drain = []
        perc_disc = []
        for i in range(1,self.nper+1):
            data = mf_list.get_data(totim=i)
            drain.append(data[7])
            perc_disc.append(data[11])

        drain_tot = []
        for i in range(self.nper):
            drain_tot.append(drain[i][1])

        drain_tot = [0] + drain_tot
        drain_val = -np.diff(drain_tot)

        plt.figure(2)
        plt.plot(drain_val)
        plt.ylabel('drain (m3/d)')
        plt.xlabel('Time (h)')
        plt.title('Flowrate a the outlet as a function of time')
        plt.grid(True)
        plt.show()

#        plt.figure(3)
#        plt.plot(head_out_val)
#        plt.ylabel('Flow through constant head (m3/d)')
#        plt.xlabel('Time (h)')
#        plt.title('Flowrate a the outlet as a function of time')
#        plt.grid(True)
#        plt.show()

        plt.figure(4)
        plt.plot((drain_val)/3600)
        plt.ylabel('Total out flow (m3/d)')
        plt.xlabel('Time (h)')
        plt.title('Flowrate a the outlet as a function of time')
        plt.grid(True)
        plt.show()

        plt.figure(5)
        plt.plot(perc_disc)
        plt.ylabel('Percent discrepancy')
        plt.xlabel('Time (h)')
        plt.title('Percent discrepancy as a function of time')
        plt.grid(True)
        plt.show()

        self.t_res = np.arange(len(drain_val))
        self.Q_r =  drain_val
        Q_reel = trapz(self.Q_r, self.t_res)
        print('\n Total volume through the outlet : ', Q_reel, ' m3 over ', self.t_res[-1]+1, ' days')
        print('Total Surface of the Watershed : ', len(watershed.xcustom)*(watershed.cell_size**2), ' m²')
        t = np.arange(len(hydro.recharge_true))
        rec = np.trapz(np.reshape(hydro.recharge_true,(1,len(hydro.recharge_true))), np.reshape(t,(1,len(t))))
        print('Total volume in : ', rec*len(watershed.xcustom)*(watershed.cell_size**2), ' m3 over 35 d')

