# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 14:27:38 2016

@author: Quentin Courtois
"""
import os
import numpy as np
import numpy.matlib as matlib
import flopy
#import LoadTest as LT
import flopy.utils.binaryfile as bf
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class ModflowModel1D(object):
    """
        Class defining a 1D Modflow Model using the library flopy. The model
        corresponds to a 1D model over an hillslope. The low topographic is a
        constant head boundary condition (corresponding to a river), and the
        high topographic is a no flow boundary condition (corresponding to a
        ridge).
        #######################################################################
        INPUT :
            WARNING ! All Inputs must have the same unit

        - nlay : number of layer (default 1)
        - ncol : number of columns (default 1)
        - nrow : number of lines (default 100)
        - ntime : maximal time value (default 35)
        - unit : time unit used (default 'day')
        - hk : hydraulic conductivity along the hillslope (default 24 m/d)
        - sy : specific yield (defautl 0.3)
        - soil_depth : difference between top and bottom layer (default 1)
        - alpha : angle of the hillslope (default 0.05 rad)
        - Lx : width of the model (default 500)
        - Ly : length of the model (default 100)
        - nper : number of time period on which the model will be calculated
            for each time unit (default is 10*24-1)
        - model_name : name of the model (default : 'BoussinesqEvaluation')
        - percentage_loaded : percentage of the cells' volume initially filled
            with water (default 0.5)
        - xcustom : list or numpy array containing size of the cells along the
            hillslope (default : -1, meaning that the spacing is constant and
            calculated by Ly/nrow)
        - recharge_custom : list or numpy array or value corresponding to the
            recharge in distance unit / time unit (default -1 : load a recharge
            chronicle)
        - top_custom : list or numpy array or value corresponding to the
            altitude of the roof of the first layer (default -1 : computed
            using angle alpha)
        - bot_custom : list or numpy array or value corresponding to the
            altitude of the bottom of the first layer (default -1 : computed
            using angle alpha)
        #######################################################################

    """


    def __init__(self, nlay=1, ncol=1, nrow=20, ntime=35, unit='day', hk=1, sy=0.3, \
                 soil_depth=2, alpha=0.05, Lx=500, Ly=100, nper=10*24, \
                 model_name='BoussinesqEvaluation', percentage_loaded=0.8, \
                 xcustom=-1, recharge_custom=-1, top_custom=-1, bot_custom=-1, strt=-1, std=False):
        self.name = model_name
        self.hk = hk
        self.sy = sy
        self.nper = ntime*nper
        perlen = np.ones(self.nper)
        nstp = np.ones(self.nper)

        steady = [std, False]
        for i in range(self.nper-2):
            steady.append(False)

        if isinstance(xcustom, int):
            self.delc = Ly/nrow
        self.delr = Lx/ncol

        Topo, Bottom = self.topo(soil_depth, alpha, nrow, top_custom, bot_custom, Ly)
        self.ztop = np.reshape(Topo, (len(Topo), 1))
        self.zbot = np.reshape(Bottom, (len(Bottom), 1))

        self.ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
#        self.ibound[0, 0, 0] = -1

        if (isinstance(strt,int)):
            self.strt = np.reshape(self.zbot, (nlay, nrow, ncol)) + (percentage_loaded*soil_depth)
            print('Using defaul init')
        else:
            self.strt = np.reshape(self.zbot, (nlay, nrow, ncol)) + (percentage_loaded*soil_depth)
            self.strt[0,:,:] = np.reshape(strt,(nrow,ncol))
            print('Using loaded init')

        self.rech = self.set_recharge(nrow, ncol, recharge_custom)

#        self.bound_sp = self.set_boundaries()

        self.lrcec = self.set_drain(nrow, ncol)

        mf = self.create_model(model_name, nlay, nrow, ncol, perlen, nstp, steady)

        success = self.run_model(mf)

        if success!= 0:
            self.stock, self.piezo = self.model_output(model_name, nrow)

    def topo(self, soil_depth, alpha, nrow, top_custom, bot_custom, Ly):
        """
            Compute elevation of the top and the bottom of the layer, based on
            either custom elevation or an angle.
        """
        if isinstance(bot_custom, int) and isinstance(top_custom, int):
            if bot_custom == -1 and top_custom == -1:
                x = np.arange(0, Ly, Ly/nrow)
                alpha = 0.05
                Topo = []
                Bottom = []
                for rownum in range(nrow):
                    Topo.append(np.sin(alpha)*x[rownum]+soil_depth)
                    Bottom.append(np.sin(alpha)*x[rownum])
            elif bot_custom < top_custom:
                Topo = np.ones((nrow, 1))*top_custom
                Bottom = np.ones((nrow, 1))*bot_custom
        else:
            Topo = np.reshape(top_custom, (len(top_custom), 1))
            Bottom = np.reshape(bot_custom, (len(bot_custom, 1)))

        return Topo, Bottom

    def set_recharge(self, nrow, ncol, recharge_custom):
        """
            Defines the recharge of the model, based on a custom recharge or a
            value or a default squared-shpae reference recharge
        """
        if isinstance(recharge_custom, int) or isinstance(recharge_custom, float):
            if recharge_custom == -1:
                A = LT.LoadTest()
                A.recharge_chronicle = A.recharge_chronicle*86400
                Recharge = np.zeros((self.nper, nrow, ncol))
                for rownum in range(self.nper):
                    Recharge[rownum, :, :] = matlib.repmat(A.recharge_chronicle[rownum], nrow, ncol)
                rech = {}
                for spi in range(self.nper):
                    rech[spi] = Recharge[spi]
            elif recharge_custom >= 0:
                recharge_chronicle = np.ones((self.nper, 1))*recharge_custom
                Recharge = np.zeros((self.nper, nrow, ncol))
                for rownum in range(self.nper):
                    Recharge[rownum, :, :] = matlib.repmat(recharge_chronicle[rownum], nrow, ncol)
                rech = {}
                for spi in range(self.nper):
                    rech[spi] = Recharge[spi]
            else:
                print('Non-recognized recharge value. Please enter a valid value')

        else:
            Recharge = np.zeros((self.nper, nrow, ncol))
            for rownum in range(self.nper):
                Recharge[rownum, :, :] = matlib.repmat(recharge_custom[rownum], nrow, ncol)
            rech = {}
            for spi in range(self.nper):
                rech[spi] = Recharge[spi]
        return rech

    def set_boundaries(self):
        """
            Set boundaries for the model : constant level in x = 0 and no flow
            up stream (x = xmax)
        """
        bound_sp = []
        conductance = self.hk * (self.ztop[0]) * self.delr/self.delc
        bound_sp.append([0, 0, 0, self.strt[0, 0, 0], conductance])
        return bound_sp

    def set_drain(self, nrow, ncol):
        """
            set drain on the topographic level to simulate overland flow
        """
        drain = []
        for ir in range(nrow):
            for ic in range(ncol):
                drain.append([0, ir, ic, self.ztop[ir][ic], self.hk*self.delc*self.delr])
        lrcec = {0:drain}
        return lrcec

    def create_model(self, model_name, nlay, nrow, ncol, perlen, nstp, steady):
        """
            Creates the model based on the previous calculations and flopy
            properties
        """
        mf = flopy.modflow.Modflow(model_name, exe_name='MODFLOW-NWT', version='mfnwt', verbose=True)
        flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr=self.delr, \
                                 delc=self.delc, top=self.ztop, botm=self.zbot, \
                                 nper=self.nper, perlen=perlen, nstp=nstp, \
                                 steady=steady, itmuni=3, lenuni=2)
        flopy.modflow.ModflowBas(mf, ibound=self.ibound, strt=self.strt, stoper=10)
        flopy.modflow.mfupw.ModflowUpw(mf, laytyp=1, hk=self.hk, sy=self.sy, ss=10**-5)
        flopy.modflow.ModflowNwt(mf, Continue=True, maxiterout=500, linmeth=2)
        flopy.modflow.ModflowRch(mf, nrchop=3, rech=self.rech)
        flopy.modflow.ModflowDrn(mf, stress_period_data=self.lrcec)
#        flopy.modflow.ModflowGhb(mf, stress_period_data=self.bound_sp)

        spd = {}
        for i in range(self.nper):
            spd[(0, i)] = ['save head', 'save budget']
        flopy.modflow.ModflowOc(mf, stress_period_data=spd)
        return mf
#        mf = flopy.modflow.Modflow(model_name, exe_name='mf2005')
#        flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr=self.delr, \
#                                 delc=self.delc, top=self.ztop, botm=self.zbot, \
#                                 nper=self.nper, perlen=perlen, nstp=nstp, \
#                                 steady=steady)
#        flopy.modflow.ModflowBas(mf, ibound=self.ibound, strt=self.strt)
#        flopy.modflow.ModflowLpf(mf, laytyp=1, hk=self.hk, sy=self.sy, laywet=1, \
#                                 wetdry=1, ss=self.sy)
#        flopy.modflow.ModflowPcg(mf)
#        flopy.modflow.ModflowRch(mf, nrchop=1, rech=self.rech)
#        flopy.modflow.ModflowDrn(mf, stress_period_data=self.lrcec)
#        flopy.modflow.ModflowGhb(mf, stress_period_data=self.bound_sp)
#
#        spd = {}
#        for i in range(self.nper):
#            spd[(0, i)] = ['save head', 'save budget', 'print head']
#        flopy.modflow.ModflowOc(mf, stress_period_data=spd)
#        return mf

    def run_model(self, mf):
        """
            Run the model using modflow 2005
        """
        mf.write_input()
        success, mfoutput = mf.run_model(silent=False, pause=False)
        if not success:
           print('MODFLOW did not terminate normally.')
        return success
    def model_output(self, model_name, nrow):
        """
            Save model ouptuts as a temporary variable
        """
        # Create the headfile object
        headobj = bf.HeadFile(model_name+'.hds')

        #Head stocked in temporary(n_time,n_lay,n_row,n_col)
        temporary = headobj.get_alldata()

        temp = np.zeros((self.nper+1, nrow))
        temp[1:, :] = np.squeeze(temporary)
        temp[0, :] = np.reshape(np.squeeze(self.strt[0]), (1, len(self.strt[0])))
        temporary_diff = (temp - np.reshape(self.zbot, (1, len(self.zbot))))*self.sy*self.delr

        self.stock = np.squeeze(temporary_diff)
        self.piezo = np.squeeze(temp)

        #Save Stock and Piezometry computed by Modflow
        name_file = os.getcwd() + "/modflow_Stock_variation"
        with open(name_file, "wb") as f:
            np.savetxt(f, self.stock, fmt='%1.8e', delimiter="\t", newline='\n')

        name_file = os.getcwd() + "/modflow_Piezometry"
        with open(name_file, "wb") as f:
            np.savetxt(f, self.piezo, fmt='%1.8e', delimiter="\t", newline='\n')

        mf_list = flopy.utils.MfListBudget(model_name + ".list", timeunit='hours')
        drain = []
        for i in range(1,self.nper+1):
            data = mf_list.get_data(totim=i)
            drain.append(data[7])

        drain_tot = []
        for i in range(self.nper):
            drain_tot.append(drain[i][1])

        drain_tot = [0] + drain_tot
        drain_val = -np.diff(drain_tot)

        name_file = os.getcwd() + "/" + model_name + "/" + model_name + "_drain"
        with open(name_file,"wb") as f:
            np.savetxt(f,drain_val, fmt='%1.8e', delimiter="\t", newline='\n')

        self.drain = drain_val

        return self.stock, self.piezo

    def model_comparison(self, nrow):
        """
            Comparison between modflow's results and Boussinesq's results
        """
        #Load the file containing Boussinesq's simulation results
        A = LT.LoadTest()
        test = A.S - self.stock
        print('Sum of differences between Modflow and Boussinesq ', np.sum(test))
        print('Mean difference between Modflow and Boussinesq ', np.mean(test))

        #Plot of the differences between Boussinesq and Modflow results
        fig6 = plt.figure(num=6)
        ax6 = fig6.gca(projection='3d')
        X = np.arange(1, nrow+1)
        Y = np.arange(0, self.nper+1)
        X, Y = np.meshgrid(X, Y)
        surf = ax6.plot_surface(X, Y, test, rstride=1, cstride=1, cmap=cm.coolwarm, \
                                linewidth=0, antialiased=False)
        fig6.colorbar(surf, shrink=0.5, aspect=5)

        #Plot of Stock computed by Modflow as a function of space and time
        fig7 = plt.figure(num=7)
        ax7 = fig7.gca(projection='3d')
        X = np.arange(1, nrow+1)
        Y = np.arange(0, self.nper+1)
        X, Y = np.meshgrid(X, Y)
        surf = ax7.plot_surface(X, Y, self.stock, rstride=1, cstride=1, cmap=cm.coolwarm, \
                                linewidth=0, antialiased=False)
        fig7.colorbar(surf, shrink=0.5, aspect=5)

