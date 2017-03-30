# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 14:49:04 2016

@author: Quentin Courtois
"""
#Adding Modules path
import sys
import os
sys.path.append(os.getcwd() + '/Hillslope1D/process/')
sys.path.append(os.getcwd() + '/modflow/')
sys.path.append(os.getcwd() + '/modflow/process/')

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import BoussinesqSimulation as BS
import numpy as np
import os
import LoadTest as LT
import ModflowModel1D as MM1D

plot_options = 0
bouss = 1
output = 1

test = 1

modflow = 0
comp = 0

if bouss != 0:
    #Initializaton
    A = BS.BoussinesqSimulation(period=None, recharge_type='square', recharge_rate=30, \
                                tmin=0, tmax=35, Nt=35*24*10, unit='days', \
                                boundary_type=['S','Q'], boundary_value=[0,0], \
                                xmin=0, xmax=100, nx=100, discretization_type='linear', x_custom=-1, \
                                angle=0, w=0, soil_depth=0, k=1/3600, f=0.3, \
                                Id=1, percentage_loaded=0.5)
    #Computation
    A.implicit_scheme_solver()

    #Tests about limits and computed values
    print('\nNumber of values computed : ', np.size(A.SR.S))
    test_sup = A.SR.S > np.reshape(A.IC.Smax,(1,len(A.IC.Smax)))
    print('\nNumber of values > Smax : ',np.sum(test_sup))
    print('% of values > Smax : ', 100*np.sum(test_sup)/np.size(A.SR.S))
    print('Maximal Gap bewteen S and Smax : ', np.max(A.SR.S -np.reshape(A.IC.Smax,(1,len(A.IC.Smax)))))
    coord_sup = np.where(test_sup == 1)
    test_inf = A.SR.S < 0
    print('\nNumber of values < 0 : ', np.sum(test_inf))
    print('% of values < 0 : ', 100*np.sum(test_inf)/np.size(A.SR.S))
    print('Minimal value < 0 : ', np.min(A.SR.S))
    coord_inf = np.where(test_inf == 1)

    #Plot of the water head for each time and each x
    S_sol = A.SR.S
    Q_sol = A.SR.Q
    QS_sol = A.SR.QS


    if plot_options !=0 :
        fig1 = plt.figure(num =1)
        ax1 = fig1.gca(projection = '3d')
        X = A.SD.compute_x_node()
        Y = A.SR.t
        X,Y = np.meshgrid(X,Y)
        surf = ax1.plot_surface(X,Y,S_sol,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0,antialiased=False)
        fig1.colorbar(surf, shrink=0.5, aspect=5)

    #    fig2 = plt.figure(num=2)
    #    ax2 = fig2.gca(projection = '3d')
    #    X = A.SD.x_edges
    #    Y = A.t_res
    #    X,Y = np.meshgrid(X,Y)
    #    surf = ax2.plot_surface(X,Y,Q_sol,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0,antialiased=False)
    #    fig2.colorbar(surf,shrink=0.5, aspect=5)
    #
    #    fig3 = plt.figure(num=3)
    #    ax3 = fig3.gca(projection = '3d')
    #    X = A.SD.compute_x_node()
    #    Y = A.t_res
    #    X,Y = np.meshgrid(X,Y)
    #    surf = ax3.plot_surface(X,Y,QS_sol,cmap=cm.coolwarm,linewidth=0,antialiased=False)
    #    fig3.colorbar(surf,shrink=0.5, aspect=5)
    #
    if output != 0:
        out_folder = os.getcwd()
        A.output_simu(out_folder)

    if test != 0:
        Test_diff = LT.LoadTest(os.getcwd())
        diff_S = A.SR.S - Test_diff.S
        diff_Q = A.SR.Q - Test_diff.Q
        diff_QS = A.SR.QS - Test_diff.QS
        t_res = np.zeros((len(A.SR.t),1))
        i = 0
        while i < len(t_res):
            t_res[i] = A.SR.t[i]
            i += 1
        diff_t = t_res - Test_diff.t_res

        print("\nMaximal gap between calculated S and reference value : ", np.max(np.absolute(diff_S)), " (stock variation)")
        print("Mean gap between calculated S and reference value : ", np.mean(diff_S), " (stock variation)\n")
        print("Maximal gap between calculated Q and reference value : ", np.max(np.absolute(diff_Q)), " (m3/s)")
        print("Mean gap between calculated Q and reference value : ", np.mean(diff_Q), " (m3/s)\n")
        print("Maximal gap between calculated QS and reference value : ", np.max(np.absolute(diff_QS)), " (m3/s)")
        print("Mean gap between calculated QS and reference value : ", np.mean(diff_QS), " (m3/s)\n")
        print("Maximal gap between integration time and reference value : ", np.max(np.absolute(diff_t)), " (s)")
        print("Mean gap between integration time and reference value : ", np.mean(diff_t), " ()\n")

if modflow !=0:
    nr = 100
    TEST = MM1D.ModflowModel1D(nlay=1, ncol=1, nrow=nr, ntime=35,unit='day', hk=24,sy=0.3, \
                               soil_depth=1, alpha=0.05, Lx=500, Ly=100, nper=10*24, \
                               model_name='BoussinesqEvaluation',percentage_loaded=0.5, \
                               xcustom=-1,recharge_custom=-1,top_custom=-1,bot_custom=-1)
    if comp !=0:
        TEST.model_comparison(nrow=nr)
#    fig4 = plt.figure(num=4)
#    ax4 = fig1.gca()
#    plot = ax4.plot(diff_t)
#
        fig8 = plt.figure(num=8)
        ax8 = fig8.gca()
        ax8.plot(A.SR.S[:,nr-1])
        ax8.plot(TEST.stock[:,nr-1])
        fig9 = plt.figure(num=9)
        ax9 = fig9.gca()
        ax9.plot(A.SR.S[:,nr-1] - TEST.stock[:,nr-1])
        fig10 = plt.figure(num=10)
        ax10 = fig10.gca()
        ax10.plot(A.SR.S[3300,90:99])
        ax10.plot(TEST.stock[3300,90:99])
#    fig5 = plt.figure(num =5)
#    ax5 = fig1.gca(projection = '3d')
#    X = A.SD.compute_x_node()
#    Y = A.SR.t
#    X,Y = np.meshgrid(X,Y)
#    surf = ax5.plot_surface(X,Y,diff_S,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0,antialiased=False)
#    fig5.colorbar(surf, shrink=0.5, aspect=5)
