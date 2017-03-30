# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 10:06:17 2017

@author: Quentin Courtois
"""

#Adding all modules path to the main path
import sys
import os
sys.path.append(os.getcwd() + '/Hillslope1D/process/model_input/')
sys.path.append(os.getcwd() + '/Hillslope1D/process/simulation/')
sys.path.append(os.getcwd() + '/Hillslope1D/process/tools')
sys.path.append(os.getcwd() + '/Hillslope1D/process/tools/boussinesq_tools')
sys.path.append(os.getcwd() + '/Hillslope1D/process/tools/file_management')
sys.path.append(os.getcwd() + '/Hillslope1D/process/tools/time')
sys.path.append(os.getcwd() + '/utils/')
sys.path.append(os.getcwd() + '/Kerbernez/')

import numpy as np
import LoadCasTest as LCT
import GeologicInputs as GI
import HydrologicInputs as HI
import MorphologicInputs as MI
import BoussinesqSimulation as BS
import matplotlib.pyplot as plt
import LoadTest as LT

"""Building of a list of hillslopes"""
###############################################################################
#Flags and input parameters definition

#Definition of the flag used to specify the test case
flag = 4
cst = False
out_put = 0
interp = 0
initial_state = 1
#Definition of folder and name for custom cases (if needed, set flag to 4)
src_path = os.getcwd()
if "\\" in src_path:
    src_path = src_path.replace("\\", "/")

watershed_path = src_path[:-10] + 'test_case/matlab/Kerbernez/'
list_dir = os.listdir(watershed_path)

cst_name = []
cst_folder = []
var_name =  []
var_folder = []
for i in range(len(list_dir)):
    temp = list_dir[i]
    if temp[-3:] == 'cst':
        cst_name.append(temp)
        cst_folder.append("Kerbernez/"+temp)
    elif temp[-3:] == 'var':
        var_name.append(temp)
        var_folder.append("Kerbernez/"+temp)
    else:
        print("Trouble reading watershed's folders")

if cst == True:
    custom_folder = cst_folder
    name_custom = cst_name
    print(len(cst_name),' hillslopes found and loaded')
else:
    custom_folder = var_folder
    name_custom = var_name
    print(len(var_name),' hillslopes found and loaded')
#Output option for BoussinesqSimulation results
output = 1

#Plot options for vizualisation of results
plot_option = 1

###############################################################################
"""Compute initial state for all hillslopes"""
#Loading and running model
if initial_state != 0:
    success = np.zeros((len(name_custom),1))
    for i in range(len(name_custom)):
        #Reading the test case and importing data
        TestCase = LCT.LoadCasTest(flag, custom_folder[i], name_custom[i], out_put, interp)
        print("Test Case ", i, " files loaded")

        perc = 0.0
        z=-1
        rec = []
        for j in range(100000):
            rec.append(np.mean(TestCase.input.recharge_chronicle[119:838]))
        rec = np.reshape(rec,(1, len(rec)))
        TestCase.input.recharge_chronicle = rec[0]
        TestCase.input.time = np.arange(100000) * 3600

        #Building Model inputs
        Morpho = MI.MorphologicInputs(nx=len(TestCase.input.x), discretization_type='custom',\
                                      x_custom = TestCase.input.x, angle=TestCase.input.i, \
                                      z_custom = z, w = TestCase.input.w)

        Geol = GI.GeologicInputs(k=TestCase.input.k, f=TestCase.input.f, soil_depth=TestCase.input.soil_depth)

        Hydro = HI.HydrologicInputs(recharge_rate=TestCase.input.recharge_chronicle, time_custom=TestCase.input.time, \
                                    recharge_type='databased',perc_loaded=perc, unit = 'hour')

        #Creating the BoussinesqSimulation model using Model Inputs
        Simulation = BS.BoussinesqSimulation(Morpho, Geol, Hydro, Id=TestCase.name)
        print("Model built")

        #Running BoussinesqSimulation model
        Simulation.implicit_scheme_solver()
    #    success = success.append(Simulation.success)
        success[i] = Simulation.success
        print("Simulation ", i, " over")

    ###############################################################################

    ###############################################################################
        #Output of the simulation

        #Building output files
        if success[i] == 1:
            if not os.path.exists(os.getcwd() + '/simulation_results'):
                os.makedirs(os.getcwd() + '/simulation_results')
                print('Folder for all simulation results saving created')

            if not os.path.exists(os.getcwd() + '/simulation_results/' + TestCase.name):
                os.makedirs(os.getcwd() + '/simulation_results/' + TestCase.name)
                print('Output folder created for the simulated Hillslope')

            out_folder = os.getcwd() + '/simulation_results/' + TestCase.name
            file = out_folder + "/" + TestCase.name + "_init"
            with open(file,"wb") as f:
                np.savetxt(f, Simulation.SR.S[-1,:], fmt='%1.8e', delimiter="\t", newline='\n')
            file = out_folder + "/" + TestCase.name + "_Smax"
            with open(file,"wb") as f:
                np.savetxt(f, Simulation.IC.Smax, fmt='%1.8e', delimiter="\t", newline='\n')
            print('Simulation ', i, ' initial state saved\n')
        else:
            print('Simulation ', i, ' failed, no state saved\n')
        Simulation = 0

print('\n Runing the model\n')
"Loading all hillslopes"
###############################################################################
#Loading and running model
success = np.zeros((len(name_custom),1))
for i in range(len(name_custom)):
    #Reading the test case and importing data
    TestCase = LCT.LoadCasTest(flag, custom_folder[i], name_custom[i], out_put, interp)
    print("Test Case ", i, " files loaded")

    z = -1
    if out_put != 0:
        if np.size(TestCase.input.x)-1 == np.size(TestCase.output.storage[:,0]):
            print("using saved initial state")
            perc = TestCase.output.perc_storage[:,0]/100
            z = np.reshape(TestCase.input.z_mod,(len(TestCase.input.z_mod),1))
        else:
            print("using default inital state")
            perc = 0.5
    else:
        folder = os.getcwd() + '/simulation_results/' + name_custom[i] + "/" + name_custom[i]
        if os.path.exists(folder + "_init"):
            temp_init = np.loadtxt(folder + "_init")
            temp_Smax = np.loadtxt(folder + "_Smax")
            perc = np.reshape(temp_init, (len(temp_init),1))/np.reshape(temp_Smax, (len(temp_Smax),1))
        else:
            perc = 0.5
    #Building Model inputs
    rec_temp = np.matlib.repmat(TestCase.input.recharge_chronicle[119:838],1,10)
    rec_temp = rec_temp[0]
    t_temp = np.squeeze(np.zeros((1, len(rec_temp))))
    for j in range(len(rec_temp)):
        t_temp[j] = j*3600
    TestCase.input.recharge_chronicle = rec_temp
    TestCase.input.time = t_temp
    Morpho = MI.MorphologicInputs(nx=len(TestCase.input.x), discretization_type='custom',\
                                  x_custom = TestCase.input.x, angle=TestCase.input.i, \
                                  z_custom = z, w = TestCase.input.w)

    Geol = GI.GeologicInputs(k=TestCase.input.k, f=TestCase.input.f, soil_depth=TestCase.input.soil_depth)

    Hydro = HI.HydrologicInputs(recharge_rate=TestCase.input.recharge_chronicle, time_custom=TestCase.input.time, \
                                recharge_type='databased',perc_loaded=perc)

    #Creating the BoussinesqSimulation model using Model Inputs
    Simulation = BS.BoussinesqSimulation(Morpho, Geol, Hydro, Id=TestCase.name)
    print("Model built")

    #Running BoussinesqSimulation model
    Simulation.implicit_scheme_solver()
#    success = success.append(Simulation.success)
    success[i] = Simulation.success
    if i == 0:
        Simu_temp = Simulation
    print("Simulation ", i, " over\n")

###############################################################################

###############################################################################
    #Output of the simulation

    #Building output files
    if output != 0:
        if success[i] == 1:
            if not os.path.exists(os.getcwd() + '/simulation_results'):
                os.makedirs(os.getcwd() + '/simulation_results')
                print('Folder for all simulation results saving created')

            if not os.path.exists(os.getcwd() + '/simulation_results/' + TestCase.name):
                os.makedirs(os.getcwd() + '/simulation_results/' + TestCase.name)
                print('Output folder created for the simulated Hillslope')

            out_folder = os.getcwd() + '/simulation_results/' + TestCase.name
            Simulation.output_simu(out_folder)
            print('Simulation ', i, ' results saved')
        else:
            print('Simulation ', i, ' failed, no results saved')
    Simulation = 0

###############################################################################
#Loading all results
Simu = []
for i in range(len(name_custom)):
    if success[i] == 1:
        folder = os.getcwd() + '/simulation_results/' + name_custom[i]
        Simul_load = LT.LoadTest(folder)
        Simu.append(Simul_load)
    else:
        Simu.append(0)


#Ploting results
if plot_option != 0:
    Q_riv = 0
    surf = []
    recharge = []
    for i in range(len(name_custom)):
        if success[i] == 1:
            temp = Simu[i]
            Q_riv = Q_riv + np.abs(temp.Q[:,1]) - np.abs(temp.QS[:,0])*5 + np.sum(temp.QS*5,1)
            surf_temp = np.trapz(np.reshape(temp.w_edges,(1,len(temp.w_edges))), np.reshape(temp.x_Q,(1,len(temp.x_Q))))
            surf.append(surf_temp)
            rec = TestCase.input.recharge_chronicle
    recharge.append(np.trapz(np.reshape(rec,(1,len(rec))), np.reshape(temp.t_res[1:],(1,len(temp.t_res)-1))))

    plt.figure(6)
    plt.plot(Simu[0].t_res, Q_riv)
    plt.ylabel('Q river (m3/s)')
    plt.xlabel('Time (s)')
    plt.title('Flowrate a the outlet as a function of time')
    plt.grid(True)
    plt.show()

    Q_r = np.trapz(np.reshape(Q_riv,(1,len(Q_riv))), np.reshape(Simu[0].t_res, (1, len(Simu[0].t_res))))
    print('Total volume through the outlet : ', float(Q_r), ' m3 over ', int((Simu[0].t_res[-1]+1)/86400), ' days')
    print('Total Surface of the watershed is : ', np.sum(surf), ' m²')
    print('Total Volume In is : ', np.sum(recharge) * np.sum(surf), ' m3 over 35 d')
###############################################################################

plt.figure(num=42)
for i in range(len(Simu)):
    plt.plot(Simu[i].x_Q,Simu[i].w_edges)
plt.show()


surf_temp = []
for i in range(len(Simu)):
    surf_temp.append(np.sum(Simu[i].w_edges)*10)
