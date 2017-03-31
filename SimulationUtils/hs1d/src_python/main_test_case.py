# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 14:32:18 2017

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

import numpy as np
import LoadCasTest as LCT
import GeologicInputs as GI
import HydrologicInputs as HI
import MorphologicInputs as MI
import BoussinesqSimulation as BS
import matplotlib.pyplot as plt

###############################################################################
#Flags and input parameters definition
plt.close("all")

#Definition of the flag used to specify the test case
flag = 0

#Definition of folder and name for custom cases (if needed, set flag to 4)
custom_case = "custom_input"
name_custom = ""

#Output option for BoussinesqSimulation results
output = 1

#Plot options for vizualisation of results
plot_option = 0
out_put = 0

###############################################################################

###############################################################################
#Loading and running model

#Reading the test case and importing data
TestCase = LCT.LoadCasTest(flag, custom_case, name_custom, output=out_put, interp=0)
print("Test Case files loaded")

z = -1

if out_put != 0:
    if np.size(TestCase.input.x)-1 == np.size(TestCase.output.storage[:,0]):
        print("using saved initial state")
        perc = TestCase.output.perc_storage[:,0]/100
        z = np.reshape(TestCase.input.z_mod,(len(TestCase.input.z_mod),1))
else:
    print("using default inital state")
    perc = 0.1

#Building Model inputs
Morpho = MI.MorphologicInputs(nx=len(TestCase.input.x)-1, discretization_type='custom',\
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
print("Simulation over")

###############################################################################

###############################################################################
#Output of the simulation

#Building output files
if output != 0:
    if not os.path.exists(os.getcwd() + '/simulation_results'):
        os.makedirs(os.getcwd() + '/simulation_results')
        print('Folder for all simulation results saving created')

    if not os.path.exists(os.getcwd() + '/simulation_results/' + TestCase.name):
        os.makedirs(os.getcwd() + '/simulation_results/' + TestCase.name)
        print('Output folder created for the simulated Hillslope')

    out_folder = os.getcwd() + '/simulation_results/' + TestCase.name
    Simulation.output_simu(out_folder)
    print('Simulation results saved')

#Ploting some results
if plot_option !=0:

#    plt.figure(num=1)
#    plt.plot(Simulation.SR.x_node, Simulation.SR.S[-1,:])
#    plt.plot(TestCase.output.x_nodes, TestCase.output.storage[:,-1])
#    plt.ylabel('Storage (m²)')
#    plt.xlabel('X (m)')
#    plt.title('Storage along the hillslope at the last time step')
#    plt.show()

    plt.figure(num=2)
    plt.plot(Simulation.SD.x_node, Simulation.SD.w_node)
    plt.plot(Simulation.SD.x_edges, Simulation.SD.Hs.w_edges)
    plt.ylabel('Width (m)')
    plt.xlabel('X (m)')
    plt.title('Width of the Hillslope (1) nodes (2) edges')
    plt.show()


###############################################################################
