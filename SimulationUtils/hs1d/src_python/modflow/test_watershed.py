# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 16:17:01 2016

@author: Quentin Courtois
"""
#Add module's path
import sys
import os
sys.path.append(os.getcwd() + '/process/')
sys.path.append(os.getcwd() + '/utils/')

import numpy as np
import ModflowModel2D as MM2D
import WatershedObj as WO
import HydroObj as HO
import WatershedReader as WR
import os
import time
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#Set Options for the Run
plot_options = 1
output = 1
initial = 1

#name = 'Laita'
#name = 'regular_watershed'
#name = 'Odet'
#name = 'Vilaine'
name = 'Kerbernez'

watershed = WR.WatershedReader(name)
#Load recharge
cwd = os.getcwd()
filename = cwd + "\\Laita\\Values_mean_drain_1"
with open(filename,'r'):
    recharge = np.loadtxt(filename)
recharge = recharge.astype(np.float)
recharge = recharge[:,1]
a = np.where(recharge != 0)
a = a[0]
rec_temp = np.loadtxt(os.getcwd()+'/Kerbernez/hydrologic.input', skiprows=3)
recharge_true = []
time_rec = []
n = 1
for i in range(int(np.ceil(len(rec_temp[:,1])/n))):
    recharge_true.append(rec_temp[n*i,1]*3600)
    time_rec.append(rec_temp[n*i,0])
#for i in range(len(recharge_true)):
#    if recharge_true[i]>0:
#        recharge_true[i]=0.00235
#recharge_true = temp[:,1]*86400
rec_temp_bis = recharge_true[119:838]
rec_temp = []
for j in range(1,11):
    for i in range(1,len(rec_temp_bis)):
        rec_temp.append(rec_temp_bis[i])

recharge_true = rec_temp
plt.figure(1)
plt.plot(recharge_true)
plt.show()
#recharge_true = -recharge[a[0]+1:130+a[0]]
recharge_constant = np.zeros((1,20000))+np.mean(recharge_true)
recharge_constant = recharge_constant[0]

start_time = time.time()

temporary = np.diff(watershed.coord_x)
temporary = np.unique(temporary)
coord_temporary = np.where(temporary>0)
temporary = temporary[coord_temporary]

hydro = HO.HydroObj(hk=1, sy=0.3, ss=0.3, soil_depth=2, percentage_loaded=0.1, \
                    recharge_initial = recharge_constant, recharge_true = recharge_true, \
                    time_select=len(recharge_true)-1)

watershed_obj = WO.WatershedObj(name=name, coord_x=watershed.coord_x, coord_y= watershed.coord_y, \
                                elevation=watershed.elevation, bottom=-1, cell_size=temporary[0], \
                                outlet=watershed.outlet)


Mod = MM2D.ModflowModel2D(hydro = hydro, watershed = watershed_obj, output = output, init = initial)
print("------- %s seconds -------" % (time.time() - start_time))



if plot_options != 0:
    Mod.plot_output(name, watershed_obj, hydro)
#    fig2 = plt.figure(num=1)
#    ax2 = fig2.gca(projection = '3d')
#    X = np.unique(watershed_obj.xcustom)
#    Y = np.unique(watershed_obj.ycustom)
#    X,Y = np.meshgrid(X,Y)
#    surf = ax2.plot_surface(X,Y,Mod.piezo,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0,antialiased=False)
#    fig2.colorbar(surf,shrink=0.5, aspect=5)
