##  file: optimizerNAnpyMeep2.py
##
##  Python Script which utilizes meep's fdtd software to simulate a
##  dipole array in a non-uniform ice profile in order to find the optimal
##  phase offset and depth of the array.  
##
##  Programmer:  Leo Deer  deer.279@osu.edu
##
##  Revision history:
##  05/10/22 - adding comments and cleaning up
##
##  Notes:  
##   * Uses single dipole as comparison metric
##   * Uses point like dipoles
##   * Plotting gives the option for 3 plots to appear
##  - Array E-field squared in dB
##  - Dipole E-field squared in dB
##  - The difference in the two E-fields
##   * This was abandoned because ParaProp works better for this project
##*********************************************************************##

##******************Importing useful packages**************************##
from __future__ import division

import meep as mp
import math
import numpy as np
from numpy import load
import matplotlib.pyplot as plt
##*********************************************************************##
name = "9_27_21_long"


#make sure these match optimizerNAnpyMeep1.py 
hice = 60     # depth of ice
hair = 10     # height of air
R = 100        # radial distance of ice

depth_array = np.arange(0,20,2)
phase_array = np.arange(10,40,3)

resolution = 10
mpp = 1/resolution    # meters per pixel
column_size = R/mpp

measure_depth = 3      # depth(m) at which area integral calculated
Rmax = 100      # max r to area integral
pixel_max = Rmax/mpp

NPYData = np.load('EFields'+str(name)+'meas'+str(measure_depth)+'mp.npy', allow_pickle=True)


### creating array for optimization maps ###
rows = len(phase_array)
columns = len(depth_array)

area_matrix_dBd = []
for _ in range(rows):
    row = []
    for _ in range(columns):
        row.append(0)
    area_matrix_dBd.append(row)

area_matrix_abs = []
for _ in range(rows):
    row = []
    for _ in range(columns):
        row.append(0)
    area_matrix_abs.append(row)

    
# set the levels for the area calculation
threshold_dBd = 11
threshold_abs = -55 

for d in range(len(depth_array)):
	print("in d for loop with d =" + str(d))
	depth_64 = depth_array[d]
	depth = float(depth_64)
	for p in range(len(phase_array)):
		print("in p for loop with p =" + str(p))
		phi = np.deg2rad(phase_array[p])
		
		area_integral_dBd = 0
		for r_out_dBd in range(round(column_size/2)-1):
			if r_out_dBd < pixel_max:
				if threshold_dBd < 10*math.log10((NPYData[0][0][p][d][round(r_out_dBd)]/NPYData[0][1][d][round(r_out_dBd)])**2):
					area_integral_dBd += r_out_dBd*mpp*mpp
		area_matrix_dBd[p][d] = 2*np.pi*area_integral_dBd

		area_integral_abs = 0   
		for r_out_abs in range(round(column_size/2)-1):
			if r_out_abs < pixel_max:
				if NPYData[0][0][p][d][round(r_out_abs)] == 0:
		    			NPYData[0][0][p][d][round(r_out_abs)] = 0.0000000001
				if threshold_abs < 10*math.log10((NPYData[0][0][p][d][round(r_out_abs)])**2):
					area_integral_abs += r_out_abs*mpp*mpp
		area_matrix_abs[p][d] = 2*np.pi*area_integral_abs                 
                        
                
fig, ax = plt.subplots()
ax.set_aspect('auto')
heat_map_dBd = area_matrix_dBd

plt.imshow(heat_map_dBd, extent=[depth_array[0],depth_array[-1],phase_array[-1],phase_array[0]], aspect = "auto", vmin=0, vmax=31420, interpolation='spline36', cmap='gist_heat', alpha=1)
ax.set_xlabel("depth (m)")
ax.set_ylabel("phase (degrees)")
plt.axis('on')
plt.colorbar()
plt.show()                 
                
fig, ax = plt.subplots()
ax.set_aspect('auto')
heat_map_abs = area_matrix_abs

plt.imshow(heat_map_abs, extent=[depth_array[0],depth_array[-1],phase_array[-1],phase_array[0]], aspect = "auto", vmin=0, vmax=31420, interpolation='spline36', cmap='gist_heat', alpha=1)
ax.set_xlabel("depth (m)")
ax.set_ylabel("phase (degrees)")
plt.axis('on')
plt.colorbar()
plt.show()
            
dBd_max = 0 
phase_dBd_max = 0
depth_dBd_max = 0
for p in range(len(phase_array)):
    for d in range(len(depth_array)):
        if dBd_max < area_matrix_dBd[p][d]:
            dBd_max = area_matrix_dBd[p][d]
            phase_dBd_max = phase_array[p]
            depth_dBd_max = depth_array[d]
                
                
print('The best parameters from dBd are \n depth:' + str(depth_dBd_max) + '\n phase:' + str(phase_dBd_max) + '\n with an area of:' + str(dBd_max))
                
abs_max = 0 
phase_abs_max = 0
depth_abs_max = 0
for p in range(len(phase_array)):
    for d in range(len(depth_array)):
        if abs_max < area_matrix_abs[p][d]:
            abs_max = area_matrix_abs[p][d]
            phase_abs_max = phase_array[p]
            depth_abs_max = depth_array[d]              

                
print('The best parameters from abs are \n depth:' + str(depth_abs_max) + '\n phase:' + str(phase_abs_max) + '\n with an area of:' + str(abs_max))


##*******************************************************************##
##****************************** fin ********************************##
##*******************************************************************##


