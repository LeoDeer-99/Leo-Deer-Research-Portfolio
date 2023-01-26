##  file: OptimizerMeep1.npy
##  phase offset and depth of the array.  
##
##  Programmer:  Leo Deer  deer.279@osu.edu
##
##  Revision history:
##  05/10/22 - Revised comments and cleaned up
##
##  Notes:  
##   * Uses single dipole as comparison metric
##   * Uses point like dipoles
##   * Plotting gives the option for 3 plots to appear
##	- Array E-field squared in dB
##	- Dipole E-field squared in dB
##	- The difference in the two E-fields
##   * This was abandoned because ParaProp works better for this project
##*********************************************************************##

##******************Importing useful packages**************************##
from __future__ import division

import meep as mp
import math
import numpy as np
from numpy import asarray
from numpy import save
import matplotlib.pyplot as plt
##*********************************************************************##

##*********************************************************************##
##************************** Set up ***********************************##
##*********************************************************************##

##****************Defining simulation parameters***********************##
name = "9_27_21_long" 

plotting = False     # if true, the 3 E-field plots appear

f = 100         # characteristic frequency in MHz 
a = 1           # characteristic length is 1m

nice = 1.3    # index of ice
nair = 1      # index of air
hice = 60     # depth of ice
hair = 10     # height of air
R = 100        # radial distance of ice
pad = 2      # padding size

depth_array = np.arange(0,20,2)
phase_array = np.arange(10,40,3)

#Fuction for a material with an index of reflection that depends of depth n(z).
def matFunc(R):
    z = R[2]
    A = 1.78
    B = 0.43
    C = 0.0132 #1/m
    return mp.Medium(index=A-B*math.exp(-C*(z + Th/2 - hair)))

NumAntennas = 8 # number of antennas in the array

resolution = 10

measure_depth = 3    # depth(m) at which area integral calculated
Rmax = 100    # max r to area integral
##********************Calculated parameters****************************##
wl = 300/f   # wavelength in meep units (300 is c in m*MHz)
meepf = 1/wl       # frequency in meep units
k = 2*np.pi/wl   # wavenumber

TR = R+pad   # radial distance of ice
Th = hice+hair     # height of set up
r = wl/10          # radius borehole
rcent = r/2        # center radial value of borehole
Rcent = r/2 + TR/2 # center radial value of ice
haircent = -hice/2 # center height of air
hicecent = hair/2  # center height of  ice
Rice = TR-r    # size of ice block

mpp = 1/resolution    # meters per pixel
column_size = R/mpp

antennaSep = wl
antennaSide = wl/4

vice = 1/nice   # phase velocity in ice
t_start = 2*(TR)/vice     # approximate time to reach steady state


##********************Creating Arrays****************************##
rows = len(phase_array)
columns = len(depth_array)
size = column_size

#array for npy file
arrayESlices = []
arrayNormESlices = []
dipoleESlices = []
dipoleNormESlices = []
for _ in range(rows):
	row = []
	for _ in range(columns):
		column = []
		for _ in range(int(size)):
			column.append(0)
		row.append(column)
	arrayNormESlices.append(row)

for _ in range(rows):
	row = []
	for _ in range(columns):
		column = []
		for _ in range(int(size)):
			column.append(0)
		row.append(column)
	arrayESlices.append(row)    

for _ in range(columns):
	column = []
	for _ in range(int(size)):
		column.append(0)
	dipoleNormESlices.append(column)

for _ in range(columns):
	column = []
	for _ in range(int(size)):
		column.append(0)
	dipoleESlices.append(column)

#array for optimization maps
rows = len(phase_array)
columns = len(depth_array)
area_matrix_dBd = []
for _ in thresholds_dBd_array:
	thresholds = []
	for _ in range(rows):
		row = []		
		for _ in range(columns):					
			row.append(0)
		thresholds.append(row)
	area_matrix_dBd.append(thresholds)

area_matrix_abs = []
for _ in thresholds_abs_array:
	thresholds = []
	for _ in range(rows):
		row = []		
		for _ in range(columns):					
			row.append(0)
		thresholds.append(row)
	area_matrix_abs.append(thresholds)


##**********************Simulation Setup*******************************##
cell = mp.Vector3(2*TR, mp.inf, Th)

dimensions = mp.CYLINDRICAL

geometry = [mp.Block(center=mp.Vector3(rcent,0),
		size=mp.Vector3(r,mp.inf,hice),
		material=mp.Medium(index=nair)),
	mp.Block(center=mp.Vector3(Rcent,0,haircent),
		size=mp.Vector3(Rice,mp.inf,hair),
		material=mp.Medium(index=nair)),
	mp.Block(center=mp.Vector3(Rcent,0,hicecent),
		size=mp.Vector3(Rice,mp.inf,hice),
		material=matFunc)]
pml_layers = [mp.PML(pad)]

for d in range(len(depth_array)):
	print("in d loop with d = " + str(d))
	depth = depth_array[d]

    #************************* Dipole Reference *************************##
    #   At each depth we need a dipole to compare the phased array to    ##
    #********************************************************************##
    # dipole ice geometry
	geometry_dipole = [mp.Block(center=mp.Vector3(rcent,0),
				size=mp.Vector3(r,mp.inf,hice),
				material=mp.Medium(index=nair)),
				mp.Block(center=mp.Vector3(Rcent,0,haircent),
				size=mp.Vector3(Rice,mp.inf,hair),
				material=mp.Medium(index=nair)),
				mp.Block(center=mp.Vector3(Rcent,0,hicecent),
				size=mp.Vector3(Rice,mp.inf,hice),
				material=matFunc)]

	# create the source
	sources_dipole = []
	sources_dipole.append(mp.Source(mp.ContinuousSource(frequency=meepf),
				component=mp.Ez,
				center=mp.Vector3(0,0,-Th/2 + hair + depth + (antennaSep/2)*(NumAntennas - 1)),
				size=mp.Vector3(0,0,0)))
	# create simulation
	sim_dipole = mp.Simulation(force_complex_fields=True,
					cell_size=cell,
					dimensions=mp.CYLINDRICAL,
					boundary_layers=pml_layers,
					geometry=geometry_dipole,
					sources=sources_dipole,
					resolution=resolution)

	# this runs the simulation until we have reached the steady state	
	sim_dipole.run(until=t_start)

	# E-field quantities	
	ez_dipole_data = sim_dipole.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
	er_dipole_data = sim_dipole.get_array(center=mp.Vector3(), size=cell, component=mp.Er)
	esquare_dipole_data = sim_dipole.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
	esquare_dipole_data_log = sim_dipole.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
	eabs_dipole_data = sim_dipole.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)

	# Array quantites
	row_size = len(ez_dipole_data)			# number of columns in E-field array
	cent_row = ((Th/2 + hicecent)/Th)*row_size	# origin column
	cent_col = column_size/2			# origin row
	boarder = round(pad/mpp)			# boarder size in pixels

	# calculating the Dipole E-field quantities
	for col in range(int(column_size)):
		for row in range(int(row_size)):
			esquare_dipole_data[row,col] = abs(ez_dipole_data[row,col])**2 + abs(er_dipole_data[row,col])**2
			eabs_dipole_data[row,col] = esquare_dipole_data[row,col]**(0.5)
			esquare_dipole_data_log[row,col] = (math.log10(esquare_dipole_data[row,col]))	# converting that to dB
    
    #***********************************************************************************#
    
	for p in range(len(phase_array)):
		print("in p loop with p = " + str(p))
		print("in d loop with d = " + str(d))

	#*************************** Phased Array ***************************##
	#  	          This creates our array at each angle 		     ##
	#********************************************************************##
		dphi = phase_array[p]
		dphir = dphi*np.pi/180

		#Creating the sources
		sources = []
		for ns in range(NumAntennas):
			sources.append(mp.Source(mp.ContinuousSource(frequency=meepf),
				component=mp.Ez,
				center=mp.Vector3(0,0,-Th/2+hair+depth+antennaSep*ns),
				size=mp.Vector3(0,0,0),
				amplitude = np.exp((-1j)*k*antennaSep*ns*dphir)))

		# create simulation
		sim = mp.Simulation(force_complex_fields=True,
			    cell_size=cell, 
			    dimensions=mp.CYLINDRICAL, 
			    boundary_layers=pml_layers, 
			    geometry=geometry, 
			    sources=sources, 
			    resolution=resolution)

		# this rus the simulation until we have reached the steady state	
		sim.run(until=t_start)

		# E-field quantities		
		eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
		ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
		er_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Er)
		esquare_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
		esquare_data_log = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
		eabs_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
		residual = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)

		# calculating the E-fields and difference in Dipole and Array E-fields
		for col2 in range(int(column_size)):
			for row2 in range(int(row_size)):
				esquare_data[row2,col2] = abs(ez_data[row2,col2])**2 + abs(er_data[row2,col2])**2	
				eabs_data[row2,col2] = esquare_data[row2,col2]**(0.5)
				esquare_data_log[row2,col2] = (math.log10(esquare_data[row2,col2])) # converting to dB
				residual[row2,col2] = eabs_data[row2,col2] - eabs_dipole_data[row2,col2]

##************************** E-field plotting *****************************##
		if plotting == True:
			# phased array E-field
			fig, ax = plt.subplots()
			ax.set_aspect('auto')
			plt.imshow(eps_data, interpolation='spline36', cmap='binary')
			plt.imshow(eabs_data, aspect='auto', cmap='hot',  vmin=1e-5, vmax=1e-2, extent=(-R, R, -hice, hair), alpha = 0.9)
			plot_title = "esquare_"+name+"_d_"+str(depth)+"_angle_"+str(dphi)+".png"
			plt.title(plot_title)
			plt.axis('on')
			ax.set_xlabel("radial direction (m)")
			ax.set_ylabel("depth (m)")
			plt.colorbar()
			#plt.show()
			fig.savefig(plot_title)

			# dipole E-field	
			fig, ax = plt.subplots()
			ax.set_aspect('auto')
			plt.imshow(eps_data, interpolation='spline36', cmap='binary')
			plt.imshow(eabs_dipole_data, aspect='auto', cmap='hot',  vmin=1e-5, vmax=1e-2, extent=(-R, R, -hice, hair), alpha=0.9)
			plot_title = "dipole_"+name+"_d_"+str(depth)+".png"
			plt.title(plot_title)
			plt.axis('on')
			ax.set_xlabel("radial direction (m)")
			ax.set_ylabel("depth (m)")
			plt.colorbar()
			#plt.show()
			fig.savefig(plot_title)

			# residual E-field	
			fig, ax = plt.subplots()
			ax.set_aspect('auto')
			plt.imshow(eps_data, interpolation='spline36', cmap='binary')
			plt.imshow(residual, aspect='auto', cmap='hot',  vmin=1e-5, vmax=1e-2, extent=(-R, R, -hice, hair), alpha=0.9)
			plot_title = "residual_"+name+"_d_"+str(depth)+"_angle_"+str(dphi)+".png"
			plt.title(plot_title)
			plt.axis('on')
			ax.set_xlabel("radial direction (m)")
			ax.set_ylabel("depth (m)")
			plt.colorbar()
			#plt.show()
			fig.savefig(plot_title)

		##*******************************************************************##
		##**************************Optimization*****************************##
		##*******************************************************************##

		row_depth = round(((measure_depth + hair)/Th)*row_size)	# measure depth in pixels

		for r_out_dBd in range(round(column_size)-1):
			arrayESlices[p][d][round(r_out_dBd)] = eabs_data[row_depth,round(r_out_dBd+cent_col)]
			dipoleESlices[d][round(r_out_dBd)] = eabs_dipole_data[row_depth,round(r_out_dBd+cent_col)]

		row_depth = round(((measure_depth + hair)/Th)*row_size)	# measure depth in pixels

		pixel_max = round(Rmax/mpp)	# max r to take area integral in pixels
		area_integral_dBd = 0
		for t_dBd in range(len(thresholds_dBd_array)):
			threshold_dBd = thresholds_dBd_array[t_dBd]
			area_integral_dBd = 0
			for r_out_dBd in range(int(round(column_size/2)-1)):
				if r_out_dBd < pixel_max:
					if threshold_dBd < 10*(esquare_data_log[row_depth, round(r_out_dBd + cent_col)] - esquare_dipole_data_log[row_depth, round(r_out_dBd + cent_col)]):		
						area_integral_dBd += r_out_dBd*mpp*mpp

			area_matrix_dBd[t_dBd][p][d] = 2*np.pi*area_integral_dBd

		area_integral_abs = 0
		for t_abs in range(len(thresholds_abs_array)):
			threshold_abs = thresholds_abs_array[t_abs]
			area_integral_abs = 0
			for r_out_abs in range(int(round(column_size/2)-1)):
				if r_out_abs < pixel_max:
					if threshold_abs < 10*math.log10(eabs_data[row_depth, round(r_out_abs + cent_col)]**2):		
						area_integral_abs += r_out_abs*mpp*mpp

			area_matrix_abs[t_abs][p][d] = 2*np.pi*area_integral_abs


# color plots
for t_dBd in range(len(thresholds_dBd_array)):
	threshold_dBd=thresholds_dBd_array[t_dBd]
	fig, ax = plt.subplots()
	ax.set_aspect('auto')
	plt.imshow(area_matrix_dBd[t_dBd], extent=[depth_array[0],depth_array[-1],phase_array[-1],phase_array[0]], aspect = "auto", vmin=0, vmax=31420, interpolation='spline36', cmap='gist_heat', alpha=1)
	ax.set_xlabel("depth (m)")
	ax.set_ylabel("phase (degrees)")
	plt.axis('on')	
	plot_title_dBd = "dBd test"
	plt.title(plot_title_dBd)
	plt.colorbar()
	fig.savefig(plot_title_dBd)

for t_abs in range(len(thresholds_abs_array)):
	threshold_abs=thresholds_abs_array[t_abs]
	fig, ax = plt.subplots()
	ax.set_aspect('auto')
	plt.imshow(area_matrix_abs[t_abs], extent=[depth_array[0],depth_array[-1],phase_array[-1],phase_array[0]], aspect = "auto", vmin=0, vmax=31420, interpolation='spline36', cmap='gist_heat', alpha=1)
	ax.set_xlabel("depth (m)")
	ax.set_ylabel("phase (degrees)")
	plt.axis('on')	
	plot_title_abs = "abs test"
	plt.title(plot_title_abs)
	plt.colorbar()
	plt.show()



# create NPY file
NPYData = asarray([[arrayESlices,dipoleESlices]],dtype=object)

save('EFields'+str(name)+'meas'+str(measure_depth)+'mp.npy', NPYData)


# Use optimizerNAnpyMeep2.py from this point on to more efficiently play with thresholds
# The rest is just here so the optimization can finish, but I suggest running optimizerNAnpyMeep2.py
# once you have your NPY file.

rows = len(phase_array)
columns = len(depth_array)

pixel_max = Rmax/mpp

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

    
threshold_dBd = 1
threshold_abs = -65 

for d in range(len(depth_array)):
	depth_64 = depth_array[d]
	depth = float(depth_64)
	for p in range(len(phase_array)): 
        #print("in p for look with p =" + str(p))
		phi = np.deg2rad(phase_array[p])
        
		area_integral_dBd = 0
		for r_out_dBd in range(round(column_size/2)-1):
			if r_out_dBd < pixel_max:
				#if threshold_dBd < 10*math.log10((arrayESlices[p][d][round(r_out_dBd)]/dipoleESlices[d][round(r_out_dBd)])**2):
				if threshold_dBd < 10*math.log10((NPYData[0][0][p][d][round(r_out_dBd)]/NPYData[0][1][d][round(r_out_dBd)])**2):
					area_integral_dBd += r_out_dBd*mpp*mpp
		area_matrix_dBd[p][d] = 2*np.pi*area_integral_dBd

		area_integral_abs = 0   
		for r_out_abs in range(round(column_size/2)-1):
			if r_out_abs < pixel_max:
				if NPYData[0][0][p][d][round(r_out_abs)] == 0:
					#arrayESlices[p][d][round(r_out_abs)] = 0.0000000001
                    			NPYData[0][0][p][d][round(r_out_abs)] = 0.0000000001
				if threshold_abs < 10*math.log10((arrayESlices[p][d][round(r_out_abs)])**2):
				#if threshold_abs < 10*math.log10((NPYData[0][0][p][d][round(r_out_abs)])**2):
					area_integral_abs += r_out_abs*mpp*mpp
		area_matrix_abs[p][d] = 2*np.pi*area_integral_abs                 
                
fig, ax = plt.subplots()
ax.set_aspect('auto')
	#area_trans_dBd = np.transpose(area_matrix_dBd)
	#heat_map_dBd = np.transpose(area_trans_dBd[0])
heat_map_dBd = area_matrix_dBd

plt.imshow(heat_map_dBd, extent=[depth_array[0],depth_array[-1],phase_array[-1],phase_array[0]], aspect = "auto", vmin=0, vmax=31420, interpolation='spline36', cmap='gist_heat', alpha=1)
ax.set_xlabel("depth (m)")
ax.set_ylabel("phase (degrees)")
plt.axis('on')
plot_title_dBd = "dbd npy"
plt.title(plot_title_dBd)
plt.colorbar()
plt.show()
	#fig.savefig(plot_title_dBd)                 
                
fig, ax = plt.subplots()
ax.set_aspect('auto')
	#area_trans_abs = np.transpose(area_matrix_abs)
	#heat_map_abs = np.transpose(area_trans_abs[0])
heat_map_abs = area_matrix_abs

plt.imshow(heat_map_abs, extent=[depth_array[0],depth_array[-1],phase_array[-1],phase_array[0]], aspect = "auto", vmin=0, vmax=31420, interpolation='spline36', cmap='gist_heat', alpha=1)
ax.set_xlabel("depth (m)")
ax.set_ylabel("phase (degrees)")
plt.axis('on')
plot_title_abs = "abs npy"
plt.title(plot_title_abs)
plt.colorbar()
plt.show()
	#fig.savefig(plot_title_dBd)  
            
dBd_max = 0 
phase_dBd_max = 0
depth_dBd_max = 0
max_threshold_dBd = 0
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
max_threshold_abs = 0
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



