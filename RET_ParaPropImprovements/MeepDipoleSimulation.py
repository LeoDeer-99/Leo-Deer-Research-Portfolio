##  file: MeepDipoleSimulation.py
##
##  Python Script which utilizes meep's fdtd software to simulate a
##  dipole in a non-uniform ice profile. Outputs a plot of the electic
##  everywhere in the grid and an npy array of the electric field arrays
##
##  Programmer:  Leo Deer  deer.279@osu.edu
##
##  Revision history:
##      05/10/22 First version
##
##  Notes:  
##   * Uses point like dipole
##
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

##****************Defining simulation parameters***********************##
name = "no_Er_comp" 

plotting = True     # if true, E-field plots appear

resolution = 12.5	#how many pixels per meter

f = 100         # characteristic frequency in MHz 
a = 1           # characteristic length is 1m

nair = 1      # index of air
hice = 75     # depth of ice
hair = 25     # height of air
R = 100        # radial distance of ice
r = 0          # radius borehole
pad = 2      # padding size for absorbing boundary conditions

depth = 10	# depth of dipole below the ice surface.

#Fuction for a material with an index of reflection that depends of depth n(z).
def matFunc(R):
    z = R[2]
    A = 1.78
    B = 0.43
    C = 0.0132 #1/m
    return mp.Medium(index=A-B*math.exp(-C*(z + Th/2 - hair)))
##*********************************************************************##

##********************Calculated parameters****************************##
wl = 300/f   # wavelength in meep units (300 is c in m*MHz)
meepf = 1/wl       # frequency in meep units
k = 2*np.pi/wl   # wavenumber

TR = R+pad   # total radial distance of grid
Th = hice+hair+2*pad    # height of grid
rcent = r/2        # center radial position of borehole
Rcent = r/2 + TR/2 # center radial position of ice
haircent = -hice/2 # center height position of air
hicecent = hair/2  # center height position of ice
Rice = TR-r	   # size of ice block

mpp = 1/resolution    # meters per pixel
column_size = R/mpp	# number of pixels in the column

vice = 1/nice   	# phase velocity in ice
t_start = 2*(TR)/vice     # approximate time to reach steady state
##*********************************************************************##

##**********************Simulation Setup*******************************##
cell = mp.Vector3(2*TR, mp.inf, Th)

dimensions = mp.CYLINDRICAL

pml_layers = [mp.PML(pad)]

# ice geometry
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
            center=mp.Vector3(0,0,-Th/2 + pad + hair + depth),
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
##*********************************************************************##

##********************** E-Field Qunatities *******************************##
eps_data = sim_dipole.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)	
ez_dipole_data = sim_dipole.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
er_dipole_data = sim_dipole.get_array(center=mp.Vector3(), size=cell, component=mp.Er)
esquare_dipole_data = sim_dipole.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
esquare_dipole_data_log = sim_dipole.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
eabs_dipole_data = sim_dipole.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)

# Array quantites
boarder = int(pad/mpp)				# boarder size in pixels
column_size = len(ez_dipole_data[0]) - 2*boarder	# number of rows in E-field array
row_size = len(ez_dipole_data) - 2*boarder	# number of columns in E-field array
cent_row = ((Th/2 + hicecent)/Th)*row_size	# origin column
cent_col = int(column_size/2)			# origin row

# calculating the E-field
for col in range(int(column_size)):
    for row in range(int(row_size)):
        esquare_dipole_data[row + boarder,col + boarder] = abs(ez_dipole_data[row + boarder,col + boarder])**2 + abs(er_dipole_data[row + boarder,col + boarder])**2	
        eabs_dipole_data[row + boarder,col + boarder] = esquare_dipole_data[row + boarder,col + boarder]**(0.5)
        esquare_dipole_data_log[row + boarder,col + boarder] = (math.log10(esquare_dipole_data[row + boarder,col + boarder]))	# converting that to dB

# Converts array to correct dimensions
comp_data = [[abs(ez_dipole_data[row + boarder,col + boarder + cent_col]) for row in range(int(row_size))] for col in range(int(column_size/2))]
plot_data = [[abs(ez_dipole_data[row + boarder,col + boarder + cent_col]) for row in range(int(row_size))] for col in range(int(column_size/2))]

comp_data = np.transpose(comp_data)

#***********************************************************************************#


##************************** E-field plotting *****************************##
if plotting == True:
    

    # dipole E-field	
    fig, ax = plt.subplots()
    ax.set_aspect('auto')
    plt.imshow(plot_data, aspect='auto', cmap='hot',  vmin=1e-5, vmax=5e-3, extent=(0, R, -hice, hair), alpha=0.9)
    plt.axis('on')
    plot_title = "dipole_"+name+"_d_"+str(depth)
    plt.title(plot_title)
    ax.set_xlabel("radial direction (m)")
    ax.set_ylabel("depth (m)")
    plt.colorbar()
    plt.show()
    fig.savefig(plot_title+".png")
##*******************************************************************##

##**************************Output File******************************##
np.save('MeepDipoleNoEr', comp_data)
print("the array is saved in the file MeepDipoleNoEr.npy")
##*******************************************************************##




