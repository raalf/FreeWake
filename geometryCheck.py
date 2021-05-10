#Plotting wake results that were generated with #FreeWake2020_dev V0.4
#Plotting wake of timestep 0

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def set_axes_equal(ax):
	x_limits = ax.get_xlim3d()
	y_limits = ax.get_ylim3d()
	z_limits = ax.get_zlim3d()
	x_range = abs(x_limits[1] - x_limits[0])
	x_middle = np.mean(x_limits)
	y_range = abs(y_limits[1] - y_limits[0])
	y_middle = np.mean(y_limits)
	z_range = abs(z_limits[1] - z_limits[0])
	z_middle = np.mean(z_limits)
	plot_radius = 0.5*max([x_range, y_range, z_range])
	ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
	ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
	ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect('equal')

ax.set_xlabel('Global X')
ax.set_ylabel('Global Y')
ax.set_zlabel('Global Z')

ax.plot([5.000000,5.175000,5.175000,5.000000,5.000000],[-1.250000,-1.250000,1.250000,1.250000,-1.250000],[-2.165064,-2.165064,2.165064,2.165064,-2.165064],'k')
ax.plot([5.175000,5.326554,5.326554,5.175000,5.175000],[-1.250000,-1.174223,1.325777,1.250000,-1.250000],[-2.165064,-2.208814,2.121314,2.165064,-2.165064],'k')

set_axes_equal(ax)
plt.show()
