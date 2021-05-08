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

ax.plot([5.696960,5.883747,5.883747,5.696960,5.696960],[-3.163866,-3.162441,-0.173857,-0.175282,-3.163866],[1.742010,1.725730,1.987198,2.003477,1.742010],'k')
ax.plot([5.883747,5.945297,5.945297,5.883747,5.883747],[-3.162441,-3.163387,-0.174803,-0.173857,-3.162441],[1.725730,1.736542,1.998009,1.987198,1.725730],'k')
ax.plot([5.696960,5.883747,5.883747,5.696960,5.696960],[-0.175282,-0.173857,2.814727,2.813303,-0.175282],[2.003477,1.987198,2.248665,2.264944,2.003477],'k')
ax.plot([5.883747,5.945297,5.945297,5.883747,5.883747],[-0.173857,-0.174803,2.813781,2.814727,-0.173857],[1.987198,1.998009,2.259477,2.248665,1.987198],'k')

set_axes_equal(ax)
plt.show()
