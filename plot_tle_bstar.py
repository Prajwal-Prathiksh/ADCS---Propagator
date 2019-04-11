import satellite as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import testcase_satellite as tc
from mpl_toolkits.mplot3d import Axes3D
import epoch_drag_analysis as eda


Radius_Earth = 6378.1363
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

def plots(ch):
	if ch == 1:
		p = ax.scatter(eda.rx_list, eda.ry_list, eda.rz_list, marker = 'o',c = eda.drag_list)

		plt.show(p)
		fig.colorbar(p)
		
		
	
	elif ch == 2:
		p = ax.scatter(eda.rx_list, eda.ry_list, eda.rz_list, marker = 'o',c = eda.drag_list) 
		plt.show(p)
		fig.colorbar(p)

		u,v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
		x = Radius_Earth*np.cos(u)*np.sin(v)
		y = Radius_Earth*np.sin(u)*np.sin(v)
		z = Radius_Earth*np.cos(v)
		ax.plot_wireframe(x,y,z)

		ax.set_aspect("equal")
		plt.show()


ch = int(input('Enter Choice:'))
plots(ch)




