import satellite as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import testcase_satellite as tc
from mpl_toolkits.mplot3d import Axes3D
import epoch_drag_analysis as eda
from mayavi import mlab


Radius_Earth = 6378.1363
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.set_xlabel(r'x-axis $\rightarrow$')
ax.set_ylabel(r'y-axis $\rightarrow$')
ax.set_zlabel(r'z-axis $\rightarrow$')

def plots(ch, i):
    '''
    :ch: choice number
    
    :i: mask_points
    '''
    
    if ch == 1:
        x, y, z = eda.rx_list[::i], eda.ry_list[::i], eda.rz_list[::i]
        b_list = eda.drag_list[::i]
        
        ax.scatter(x,y,z, marker = 'o',c = b_list)
        fig.colorbar(ax.scatter(x,y,z, marker = 'o',c = b_list))
        plt.show()
        ax.set_aspect("equal")
        
        
    elif ch == 2:
        x, y, z = eda.rx_list[::i], eda.ry_list[::i], eda.rz_list[::i]
        b_list = eda.drag_list[::i]
        
        ax.scatter(x,y,z, marker = 'o',c = b_list)
        fig.colorbar(ax.scatter(x,y,z, marker = 'o',c = b_list))
        plt.show()      
        
        u,v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = Radius_Earth*np.cos(u)*np.sin(v)
        y = Radius_Earth*np.sin(u)*np.sin(v)
        z = Radius_Earth*np.cos(v)
        ax.plot_wireframe(x,y,z)        
        ax.set_aspect("equal")
	

'''
pts = mlab.points3d(eda.rx_list, eda.ry_list, eda.rz_list,eda.drag_list, mask_points = 2)
phi, theta = np.mgrid [0:np.pi:30j, 0:2*np.pi:20j]
x = Radius_Earth*np.sin(phi)*np.cos(theta)
y = Radius_Earth*np.sin(phi)*np.sin(theta)
z = Radius_Earth*np.cos(phi)

mlab.mesh(x,y,z,representation = 'wireframe')
mlab.show()
'''
