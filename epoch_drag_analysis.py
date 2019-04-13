import satellite as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import testcase_satellite as tc
from mayavi import mlab

    
gregdate_list, gregtime_list, drag_list = [], [], []
rx_list, ry_list, rz_list = [], [], []
vx_list, vy_list, vz_list = [], [], []

analysis_tle = tc.niusat_tle_complete

for i in range( 0, len( analysis_tle) - 2, 2 ):
    l1, l2 = analysis_tle[i], analysis_tle[i+1]
    t_temp = st.satellite(l1, l2)
    
    x,y,z = t_temp.r_init
    rx_list.append(x)
    ry_list.append(y)
    rz_list.append(z) 
    
    x,y,z = t_temp.v_init
    vx_list.append(x)
    vy_list.append(y)
    vz_list.append(z)
    
    gregdate_list.append(t_temp.gregdate)
    gregtime_list.append(t_temp.gregtime)
    drag_list.append(t_temp.b_star)
    
sat_df  = pd.DataFrame( {'Date': gregdate_list, 'Time': gregtime_list, 'R_x': rx_list,
                       'R_y': ry_list, 'R_z': rz_list, 'V_x': vx_list, 'V_y': vy_list,
                        'V_z':vz_list, 'BSTAR':drag_list} )
                        
def E_anomaly(row):
    rx, ry, rz = row['R_x'], row['R_y'], row['R_z']
    vx ,vy, vz = row['V_x'], row['V_y'], row['V_z']
    
    r = np.array( [rx, ry, rz], dtype = np.float64 )
    v = np.array( [vx, vy, vz], dtype = np.float64 )
    h = np.cross( r, v )
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)
    

    mu = 398600.4415
    
    e = np.cross(v,h)/mu - r/r_norm
    e_norm = np.linalg.norm(e)
    
    theta_temp = np.arccos( np.dot(e,r)/( e_norm*r_norm) )
    if np.dot(r,v) >= 0:
        theta = theta_temp
    else:
        theta = 2*np.pi - theta_temp
        
    E_anly = 2* np.arctan ( np.tan(theta/2) * ( (1 - e_norm)/(1 + e_norm) )**0.5 )
    E_anly = E_anly*180/np.pi
    
    return E_anly


def eccen(row):
    rx, ry, rz = row['R_x'], row['R_y'], row['R_z']
    vx ,vy, vz = row['V_x'], row['V_y'], row['V_z']
    
    r = np.array( [rx, ry, rz], dtype = np.float64 )
    v = np.array( [vx, vy, vz], dtype = np.float64 )
    h = np.cross( r, v )
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)    

    mu = 398600.4415
    
    e = np.cross(v,h)/mu - r/r_norm
    e_norm = np.linalg.norm(e)
    return e_norm


def radius_norm(row):
    rx, ry, rz = row['R_x'], row['R_y'], row['R_z']
    r = np.array( [rx, ry, rz], dtype = np.float64 )
    r_norm = np.linalg.norm(r)
    return r_norm

    
sat_df['Eccentric_Anomaly'] = sat_df.apply( E_anomaly, axis = 1)
sat_df['Radius_Norm'] = sat_df.apply( radius_norm, axis = 1)
sat_df['Eccentricity'] = sat_df.apply( eccen, axis = 1)
                        
columnTitles = ['Date', 'Time', 'R_x', 'R_y', 'R_z', 'V_x', 'V_y', 'V_z', 'BSTAR', 'Eccentric_Anomaly', 'Radius_Norm', 'Eccentricity']
sat_df = sat_df.reindex(columns = columnTitles)
sat_df.sort_values (by = ['Eccentric_Anomaly'], inplace = True)


csv_title = 'pratham' + '_Epoch_Analysis'
#------Exporting the Panda Dataframe as a CSV------
# NOTE: Change the Title of the CSV File before saving.

'''sat_df.to_csv(csv_title + '.csv', index = True)
sat_df.describe().to_csv(csv_title + '_Describe.csv', index = True)'''


#print(sat_df['Date'])

'''Radius_Earth = 6378.1363
pts = mlab.points3d(rx_list, ry_list, rz_list,drag_list, mask_points = 1)
phi, theta = np.mgrid [0:np.pi:30j, 0:2*np.pi:20j]
x = Radius_Earth*np.sin(phi)*np.cos(theta)
y = Radius_Earth*np.sin(phi)*np.sin(theta)
z = Radius_Earth*np.cos(phi)

mlab.mesh(x,y,z,representation = 'wireframe')
mlab.show()'''
