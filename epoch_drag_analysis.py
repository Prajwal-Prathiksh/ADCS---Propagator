import satellite as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import testcase_satellite as tc
    
gregdate_list, gregtime_list, drag_list = [], [], []
rx_list, ry_list, rz_list = [], [], []
vx_list, vy_list, vz_list = [], [], []

    
for i in range( 0, len( tc.tleINS_full), 2 ):
    l1, l2 = tc.tleINS_full[i], tc.tleINS_full[i+1]
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
    
sat_df  =pd.DataFrame( {'Date': gregdate_list, 'Time': gregtime_list, 'R_x': rx_list,
                       'R_y': ry_list, 'R_z': rz_list, 'V_x': vx_list, 'V_y': vy_list,
                        'V_z':vz_list, 'BSTAR':drag_list} )
                        
columnTitles = ['Date', 'Time', 'R_x', 'R_y', 'R_z', 'V_x', 'V_y', 'V_z', 'BSTAR']
sat_df = sat_df.reindex(columns = columnTitles)

#------Exporting the Panda Dataframe as a CSV------
#sat_df.to_csv('hell.csv', index = False)