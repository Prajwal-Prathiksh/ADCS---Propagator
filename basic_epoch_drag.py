import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('ins_1c_tle.txt', sep = ' ',header = None)

df.columns = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15']

drag_df = df[ ['1','6','12'] ][0::2]

def str_to_float(x):
    """ Convert the B-star term given by the TLE into a 
        useful floating-point number.
               
        Parameters:
        -----------
        x: string
           Takes the B-star value given according to TLE
           format.
        """ 
    y = x.split('-')
    z = (float(y[0])*10**-5) * (10**(-float(y[1])))
    return z

drag_df['12'] = drag_df['12'].apply(str_to_float)
drag_df.columns = ['Line Number', 'Satellite - Epoch', 'B-star']
print(drag_df.describe())
