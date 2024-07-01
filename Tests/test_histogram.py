# -*- coding: utf-8 -*-
"""
generate a histogram of steady state heights.
"""

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import Functions as f
from BacPlot import plot
from histogram import plot_histogram

# def steady_state(data_path, config_path):
#     '''

#     Parameters
#     ----------
#     data_path : string
#         path to folder containing subfolders (run folders) of raw output data from BacStroke.
#     no_data_points : integer
#         number of data points provided in each data file (must be the same for all runs)

#     Returns
#     -------
#     None.

#     '''
    
#     # number of data sets from this config
#     runs = os.listdir(data_path)
    
#     # config constants
#     constants = f.read_config(config_path)
#     dt = int(constants[1]) # timestep
#     T = int(constants[2]) # total sim time
#     output_interval = int(constants[-1]) # interval between positional readout
    
#     # how many data points were output in each sim? (avoids prematurely opening a file and reading no. lines)
#     no_data_points = int((T/dt)/output_interval)
        
#     # array to store all y values
#     ys = np.zeros((no_data_points, len(runs)))
    
#     for i in range(len(runs)):
        
#         # read data file
#         path = f'{data_path}/{runs[i]}/output.csv'
        
#         # read in output from bacstroke
#         df = np.array(pd.read_csv(path, sep=',', header = None))
        
#         # grab y vals - store
#         ys[:, i] = df[:, 1]
        
    
        
        
#         #plot(config_path, path, 'na')
        
#     # # find combined standard deviation at each time point
#     # stdev = np.std(ys, axis = 1)
    
#     # # expect back end of the data to be in steady state (MUST HAVE BEEN RUN FOR LONG TIME)
#     back_90 = int((no_data_points/100)*90) # index marking 90% of data points
    
#     # mean and standard deviation of backend of the data set ()
    
#     mean_90 = np.mean(back_90)
    
def steady_state_single(data_path, config_path):
    '''
    Gets steady state of a data set

    Parameters
    ----------
    data_path : TYPE
        DESCRIPTION.
    config_path : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    # config constants
    constants = f.read_config(config_path)
    dt = int(constants[1]) # timestep
    T = int(constants[2]) # total sim time
    output_interval = int(constants[-1]) # interval between positional readout
    
    # how many data points were output in each sim? (avoids prematurely opening a file and reading no. lines)
    no_data_points = int((T/dt)/output_interval)
    
    # read in output from bacstroke
    df = np.array(pd.read_csv(data_path, sep=',', header = None))
    
    # grab y vals - store
    y = df[:, 1]
    
    # index marking back 90% of data points
    ind90 = int((no_data_points/100)*90)
    
    # mean and standard deviation of the back end of the data
    mean90 = np.mean(y[ind90:])
    stddev90 = np.abs(np.std(y[ind90:]))
    
    # find which points are close to the cut off line within a standard deviation
    mask = np.isclose(mean90, y, atol = stddev90)
    
    # find first instance where data comes close to ss
    ind_ss = np.argmax(mask)
    
    # plt.hist(y[ind_ss:])
    # plt.yscale('log')
    
    return y[ind_ss:] # y values in the steady state
    
# raw data in steady state
y_ss = steady_state_single('C:/Users/kenzi/Documents/Masters/Summer/Studies/Tests/RawData/g=0.098,vs=0.0/run1/output.csv', 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Tests/Configs/g=0.098')  

# read in the config file
constants = f.read_config('C:/Users/kenzi/Documents/Masters/Summer/Studies/Tests/Configs/g=0.098')

R = float(constants[4]) # clinostat radius
r = float(constants[5]) # inner clinostat radius

g = float(constants[9])

bm = 1.1309733552923218E-16 # kg

kT = 4.11E-21 # J

H = kT/(bm*g)



# reference point is -R 
h = y_ss + R # changes from x = 0 to x = R as reference point, i.e height is main measurement

# plot the histogram and save output
plot_histogram(h, R = R, r = r, ylog = True, xlog = False, xlabel = 'Height, h (m)', ylabel = 'PDF(h)', output = 'histogram.png')

