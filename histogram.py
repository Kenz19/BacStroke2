# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:12:42 2024

Plot histogram of y positions from BacStroke output. Raw data only.

@author: kenzi
"""

import numpy as np
import pandas as pd
import Functions as f 
import matplotlib.pyplot as plt

plt.rcParams.update({
    'text.usetex': False,
    'font.family': 'roman',
})

# #handle graph formatting and style
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed

def histogram(data, no_bins = 10, ylog = True, xlog = True, xlabel = 'Height, h (m)', ylabel = 'PDF(h)', output = 'histogram.png'):
    '''
    Display histogram of selected data set
    
    :param data: array
            data of which histogram should be created from
    :param no_bins: float
            number of bins that data should be split into.
    :param ylog: boolean
            should the yscale be log?
    :param xlog: boolean
            should the xscale be log?
    :param xlabel: string
            x axis label
    :param ylabel: string
            y axis label
    :param output: string
            where the histrogram should be saved and output to
    '''
    
    # g = 0.0

    # bm = 1.1309733552923218e-16 # kg

    # kT = 4.11e-21 # J

    # H = kT/(bm*g)
    
    hist, bin_edges = np.histogram(data, bins = no_bins, density = True)  
    
    bin_width = np.abs(bin_edges[0] - bin_edges[1])
    
    bin_centres = bin_edges[:-1] + 0.5*bin_width
    
    # area of each bin (need to take away inner circle)
    #areas = Ah(bin_edges, R, r)
    
    return hist, bin_centres


def plot_histogram(hist, bin_centres, output, ylog = True, xlog = False, xlabel = 'Height, h (m)', ylabel = 'PDF(h)'):
    '''
    Display histogram of selected data set
    
    :param hist: array
            PDF values of data set
    :param bin_centres: float
            centres of bins corresponding to PDF values
    :param ylog: boolean
            should the yscale be log?
    :param xlog: boolean
            should the xscale be log?
    :param xlabel: string
            x axis label
    :param ylabel: string
            y axis label
    :param output: string
            where the histrogram should be saved and output to

    '''
    
    # plot the histogram as a scatter plot
    plt.scatter(bin_centres, hist)
    
    if ylog == True:
        plt.yscale('log')
        
    if xlog == True:
        plt.xscale('log')
        
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    plt.savefig(output, dpi = 300)
    

def steady_state_single_Y(data_path, config_path):
    '''
    Intake an output data set from BacStroke and its corresponding config file,
    the y data is isolated and steady state y positions are output.

    Parameters
    ----------
    data_path : string
        Absolute path to a BacStroke output file.
    config_path : string
        Absolute path to the corresponding BacStroke configuration file.

    Returns
    -------
    y values of the BacStroke output file in the steady state.

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
    
    return y[ind_ss:] # y values in the steady state
    

# def Ah(bin_edges, R, r):
#     '''
#     get area between bin edges 
#     '''
#     # area between h = 0 and h = h' - are negative???
    
#     print(bin_edges, R)
#     areas = ((R**2)*np.arcsin((bin_edges - R)/R)) + ((R*(bin_edges - R))*np.sqrt(1 - (bin_edges - R)/R)**2)
#     print(areas)
    
#     #print(((r**2)*np.arcsin((bin_edges - r)/r)) + ((r*(bin_edges - r))*np.sqrt(1 - (bin_edges - r)/r)**2))
    
#     # area encompassed by each bin (A(h1) - A(h2))
#     bin_areas = np.diff(areas)
#     print(bin_areas)
    
#     return bin_areas
    
# # read in output file
# df = np.array(pd.read_csv('output.csv', sep = ',', header = None))
# #print(df)

# # y positions
# y = df[:, 1]

# # read in the config file
# constants = f.read_config('test_config.txt')

# R = float(constants[4]) # clinostat radius
# r = float(constants[5]) # inner clinostat radius

# # reference point is -R 
# h = y + R # changes from x = 0 to x = R as reference point, i.e height is main measurement

# # plot the histogram and save output
# plot_histogram(h, R = R, r = r, ylog = True, xlog = False, xlabel = 'Height, h (m)', ylabel = 'PDF(h)', output = 'histogram.png')


