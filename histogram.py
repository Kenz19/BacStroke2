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

def plot_histogram(data, no_bins = 10, ylog = True, xlog = True, xlabel = 'Height, h (m)', ylabel = 'PDF(h)', output = 'histogram.png'):
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
    
    hist, bin_edges = np.histogram(h, bins = no_bins, density = True)
    
    bin_width = np.abs(bin_edges[0] - bin_edges[1])
    
    bin_centres = bin_edges[:-1] + 0.5*bin_width
    
    plt.scatter(bin_centres, hist)
    
    if ylog == True:
        plt.yscale('log')
        
    if xlog == True:
        plt.xscale('log')
        
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    plt.savefig(output, dpi = 300)

# read in output file
df = np.array(pd.read_csv('output.csv', sep = ',', header = None))
#print(df)

# y positions
y = df[:, 1]

# read in the config file
constants = f.read_config('test_config.txt')

R = float(constants[4]) # clinostat radius

# reference point is -R 
h = y + R # changes from x = 0 to x = R as reference point, i.e height is main measurement

# plot the histogram and save output
plot_histogram(h, ylog = True, xlog = True, xlabel = 'Height, h (m)', ylabel = 'PDF(h)', output = 'histogram.png')


