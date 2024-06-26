# -*- coding: utf-8 -*-
"""
Visualise output.csv
"""

import numpy as np
import pandas as pd
import Functions as f 
import matplotlib.pyplot as plt

# style sheet

plt.rcParams.update({
    'text.usetex': False,
    'font.family': 'roman',
})

# #handle graph formatting and style
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed

def plot(config_path, data, output_path, save = False):


    # get constants related to the data file
    constants = f.read_config(config_path)
    
    R = float(constants[4])
    r = float(constants[5])
    H = float(constants[6])
    
    #print(R, r, H)
    
    # read in output from bacstroke
    df = np.array(pd.read_csv(data, sep=',', header = None))
    
    x_coords = df[:, 0]
    y_coords = df[:, 1]
    z_coords = df[:, 2]
    time = df[:, 3]
    
    no_plot = 3
    
    fig, ax = plt.subplots(1, 3, figsize = (18, 6.5), sharey = True) # 18 6.5 gives square look for the panels
    
    # first panel y v z positional coordinates
    ax[0].scatter(z_coords[::no_plot], y_coords[::no_plot])
    ax[0].set_xlabel('z')
    ax[0].set_ylabel('y')
    
    # second panel y v x coordinates
    ax[1].scatter(x_coords[::no_plot], y_coords[::no_plot], zorder = 100)
    ax[1].set_xlabel('x')
    
    # outer circle patch
    cir = plt.Circle((0, 0), R, facecolor='#c7c7c7', alpha=1, linewidth=3, linestyle='--', edgecolor='black')#color='darkorange',fill=False)
    ax[1].add_patch(cir)
    
    # inner circle patch
    cir2 = plt.Circle((0, 0), r, facecolor='white', alpha=1, linewidth=3, linestyle='--', edgecolor='black')#color='darkorange',fill=False)
    ax[1].add_patch(cir2)
    
    # third panel y v time
    ax[2].scatter(time[::no_plot], y_coords[::no_plot])
    ax[2].set_xlabel('t  (s)')
    
    if save == True:
        plt.savefig(output_path)

#plot('test_config.txt', 'output.csv', 'PositionPanel.png')