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
from histogram import *

    
# raw data in steady state
y_ss = steady_state_single_Y('C:/Users/kenzi/Documents/Masters/Summer/Studies/Tests/RawData/g=0.098,vs=0.0/run1/output.csv', 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Tests/Configs/g=0.098')  

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
hist, bin_centres = histogram(h, ylog = True, xlog = False, xlabel = 'Height, h (m)', ylabel = 'PDF(h)', output = 'histogram.png')

plt.scatter(bin_centres, hist, label = 'Simulation')#/areas)
plt.plot(bin_centres, np.exp(-(bin_centres/H)), label = 'Exponential')
#plt.scatter(bin_centres, np.zeros_like(bin_centres)-(1/H), label = 'Exponential')
#plt.scatter(bin_centres, np.exp(-bin_centres/H))

plt.xlabel('Height, h (m)')
plt.ylabel('PDF(h)')
plt.yscale('log')
plt.legend()

plt.savefig('histogram.png', dpi = 300)