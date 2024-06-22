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

# read in output file
df = np.array(pd.read_csv('output.csv', sep=',', header = None))
#print(df)

# y positions
y = df[:, 1]

# read in the config file
constants = f.read_config('test_config.txt')

R = float(constants[4]) # clinostat radius

# reference point is -R 
h = y + R # changes from x = 0 to x = R as reference point, i.e height is main measurement

# counts and bin edges, normalised
hist, bin_edges = np.histogram(h, bins = 10, density = True)

# width of bin, all are equal in width
bin_width = np.abs(bin_edges[0] - bin_edges[1])

# bin centres
bin_centres = bin_edges[:-1] + 0.5*bin_width

plt.scatter(bin_centres, hist)
plt.xlabel('Height, h (m)')
plt.ylabel('PDF(h)')

