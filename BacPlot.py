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


# get constants related to the data file
constants = f.read_config('test_config.txt')

R = float(constants[4])
r = float(constants[5])
H = float(constants[6])

print(R, r, H)

# read in output from bacstroke
df = np.array(pd.read_csv('output.csv', sep=',', header = None))

x_coords = df[:, 0]
y_coords = df[:, 1]
z_coords = df[:, 2]
time = df[:, 3]

fig, ax = plt.subplots(1, 3, figsize = (18, 6), sharey = True)

# first panel y v z positional coordinates
ax[0].scatter(z_coords, y_coords)
ax[0].set_xlabel('z')
ax[0].set_ylabel('y')

# second panel y v x coordinates
ax[1].scatter(x_coords, y_coords, zorder = 100)
ax[1].set_xlabel('x')

# outer circle patch
cir = plt.Circle((0, 0), R, facecolor='#c7c7c7', alpha=1, linewidth=3, linestyle='--', edgecolor='black')#color='darkorange',fill=False)
ax[1].add_patch(cir)

# inner circle patch
cir2 = plt.Circle((0, 0), r, facecolor='white', alpha=1, linewidth=3, linestyle='--', edgecolor='black')#color='darkorange',fill=False)
ax[1].add_patch(cir2)

# third panel y v time
ax[2].scatter(time, y_coords)
ax[2].set_xlabel('t  (s)')

plt.savefig('PositionPanel.png')
