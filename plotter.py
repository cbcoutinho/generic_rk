#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# filename = 'data.out'
filename = 'raw.out'

fig = plt.figure()

df = pd.read_csv(filename,
                 delim_whitespace=True,
                 names=['t', 'y1', 'y2'])
                #  names=['t', 'x1', 'x2', 'x3'])
                #  names=['t', 'x', 'y', 'z'])

df.sort_values(by='t', inplace=True)

t   = df.t.values
x1  = df.y1.values
x2  = df.y2.values
# x3  = df.x3.values

ax = fig.add_subplot(111)
ax.plot(t, x1, marker='o', linestyle='-', label='x1')
ax.plot(t, x2, marker='o', linestyle='-', label='x2')
# ax.plot(t, x3, marker='o', linestyle='None', label='x3')

# tt = np.linspace(t[0], t[-1], 100)
# ax.plot(tt, np.cos(tt), linestyle='-', label='cos(t)')
# ax.plot(tt, np.sin(tt), linestyle='-', label='sin(t)')
# ax.plot(tt, -np.sin(tt), linestyle='-', label='-cos(t)')

# ax = fig.add_subplot(111, projection='3d')
# ax.plot(x, y, z, label='lorenz attractor')

plt.savefig('plot.png')
plt.show()
