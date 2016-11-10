#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data_file = 'data.out'
raw_file = 'raw.out'

fig = plt.figure()

df = pd.read_csv(data_file,
                 delim_whitespace=True,
                 names=['t', 'y1', 'y2'])
                #  names=['t', 'x', 'y', 'z'])
                #  names=['t', 'x1', 'x2', 'x3'])

df_raw = pd.read_csv(raw_file,
                     delim_whitespace=True,
                     names=['t', 'dt', 'y1', 'y2'])

df.sort_values(by='t', inplace=True)
df_raw.sort_values(by='t', inplace=True)

t   = df_raw.t.values
y1  = df_raw.y1.values
y2  = df_raw.y2.values
# y3  = df.y3.values

ax = fig.add_subplot(211)
ax.plot(t, y1, marker='o', linestyle='-', label='y1')
ax.plot(t, y2, marker='o', linestyle='-', label='y2')
# ax.plot(t, x3, marker='o', linestyle='None', label='x3')

ax.set_xlabel('t')
ax.set_ylabel('y')
ax.legend()

# tt = np.linspace(t[0], t[-1], 100)
# ax.plot(tt, np.cos(tt), linestyle='-', label='cos(t)')
# ax.plot(tt, np.sin(tt), linestyle='-', label='sin(t)')
# ax.plot(tt, -np.sin(tt), linestyle='-', label='-cos(t)')

# x  = df.x.values
# y  = df.y.values
# z  = df.z.values
#
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(x, y, z, label='lorenz attractor')

ax = fig.add_subplot(212)
ax.plot(df_raw.t.values, df_raw.dt.values,
        marker='+',
        linestyle='None',
        label='Acceptable dt')

ax.set_xlabel('t')
ax.set_ylabel('dt')
ax.legend()

plt.savefig('plot.png')
plt.show()
