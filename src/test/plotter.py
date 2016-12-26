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
                 names=['t', 'y1', 'y2', 'y3'])
                #  names=['t', 'y1', 'y2'])
                #  names=['t', 'x', 'y', 'z'])
                #  names=['t', 'x1', 'x2', 'x3'])

df_raw = pd.read_csv(raw_file,
                     delim_whitespace=True,
                     names=['t', 'dt', 'y1', 'y2', 'y3'])

df.sort_values(by='t', inplace=True)
df_raw.sort_values(by='t', inplace=True)

t   = df.t.values
y1  = df.y1.values
y2  = df.y2.values
y3  = df.y3.values

ax1 = fig.add_subplot(211)
ax1.plot(t, y1, marker='o', linestyle='None', label='y1', color='b')
ax1.plot(t, y2, marker='o', linestyle='None', label='y2', color='g')
ax1.plot(t, y3, marker='o', linestyle='None', label='y3', color='r')

t   = df_raw.t.values
y1  = df_raw.y1.values
y2  = df_raw.y2.values
y3  = df_raw.y3.values

ax1.plot(t, y1, marker='None', linestyle='-', color='b')
ax1.plot(t, y2, marker='None', linestyle='-', color='g')
ax1.plot(t, y3, marker='None', linestyle='-', color='r')



ax1.set_xlabel('t')
ax1.set_ylabel('y')
ax1.legend()

# tt = np.linspace(t[0], t[-1], 100)
# ax1.plot(tt, np.cos(tt), linestyle='-', label='cos(t)')
# ax1.plot(tt, np.sin(tt), linestyle='-', label='sin(t)')
# ax1.plot(tt, -np.sin(tt), linestyle='-', label='-cos(t)')

# x  = df.x.values
# y  = df.y.values
# z  = df.z.values
#
# ax1 = fig.add_subplot(111, projection='3d')
# ax1.plot(x, y, z, label='lorenz attractor')

ax2 = fig.add_subplot(212, sharex=ax1)
ax2.semilogy(df_raw.t.values, df_raw.dt.values,
             marker='+',
             linestyle='None',
             label='Acceptable dt')

ax2.set_xlabel('t')
ax2.set_ylabel('dt')
ax2.legend()

plt.savefig('plot.png')
plt.show()
