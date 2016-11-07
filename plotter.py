#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

filename = 'data.out'

df = pd.read_csv(filename,
                 names=['t', 'x', 'y', 'z'],
                 delim_whitespace=True)

# df.sort_values(by='t', inplace=True)
# df.plot('t', ['x(t)', 'y(t)', 'z(t)'], linestyle='-', marker='o')

x = df.x.values
y = df.y.values
z = df.z.values

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, label='lorenz attractor')

plt.show()
