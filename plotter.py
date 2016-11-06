#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

filename = 'data.out'

df = pd.read_csv(filename,
                 names=['t', 'x(t)', 'x`(t)'],
                 delim_whitespace=True)

df.sort_values(by='t', inplace=True)

df.plot('t', ['x(t)', 'x`(t)'], marker='o')

plt.show()
