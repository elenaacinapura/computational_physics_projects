import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

data = pd.read_csv('parameters.csv', delimiter='\t')
N = data['N'][0]
L = data['L'][0]

data = pd.read_csv('func.csv', delimiter='\t')
x = data['x'].to_numpy()
f = data['f'].to_numpy()

data = pd.read_csv('derivative.csv', delimiter='\t')
d = data['f'].to_numpy()

plt.figure('Derivative')
plt.ylim(0, 10)
plt.plot(x, f, c = 'r', label='f')
plt.plot(x, d, c='b', label='derivative')
plt.xlabel('x')
plt.legend()
plt.show()