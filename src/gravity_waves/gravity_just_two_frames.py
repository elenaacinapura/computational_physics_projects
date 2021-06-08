import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

data = pd.read_csv('parameters.csv', delimiter='\t')
N = data['N'][0]
T = data['T'][0]
dt = data['dt'][0]
L = data['L'][0]

data = pd.read_csv('f.csv', delimiter='\t')
x = data['x'].to_numpy()
f = data['f'].to_numpy()

l = len(x)
x1 = x[0:N]
x2 = x[l-N-1:l]
f1 = f[0:N]
f2 = f[l-N-1:l]

fig = plt.figure("Gravitational wave")
plt.plot(x1, f1, label = 'Start', c = 'k')
plt.plot(x2, f2, label = 'End', c = 'r')
plt.legend()
plt.show()
