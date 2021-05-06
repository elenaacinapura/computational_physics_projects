import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

data = pd.read_csv('parameters.csv', delimiter='\t')
N = data['N'][0]
T = data['T'][0]
dt = data['dt'][0]
L = data['L'][0]

fig = plt.figure('Quantum particle')
ax1 = fig.add_subplot(1,1,1)

data = pd.read_csv('schrodinger.csv', delimiter='\t')
x = data['x'].to_numpy()
f = data['f'].to_numpy()

V = [0 if abs(x[i] - L/2) > 1.0 else 1.0 for i in range(len(x))]

def animate (i):
    xframe = x[i*N: (i+1)*N]
    fframe = f[i*N: (i+1)*N]
    ax1.clear()
    ax1.set_xlim(0, L)
    ax1.set_ylim(-1.5, 1.5)
    ax1.set_xlabel('x')
    ax1.set_ylabel('f(x)')
    ax1.plot(xframe, fframe)
    ax1.plot(x, V)
    return ax1

a = animation.FuncAnimation(fig, animate, frames=np.arange(T), interval=1)
plt.show()