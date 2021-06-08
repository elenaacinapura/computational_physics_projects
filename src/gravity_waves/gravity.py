import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)

data = pd.read_csv('parameters.csv', delimiter='\t')
N = data['N'][0]
T = data['T'][0]
dt = data['dt'][0]
L = data['L'][0]

fig = plt.figure('Gravity waves')
ax1 = fig.add_subplot(1,1,1)

data = pd.read_csv('f.csv', delimiter='\t')
x = data['x'].to_numpy()
f = data['f'].to_numpy()


def animate (i):
    xframe = x[i*N: (i+1)*N]
    fframe = f[i*N: (i+1)*N]
    ax1.clear()
    ax1.set_title("Dispersion of a wave packet")
    ax1.set_xlim(0, L)
    ax1.set_ylim(-1.5, 1.5)
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$f(x)$', rotation=0)
    ax1.plot(xframe, fframe)
    return ax1

a = animation.FuncAnimation(fig, animate, frames=np.arange(int(T/dt)), interval=10)
plt.show()