import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', size=12)

ANIMATION = 1
START_END = 0
T_PLOT = 0
EVOLUTION = 1

data = pd.read_csv('parameters.csv', delimiter='\t')
N = data['N'][0]
T = data['T'][0]
L = data['L'][0]
skip = int(10)

#========================== EVOLUTION ==========================
if (EVOLUTION):
    E = data['E'][0]

    data = pd.read_csv('schrodinger.csv', delimiter='\t')
    x = data['x'].to_numpy()
    f = data['f'].to_numpy()


    V = [0 if abs(x[i] - L/2) > 1.0 else 1.0 for i in range(N)]

    def animate (i):
        xframe = x[i*skip*N: (i*skip+1)*N]
        fframe = f[i*skip*N: (i*skip+1)*N]
        ax1.clear()
        ax1.set_title('Quantum wave packet with E = {}'.format(E))
        ax1.set_xlim(0, L)
        ax1.set_ylim(-1.5, 1.5)
        ax1.set_xlabel('x')
        ax1.set_ylabel('f(x)')
        ax1.plot(xframe, fframe, c='r', label='Wave packet')
        ax1.plot(x[0:N], V, c='k', label='Potential')
        ax1.legend()
        return ax1

    if (ANIMATION):
        fig = plt.figure('Quantum particle')
        ax1 = fig.add_subplot(1,1,1)    
        a = animation.FuncAnimation(fig, animate, frames=np.arange(int(T/skip)), interval=0.001)


    l = len(x)
    x1 = x[0:N]
    x2 = x[l-N-1:l]
    f1 = f[0:N]
    f2 = f[l-N-1:l]

    if (START_END):
        plt.figure('Start and end')
        plt.title('Quantum wave packet with E = {}'.format(E))
        plt.plot(x[0:N], V, c='k', label='Potential')
        plt.plot(x1, f1, label = 'Start', c = 'royalblue')
        plt.plot(x2, f2, label = 'End', c = 'firebrick')
        plt.legend()


#=================================== T(E) =======================================
if (T_PLOT):
    data = pd.read_csv('T.csv', delimiter='\t')
    E = data['E'].to_numpy()
    T = data['T'].to_numpy()

    plt.figure("Quantum wave packet")
    plt.title("Transmission coefficient")
    plt.plot(E, T, c='k')
    plt.xlabel('E')
    plt.ylabel('T', rotation=0)

plt.show()

print(" ")