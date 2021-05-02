import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

N = 256

fig = plt.figure('Gravity waves')
ax1 = fig.add_subplot(1,1,1)

data = pd.read_csv('f.csv', delimiter='\t')
x = data['x'].to_numpy()
f = data['f'].to_numpy()


def animate (i):
    xframe = x[i*N: (i+1)*N]
    fframe = f[i*N: (i+1)*N]
    ax1.clear()
    ax1.set_xlim(0, 30)
    ax1.set_ylim(-1.0, 1.0)
    ax1.plot(xframe, fframe)
    return ax1

a = animation.FuncAnimation(fig, animate, frames=np.arange(100), interval=100)
plt.show()