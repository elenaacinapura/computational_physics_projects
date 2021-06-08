import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=18)

data = pd.read_csv('sir.csv', delimiter='\t')
t = data['t'].to_numpy()
S = data['S'].to_numpy()
I = data['I'].to_numpy()
R = data['R'].to_numpy()

plt.figure('Sir model')
plt.title('Modello SIR')
plt.plot(t, S, label='s', c='royalblue')
plt.plot(t, I, label='i', c='firebrick')
plt.plot(t, R, label='r', c='yellowgreen')
plt.xlabel('t')
plt.legend()
plt.show()

print(" ")