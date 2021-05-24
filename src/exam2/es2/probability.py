import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

data = pd.read_csv('probability.csv', delimiter='\t')
t = data.iloc[:, 0]
prob = data.iloc[:, 1]

plt.figure('Probability')
plt.title(r'$p(t)$')
plt.plot(t, prob, c = 'royalblue')
plt.xlabel(r'$t$')
plt.ylabel(r'$p(t)$')

data = pd.read_csv('spectrum.csv', delimiter='\t')
w = data.iloc[:, 0]
f = data.iloc[:, 1]

plt.figure('Spectrum')
plt.title(r'Spectrum')
plt.plot(w, f, c = 'royalblue')
plt.xlim(-1, 1)
plt.xlabel(r'$\omega$')
plt.ylabel(r'Spectrum')
plt.show()
