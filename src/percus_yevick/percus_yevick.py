import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

#============= PARAMETERS =============
data = pd.read_csv('parameters.csv', delimiter='\t')
rho = data['rho'][0]
T = data['T'][0]

#============= DATA =============
data = pd.read_csv('percus.csv', delimiter='\t')
r = data['r'].to_numpy()
g = data['g'].to_numpy()

#============= PLOT =============
plt.figure('Percus-Yevick')
plt.plot(r, g, '.', c='r')
plt.title('Pair distribution function for rho = {} and T = {}'.format(rho, T))
plt.xlabel('r')
plt.ylabel('g', rotation=0)

#============= TEST =============
data = pd.read_csv('test.csv', delimiter='\t')
f1 = data['F1'].to_numpy()
f2 = data['F2'].to_numpy()

plt.figure('Test')
plt.title('Convolutions with rho = {}, T = {}'. format(rho, T))
plt.plot(r, f1, label='FFT')
plt.plot(r, f2, label='Integral')
plt.legend()
plt.show()
