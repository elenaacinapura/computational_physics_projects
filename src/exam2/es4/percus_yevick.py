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
plt.show()

print(' ')