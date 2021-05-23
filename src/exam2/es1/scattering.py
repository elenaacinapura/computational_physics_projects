import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', size=16)

data = pd.read_csv('sigma_tot.csv', delimiter='\t')
E = data.iloc[:, 0]
S = data.iloc[:, 1]

plt.figure('scattering')
plt.title('Total cross section')
plt.plot(E, S, c='lightgreen')
plt.xlabel(r'$E$ [meV]')
plt.ylabel(r'$\Sigma _{TOT}~[A^2]$')
plt.show()

print(' ')