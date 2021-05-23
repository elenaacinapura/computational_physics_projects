import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

data = pd.read_csv('eigenfunction.csv', delimiter='\t')
x = data['x'].to_numpy()
psi = data['psi'].to_numpy()

dim = len(x)
dx = x[1]-x[0]
print(dim)

for i in range(dim):
    x = np.append(x, x[dim - 1 + i]+dx)
    psi = np.append(psi, 0.0)

plt.figure('Ground state')
plt.title('Ground state')
plt.plot(x, psi, c='firebrick', label=r'$\psi_0$')
plt.xlabel(r'$x$')
plt.legend()
plt.show()