import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

data = pd.read_csv('eigenfunction_diag.csv', delimiter='\t')
x = data.iloc[:, 0]
psi = data.iloc[:, 1]

dim = len(x)
dx = x[1]-x[0]


for i in range(dim):
    x = np.append(x, x[dim - 1 + i]+dx)
    psi = np.append(psi, 0.0)

plt.figure('Ground state diag')
plt.title('Ground state with diagonalization')
plt.plot(x, psi, c='firebrick', label=r'$\psi_0$')
plt.xlabel(r'$x$')
plt.legend()
plt.show()

print(' ')