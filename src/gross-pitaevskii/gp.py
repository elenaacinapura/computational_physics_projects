import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#============= DATA =============
data = pd.read_csv('solution.csv', delimiter='\t')
r = data.iloc[:, 0]
u = data.iloc[:, 1]


#============= PLOT =============
plt.figure('Gross-Pitaevskii')
plt.plot(r, u, '.', c='r')
plt.title('Gross-Pitaevskii')
plt.xlabel('r')
plt.ylabel('u', rotation=0)
plt.show()

print(' ')
