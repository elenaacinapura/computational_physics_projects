import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

data = pd.read_csv('eigenfunction.csv', delimiter='\t')
x = data['x'].to_numpy()
psi = data['psi'].to_numpy()

plt.figure('Ground state')
plt.title('Ground state')
plt.plot(x, psi)
plt.xlabel('x')
plt.ylabel('psi')
plt.show()