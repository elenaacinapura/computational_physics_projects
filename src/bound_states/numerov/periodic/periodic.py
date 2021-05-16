import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

data = pd.read_csv("periodic.csv", delimiter='\t')
K = data['K'].to_numpy()
E0 = data['E0'].to_numpy()
E1 = data['E1'].to_numpy()
E2 = data['E2'].to_numpy()

plt.figure('Dispersion relation')
plt.title ('Dispersion relation')
plt.plot(K, E0, '.', label='E_0')
plt.plot(K, E1, '.', label='E_1')
plt.plot(K, E2, '.', label='E_2')
plt.xlabel('K')
plt.ylabel('E', rotation=0)
plt.legend()
plt.show()