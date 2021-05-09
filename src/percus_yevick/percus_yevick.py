import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

data = pd.read_csv('percus.csv', delimiter='\t')
r = data['r'].to_numpy()
g = data['g'].to_numpy()

plt.figure('Percus-Yevick')
plt.plot(r, g, '.', c='r')
plt.title('g(r)')
plt.xlabel('r')
plt.ylabel('g', rotation=0)
plt.show()