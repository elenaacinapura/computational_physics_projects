import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

data = pd.read_csv('structure.csv', delimiter='\t')
r = data['r'].to_numpy()
S = data['S'].to_numpy()

plt.figure('Structure Factor')
plt.plot(r, S, '.', c='r')
plt.title('Structure factor')
plt.xlabel('q')
plt.ylabel('S', rotation=0)
plt.show()