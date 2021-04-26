import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

data = pd.read_csv('universe.csv', delimiter='\t')
t = data['t'].to_numpy()
a = data['a'].to_numpy()

plt.figure('Universe')
plt.plot(t, a)

plt.title('Universal scale factor')
plt.xlabel('t')
plt.ylabel('a(t)')
plt.show()