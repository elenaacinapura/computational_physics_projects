import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

#============= DATA =============
data = pd.read_csv('virial.csv', delimiter='\t')
T = data['T']
B = data['B']


#============= PLOT =============
plt.figure('Second virial coefficient')
plt.plot(T, B, '.', c='r')
plt.title('Second virial coefficient')
plt.xlabel('T')
plt.ylabel('B', rotation=0)
plt.show()

print(' ')
