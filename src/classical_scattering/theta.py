import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)

#============= THETA(b) =============
data = pd.read_csv('theta.csv', delimiter='\t')
b = data['b']
theta = data['theta']
theta_theo = data['theta_theo']

plt.figure('Classical scattering - theta')
plt.plot(b, theta, '.', c='r', label='Verlet')
plt.plot(b, theta_theo, c = 'b', label='Analytically')
plt.title(r'$\theta$ as a function of the impact parameter $b$')
plt.xlabel(r'$b$')
plt.ylabel(r'$\theta$', rotation=0)
plt.legend()

#============= CROSS SECTION =============
data = pd.read_csv('cross_section.csv', delimiter='\t')
theta = data['theta']
cnt = data['cnt']

plt.figure('Classical scattering - cross section')
plt.plot(theta, cnt, c='lightgreen')
plt.title(r'Reduced cross section')
plt.xlabel(r'$\theta$')
plt.ylabel(r'$\bar{\sigma}$', rotation=0)

plt.show()

print(' ')
