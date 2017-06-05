#!/Users/Lenardo/anaconda/envs/py36/bin/python
import numpy as np
from matplotlib import pyplot as plt


er = np.genfromtxt('sim_er_4-10.txt')
nr = np.genfromtxt('sim_nr_4-10.txt')
print(er)
print(nr)

plt.plot(er[:,0],er[:,1],'--',linewidth=3,color=(0.,0.,0.9))
plt.plot(er[:,0],er[:,2],'-',linewidth=3,color=(0.,0.,0.9))
plt.plot(er[:,0],er[:,3],'--',linewidth=3,color=(0.,0.,0.9))
plt.plot(nr[:,0],nr[:,1],'--',linewidth=3,color=(0.9,0.,0.))
plt.plot(nr[:,0],nr[:,2],'-',linewidth=3,color=(0.9,0.,0.))
plt.plot(nr[:,0],nr[:,3],'--',linewidth=3,color=(0.9,0.,0.))


plt.xlabel('S1 size')
plt.ylabel('Prompt Fraction')

plt.show()

