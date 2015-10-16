from numpy import *
import matplotlib.pyplot as plt


data = loadtxt("force_hist.txt")

drag = data[:,0]
lift = data[:,1]

iiter = 5*arange(data.shape[0])

plt.subplot(1,2,1)
plt.plot(iiter,lift)
plt.xlabel('time step')
plt.ylabel('lift')

plt.subplot(1,2,2)
plt.plot(iiter,drag)
plt.xlabel('time step')
plt.ylabel('drag')


plt.show()
