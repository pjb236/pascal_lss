from numpy import *
import matplotlib.pyplot as plt

data0 = loadtxt("force_hist.txt")
N = 5*data0.shape[0]

data = loadtxt("adj_hist.txt")

aru = data[:,0]
arv = data[:,1]

iiter = arange(aru.shape[0])

iiter = N * ones(iiter.shape) - iiter

plt.semilogy(iiter,aru)
plt.semilogy(iiter,arv)

plt.ylabel("adjoint momentum norm")
plt.xlabel("time step")

plt.legend(('x','y'))
plt.show()
