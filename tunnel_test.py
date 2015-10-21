################################################################################
#                                                                              #
#       tunnel_test.py copyright 2015 Patrick Blonigan (pjblonigan@gmail.com)  #
#                                                                              #
################################################################################

import matplotlib
matplotlib.use('agg')
matplotlib.interactive(False)

import pdb
import sys
import time
import argparse
from pylab import *

from tunnel import *


parser = argparse.ArgumentParser()
parser.add_argument('--restart', type=str, default='')
args = parser.parse_args()

if len(args.restart) == 0:
    w = grid.zeros([4]) + w0
    w[:] *= 1 + 0.01 * (grid.random() - 0.5)
    # tangent initial condition
    v0 = grid.random([4])
    # adjoint terminal condition
    wh1 = grid.random([4])

else:
    print('restarting from ', args.restart)
    w = load(args.restart)
    assert w.shape == (Nx, Ny, 4)

s = 0.1 * ones(2)



df = step.tan_force(w,s)

print df._data.shape

'''
# test discrete consitency of tangent and adjoint
t0 = time.time()
v1,w1 = step.tangent(v0,w,s)
t1 = time.time()
wh0 = step.adjoint(wh1,w,s)
t2 = time.time()


a0 = (v0._data * wh0._data).sum()
a1 = (v1._data * wh1._data).sum()

print "tangent:", t1-t0
print "adjoint:", t2-t1

print a0, a1, (a0 - a1)/a0
'''
################################################################################
################################################################################
################################################################################
