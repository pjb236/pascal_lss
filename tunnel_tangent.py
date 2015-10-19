################################################################################
#                                                                              #
#       tunnel_adjoint.py copyright 2015 Qiqi Wang (qiqi.wang@gmail.com)       #
#                                                                              #
################################################################################

import matplotlib
matplotlib.use('agg')
matplotlib.interactive(False)

import os
import pdb
import sys
import time
import argparse

from pylab import *

from tunnel import *

parser = argparse.ArgumentParser()
parser.add_argument('nStart', type=int)
parser.add_argument('nEnd', type=int)
args = parser.parse_args()

assert args.nStart < args.nEnd
assert os.path.exists('w{0:06d}.npy'.format(args.nStart))
w = grid.load('w{0:06d}.npy'.format(args.nStart))


if os.path.exists('tan{0:06d}.npy'.format(args.nStart)):
    tan = grid.load('tan{0:06d}.npy'.format(args.nStart))
else:
    tan = grid.zeros([4])


fout = open('tan_hist.txt','w')

figure(figsize=(28,10))

for iplot in range(args.nStart, args.nEnd):
    for iprint in range(nPrintsPerPlot):
        for istep in range(nStepPerPrint):
            # tangent forcing
            tan[1:3] += c0 * obstacle * w[1:3] # sensitivity to obstacle weight
            ####################################################################
            tan,w = step.tangent(tan, w)
        tru_mag = norm(tan[1]._data) 
        trv_mag = norm(tan[2]._data)
        print('%f %f' % (tru_mag, trv_mag))
        fout.write('%f %f\n' % (tru_mag, trv_mag))

    tan.save('tan{0:06d}.npy'.format(iplot+1))
    clf()
    subplot(2,1,1); contourf(x._data.T, y._data.T, tan[1]._data.T, 200);
    axis('scaled'); colorbar()
    subplot(2,1,2); contourf(x._data.T, y._data.T, tan[2]._data.T, 200);
    axis('scaled'); colorbar()
    savefig('tan{0:06d}.png'.format(iplot+1))

fout.close()
################################################################################
################################################################################
################################################################################
