################################################################################
#                                                                              #
#   lsspde.py copyright 2015 Patrick Blonigan (pjblonigan@gmail.com)           #
#                                                                              #
################################################################################

# 

import sys
import os
import histstack
import time
################################################################################
# Import primal function 
# Must include primal time step, tangent step step.tangent, adjoint step 
# step.adjoint, objective function obj, tangent and adjoint forcing tan_force 
# and adj_force
from tunnel import *
# time step size dt is a global variable
################################################################################


class primal(object):
    def __init__(self, m, K, u0=None, T0=0, obj=obj, restart=''):
        # inputs:
        # m: number of time steps per segment
        # K: number of time segments
        # obj: objective function, takes state u, parameters s as inputs

        self.m, self.K = int(m), int(K)
        self._obj = obj
        self.restart = restart


        if len(self.restart) == 0:
            if u0 is None:
                u0 = grid.zeros([4]) + w0
                u0[:] *= 1 + 0.01 * (grid.random() - 0.5)
            
            self.spin_up(u0,T0)
            self.restart = 'u000000.npy'

        self.primal_horiz()

        # print file with information on primal solution
        fout = open('primal_info.txt','w')
        fout.write('Information on current primal solution:\n')
        fout.write('Started from: ' + self.restart + '\n')
        fout.write('Spin up time: ' + str(T0) + '\n')
        fout.write('# Time segments: ' + str(self.K) + '\n')
        fout.write('# steps per segment: ' + str(self.m) + '\n')



    def spin_up(self, u0, T0):
        # solve primal for spin up time, save final step
        for i in range(int(T0/dt)):
            u0 = step(u0)
        
        print type(u0)
        u0.save('u000000.npy') 

    def primal_horiz(self):
        # solve for primal on time horizon
        # save solution to disk at each check point (maybe save to memory?)
        print('restarting from ', self.restart)
        u = grid.load(self.restart)
 
        # initialize objective function computation
        Jbar = obj(u) / (self.m * self.K)

        for i in range(self.K):
            for j in range(self.m):
                u = step(u)
                J = obj(u)
                print J
                Jbar = Jbar + J / (self.m * self.K)

            u.save('u{0:06d}.npy'.format(i+1))
            print 'segment ' + str(i+1) + '/' + str(self.K) + ' complete'

        print 'Time average: ' + str(Jbar)
        return Jbar

class lss(object):
     # class containing all functions needed for multiple shooting LSS
     # inputs:
     # m: time steps per segment
     # K: number of time segments
     # eps: filtering parameter (default 0)
     # obj: objective function 
     # tan_force
     # adj_force


     def __init__(self,m,K,obj,tan_force,adj_force,eps=0.0,proj=True):
         self.eps = float(eps) 
         self._obj = obj
         self._tan_force = tan_force
         self._adj_force = adj_force
         self.m, self.K = m, K
         self.proj = proj
         print self.m, self.K 

     def project_ddt(self,u,v):
         # remove component of v in the ddt(u) direction
         if self.proj:
             du = ddt(u)
             du2 = (du._data ** 2).sum()
             vu = (v._data * du._data).sum()
             v = v - (vu/du2) * du
         
         return v

     def forward(self, iSeg, vt0, inhomo=0.0):
         # solve tangent for a time segment

         # load primal
         u = grid.load('u{0:06d}.npy'.format(iSeg))  
         vt = vt0.copy()

         # iterate m times
         for i in range(self.m):
             vt += inhomo * self._tan_force(u)
             vt,u = step.tangent(vt,u)

         return vt

     def backward(self, iSeg, vaT, inhomo=0.0):
         # solve adjoint for a time segment

         history = histstack.HistoryStack(self.m, step)
         history.populate(grid.load('u{0:06d}.npy'.format(iSeg)))
         va = vaT.copy()
         
         for i in range(self.m-1,-1,-1):
             u = history.pop()
             va = step.adjoint(va, u)
             va += inhomo * self.adj_force(u)

         return va


# TODO: Tangent Gradient
# - compute gradient(s) from tangent solution
# - DO NOT forget dJds terms!

# TODO: Adjoint Gradient
# - compute gradient(s) from adjoint solution
# - DO NOT forget dJds terms!

# TODO: grid to 1D array (and vice-versa)


# TODO: matvec
# - matrix vector opertion to:
# -- generate right hand side
# -- perform operation of KKT matrix

# TODO: matvec_pre
# - matvec for preconditioner

# TODO: tests
# - test tangent and adjoint consistency
# - tests from NASA?


# SOLVER CLASS

# TODO: precon
# - preconditioner: fixed number of krylov iterations with some sort of simplified matrix?

# TODO: solver
# - solve KKT system -> GMRES! -> careful, need parallel implementation, note that minres restart capability is bad for python
