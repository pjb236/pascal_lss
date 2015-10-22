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
import pdb

from scipy import sparse
import scipy.sparse.linalg as splinalg

################################################################################
# Import primal function 
# Must include primal time step, tangent step step.tangent, adjoint step 
# step.adjoint, objective function obj, tangent and adjoint forcing tan_force 
# and adj_force
from tunnel import *
# time step size dt is a global variable
################################################################################


class primal(object):
    def __init__(self, s, m, K, u0=None, T0=0, obj=obj, restart=''):
        # inputs:
        # m: number of time steps per segment
        # K: number of time segments
        # obj: objective function, takes state u, parameters s as inputs
        self.s = s
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
        fout.write('Parameter values:' + str(s) + '\n')
        fout.write('Started from: ' + self.restart + '\n')
        fout.write('Spin up time: ' + str(T0) + '\n')
        fout.write('# Time segments: ' + str(self.K) + '\n')
        fout.write('# steps per segment: ' + str(self.m) + '\n')

        fout.close()


    def spin_up(self, u0, T0):
        # solve primal for spin up time, save final step
        for i in range(int(T0/dt)):
            u0 = step(u0,self.s)
        
        print type(u0)
        u0.save('u000000.npy') 

    def primal_horiz(self):
        # solve for primal on time horizon
        # save solution to disk at each check point (maybe save to memory?)
        print('restarting from ', self.restart)
        u = grid.load(self.restart)
 
        # initialize objective function computation
        Jbar = obj(u,self.s) / (self.m * self.K)

        for i in range(self.K):
            for j in range(self.m):
                u = step(u,self.s)
                J = obj(u,self.s)
                print J
                Jbar = Jbar + J / (self.m * self.K)

            u.save('u{0:06d}.npy'.format(i+1))
            print 'segment ' + str(i+1) + '/' + str(self.K) + ' complete'

        print 'Time average: ' + str(Jbar)
        return Jbar

class lss(object):
     # class containing all functions needed for multiple shooting LSS
     # inputs:
     # s: float(s), array of design parameters
     # m: int, time steps per segment
     # K: int, number of time segments
     # obj: objective function 
     # tan_force: tangent forcing function
     # adj_force: adjoint forcing function
     # eps: float, filtering parameter
     # proj: boolean, use projection functions or not

     def __init__(self,s,m,K,obj,tan_force,adj_force,eps=0.0,proj=True):
         self.s = s
         self._obj = obj
         self._tan_force = tan_force
         self._adj_force = adj_force
         self.m, self.K = m, K
         self.n = int( (Lx/dx) * (Ly/dy) * 4) # number of degrees of freedom
         self.eps = float(eps) 
         self.proj = proj

     def grid_dot(self,u,v):
         # dot product on all axes for grid2d objects u and v
         return (u._data * v._data).sum()

     def project_ddt(self,u,v):
         # remove component of v in the ddt(u) direction
         if self.proj:
             du = ddt(u,self.s)
             du2 = self.grid_dot(du,du)
             vu = self.grid_dot(v,du)
             v = v - (vu/du2) * du
         
         return v

     def forward(self, iSeg, vt0, inhomo=False):
         # solve tangent for a time segment

         # load primal
         u = grid.load('u{0:06d}.npy'.format(iSeg))  
         # projection for time dilation
         vt = self.project_ddt(u,vt0.copy())
         grad = 0.0

         # iterate m times
         for i in range(self.m):
             if inhomo:
                  grad += self.grid_dot(self._adj_force(u,self.s),vt) / (self.m*self.K)
                  vt += self._tan_force(u,self.s) 
             vt,u = step.tangent(vt,u,self.s)
         # projection for time dilation
         vt = self.project_ddt(u,vt) 
         return vt, grad

     def backward(self, iSeg, vaT, inhomo=False):
         # solve adjoint for a time segment
        
         # projection for time dilation
         u = grid.load('u{0:06d}.npy'.format(iSeg+1))
         va = self.project_ddt(u,vaT.copy())
         
         step_fx = lambda x : step(x,self.s) 
         history = histstack.HistoryStack(self.m, step_fx)
         history.populate(grid.load('u{0:06d}.npy'.format(iSeg)))
         grad = 0.0

         for i in range(self.m-1,-1,-1):
             u = history.pop()
             va = step.adjoint(va, u, self.s)
              
             if inhomo:
                 va += self._adj_force(u,self.s) / (self.m * self.K)
                 grad += self.grid_dot(self._tan_force(u,self.s),va)
        
         # projection for time dilation
         va = self.project_ddt(u,va)
         return va, grad


# TODO: Tangent Gradient
# - compute gradient(s) from tangent solution
# - DO NOT forget dJds terms!

# TODO: Adjoint Gradient
# - compute gradient(s) from adjoint solution
# - DO NOT forget dJds terms!


     def grid2array(self,v_grid):
         # convert grid object to 1D array
         v_arr = v_grid._data
         v_arr = np.reshape(v_arr, (self.n,))
         return v_arr

     def array2grid(self,v_arr):
         # convert 1D array to grid object
         Nx,Ny = int(Lx / dx), int(Ly / dy)
         assert v_arr.shape[0] == self.n
         v_arr = np.reshape(v_arr, (Nx,Ny,4))
         v_grid = grid.array(v_arr, v_arr.shape[2:])
         return v_grid


     def matvec(self, x, inhomo=False, adjoint=True):
         # matrix vector multiplication Ax + inhomo*b
         # A is the KKT matrix, b is the right hand side
         assert x.shape == (self.n * self.K,)
         w = x.reshape([-1,self.n])
         
         if adjoint:
             inhomo_a = True
             inhomo_t = False
         else:
             inhomo_a = False
             inhomo_t = True

         R_w = np.zeros([self.K, self.n])
         
         v = [] # empty list for tangent initial conditions
         for i in range(self.K):
             # solve adjoint
             wim,g = self.backward(i,self.array2grid(w[i]),inhomo=inhomo_a)
             if i > 0:
                 v.append(wim - self.array2grid(w[i-1]))
             else:
                 v.append(wim)

             print "adjoint for ", i, " complete"
         v.append(self.array2grid(-w[-1]))

         for i in range(self.K):
             # solve tangent
             vip,g = self.forward(i,v[i],inhomo=inhomo_t)
             R_w[i] = self.grid2array(vip - v[i+1]) + self.eps * w[i]
             print "tangent for ", i, " complete"
         return np.ravel(R_w)


# TODO: matvec_pre
# - matvec for preconditioner

# TODO: tests
# - test tangent and adjoint consistency
# - tests from NASA?


# SOLVER CLASS
class lss_solver(object):
     # Class containing solvers, using scipy's sparse linear algebra package
     def __init__(self,lss):
         # inputs:
         # lss: lss object
         # ...
         # TODO: all solver parameters and options go here!
         
         self.m, self.K, self.n = lss.m, lss.K, lss.n
         self.matvec = lss.matvec

     def callback(self,rk):
         # TODO: fix callback!
         print sum(x*x) ** 0.5 


     def tangent(self):
         # tangent LSS
         return 0.0


     def adjoint(self):
         # adjoint LSS
         
         # compute right hand side
         x = np.zeros(self.K * self.n)
         b = -self.matvec(x)

         # build linear operator
         oper = splinalg.LinearOperator((x.size,x.size), matvec=self.matvec, dtype=float)

         # set up call back
         callback = lambda x : self.callback(x)
         callback(x)
         

         # solve system
         w, info = splinalg.gmres(oper, b, x0=x, maxiter = 20, callback = callback)
         
         print w.shape, np.linalg.norm(w)

         # compute gradient(s)

# TODO: restart capability: restart solution given current adjoints, 

# TODO: precon
# - preconditioner: fixed number of krylov iterations with some sort of simplified matrix?

# TODO: solver
# - solve KKT system -> GMRES! -> careful, need parallel implementation, note that minres restart capability is bad for python


# TODO: Post Process Class: generate plots of primal and adjoint fields, norms, etc from solution.
