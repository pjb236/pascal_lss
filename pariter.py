#!/usr/bin/env python
import numpy, math
from numpy import *
from scipy.sparse import linalg as splinalg

__all__ = ['par_cg', 'par_minres']

def par_cg(A, b, x, dot, tol=1e-5, maxiter=None, callback=None, M=None):
    """
    Implementation of conjugate gradient for solving symmetric positive
    definite linear systems, mimics api of scipy.sparse.linalg.cg with
    addition of user supplied inner product. --- Also uses Polak-Ribiere
    formula for beta to improve performance for non fixed preconditioners
    """
    sqrt = math.sqrt
    n = len(b)

    if M is None:
        M = splinalg.LinearOperator( (n,n), matvec=lambda x: x, dtype=float )

    if maxiter is None:
        maxiter = n*10

    r = b - A*x
    z = M*r
    p = z.copy()

    rr0 = dot(r, r)


    for k in xrange(maxiter):

        Ap = A*p
        pAp = dot(p, Ap)
        zr = dot(z, r)


        alpha = zr / pAp
        x += alpha*p
        Ap *= -alpha; r += Ap; # r = r - alpha*Ap


        rr = dot(r,r)

        if callback is not None:
            #callback(x, sqrt(rr) )
            callback(x)

        if rr/rr0 < tol*tol:
            err = (0, k+1)
            break

        z = M*r
        beta = dot(z, Ap)/zr # at this point Ap = r_{k+1} - r_k
        
        p*=beta; p+=z #p = z + beta*p

    else: # if not break
        err = (-1, maxiter)

    return x, err

def par_minres(A, b, x, dot, tol=1e-5, maxiter=None, callback=None):
    """
    Implementation of unpreconditioned MINRES for solving symmetric
    linear systems mimics api of scipy.sparse.linalg.minres with addition
    of user supplied inner product.
    """
    n = len(b)

    sqrt = math.sqrt # local lookup

    if maxiter is None:
        maxiter = n*10

    vi = b - A*x
    vim = zeros_like(vi)

    betai = sqrt( dot(vi,vi) )

    r0 = betai

    eta = betai
    gammai = gammaim = 1; sigmai = sigmaim = 0;

    wi = zeros_like(vi); wim = zeros_like(vi)

    for k in xrange(maxiter):

        #Lancczos
        vi /= betai;

        Avi = A * vi;
        alphai = dot(vi, Avi)
        vip = Avi - alphai*vi - betai*vim
        betaip = sqrt( dot(vip,vip) )

        # QR Factor
        delta = gammai*alphai - gammaim*sigmai*betai
        p1 = sqrt( delta*delta + betaip*betaip )
        p2 = sigmai*alphai + gammaim*gammai*betai
        p3 = sigmaim*betai
        # Givens Rotation
        gammaip = delta/p1; sigmaip = betaip/p1
        #Update soln
        wip = (vi - p3*wim - p2*wi )/p1
        x += gammaip*eta*wip
        eta *= -sigmaip

        # roll over variables
        betai = betaip;
        sigmaim = sigmai; sigmai = sigmaip;
        gammaim = gammai; gammai = gammaip;
        wim = wi; wi = wip;
        vim = vi; vi = vip;

        if callback is not None:
            #callback(x, abs(eta) )
            callback(x)

        if abs(eta)/r0 < tol:
            return x, (0, k+1)

    else: # if not break
        return x, (1, maxiter)
