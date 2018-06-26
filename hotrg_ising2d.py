#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import scipy as sc
import scipy.linalg as scl
from tensorfunc import tensorsvd
import matplotlib.pyplot as plt
from dncon import dncon as ncon

#to compute Ising 2D tensor partition function
def isingtpf2d(beta,h):
    H_local = np.array([[-1.-h/2,1.],[1.,-1.+h/2]]) #h is intended to be a small magnetic field to break symmetry
    Q = np.exp(-beta*H_local)
    
    delta = np.zeros((2,2,2,2))
    delta[0,0,0,0] = 1.
    delta[1,1,1,1] = 1.

    Qsr = scl.sqrtm(Q)

    a = ncon([delta,Qsr,Qsr,Qsr,Qsr],[[1,2,3,4],[-1,1],[-2,2],[3,-3],[4,-4]])
    
    return a

#To compute the tensor that is used to measure one-site magnetization
def isings2d(beta,h):
    H_local = np.array([[-1.-h/2,1.],[1.,-1.+h/2]])
    Q = np.exp(-beta*H_local)
    
    g = np.zeros((2,2,2,2))
    g[0,0,0,0] = 1.
    g[1,1,1,1] = -1.

    Qsr = scl.sqrtm(Q)

    b = ncon([g,Qsr,Qsr,Qsr,Qsr],[[1,2,3,4],[-1,1],[-2,2],[3,-3],[4,-4]])
    
    return b

#coarse-grain a sub-network with the form
# b - a
# a - a
def coarse_graining_step(b,a,D='infinity'):
    
    A = ncon([a,a],[[-2,-3,-4,1],[-1,1,-5,-6]])
    U, s, V = tensorsvd(A,[0,1],[2,3,4,5],D) 
    A = ncon([U,A,U],[[1,2,-1],[1,2,-2,3,4,-4],[4,3,-3]])
    
    B = ncon([b,a],[[-2,-3,-4,1],[-1,1,-5,-6]])    
    B = ncon([U,B,U],[[1,2,-1],[1,2,-2,3,4,-4],[4,3,-3]])
    
    AA = ncon([A,A],[[-1,-2,1,-6],[1,-3,-4,-5]])
    U, s, V = tensorsvd(AA,[1,2],[0,3,4,5],D)  
    AA = ncon([U,AA,U],[[1,2,-2],[-1,1,2,-3,4,3],[3,4,-4]])  
    
    BA = ncon([B,A],[[-1,-2,1,-6],[1,-3,-4,-5]])
    BA = ncon([U,BA,U],[[1,2,-2],[-1,1,2,-3,4,3],[3,4,-4]])  
    
    maxAA = np.max(AA)
    
    BA = BA/maxAA #divides over largest value in the tensor
    AA = AA/maxAA
    
        
    return BA, AA, s, maxAA

