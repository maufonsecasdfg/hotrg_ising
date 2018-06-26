#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import scipy as sc
import scipy.linalg as scl
from tensorfunc import tensorsvd, tensoreig
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

def isingtpf3d(beta,h):
    H_local = np.array([[-1.-h/2,1.],[1.,-1.+h/2]]) #h is intended to be a small magnetic field to break symmetry
    Q = np.exp(-beta*H_local)
    
    delta = np.zeros((2,2,2,2,2,2))
    delta[0,0,0,0,0,0] = 1.
    delta[1,1,1,1,1,1] = 1.

    Qsr = scl.sqrtm(Q)

    a = ncon([delta,Qsr,Qsr,Qsr,Qsr,Qsr,Qsr],[[1,2,3,4,5,6],[-1,1],[-2,2],[3,-3],[4,-4],[-5,5],[6,-6]])
    
    return a
   

def isings3d(beta,h):
    H_local = np.array([[-1.-h,1.],[1.,-1.+h]])
    Q = np.exp(-beta*H_local)
    
    g = np.zeros((2,2,2,2,2,2))
    g[0,0,0,0,0,0] = 1.
    g[1,1,1,1,1,1] = -1.

    Qsr = scl.sqrtm(Q)

    b = ncon([g,Qsr,Qsr,Qsr,Qsr,Qsr,Qsr],[[1,2,3,4,5,6],[-1,1],[-2,2],[3,-3],[4,-4],[-5,5],[6,-6]])
    
    return b

def isingtpf(dim,beta,h):
    if dim == 2:
        return isingtpf2d(beta,h)
    elif dim == 3:
        return isingtpf3d(beta,h)
    
def isings(dim,beta,h):
    if dim == 2:
        return isings2d(beta,h)
    elif dim == 3:
        return isings3d(beta,h)
        

#coarse-grain a sub-network with the form
# b - a
# a - a
def coarse_graining_step_2d(a,b,D='infinity'):
    
    A = ncon([a,a],[[-2,-3,-4,1],[-1,1,-5,-6]])
    U, s, V = tensorsvd(A,[0,1],[2,3,4,5],D) 
    A = ncon([U,A,U],[[1,2,-1],[1,2,-2,3,4,-4],[4,3,-3]])
    if b != None:
        B = ncon([b,a],[[-2,-3,-4,1],[-1,1,-5,-6]])    
        B = ncon([U,B,U],[[1,2,-1],[1,2,-2,3,4,-4],[4,3,-3]])
    
    AA = ncon([A,A],[[-1,-2,1,-6],[1,-3,-4,-5]])
    U, s, V = tensorsvd(AA,[1,2],[0,3,4,5],D)  
    AA = ncon([U,AA,U],[[1,2,-2],[-1,1,2,-3,4,3],[3,4,-4]])  
    if b != None:
        BA = ncon([B,A],[[-1,-2,1,-6],[1,-3,-4,-5]])
        BA = ncon([U,BA,U],[[1,2,-2],[-1,1,2,-3,4,3],[3,4,-4]])  
    
    maxAA = np.max(AA)
     
    AA = AA/maxAA #divides over largest value in the tensor
    if b != None:
        BA = BA/maxAA
    
        
    return AA, BA, maxAA

def coarse_graining_step_3d(a,b=None,D='infinity'):
    
    #z direction
    A = ncon([a,a],[[-1,-3,-5,-7,1,-10],[-2,-4,-6,-8,-9,1]])
    Ux, s, V = tensorsvd(A,[2,3],[0,1,4,5,6,7,8,9],D)
    Uy, s, V = tensorsvd(A,[0,1],[2,3,4,5,6,7,8,9],D)
    A = ncon([Ux,Uy,A,Uy,Ux],[[3,4,-2],[1,2,-1],[1,2,3,4,5,6,7,8,-5,-6],[5,6,-3],[7,8,-4]])
    if b != None:
        B = ncon([b,a],[[-1,-3,-5,-7,1,-10],[-2,-4,-6,-8,-9,1]])    
        B = ncon([Ux,Uy,B,Uy,Ux],[[3,4,-2],[1,2,-1],[1,2,3,4,5,6,7,8,-5,-6],[5,6,-3],[7,8,-4]])
    
    #y direction
    AA = ncon([A,A],[[-1,-2,1,-6,-7,-8],[1,-3,-4,-5,-9,-10]])
    Uz, s, V  = tensorsvd(AA,[6,8],[0,1,2,3,4,5,7,9],D)
    Ux, s, V  = tensorsvd(AA,[1,2],[0,3,4,5,6,7,8,9],D) 
    AA = ncon([Uz,Ux,AA,Ux,Uz],[[5,7,-5],[1,2,-2],[-1,1,2,-3,3,4,5,6,7,8],[4,3,-4],[6,8,-6]])  
    if b != None:
        BA = ncon([B,A],[[-1,-2,1,-6,-7,-8],[1,-3,-4,-5,-9,-10]])
        BA = ncon([Uz,Ux,BA,Ux,Uz],[[5,7,-5],[1,2,-2],[-1,1,2,-3,3,4,5,6,7,8],[4,3,-4],[6,8,-6]])
    
    #x direction
    
    AAAA = ncon([AA,AA],[[-2,-3,-4,1,-7,-8],[-1,1,-5,-6,-9,-10]])
    Uz, s, V  = tensorsvd(AAAA,[6,8],[0,1,2,3,4,5,7,9],D)
    Uy, s, V  = tensorsvd(AAAA,[3,4],[0,1,2,5,6,7,8,9],D)
    AAAA = ncon([Uz,Uy,AAAA,Uy,Uz],[[5,7,-5],[3,4,-3],[1,2,-2,3,4,-4,5,6,7,8],[2,1,-1],[6,8,-6]])
    if b != None:
        BAAA = ncon([BA,AA],[[-2,-3,-4,1,-7,-8],[-1,1,-5,-6,-9,-10]])
        BAAA = ncon([Uz,Uy,BAAA,Uy,Uz],[[5,7,-5],[3,4,-3],[1,2,-2,3,4,-4,5,6,7,8],[2,1,-1],[6,8,-6]])
    
    maxAAAA = np.max(AAAA)

    AAAA = AAAA/maxAAAA #divides over largest value in the tensor
    if b != None:
        BAAA = BAAA/maxAAAA
    
        
    return AAAA, BAAA, maxAAAA

def coarse_graining_step(dim,a,b=None,D='infinity'):
    if dim == 2:
        return coarse_graining_step_2d(a,b=b,D=D)
    elif dim == 3:
        return coarse_graining_step_3d(a,b=b,D=D)
    
def final_contraction(dim,a,b=None):
    if dim == 2:
        ap = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
        if b != None:
            bp = ncon([b,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
    elif dim == 3:
        ap = ncon([a,a,a,a,a,a,a,a],[[1,2,4,5,6,3],[4,7,1,8,9,10],[11,12,13,14,3,6],[13,15,11,16,10,9],[20,14,19,12,17,18],[19,16,20,15,22,21],[24,5,23,2,18,17],[23,8,24,7,21,22]])
        if b != None:
            bp = ncon([b,a,a,a,a,a,a,a],[[1,2,4,5,6,3],[4,7,1,8,9,10],[11,12,13,14,3,6],[13,15,11,16,10,9],[20,14,19,12,17,18],[19,16,20,15,22,21],[24,5,23,2,18,17],[23,8,24,7,21,22]])
    return ap, bp
    
