#!/usr/bin/env python2
# -*- coding: utf-8 -*-
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

def coarse_graining_step_3D(b,a,D='infinity'):
    
    #z direction
    A = ncon([a,a],[[-1,-3,-5,-7,1,-10],[-2,-4,-6,-8,-9,1]])
    Uy, s, V = tensorsvd(A,[6,7],[0,1,2,3,4,5,8,9],D)
    Ux, s, V = tensorsvd(A,[0,1],[2,3,4,5,6,7,8,9],D)
    A = dncon([Ux,Uy,A,Uy,Ux],[[1,2,-1],[7,8,-4],[1,2,3,4,5,6,7,8,-5,-6],[3,4,-2],[5,6,-3]])
    
    B = dncon([b,a],[[-1,-3,-5,-7,1,-10],[-2,-4,-6,-8,-9,1]])    
    B = dncon([Ux,Uy,B,Uy,Ux],[[1,2,-1],[7,8,-4],[1,2,3,4,5,6,7,8,-5,-6],[3,4,-2],[5,6,-3]])
    
    #y direction
    AA = dncon([A,A],[[-1,-2,1,-6,-7,-8],[1,-3,-4,-5,-9,-10]])
    Uz, s, V = tensorsvd(AA,[6,8],[0,1,2,3,4,5,7,9],D)
    Ux, s, V = tensorsvd(AA,[1,2],[0,3,4,5,6,7,8,9],D) 
    AA = dncon([Uz,Ux,AA,Ux,Uz],[[5,7,-5],[1,2,-2],[-1,1,2,-3,3,4,5,6,7,8],[4,3,-4],[6,8,-6]])  
    
    BA = dncon([B,A],[[-1,-2,1,-6,-7,-8],[1,-3,-4,-5,-9,-10]])
    BA = dncon([Uz,Ux,BA,Ux,Uz],[[5,7,-5],[1,2,-2],[-1,1,2,-3,3,4,5,6,7,8],[4,3,-4],[6,8,-6]])
    
    #x direction
    
    AAAA = dncon([AA,AA],[[-2,-3,-4,1,-7,-8],[-1,1,-5,-6,-9,-10]])
    Uz, s, V = tensorsvd(AAAA,[6,8],[0,1,2,3,4,5,7,9],D)
    Uy, s, V = tensorsvd(AAAA,[3,4],[0,1,2,5,6,7,8,9],D)
    AAAA = dncon([Uz,Uy,AAAA,Uy,Uz],[[5,7,-5],[3,4,-3],[1,2,-2,3,4,-4,5,6,7,8],[2,1,-1],[6,8,-6]])
                        
    BAAA = dncon([BA,AA],[[-2,-3,-4,1,-7,-8],[-1,1,-5,-6,-9,-10]])
    BAAA = dncon([Uz,Uy,BAAA,Uy,Uz],[[5,7,-5],[3,4,-3],[1,2,-2,3,4,-4,5,6,7,8],[2,1,-1],[6,8,-6]])
    
    maxAAAA = np.max(AAAA)
    
    BAAA = BAAA/maxAAAA #divides over largest value in the tensor
    AAAA = AAAA/maxAAAA
    
        
    return BAAA, AAAA, s, maxAAAA

