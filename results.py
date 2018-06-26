#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from hotrg_ising2d import *

h = 10**(-8)
temperatures = np.linspace(0.1,4,20)
results = {}
Ds = [2,10,15]

for D in Ds:
    results[D] = []
    for temp in temperatures:
        beta = 1/temp
        a = isingtpf2d(beta,h)
        b = isings2d(beta,h)
        converg_criteria = 10**(-2)
        
        i = 0
        s_i = np.eye(D)
        while True:
            b, a, s = coarse_graining_step(b,a,D=D)
            if i >= 5:
                delta_s = np.abs(s - s_i)
                if np.max(delta_s) <= converg_criteria:
                    print i
                    break
                else:
                    s_i = s
            i += 1
        '''
        for i in range(20):
            b, a, s = coarse_graining_step(b,a,D=D)
        '''
        ap = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
        bp = ncon([b,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
        
        results[D].append(bp/ap)

for D in Ds:
    plt.plot(temperatures,results[D])
#plt.axvline(x=0.44069,linestyle='--')
plt.show()
