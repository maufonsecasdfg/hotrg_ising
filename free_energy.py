#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from hotrg_ising import *

def compute_free_energy_temps(dim,Ds,h,temperatures):

    results = {}
    
    for D in Ds:
        results[D] = []
        for temp in temperatures:
            print 'temp: %s' %(temp)
            beta = 1/temp
            a = isingtpf(dim,beta,h)
            converg_criteria = 10**(-10)
            
            i = 0
            Z = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
            f_i = -temp*(np.log(Z))/(4)
            C = 0
            N = 1
    
            while True:
                a, a, s, maxAA = coarse_graining_step(dim,a,a,D=D)
                C = np.log(maxAA)+4*C
                N *= 4.
                if i >= 10:
                    Z = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
                    f = -temp*(np.log(Z)+4*C)/(4*N)
                    delta_f = np.abs((f - f_i)/f)
                    #print delta_f
                    if delta_f <= converg_criteria:
                        print i
                        break
                    else:
                        f_i = f
                i += 1
            '''
            for i in range(20):
                b, a, s, maxAA = coarse_graining_step(b,a,D=D)
                C *= maxAA
                N += 4
            '''
            Z = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
            f = -temp*(np.log(Z)+4*C)/(4*N)
            
            results[D].append(f)
            
        return results

def compute_free_energy_hs(dim,Ds,hs,temp):

    results = {}
    beta = 1/temp
    for D in Ds:
        results[D] = []
        for h in hs:
            print 'h: %s' %(h)
            a = isingtpf(dim,beta,h)
            converg_criteria = 10**(-10)
            
            i = 0
            Z = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
            f_i = -temp*(np.log(Z))/(4)
            C = 0
            N = 1
    
            while True:
                a, a, s, maxAA = coarse_graining_step(dim,a,a,D=D)
                C = np.log(maxAA)+4*C
                N *= 4.
                if i >= 10:
                    Z = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
                    f = -temp*(np.log(Z)+4*C)/(4*N)
                    delta_f = np.abs((f - f_i)/f)
                    #print delta_f
                    if delta_f <= converg_criteria:
                        print i
                        break
                    else:
                        f_i = f
                i += 1
            '''
            for i in range(20):
                b, a, s, maxAA = coarse_graining_step(b,a,D=D)
                C *= maxAA
                N += 4
            '''
            Z = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
            f = -temp*(np.log(Z)+4*C)/(4*N)
            
            results[D].append(f)
            
        return results
'''
temp = 0.5
hs = np.linspace(-1,1,50)
Ds = [5]

results_f = compute_free_energy_hs(Ds,hs,temp)


magnet_results = compute_magnetization_hs(Ds,hs,temp)

for D in Ds:
    plt.plot(hs,results_f[D],label=str(D))

plt.legend()
plt.show()

for D in Ds:
    plt.plot(hs,magnet_results[D],label=str(D))

plt.legend()
plt.show()
'''
'''
hs = np.linspace(-2,2,50)
temp = 1.0
Ds = [4]

results[D] = []
for h in hs:
    print 'h: %s' %(h)
    beta = 1/temp
    a = isingtpf2d(beta,h)
    converg_criteria = 10**(-10)
            
    i = 0
    Z = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
    f_i = -temp*(np.log(Z))/(4)
    C = 0
    N = 1
    
    while True:
        a, a, s, maxAA = coarse_graining_step(a,a,D=D)
        C = np.log(maxAA)+4*C
        N *= 4.
        if i >= 1:
            Z = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
            f = -temp*(np.log(Z)+4*C)/(4*N)
            delta_f = np.abs((f - f_i)/f)
            print delta_f
            if delta_f <= converg_criteria:
                print i
                break
            else:
                f_i = f
        i += 1

    Z = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
    f = -temp*(np.log(Z)+4*C)/(4*N)
            
    results[D].append(f)

for D in Ds:
    plt.plot(hs,results[D])
plt.axhline(y=0,linestyle='--')
plt.legend()
plt.show()
'''


'''
#Susceptibility:
h = 10**-8
temperatures = np.linspace(0.1,4.0,50)
Ds = [5]
delta_h = 10**(-7)
  
results = compute_free_energy(Ds,h,temperatures)
resultsplus = compute_free_energy(Ds,h+delta_h,temperatures)
resultsminus = compute_free_energy(Ds,h-delta_h,temperatures)

susc = {}
for D in Ds:
    susc[D] = np.array(resultsplus[D]) - 2*np.array(results[D]) + np.array(resultsminus[D])
    susc[D] = susc[D]/(delta_h)**2

for D in Ds:
    plt.plot(temperatures,susc[D],label=str(D))
plt.axvline(x=2.269185,linestyle='--')
plt.legend()
plt.show()
'''

h = 10**-8
temperatures = np.linspace(0.1,4.0,50)
Ds = [6]
delta_h = 10**(-8)/4
  
results = compute_free_energy_temps(Ds,h,temperatures)
resultsplus = compute_free_energy_temps(Ds,h+delta_h,temperatures)

mag = {}
for D in Ds:
    mag[D] = np.array(resultsplus[D]) - np.array(results[D])
    mag[D] = mag[D]/(delta_h)**2

for D in Ds:
    plt.plot(temperatures,mag[D],label=str(D))
plt.axvline(x=2.269185,linestyle='--')
plt.legend()
plt.show()
