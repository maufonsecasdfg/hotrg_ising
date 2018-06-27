#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from hotrg_ising import *
import scipy.integrate as integrate

def compute_free_energy_temps(dim,Ds,h,temperatures,compute_error=True):

    results = {}
    error ={}
    
    for D in Ds:
        results[D] = []
        for temp in temperatures:
            print 'temp: %s' %(temp)
            beta = 1/temp
            a = isingtpf(dim,beta,h)
            converg_criteria = 10**(-6)
            if dim == 2:
                Nplus = 4.
            elif dim == 3:
                Nplus = 6.
            
            i = 0
            #Z = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
            Z = final_contraction(dim,a)[0]
            f_i = -temp*(np.log(Z))/(Nplus)
            C = 0
            N = 1
    
            while True:
                a, a, maxAA = coarse_graining_step(dim,a,D=D)
                C = np.log(maxAA)+Nplus*C #arreglar para 3D
                N *= Nplus
                if i >= 10:
                    #Z = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
                    Z = final_contraction(dim,a)[0]
                    f = -temp*(np.log(Z)+Nplus*C)/(Nplus*N)
                    delta_f = np.abs((f - f_i)/f)
                    #print delta_f
                    if delta_f <= converg_criteria:
                        print i
                        break
                    else:
                        f_i = f
                i += 1
                
            #Z = ncon([a,a,a,a],[[7,5,3,1],[3,6,7,2],[8,1,4,5],[4,2,8,6]])
            Z = final_contraction(dim,a)[0]
            f = -temp*(np.log(Z)+4*C)/(4*N)
            
            results[D].append(f)
        if dim == 2 and compute_error:
            exact_sol = np.zeros(temperatures.size)
            def funct_to_integrate(theta1,theta2,beta):
                return np.log((np.cosh(2*beta))**2-np.sinh(2*beta)*(np.cos(theta1)+np.cos(theta2)))
            for i in range(temperatures.size):
                beta = 1/temperatures[i]
                integ = integrate.dblquad(funct_to_integrate,0,np.pi,lambda x: 0, lambda x: np.pi,args=([beta]))[0]
                exact_sol[i] = (-1/beta)*((np.log(2)+(1/(2*np.pi**2))*integ))
                
            error[D] = np.abs(exact_sol - results[D])/np.abs(exact_sol)            
            
    return results, error

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
    
def compute_internal_energy_and_heatc(dim,Ds,h,temperatures,results=False):
    if results == False:
        f, e = compute_free_energy_temps(dim,Ds,h,temperatures,compute_error=False)
    else:
        f = results
#    for D in Ds:
#        plt.plot(temperatures,f[D],label=str(D))
#    plt.axvline(x=2.269185,linestyle='--')
#    plt.legend()
#    plt.show()
    u = {}
    c = {}
    delta_temp = temperatures[1]-temperatures[0]
    for D in Ds:
        fot = f[D]/temperatures
        fot_dif = (fot[1:] - fot[:-1])/delta_temp
        u[D] = (-temperatures[1:]**2)*fot_dif        
        c[D] = (u[D][1:]-u[D][:-1])/delta_temp
    
    return u, c

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
#h = 10**-10
#temperatures = np.linspace(0.7,5.5,100)
#Ds = [2,3,4]
#
##results, error = compute_free_energy_temps(3,Ds,h,temperatures)
#
##for D in Ds:
##	plt.semilogy(temperatures[50:85],np.abs(error[D][50:85]),label=str(D))	
##plt.axvline(x=2./np.log(1+np.sqrt(2)),linestyle='--',color='grey')
##plt.xlabel('Temperature')
##plt.ylabel('Error')
##plt.legend()
##plt.savefig('freeenergy2derror.png')
##plt.show()
#
#for D in Ds:
#	plt.plot(temperatures,results[D],label=str(D))	
#
#plt.axvline(x=4.511544,linestyle='--',color='grey')
#plt.xlabel('Temperature')
#plt.ylabel('Free energy per site')
#plt.legend()
#plt.savefig('freeenergy3d.png')
#plt.show()

#h = 10**-10
#temperatures = np.linspace(0.6,3.5,7000)
#Ds = [2,4,6,8,10]
#us, cs = compute_internal_energy_and_heatc(dim,Ds,h,temperatures)
#for D in Ds:
#	plt.plot(temperatures,us[D],label=str(D))	
#plt.axvline(x=2./np.log(1+np.sqrt(2)),linestyle='--',color='grey')
##plt.axvline(x=4.5,linestyle='--',color='grey')
#plt.xlabel('Temperature')
#plt.ylabel('Internal Energy')
#plt.legend()
#plt.savefig('freeenergy2d.png')
#plt.show()


  
#results = compute_free_energy_temps(Ds,h,temperatures)
#resultsplus = compute_free_energy_temps(Ds,h+delta_h,temperatures)

#mag = {}
#for D in Ds:
#    mag[D] = np.array(resultsplus[D]) - np.array(results[D])
#    mag[D] = mag[D]/(delta_h)**2





h = 10**-10
temperatures = np.linspace(0.7,5.5,500)
Ds = [2,3]

us, cs = compute_internal_energy_and_heatc(3,Ds,h,temperatures)

for D in Ds:
    plt.plot(temperatures[1:],us[D],label=str(D))
plt.axvline(x=4.511544,linestyle='--')
plt.xlabel('Temperature')
plt.ylabel('Internal energy per site')
plt.legend()
plt.savefig('internalenergy3d.png')
plt.show()


for D in Ds:
    plt.plot(temperatures[2:],cs[D],label=str(D))
plt.axvline(x=4.511544,linestyle='--')
plt.xlabel('Temperature')
plt.ylabel('Heat capacity')
plt.legend()
plt.savefig('heatcapacity3d.png')
plt.show()