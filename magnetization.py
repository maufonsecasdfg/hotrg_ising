#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from hotrg_ising import *

def compute_magnetization_temps(dim,Ds,h,temperatures,compare=True):
    results = {}
    error = {}
    for D in Ds:
        results[D] = []
        for temp in temperatures:
            print 'temp: %s' %(temp)
            beta = 1/temp
            a = isingtpf(dim,beta,h)
            b = isings(dim,beta,h)
            converg_criteria = 10**(-6)
            
            i = 0
            s_i = np.eye(D)
            r_i = 2
            while True:
                
                a, b, maxAA = coarse_graining_step(dim,a,b=b,D=D)
                if i >= 10:
                    ap, bp = final_contraction(dim,a,b)
            
                    r = (bp/ap)
                    delta_r = np.abs((r - r_i)/r)
                    if delta_r <= converg_criteria:
                        print i
                        break
                    else:
                        r_i = r
                i += 1
    
            ap, bp = final_contraction(dim,a,b)
            
            results[D].append(bp/ap)
    if compare and dim ==2:
        exact_result = np.zeros(temperatures.size)
        for i in range(temperatures.size):
            if temperatures[i] < 2./np.log(1+np.sqrt(2)):
                try:
                    exact_result[i] = (1-(np.sinh(2/temperatures[i]))**-4)**(1/8.)
                except:
                    exact_result[i] = 0
            else:
                exact_result[i] = 0
        for D in Ds:                
            error[D] = np.abs(exact_result - results[D])/(np.ones(temperatures.size)+np.abs(exact_result))
            
    return results, error
    
def compute_magnetization_hs(dim,Ds,hs,temp):
    results = {}
    beta = 1/temp
    for D in Ds:
        results[D] = []
        for h in hs:
            print 'h: %s' %(h)
            a = isingtpf(dim,beta,h)
            b = isings(dim,beta,h)
            cconverg_criteria = 10**(-6)
            
            i = 0
            r_i = 2
            while True:
                a, b, maxAA = coarse_graining_step(dim,a,b=b,D=D)
                if i >= 10:
                    ap, bp = final_contraction(dim,a,b)
            
                    r = (bp/ap)
                    delta_r = np.abs((r - r_i)/r)
                    if delta_r <= converg_criteria:
                        print i
                        break
                    else:
                        r_i = r
                i += 1
    
            ap, bp = final_contraction(dim,a,b)
            
            results[D].append(bp/ap)
            
    return results

def compute_susceptibility(dim,Ds,h,temperatures):
    delta_h = np.abs(h)/8.
    resultsminus, e = compute_magnetization_temps(dim,Ds,h-delta_h/2,temperatures,compare=False)
    resultsplus,e  = compute_magnetization_temps(dim,Ds,h+delta_h/2,temperatures,compare=False)

    
    susc = {}
    for D in Ds:
        susc[D] = np.array(resultsplus[D]) - np.array(resultsminus[D])
        susc[D] = -susc[D]/(delta_h)
        
    return susc

def write_results(results,xaxis,docu_name):
	f = open('./%s.csv' %(docu_name),'w')
	for x in xaxis:
		f.write('%s,' %(x))
	f.write('/n')
	for D in results:
		f.write(str(D))
		f.write('/n')
		for y in results[D]:
			f.write('%s,' %(y))
		f.write('/n')
	f.close()
	


h = 10**(-10)
temperatures = np.linspace(0.7,5.5,100)
Ds = [2,3,4]

#results, error = compute_magnetization_temps(3,Ds,h,temperatures)
#write_results(results,temperatures,'magnetization2d')
#write_results(error,temperatures,'magnetization2d_errors')

#for D in Ds:
#	plt.semilogy(temperatures,np.abs(error[D]),label=str(D))	
#plt.axvline(x=2./np.log(1+np.sqrt(2)),linestyle='--',color='grey')
#plt.xlabel('Temperature')
#plt.ylabel('Error')
#plt.legend()
#plt.savefig('magnetization2derrorcomplete.png')
#plt.show()

for D in Ds:
	plt.plot(temperatures,results[D],label=str(D))	
#plt.axvline(x=2./np.log(1+np.sqrt(2)),linestyle='--',color='grey')
plt.axvline(x=4.511544,linestyle='--',color='grey')
plt.xlabel('Temperature')
plt.ylabel('Magnetization per site')
plt.legend()
plt.savefig('magnetization3d.png')
plt.show()


#h = 10**(-10)
#temperatures = np.linspace(3.0,6,25)
#Ds = [5]
#
##results, error = compute_magnetization_temps(3,Ds,h,temperatures)
##results2, error2 = compute_magnetization_temps(3,Ds,h,temperatures)
#results3, error3 = compute_magnetization_temps(3,Ds,h,temperatures)
#
#for D in Ds:
#	plt.plot(temperatures,results[4],label=str(D))
#	plt.plot(temperatures,results2[3],label=str(D))
#	plt.plot(temperatures,results3[3],label=str(D))
#	
#plt.axvline(x=4.5,linestyle='--')
#plt.legend()
#plt.show()



#
#h = 10**(-10)
#temperatures = np.linspace(0.1,3.5,100)
#Ds = [2,3,4]
##susc = compute_susceptibility(2,Ds,h,temperatures)
#
#for D in Ds:
#    plt.plot(temperatures,susc[D],label=str(D))
##plt.axvline(x=2./np.log(1+np.sqrt(2)),linestyle='--',color='grey')
#plt.xlabel('Temperature')
#plt.ylabel('Magnetic Susceptibility')
#plt.legend()
#plt.savefig('susc2d.png')
#plt.show()


