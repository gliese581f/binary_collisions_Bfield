from types import FunctionType
import marshal

import numpy as np
import itertools
import math
import random
import math
from pylab import *
from scipy import integrate
global Bzs,q,m, vmag, wps

Bzs = 1.0e-1
q  =-1.602e-19    
m = 9.109e-31  
vmag = 1
wps = 1./np.abs( q*Bzs/m )	
r2 = np.abs(np.sqrt(1)  * m/(q*Bzs))
bbar = 2*q**2/(0.5*m*vmag**2)
kp = np.abs(bbar/r2)*(1.0/np.sqrt(2))  
#global Bzs,q,m, vmag, wps

def gen(n):
    bob = []
    for xx in itertools.islice(sampler(), n):
        bob.append(xx)
    return bob


def bobby(xx):
    Bz = 1.0e1
    wp = 1./np.abs( q*Bz/m )	
    vmag = np.sqrt(0**2+1**2+5**2)
    r2 = np.abs(1*np.sqrt(2) * m/(q*Bz))
    bbar = 2*q**2/(0.5*m*vmag**2)
    kp = np.abs(bbar/r2)*(1.0/np.sqrt(2))
    width = 2*pv0[2]/(5*wp)
    t = np.arange(0*wp, width*wp, 0.1*wp)		
    pv0 = [ xx*max([bbar,r2])/20, xx*max([bbar,r2])/20, 60*max([bbar,r2]), 0,1,-5 ]
    bz1.append( xx*max([bbar,r2])/20 )
    pv = integrate.odeint(solver, pv0, t) 
    return ( ((np.sqrt(pv[:,3]**2+pv[:,4]**2))[pv[:,4].size-1]-np.sqrt(pv0[4]**2+pv0[3]**2) ) )





def bandpass(dic,freqL,freqH,sample):
	import scipy
	rate = 1./(2*sample)
	band = []
	bL = (freqL/rate)    #split into smaller bands and compare/
	bH = (freqH/rate)
	b,a = scipy.signal.butter(8,bH,btype='lowpass')
	d,c = scipy.signal.butter(8,bL,btype='highpass')

	band = scipy.signal.filtfilt(d,c,(scipy.signal.filtfilt(b,a,dic)))

	return band


def integrates_multi(n):
    # Sum elements and elements squared
    from mc import make_applicable, make_mappable, sampler, stuff
    from multiprocessing import Pool

 
    if __name__ == "__main__":
        pool    = Pool(processes=4)
        results = [pool.apply_async(*make_applicable(stuff,n)) for x in range(0,n)]
        kk = [result.get(timeout=100000000000) for result in results]
        pool.close()
        avg,var = 0.0,0.0
        for values in xrange(len(kk)):
            avg += kk[i][0]
            var    += kk[i][1]
        avg1 = avg/n
        var1 = var/n
        sample_mean = avg1/n
        sample_var = (var1 - ((avg1/n)**2)/n)/(n-1.0)

    return (1*sample_mean, 1*np.sqrt(sample_var/n))


def integrates(integrand, sampler, measure=1.0, n=100):
    import math
    import itertools
    import random

    r2 = np.abs(np.sqrt(vmag**2)  * m/(q*Bzs))
    bbar = 2*q**2/(0.5*m*vmag**2)
    # Sum elements and elements squared
    total = 0.0
    total_sq = 0.0
    values = []
    for xx in itertools.islice(sampler, n):
        f = integrand(xx)
        total += f
        total_sq += (f**2)
        values.append([xx[0],xx[1],xx[2]])

    # Return answer
    sample_mean = total/n
    sample_var = (total_sq - ((total/n)**2)/n)/(n-1.0)
    
    vpar = [min( np.array(values).T[0] ),max( np.array(values).T[0] )]
    vperp1 = [min( np.array(values).T[1] ),max( np.array(values).T[1] )]
    vperp2 = [min( np.array(values).T[2] ),max( np.array(values).T[2] )]
    
    pvolume = 10*max([r2,bbar])*(vpar[1]-vpar[0])*(vperp1[1]-vperp1[0])*(vperp2[1]-vperp2[0])
    
    return (pvolume*sample_mean, pvolume*math.sqrt(sample_var/n))


def integrand(x):
    import numpy as np
    from scipy import integrate
    import random
    Ex = 0.0	# Electric Field vector
    Ey = 0.0
    Ez = 0.0
    Bx = 0.0	# Magnetic field vector
    By = 0.0

    Bzs = 1.0e-2
    #q  =-1.602e-19    
    #m = 9.109e-31  
    #vmag = 1
    wps = 1./np.abs( q*Bzs/m )	
    
    vpar    = (x[0])
    vperp1 = (x[1])
    vperp2 = (x[2])

    r2 = np.abs(np.sqrt(vperp1**2+vperp2**2)  * m/(q*Bzs))
    bbar = 2*q**2/(0.5*m*vmag**2)
    kp = np.abs(bbar/r2)*(1.0/np.sqrt(2))   
   
    rhox = random.uniform(  0.0000*max([r2,bbar]),  10* max([r2,bbar])   )
    rhoy = random.uniform(  0.0000*max([r2,bbar]),  10* max([r2,bbar])   )
    
    
    #rhox = random.uniform(  0.0 *bbar,  10* bbar   )
    #rhoy = random.uniform(  0.0 *bbar,  10* bbar   )
    
    distr = np.exp((-vpar**2 - vperp1**2- vperp2**2)/(2.0*vmag**2))

    pv0 = [ rhox,rhoy,30*max([r2,bbar]), x[1],x[2], -np.abs(vpar) ]

    rhoxy = np.sqrt( rhox**2 + rhoy**2 )

    def solver(X, t0): # X contains x,y,z and dx,dy,dz , 6 elements
            bb=1.0
            x = X[0]	
            y = X[1]	
            z = X[2]	
            vx = X[3]
            vy = X[4]
            vz = X[5]
            r = np.sqrt(x**2+y**2+z**2)
            ax = q * (Ex + (vy * Bzs) - (vz * By) ) /m + 2*((q**2)/m)*x/r**3	# Lorentz force / mass
            ay = q * (Ey - (vx * Bzs) + (vz * Bx) ) /m + 2*((q**2)/m)*y/r**3
            az = q * (Ez + (vx * By) - (vy * Bx) ) /m + 2*((q**2)/m)*z/r**3
            return [vx, vy, vz, ax, ay, az]

    width = 2*pv0[2]/(np.abs(vpar)*wps)
    t = np.arange(0*wps, width*wps, 0.1*wps)		
    pv = integrate.odeint(solver, pv0, t)	
    
    #ddz = ( pv[:,5][ pv[:,4].size-1 ]**2 - pv0[5]**2 )
    #ddx = ( pv[:,3][ pv[:,4].size-1 ]**2 - pv0[3]**2 )
    #ddy = ( pv[:,4][ pv[:,4].size-1 ]**2 - pv0[4]**2 )
    
    #dd = ( ddx, ddy, ddz )
    
    delvperp = np.sqrt( pv[:,3][ pv[:,4].size-1 ]**2 + pv[:,4][ pv[:,4].size-1 ]**2)**2 - np.sqrt( pv0[3]**2 + pv0[4]**2)**2 
     
    return   rhoxy * np.abs(vpar) * delvperp**2 * distr



def sampler():
  
    import numpy as np
    import random
    while True:
        vpar1   = np.round(random.uniform(-5.,5.),8)
        vperpx1 = np.round(random.uniform(0.,5.),8)
        vperpy1 = np.round(random.uniform(0.,5.),8)
        vpar2   = np.round(random.uniform(-5.,5.),8)
        vperpx2 = np.round(random.uniform(0.,5.),8)
        vperpy2 = np.round(random.uniform(0.,5.),8)
        
        vthermal = np.sqrt(0.5*( vpar1**2 + vpar2**2 + vperpx1**2 + vperpx2**2 + vperpy1**2 + vperpy2**2 ))
        
        vperp1 = vperpx2 - vperpx1
        vperp2 = vperpy2 - vperpy1
        vpar   = vpar2   - vpar1

              
        if vthermal <=1.1 and vthermal >= 0.9 and np.abs(vpar) > 0.5 and vperp1>0.0 and vperp2>0.0 and np.sqrt(vperp1**2+vperp2**2) <=1.1 and np.sqrt(vperp1**2+vperp2**2) >=0.9:
            yield ( vpar, vperp1, vperp2 )   
