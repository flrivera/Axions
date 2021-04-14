from __future__ import unicode_literals
import sympy as sym
from sympy import summation, var, solve,symbol
sym.init_printing()
from matplotlib import pyplot as plt
s0 = sym.Symbol('s0')
v0 = sym.Symbol('v0')
sigma=sym.Symbol('sigma')
gamma_c=sym.Symbol('gamma_c')
z=sym.Symbol('r')
rad=sym.Symbol('rad')
import numpy as np

from IPython.display import display, Markdown, Latex
import sympy as sym
from sympy import summation, var, solve,symbol
sym.init_printing()
from matplotlib import pyplot as plt


import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
import matplotlib.pyplot as plt

r = sym.Symbol('r',integer=True)
S0 = sym.Symbol('S0',integer=True)
sigma=sym.Symbol('sigma',integer=True)
beta= sym.Symbol('beta',integer=True)
alpha= sym.Symbol('alpha',integer=True)
mu= sym.Symbol('mu',integer=True)

import numpy as np


def Coeffiecients(gamma_c,N):
    s=np.zeros((N+1),dtype=object)
    v=np.zeros((N+1),dtype=object)
   
    s[1]=0
    v[1]=0
    
    
    v[0]=0.93832284019+ 0.102743*gamma_c- 0.080310*gamma_c**2+0.058708*gamma_c**3- 0.037703*gamma_c**4-0.002557*gamma_c**5 +0.013512*gamma_c**6
    s[0]=1.02149303631- 0.390946*gamma_c+ 0.171489*gamma_c**2-0.064820*gamma_c**3+ 0.004328*gamma_c**4+0.028849*gamma_c**5 -0.017732*gamma_c**6
    #v[0]=v0
    #s[0]=s0
    for n in range(0, N-1):
        coeff= 1/((n+2)*(n+3))

        s[n+2]=-1*sum(s[m]*v[n-m] for m in range (0,n+1) ) + gamma_c* sum(sum(s[l]*s[m-l]*s[n-m] for l in range (0,m+1))for m in range (0,n+1))
            #print(s[n+2])
        s[n+2]=s[n+2]*coeff
        v[n+2]=-1*sum(s[m]*s[n-m] for m in range(0,n+1))
        v[n+2]=v[n+2]*coeff
            #print(s[n+2],v[n+2])
      
                #print(s[0])
    return s,v
            

def phi_small(r,gamma_c,nlimit):
    
   # gamma_c=5
    s0=1.02149303631- 0.390946*gamma_c+ 0.171489*gamma_c**2-0.064820*gamma_c**3+ 0.004328*gamma_c**4+0.028849*gamma_c**5 -0.017732*gamma_c**6
    
    sumS =0.0
    sumV= 0.0
    for n in range(0,nlimit+1):
        s_n,v_n=Coeffiecients(gamma_c,nlimit)
        sumS += s_n[n] * r**n
        sumV += v_n[n] * r**n
    return sumS , sumV
    
    
    
    
def Coeffiecientslarge(gamma_c,sigma,alpha,beta,N,M):
    

    s=np.zeros((N+1,M+1),dtype=object)
    
    v=np.zeros((N+1,M+1),dtype=object)
   
    sigma=1-beta
    v[0,0]=-1
    v[0,1]= 2*beta
    s[1,0]=alpha
    s[0,0]=0
    alpha=3.4951309897- 0.117682*gamma_c-0.391600*gamma_c**2+0.191882*gamma_c**3-0.041828*gamma_c**4-0.041507*gamma_c**5 +0.033020*gamma_c**6
    beta= 1.7526648513+ 0.703934*gamma_c-0.109101*gamma_c**2+0.013436*gamma_c**3+0.017778*gamma_c**4-0.018281*gamma_c**5 +0.005129*gamma_c**6
    sigma=1-beta
    
   
    for n in range(0, N+1):
       
        if n==0:
            #print('n=',n)
            for m in range (4,M+1):
                
 
                #print('m=',m)
                coeff= 1/((m-2)*(m-3))
                #s[n,m-2]=summation(s[n,q]*v[n,m-q] for q in range (0,m))
                sums=-1* sum (s[n,q]*v[n,m-q] for q in range (0,m+1)) #n==0 no need to sum over it
                gamma_term= gamma_c* (sum(sum(s[n,t]*s[n,q-t]*s[n,m-q] for t in range (0,q+1))for q in range (0,m+1)))
                
                s[n,m-2]=(sums+ gamma_term)*coeff
                v[n,m-2]= -1*sum((s[n,q]*s[n,m-q] for q in range (0,m+1)))*coeff
                #print(s[0,m-2])
                

        elif n==1:
            
       
        
           # print('n=',n)
            
            for m in range (1,M+1):
                #print('m=',m)
                #print(m)
                v[n,m]= -1*((sigma+ m-1)*(sigma+m -2)/(2*m) ) * v[n,m-1]
                s[n,m]=  -1*((sigma+ m-1)*(sigma+m -2)/(2*m))* s[n,m-1]
                #print('s',s[n,m])
                
                
        elif n>=2:
            #print('monkeys')    

            #print('n=',n)
            for m in range (0,M+1):
                
                if m==0:
                    #print('m=0',m)
                    sums=-1* sum(s[p,0]*v[n-p,0] for p in range (0,n+1) ) #m==0 no need to sum over it
                    
                    gamma_term= gamma_c* (sum(sum(s[r,0]*s[p-r,0]*s[n-p,0] for r in range (0,p+1))for p in range (0,n+1)))
                    #print(s[n+2])
                    s[n,m]=(sums+ gamma_term)/(n*n-1)
                    #print(s[n,m])
                    sumv=-1* sum(s[p,0]*s[n-p,0] for p in range (0,n+1) )
                    #print(s[n+2])
                    v[n,m]=sumv/(n*n)
                    #print(v[2,0])
                elif m==1:
 
                    #print('m=',m)
                    coeff_1= 2.0*n*(sigma*n+m-2)
                    sums=-1* sum(sum(s[p,q]*v[n-p,m-q] for q in range (0,m+1))for p in range (0,n+1) )
                    gamma_term= gamma_c* (sum(sum(sum(sum(s[r,t]*s[p-r,q-t]*s[n-p,m-q] for t in range (0,q+1))for r in range (0,p+1))for p in range (0,n+1))for q in range (0,m+1)))
                   
                    #print(s[n+2])
                
                    s[n,m]=(sums + gamma_term-coeff_1*s[n,m-1])/(n*n-1)
                    #print(s[n,m])
                    sumv=-1*sum(sum(s[p,q]*s[n-p,m-q]for q in range (0,m+1)) for p in range (0,n+1) )
                    #print(s[n+2])
                    
                    v[n,m]=(sumv-coeff_1*v[n,m-1])/(n*n)
                    
                    #if n==2 and m==1:
                      #  print(sumv)
                    
                elif m>=2:

                    #print('m=',m)
                    coeff1= 2*n*(sigma*n+m-2)
                    coeff2= (sigma*n +m-2)*(sigma*n +m-3)
                    sums=-1* sum(sum(s[p,q]*v[n-p,m-q] for p in range (0,n+1))for q in range (0,m+1))
                    gamma_term= gamma_c* (sum(sum(sum(sum(s[r,t]*s[p-r,q-t]*s[n-p,m-q] for t in range (0,q+1))for r in range (0,p+1))for p in range (0,n+1))for q in range (0,m+1)))
                     
                    #print(s[n+2])
                    sumv=-1*sum(sum(s[p,q]*s[n-p,m-q] for p in range (0,n+1) )for q in range (0,m+1))
                    s[n,m]=(sums+ gamma_term-coeff1*s[n,m-1]- coeff2*s[n,m-2])/(n*n-1)
                    v[n,m]=(sumv-coeff1*v[n,m-1]- coeff2*v[n,m-2])/(n*n)
                    #print(s[n,m])
                    #print(v[2,1])
                    #print(s[0])
                    
                    
    return  s,v
            
            
            
def phi_large(rad,gamma_c,nlimit,mlimit):
    #S0=.01
    #mu=.01
    #gamma_c=5
    alpha=3.4951309897- 0.117682*gamma_c-0.391600*gamma_c**2+0.191882*gamma_c**3-0.041828*gamma_c**4-0.041507*gamma_c**5 +0.033020*gamma_c**6
    beta= 1.7526648513+ 0.703934*gamma_c-0.109101*gamma_c**2+0.013436*gamma_c**3+0.017778*gamma_c**4-0.018281*gamma_c**5 +0.005129*gamma_c**6
    sigma=1-beta
    #mu=1
    sumS =0.0
    sumV= 0.0
    for n in range(0,nlimit+1):
        for m in range (0,mlimit+1):
            
            s_n,v_n=Coeffiecientslarge(gamma_c,sigma,alpha,beta,nlimit,mlimit)
            sumS += s_n[n,m] * (np.exp(-rad)/rad**sigma)**n * rad**(-m)
            sumV += v_n[n,m] * (np.exp(-rad)/rad**sigma)**n * rad**(-m)
    return sumS , sumV