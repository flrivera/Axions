import sympy as sym
#from sympy import summation, var, solve,symbol
sym.init_printing(use_latex=True)
from IPython.display import display, Math, Latex
from matplotlib import pyplot as plt

import numpy as np



############# This is the old recursion relations, from your papers#############

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

        s[n+2]=-1*sum(s[m]*v[n-m] for m in range (0,n+1) ) #+ gamma_c* sum(sum(s[l]*s[m-l]*s[n-m] for l in range (0,m+1))for m in range (0,n+1))
            #print(s[n+2])
        s[n+2]=s[n+2]*coeff
        v[n+2]=-1*sum(s[m]*s[n-m] for m in range(0,n+1))
        v[n+2]=v[n+2]*coeff
            #print(s[n+2],v[n+2])
      
                #print(s[0])
    return s,v
            
            
def phi_small(r,s0,gamma_c,nlimit):
    
   # gamma_c=5
    s0=1.02149303631- 0.390946*gamma_c+ 0.171489*gamma_c**2-0.064820*gamma_c**3+ 0.004328*gamma_c**4+0.028849*gamma_c**5 -0.017732*gamma_c**6
    
    sumS =0.0
    sumV= 0.0
    for n in range(0,nlimit+1):
        s_n,v_n=Coeffiecients(gamma_c,nlimit)
        sumS += s_n[n] * r**n
        sumV += v_n[n] * r**n
    return sumS , sumV
    
    
    

############# This is the modified recursion relations###################
    
def Coeffiecients_1(gamma_c,ds0,dv0,mu,N):
    S=np.zeros((N+1),dtype=object)
    V=np.zeros((N+1),dtype=object)
   

    S[1]=ds0
    V[1]=dv0
  
    m_l=1
    l=1
    V[0]=0#V0
    S[0]=0#S0
    
    s,v=Coeffiecients(gamma_c,N)
    for n in range(0, N-1):
        coeff= 1/((n+2)*(n+3)-l*(l+1))

        S[n+2]=-1*sum((v[m]*S[n-m]+s[m]*V[n-m]) for m in range (0,n+1)) +mu*m_l*S[n] #+ 3*gamma_c*sum(sum((s[o]*s[m-o]*S[n-m]) for o in range(0,m+1))for m in range(0,n+1))
            #print(s[n+2])
        S[n+2]=(S[n+2])*coeff
        V[n+2]=-2*sum(s[m]*S[n-m] for m in range(0,n+1))
        V[n+2]=V[n+2]*coeff

    return S,V

def phi_small_1(r,dS0,dV0,mu,gamma_c,nlimit):

  
    sumS =0.0
    sumV= 0.0
    for n in range(0,nlimit+1):
        s_n,v_n=Coeffiecients_1(gamma_c,dS0,dV0,mu,nlimit)
        sumS += s_n[n] * r**n
        sumV += v_n[n] * r**n
    return sumS , sumV