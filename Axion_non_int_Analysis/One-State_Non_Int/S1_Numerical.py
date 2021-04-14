import numpy as np
import pylab as pl
import matplotlib
import os
from mpmath import*
from IPython import display
from datetime import datetime
from scipy.optimize import minimize, curve_fit
from matplotlib import pyplot as plt
    
from BosonStar import BosonStar 




gamma=0
ground_state=BosonStar(gamma=gamma)

z=np.linspace(0,50,200)
s,v=ground_state.get_profile(z,z_match=1.3,n_near=10,n_far=4,m_far=4)










def func_s(z,s1,s1p,v1,v1p,l,gamma,ggamma):
    return s1p
def func_v(z,s1,s1p,v1,v1p,l,gamma,ggamma):
    return v1p
def func_ds(z,s1,s1p,v1,v1p,l,gamma,ggamma):
    profile = ground_state.get_profile([z])
    s0, v0 = profile[0][0], profile[1][0]
    return -2./z*s1p + l*(l+1)/z**2*s1 - v1*s0 - v0*s1 + ggamma*s1 + 3.*gamma*s0**2*s1 

def func_dv(z,s1,s1p,v1,v1p,l,gamma,ggamma):
    profile = ground_state.get_profile([z])
    s0, v0 = profile[0][0], profile[1][0]
    return -2./z*v1p + l*(l+1)/z**2*v1 - 2.*s1*s0 
    
    
    
def rungekutta_step(z, s1, ds1, v1, dv1 , h, l, gamma, ggamma): 
    
  #  print(ds1)
    k1_s  = func_s (z, s1, ds1, v1, dv1, l, gamma, ggamma)
    k1_ds = func_ds(z, s1, ds1, v1, dv1, l, gamma, ggamma)
    k1_v  = func_v (z, s1, ds1, v1, dv1, l, gamma, ggamma)
    k1_dv = func_dv(z, s1, ds1, v1, dv1, l, gamma, ggamma)
    
    
    
    
  #  print(z,k1_ds)
        
    k2_s  = func_s (z+h/2, s1+h/2*k1_s, ds1+h/2*k1_ds, v1+h/2*k1_v, dv1+h/2*k1_dv, l, gamma, ggamma)
    k2_ds = func_ds(z+h/2, s1+h/2*k1_s, ds1+h/2*k1_ds, v1+h/2*k1_v, dv1+h/2*k1_dv, l, gamma, ggamma)
    k2_v  = func_v (z+h/2, s1+h/2*k1_s, ds1+h/2*k1_ds, v1+h/2*k1_v, dv1+h/2*k1_dv, l, gamma, ggamma)
    k2_dv = func_dv(z+h/2, s1+h/2*k1_s, ds1+h/2*k1_ds, v1+h/2*k1_v, dv1+h/2*k1_dv, l, gamma, ggamma)
      
   # print(z,k2_ds)
    
    k3_s  = func_s (z+h/2, s1+h/2*k2_s, ds1+h/2*k2_ds, v1+h/2*k2_v, dv1+h/2*k2_dv, l, gamma, ggamma)
    k3_ds = func_ds(z+h/2, s1+h/2*k2_s, ds1+h/2*k2_ds, v1+h/2*k2_v, dv1+h/2*k2_dv, l, gamma, ggamma)
    k3_v  = func_v (z+h/2, s1+h/2*k2_s, ds1+h/2*k2_ds, v1+h/2*k2_v, dv1+h/2*k2_dv, l, gamma, ggamma)
    k3_dv = func_dv(z+h/2, s1+h/2*k2_s, ds1+h/2*k2_ds, v1+h/2*k2_v, dv1+h/2*k2_dv, l, gamma, ggamma)
        
   # print(k3_s)
        
    k4_s  = func_s (z+h, s1+h*k3_s, ds1+h*k3_ds, v1+h*k3_v, dv1+h*k3_dv, l, gamma, ggamma)
    k4_ds = func_ds(z+h, s1+h*k3_s, ds1+h*k3_ds, v1+h*k3_v, dv1+h*k3_dv, l, gamma, ggamma)
    k4_v  = func_v (z+h, s1+h*k3_s, ds1+h*k3_ds, v1+h*k3_v, dv1+h*k3_dv, l, gamma, ggamma)
    k4_dv = func_dv(z+h, s1+h*k3_s, ds1+h*k3_ds, v1+h*k3_v, dv1+h*k3_dv, l, gamma, ggamma)
        
    #print(z,k1_ds,k2_ds,k3_ds,k4_ds)    
        
    ksum_s  = k1_s + 2*k2_s + 2*k3_s + k4_s
    ksum_ds = k1_ds+ 2*k2_ds+ 2*k3_ds+ k4_ds
    ksum_v  = k1_v + 2*k2_v + 2*k3_v + k4_v
    ksum_dv = k1_dv+ 2*k2_dv+ 2*k3_dv+ k4_dv
        
    step_s  = h/6.*ksum_s
    step_ds = h/6.*ksum_ds
    step_v  = h/6.*ksum_v
    step_dv = h/6.*ksum_dv
    
    return step_s, step_ds, step_v, step_dv
    
    
    
    
def onestatesolution(l,gamma,ggamma,s11,do_write=False):
   
    #step size 
    h=.01 
    
    #initial consitions
    z=h/2.
    dS_1, dV_1 = s11, 1
   
    S_1, V_1 = dS_1*z, dV_1*z
    
    
    #prepare output
    output_z, output_S1, output_dS1, output_V1, output_dV1 = [], [], [], [], []
    
    if do_write:
        with open(f'Type Solutions{ggamma+1} l={l}.csv', 'w') as f:
            print  ('r', ',','\t', 'S',',','\t','dS',',','\t','V',',','\t','dV',file=f)
        f.close()
  
    #Condition to ensure I see full behaivior up to when the functions starts diverging.
    while(abs(S_1)<3) :
#    while(z<20) :
        if do_write:
            with open(f'Type Solutions{ggamma+1} l={l}.csv', 'a') as f:
                print  (z, ',','\t', S_1,',','\t',dS_1,',','\t',V_1,',','\t',dV_1,file=f)
            f.close()

        output_z.append(z)
        output_S1.append(S_1)
        output_dS1.append(dS_1)
        output_V1.append(V_1)
        output_dV1.append(dV_1)
 
        step_s, step_ds, step_v, step_dv= rungekutta_step(z, S_1, dS_1, V_1, dV_1, h, l, gamma, ggamma)
        z    = z    + h
        S_1  = S_1  + step_s
        V_1  = V_1  + step_v
        dS_1 = dS_1 + step_ds
        dV_1 = dV_1 + step_dv
        
    return np.array(output_z),np.array(output_S1),np.array(output_dS1),np.array(output_V1),np.array(output_dV1)