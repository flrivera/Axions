import numpy as np


from matplotlib import pyplot as plt

from decimal import *
setcontext(ExtendedContext)
inf = Decimal(1.) / Decimal(0)


neginf = Decimal(-1.) / Decimal(0)



   


def rungekutta(S,do_write=False):

    b=0
    t=0
    V=0.93832284019
    r=.01
    h=.01
    
    
    if do_write:
        with open('Ground_state_Numerical.txt', 'w') as f:
            print  ('r', ',','\t', 'S',',','V',file=f)
        f.close()
  
    while((abs(S)<3)):
    
        l1=b
        m1=t
        k1=np.float64(-2*b)/r -S*V
        n1=np.float64(-2*t)/r -S**2
    
        fr2=r + h/2.
        fu2=S + h/2. * l1
        fs2=V + h/2. * m1
        fb2=b + h/2. * k1
        ft2=t + h/2. * n1
    
        l2=fb2
        m2=ft2
        k2=np.float64(-2*fb2)/fr2 -fu2*fs2
        n2=np.float64(-2*ft2)/fr2 -fu2**2
        fr3=r + h/2.
        fu3=S + h/2. * l2
        fs3=V + h/2. * m2
        fb3=b + h/2. * k2
        ft3=t + h/2. * n2
    
        l3=fb3
        m3=ft3
        k3=np.float64(-2*fb3)/fr3 -fu3*fs3
        n3=np.float64(-2*ft3)/fr3 -fu3**2
        
        fr4=r + h
        fu4=S + h * l3
        fs4=V + h * m3
        fb4=b + h * k3
        ft4=t + h* n3
    
        l4=fb4
        m4=ft4
        k4 =np.float64(-2*fb4)/fr4 -fu4*fs4
        n4=np.float64(-2*ft4)/fr4 -fu4**2
        
        r=r + h
        S= S + (h/6)*(l1 + 2*l2 + 2*l3 + l4)
        V= V + (h/6)*(m1 + 2*m2 + 2*m3 + m4)
        b= b + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        t= t + (h/6)*(n1 + 2*n2 + 2*n3 + n4)
        
        
        if do_write:
            with open('Ground_state_Numerical.txt', 'a') as f:
                print  (r, ',', S, ',', V,file=f)
            f.close

    return ()
    
    
    
    
    
