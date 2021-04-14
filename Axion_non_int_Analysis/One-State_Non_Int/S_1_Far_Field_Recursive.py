     
import numpy as np
import matplotlib
    
    
    
    
    
    
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
    
    alpha=3.4951309897
    beta= 1.7526648513
    sigma=1-1.7526648513
    mu=1
    sumS =0.0
    sumV= 0.0
    for n in range(0,nlimit+1):
        for m in range (0,mlimit+1):
            
            s_n,v_n=Coeffiecientslarge(gamma_c,sigma,alpha,beta,nlimit,mlimit)
            sumS += s_n[n,m] * (np.exp(-rad)/rad**sigma)**n * rad**-m
            sumV += v_n[n,m] * (np.exp(-rad)/rad**sigma)**n * rad**-m
    return sumS , sumV    
    
    
    
    
    
    
    
    
    
    
    
def Coeffiecientslarge1(S1_10,V1_02,gamma_c,mu,sigma,alpha,beta,N,M):
    
    s,v= Coeffiecientslarge(gamma_c,sigma,alpha,beta,N,M)
    m_l=1
    l=1

    S_new=np.zeros((N+1,M+1),dtype=object)
    
    V_new=np.zeros((N+1,M+1),dtype=object)
   
   # alpha=3.4951309897- 0.117682*gamma_c-0.391600*gamma_c**2+0.191882*gamma_c**3-0.041828*gamma_c**4-0.041507*gamma_c**5 +0.033020*gamma_c**6
  #  beta= 1.7526648513+ 0.703934*gamma_c-0.109101*gamma_c**2+0.013436*gamma_c**3+0.017778*gamma_c**4-0.018281*gamma_c**5 +0.005129*gamma_c**6
  #  sigma=1-beta



        
       # print('monkeys')
    S_new[0,0]=0
    S_new[1,0]=S1_10
       # V[0,0]=0
       # V[0,1]=0
    V_new[1,0]=0
    V_new[0,2]=V1_02
        
 

   
    for n in range(0, N+1):
       
        if n==0:
            continue 
            #print('n=',n)
            
         #   for m in range (0,M+1):
                
              #  if m<2:
                #    coeff= 1/(1+mu)
                       
                 #   S_new[n,m]= 2*beta*S_new[0,m-1]*coeff

               # elif m>=2:
                        
                  #  if (m*(m-1)!=l*(l+1)):
                            
                    
                     #   coeff1= 1/(1+mu)
                     #   coeff2= ((m-2)*(m-3)-l*(l+1))
                     #   S_new[n,m]= (2*beta*S_new[0,m-1] + coeff2* S_new[0,m-2])*coeff1

                            
                   # else:
                           # print('monkeys2')
                        #coeff1= 1/(1+mu)
                      #  coeff2= ((m-2)*(m-3)-l*(l+1))
                       # S_new[n,m]= (2*beta*S_new[0,m-1] + coeff2* S_new[0,m-2])*coeff1
                           
                       # V_new[n,m]=V1_02
                            
                            #print(V)
                        #continue
                            
                   
                   # print(S,V)
                
        elif n==1:
            

            
            for m in range (1,M+1):
                coeff1=2*(beta + np.sqrt(1+ mu)* (sigma+ m -1))
                    
                coeff2= (l*(l+1)-(sigma +m-1)*(sigma+m-2))

                V_new[n,m]= 0
                S_new[n,m]= (S_new[n,m-1]*coeff2 - V_new[0,2]*s[1,m-1])*(1/coeff1)
                    
                
                
                
                #print('s',s[n,m])
        elif n>=2:
            #print('monkeys')    

            #print('n=',n)
            for m in range (0,M+1):
                coeff_s=((1+ mu) *n**2 - mu)
                coeff_v=1/((1+mu)*n**2)
                
                if m==0:
                    

                    #print(mu)
              
                    V_new[n,m]= -2*sum(s[p,m]*S_new[n-p,m] for p in range (0,n+1) )*coeff_v
                        
                    sums=  -1*sum((v[p,m]*S_new[n-p,m]+s[p,m]*V_new[n-p,m]) for p in range (0,n+1) )
                    gamma_term= 3*gamma_c* (sum(sum(sum(sum(s[r,t]*s[p-r,q-t]*S_new[n-p,m-q] for t in range (0,q+1))for r in range (0,p+1))for p in range (0,n+1))for q in range (0,m+1)))
                        
                    S_new[n,m]=(sums+ gamma_term)*(1/coeff_s)
                    #print('m=0',m)

                elif m==1:
 
                    #print('m=',m)
                    coeff3= 2*n*np.sqrt(1+mu)*(sigma*n+m-2)
                    sums=sum(sum((v[p,q]*S_new[n-p,m-q]+s[p,q]*V_new[n-p,m-q]) for p in range (0,n+1) )for q in range (0,m+1))
                    gamma_term=3* gamma_c* (sum(sum(sum(sum(s[r,t]*s[p-r,q-t]*S_new[n-p,m-q] for t in range (0,q+1))for r in range (0,p+1))for p in range (0,n+1))for q in range (0,m+1)))
                        
                
                    S_new[n,m]=(-1*sums+ gamma_term-coeff3*S_new[n,m-1])*(1/coeff_s)
                    #print(s[n,m])
                    sumv=sum(sum(s[p,q]*S_new[n-p,m-q] for p in range (0,n+1) )for q in range (0,m+1))
                    #print(s[n+2])
                    
                    V_new[n,m]=(-2*sumv-coeff3*V_new[n,m-1])*coeff_v
                
                elif m>=2:

                    #print('m=',m)
                    coeff1= 2*n*np.sqrt(1+mu)*(sigma*n+m-2)
                    coeff2= ((sigma*n +m-2)*(sigma*n +m-3)-l*(l+1))
                    sums=sum(sum((v[p,q]*S_new[n-p,m-q]+s[p,q]*V_new[n-p,m-q]) for p in range (0,n+1) )for q in range (0,m+1))
                    gamma_term= 3*gamma_c* (sum(sum(sum(sum(s[r,t]*s[p-r,q-t]*S_new[n-p,m-q] for t in range (0,q+1))for r in range (0,p+1))for p in range (0,n+1))for q in range (0,m+1)))
                        
                    sumv=sum(sum(s[p,q]*S_new[n-p,m-q] for p in range (0,n+1) )for q in range (0,m+1))
                        
                    S_new[n,m]=(-1*sums + gamma_term-coeff1*S_new[n,m-1]- coeff2*S_new[n,m-2])*(1/coeff_s)
                    V_new[n,m]=(-2*sumv-coeff1*V_new[n,m-1]- coeff2*V_new[n,m-2])*coeff_v

                

        

    return  S_new,V_new  

def phi_large_1(rad,S1_10,V1_02,gamma_c,mu,nlimit,mlimit):
    #S0=.01
    #mu=.01
    
    alpha=3.4951309897
    beta= 1.7526648513
    sigma=1+ beta/(-(1+mu)**(1/2))
   
    sumS =0.0
    sumV= 0.0
    for n in range(0,nlimit+1):
        for m in range (0,mlimit+1):
            
            s_n,v_n=Coeffiecientslarge1(S1_10,V1_02,gamma_c,mu,sigma,alpha,beta,nlimit,mlimit)
            sumS += s_n[n,m] * (np.exp(-np.sqrt(1+mu)*rad)/rad**sigma)**n * rad**-m
            sumV += v_n[n,m] * (np.exp(-np.sqrt(1+mu)*rad)/rad**sigma)**n * rad**-m
    return sumS , sumV