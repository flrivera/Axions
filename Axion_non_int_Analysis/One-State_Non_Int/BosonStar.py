import math

class BosonStar:
    
    def __init__(self,gamma=0):
        
        # Expansion Parameters
        self.gamma = gamma
        self.alpha = 3.4951309897
        self.beta  = 1.7526648513
        self.szero = 1.02149303631
        self.vzero = 0.93832284019
        if (self.gamma!=0):
            self.get_expansion_parameter()
        self.sigma = 1.-self.beta
        
        # Near expansion
        self.sn_near = [self.szero,0]
        self.vn_near = [self.vzero,0]
        
        # Far expansion
        self.snm_far = [[0.,0.],[self.alpha,-self.alpha*self.sigma*(self.sigma-1.)/2.]]
        self.vnm_far = [[-1.,2.*self.beta],[0.,0.]]
    
    def get_expansion_parameter(self,):
        """
            Calculates the expansion parameters alpha, beta, s0, v0 as function of gamma
            use Eq 29 and Eq 38 in 1712.06539
            """
        g=self.gamma
        if g<1:
            self.alpha = 3.495059 - 0.117682*g - 0.391600*g**2 + 0.191882*g**3 - 0.041828*g**4 - 0.041507*g**5 + 0.033020*g**6
            self.beta  = 1.752717 + 0.703934*g - 0.109101*g**2 + 0.013436*g**3 + 0.017778*g**4 - 0.018281*g**5 + 0.005129*g**6
            self.szero = 1.021494 - 0.390946*g + 0.171489*g**2 - 0.064820*g**3 + 0.004328*g**4 + 0.028849*g**5 - 0.017732*g**6
            self.vzero = 0.938204 + 0.102743*g - 0.080310*g**2 + 0.058708*g**3 - 0.037703*g**4 - 0.002557*g**5 + 0.013512*g**6
        else:
            exponent = - math.pi/2.*g**(1./2.)*math.log(math.pi/(4.*math.e)*g**(1./2.)) - 1./(6.*math.pi)*g**(-1./2.)
            alhpa_tf = (3./4.)**(1./6.) / math.gamma(1./3.) * math.pi * g**(-1./4.) * math.exp(exponent)
            beta_tf  = math.pi/2 * g**(1./2.)
            szero_tf = g**(-1./2.)
            vzero_tf = 1
            self.alpha = alhpa_tf*(0.603380 + 0.485970*g**(-1./6.) - 4.422475*g**(-2./6.) + 8.719758*g**(-3./6.) - 8.363927*g**(-4./6.) + 4.397913*g**(-5./6.) - 1.001027*g**(-1.) )
            self.beta  = beta_tf *(1.       - 0.001478*g**(-1./3.) + 0.045642*g**(-2./3.) + 0.823049*g**(-1.)    - 0.590994*g**(-4./3.) + 0.347840*g**(-5./3.) - 0.118132*g**(-2.) )
            self.szero = szero_tf*(1.       + 0.003712*g**(-1./2.) - 0.067139*g**(-1.)    - 0.436976*g**(-3./2.) - 0.107433*g**(-2.)    + 0.687868*g**(-5./2.) - 0.327405*g**(-3.) )
            self.vzero = vzero_tf*(1.       + 0.008062*g**(-1./2.) + 0.388054*g**(-1.)    - 1.245466*g**(-3./2.) + 1.486280*g**(-2.)    - 0.823199*g**(-5./2.) + 0.178605*g**(-3.) )

    def calculate_near_coefficients(self,n_near):
        """
        Calculates the coefficients s_n and v_n for the near solution, using
        Eq 22 in 1712.06539
        """
        n_current = len(self.sn_near)-1
        for n in range (n_current-1,n_near-1):
            sn1= - sum(self.sn_near[m]*self.vn_near[n-m] for m in range(0,n+1) )
            sn2=   sum(sum(self.sn_near[l]*self.sn_near[m-l]*self.sn_near[n-m] for l in range (0,m+1)) for m in range (0,n+1))
            sn = (sn1+self.gamma*sn2)/(n+2)/(n+3)
            vn = - sum(self.sn_near[m]*self.sn_near[n-m] for m in range(0,n+1) )/(n+2)/(n+3)
                
            self.sn_near.append(sn)
            self.vn_near.append(vn)
        
    def calculate_far_coefficients(self,n_far, m_far):
        """
        Calculates the coefficients s_nm and v_nm for the far solution, using
        Eq 24 in 1712.06539
        """
        
        #How much did we already calculate?
        n_current = len(self.snm_far)
        m_current = len(self.snm_far[0])
        
        #coefficients
        sigma = self.sigma
        
        #loop over n
        for n in range (0,n_far+1):
            if n==0:
                for m in range(m_current,m_far+1):
                    self.snm_far[n].append(0)
                    self.vnm_far[n].append(0)
            elif n==1:
                for m in range(m_current,m_far+1):
                    news = - (sigma+m-1.)*(sigma+m-2.)/(2.*m) * self.snm_far[n][m-1]
                    self.snm_far[n].append(news)
                    self.vnm_far[n].append(0)
            else:
                if n<n_current:
                    m_min=m_current
                else:
                    m_min=0
                    self.snm_far.append([])
                    self.vnm_far.append([])
                for m in range(m_min,m_far+1):
                    #term1
                    term1s=0
                    term1v=0
                    if m>0:
                        term1s = 2.*n*(n*sigma+m-2)*self.snm_far[n][m-1]
                        term1v = 2.*n*(n*sigma+m-2)*self.vnm_far[n][m-1]
                    #term2
                    term2s=0
                    term2v=0
                    if m>1:
                        term2s = (n*sigma+m-2)*(n*sigma+m-3)*self.snm_far[n][m-2]
                        term2v = (n*sigma+m-2)*(n*sigma+m-3)*self.vnm_far[n][m-2]
                    #term3
                    sums=0
                    sumv=0
                    for i in range (1,n):
                        for j in range(0,m+1):
                            sums += self.snm_far[i][j]*self.vnm_far[n-i][m-j]
                            sumv += self.snm_far[i][j]*self.snm_far[n-i][m-j]
                    for j in range (0,m):
                        sums += self.snm_far[n][j]*self.vnm_far[0][m-j]
                    #term4 (self coupling term)
                    sums3=0
                    for i in range(1,n):
                        for j in range(0,m+1):
                            for k in range(1,i+1):
                                for l in range(0,j+1):
                                    sums3 += self.snm_far[k][l]*self.snm_far[i-k][j-l]*self.snm_far[n-i][m-j]
                    
                    #combine
                    news = ( -term1s - term2s - sums +self.gamma*sums3 ) / (n*n-1.)
                    newv = ( -term1v - term2v - sumv ) / (n*n)
                    self.snm_far[n].append(news)
                    self.vnm_far[n].append(newv)


    def get_near_solution(self,z,n_near=10):
        """
        Calculates the near solution s(z), v(z) using Eq 21 in 1712.06539
        """
        self.calculate_near_coefficients(n_near=n_near)
        s=sum(self.sn_near[n]*z**n for n in range (0,n_near+1))
        v=sum(self.vn_near[n]*z**n for n in range (0,n_near+1))
        return s,v
                                        
    def get_far_solution(self,z,n_far=4,m_far=4):
        """
        Calculates the far solution s(z), v(z) using Eq 23 in 1712.06539
        """
        self.calculate_far_coefficients(n_far=n_far,m_far=m_far)
        s=sum(sum(self.snm_far[n][m]*(math.exp(-z)*z**(-self.sigma))**n * z**(-m) for m in range (0,m_far+1)) for n in range(0,n_far+1))
        v=sum(sum(self.vnm_far[n][m]*(math.exp(-z)*z**(-self.sigma))**n * z**(-m) for m in range (0,m_far+1)) for n in range(0,n_far+1))
        return s,v
                                                        
    def get_solution(self,z,z_match=2,n_near=10,n_far=4,m_far=4):
        """
        Calculates the solution s(z), v(z), matched at a matching points z_match
        """
        if (z<z_match):
            return self.get_near_solution(z,n_near)
        else:
                return self.get_far_solution(z,n_far,m_far)
                                                                                    
    def get_profile(self,z,z_match=2,n_near=10,n_far=4,m_far=4):
        """
        Returns the profiles s(z), v(z) for a list of z
        """
        list_s=[]
        list_v=[]
        for this_z in z:
            s,v=self.get_solution(this_z,z_match,n_near,n_far,m_far)
            list_s.append(s)
            list_v.append(v)
        return list_s,list_v
