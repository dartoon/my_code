#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 10:24:09 2018

@author: Dartoon

To random draw the BH mass function
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import fsolve

class BH_mass_function:
    '''
    A class to generate the BH-BH mass and fast calcualte the chirp mass.
    
    Parameters
    --------
        a: The power law slope alpha 
        mbh_max: The blash of blash
        mbh_min: usually define as 5 Solar mass
        vol: The number of simulated data
    '''
    def __init__(self, a=2.35, mbh_max=80., mbh_min=5., vol=1):
        self.a = a
        self.mbh_max = mbh_max
        self.mbh_min = mbh_min
        self.vol = vol
    
    def dm(self,m):
        '''
        To generate the norminazed dm power-law model.
        
        Parameter
        --------
            N: is the normalized parameter.
            
        Return
        --------
            1. A function of normalized distritbuion of m between mbh_min and mbh_max
            2. The intergral function of 1
        '''
        N = (self.mbh_max)**(-self.a+1)/(1-self.a) - (self.mbh_min)**(-self.a+1)/(1-self.a)             # The intergral function of m**a
        return m**(-self.a)/N, ((m)**(-self.a+1)/(1-self.a) - (self.mbh_min)**(-self.a+1)/(1-self.a))/N
    
    
    def gen_dm(self):
        '''
        Random draw the mass of m1, the m2 is smaller than m1 and uniform between Mmin and m1.
        In the unit of Solar mass
        Parameter
        --------
            a: The power law slope alpha 
            mbh_max: The blash of blash
            mbh_min: usually define as 5 Solar mass
            
        Return
        --------
            A sth sth
        '''
        m1 = np.zeros(self.vol)
        func = lambda m : (seed-self.dm(m)[1])
        for i in range(self.vol):
            seed = np.random.random() 
            m1[i] = fsolve(func,6)
        return m1
    
    def gen_m1m2(self):
        '''
        To random the m1 and m2 for the BH-BH. 
            1. first generate the m1 based on the power law model
            2. m2 is flat distribued between mbh_min and m1
        
        Parameter
        --------
        Return
        --------
            One list(or tuple) as (m1, m2)
        '''
        m1 = self.gen_dm()
        m2 = np.random.uniform(self.mbh_min, m1)
        return (m1, m2)
    
    def chirp_mass(self):
        m1, m2 = self.gen_m1m2()
        m_chirp = cal_chirpmass(m1, m2)
        masses = np.asarray([m_chirp, m1, m2])
        return masses

def cal_chirpmass(m1,m2):
    m_tot = m1 + m2
    return (m1*m2)**(3/5.)/m_tot**(1/5.)
    

#vol = 2000   
#BHBH = BH_mass_function(vol = vol)
#m = np.linspace(5,80,10000)
#dn, dn_int = BHBH.dm(m)
#plt.plot(m,dn)
#plt.plot(m,dn_int)
#plt.show()
#
#gendm = BHBH.gen_dm()
#plt.hist(gendm)
#plt.show()
#
#m_chirp = BHBH.chirp_mass()
#plt.hist(m_chirp)
#plt.show()

class Theta:
    '''
    A clase random generate the (an list of) Theta
    '''
    def __init__(self, vol=1):
        self.vol = vol
        
    def theta(self, cos_t, phi_p, psi_p, cos_l):
        '''
        Calculate the value of Theta given the four parameters
        
        Parameter
        --------
            Note that the phi and psi are in the unit of np.pi
        Return
        --------
            A sth sth
        '''
        phi = phi_p * np.pi
        psi = psi_p * np.pi
        F_p = 0.5*(1+cos_t**2.)*np.cos(2.*phi)*np.cos(2.*psi) - cos_t*np.sin(2.*phi)*np.sin(2.*psi)
        F_c = 0.5*(1+cos_t**2.)*np.cos(2.*phi)*np.sin(2.*psi) + cos_t*np.sin(2.*phi)*np.cos(2.*psi)
        return  2 * (F_p**2.*(1+cos_l**2)**2. + 4.*F_c**2. * cos_l**2.)**0.5
    
    def gene_theta(self):
        cos_t, phi_p, psi_p, cos_l = np.random.uniform(-1,1,size=(4,self.vol))
        value = self.theta(cos_t=cos_t, phi_p=phi_p, psi_p=psi_p, cos_l=cos_l)
        return value
    
    def gene_theta_func(self):
        '''
        Generate the Theta based on the function.
        '''
        thetas = np.zeros(self.vol)
        x = np.linspace(0,4,1000)
        Rx = np.zeros([len(x),3])      
        Rx[:,0]=x
        Rx[:,1]= 5.*x *(4-x) **3/256.
        Rx[0,2]=0
        for i in range(len(x)-1):
            Rx[i+1,2]=Rx[i,2]+Rx[i,1]*(Rx[i+1,0]-Rx[i,0])
        for i in range(self.vol):
            idx = int(np.sum(np.random.random()>Rx[:,2]))-1
            thetas[i] = Rx[idx, 0]   
        return thetas
    
#vol = 2000
#theta_class = Theta(vol = vol)
#theta = theta_class.gene_theta()
#plt.hist(theta,bins=50,normed = True)
#x=np.linspace(0, 4)
#y=5.*x*(4.-x)**3/256.
#plt.plot(x,y)
#plt.xlabel("$\Theta$")
#plt.ylabel("P($\Theta$)")
#plt.show()
#theta_func = theta_class.gene_theta()


def Ez(z,om):
    '''
    Actually, this is the reverse of the E(z)
    '''
    w = -1
    return   1/np.sqrt(om*(1+z)**3+(1-om)*(1+z)**(3*(1+w)))   
def r(z,om):
    return integrate.quad(Ez, 0, z, args=(om))[0]
vec_r=np.vectorize(r)
def dl(z,om=0.3,h0=70):
    """
    Calculate the luminosity distance. In Mpc.
    """
    c=299790.
    dl_distance = (1+z) * c/h0 * vec_r(z,om=om)
    return dl_distance
#print dl(np.array([1,2]))

def solve_z(lum_dis, om=0.3, h0=70):
    func = lambda z : (lum_dis-dl(z,om=om, h0=h0))
    zs = fsolve(func,2)
    return zs[0]
#print solve_z(15500)

def random_z(itera = 1, fname='../data/BHBHrates.dat'):
    h0=70.
    om = 0.3
    c=299790.
    dH=c/h0
    f=open(fname)
    bhbh = np.loadtxt(f)
    bhbh=bhbh[bhbh[:,0].argsort()]
    z = bhbh[:,0]
    scenario = 2        # 2 is standard low
    n0 = bhbh[:,scenario]*10**(-9) 
    Dist = vec_r(z,om)
    dotN = 4 * np.pi * dH**3. * (n0/(1+z)) * Dist**2* Ez(z,om)
    p_dotN = dotN[1:] * (z[1:]-z[:-1])
    norm_n = p_dotN.sum()
    dotN /= norm_n              # The total BHBH events, in order to normalize the total BHBH numbers
    R = np.zeros([len(z),3])
    R[:,0]=z
    R[:,1]=dotN     # The possibility density of pointing to one GW event.
    R[0,2]=0
    for i in range(len(z)-1):
        R[i+1,2]=R[i,2]+R[i,1]*(R[i+1,0]-R[i,0])
    R[:,1]=R[:,1]/R[-1,2]
    R[:,2]=R[:,2]/R[-1,2]
    zs = np.zeros(itera)
    for i in range(itera):
        idx = int(np.sum(np.random.random()>R[:,2]))    # minus one so that the idx can start from zero. (i.e. [:-1])
        zs[i] = R[idx, 0] #np.random.uniform(R[idx, 0],R[idx+1, 0])
    return zs

class gene_BHBH:
    def __init__(self, h0=70., m_type = "BHBH", rho0 = 8., scenario = 2 ):
        self.h0 = h0
#        global n0, om, rho0, r0, c
        self.c=299790. # speed of light [km/s]
        self.r0 = 1527.;	# in [Mpc] value for ET - polynomial Noise curve
#        self.r0 = 1591.;	# ET-D sensitivity
#       self. r0 = 1918.;    # ET - xylophone
        self.m_type = m_type
        self.rho0 = rho0
        self.scenario = scenario
        filename = '../data/{0}rates.dat'.format(self.m_type)
        f=open(filename)
        bhbh = np.loadtxt(f)
        self.bhbh=bhbh[bhbh[:,0].argsort()]
        self.z = self.bhbh[:,0]
        self.om = 0.3
        self.n0 = self.bhbh[:,self.scenario]*10**(-9)
        print filename, 'scenario', self.scenario, self.rho0
    def num_year_rate(self, ave_chirp_mass = 6.7):
        '''
        Numerically calculate the rate as show in Arxiv: 1409.8360
        '''
        dH=self.c/self.h0
        z = self.z
        Dist = vec_r(z,self.om)
        Mchirp = ave_chirp_mass
        N=len(z)
        x=(self.rho0/8.)*(1+z)**(1/6.)*dH*(Dist/self.r0)*(1.2/Mchirp)**(5/6.)
        Ctheta2 = np.zeros(N)
        for i in range (0,N):
            if x[i]>=0 and x[i]<=4:
                Ctheta2[i]=(1+x[i])*(4-x[i])**4/256.
            else:
                Ctheta2[i]=0
        dotN = 4 * np.pi * dH**3. * (self.n0/(1+z)) * Dist**2*Ez(z,self.om)
        dis_dotN = dotN[:-1]*Ctheta2[:-1]
        year_rate = np.sum(dis_dotN* (z[1:]-z[:-1]))
        return dis_dotN, year_rate, Ctheta2
    def mc_year_rate(self, seed_vol=10000, itera=40, a=2.35, mbh_max=80., mbh_min=5., ev_type = 'const'):
        '''
        Use MCMC to calculate the rate, with steps equals to seed_vol.
            1. the BH mass are randomly given with BH_mass_function.
            2. the Theta are randovec_rmly given using the Theta_class
        '''
        z = self.z
        dH=self.c/self.h0
        Dist = vec_r(z,self.om)
        seed_vol = seed_vol
        #==============================================================================
        #    Randomly assign the values for Chirpmass and Thetas     
        #==============================================================================
        theta_class = Theta(vol = seed_vol)
#        thetas = theta_class.gene_theta_func()
#        mass_Chirp = 6.7 
        #==============================================================================
        #    Randomly generate a list of redshift
        #==============================================================================  
        dotN = 4 * np.pi * dH**3. * (self.n0/(1+z)) * Dist**2* Ez(z,self.om)
        p_dotN = dotN[1:] * (z[1:]-z[:-1])
        norm_n = p_dotN.sum()
        dotN /= norm_n              # The total BHBH events, in order to normalize the total BHBH numbers
        R = np.zeros([len(z),3])
        R[:,0]=z
        R[:,1]=dotN     # The possibility density of pointing to one GW event.
        R[0,2]=0
        for i in range(len(z)-1):
            R[i+1,2]=R[i,2]+R[i,1]*(R[i+1,0]-R[i,0])    
        R[:,1]=R[:,1]/R[-1,2]
        R[:,2]=R[:,2]/R[-1,2]
        over_rate = np.zeros(itera)
        zs_detected, rhos_detected = np.array([]), np.array([])
        for j in range(itera):
            thetas = theta_class.gene_theta()
            zs = np.zeros(seed_vol)
            dist_zs = np.zeros(seed_vol)
            for i in range(seed_vol):
                idx = int(np.sum(np.random.random()>R[:,2]))    # minus one so that the idx can start from zero. (i.e. [:-1])
                zs[i] = R[idx, 0] #np.random.uniform(R[idx, 0],R[idx+1, 0])
                dist_zs[i] = Dist[zs[i]==z]
                
            if ev_type == 'const':
                bhmass_class = BH_mass_function(vol = seed_vol, a=a, mbh_max=mbh_max, mbh_min=mbh_min)
                mass_Chirp = bhmass_class.chirp_mass()                
            elif ev_type != 'const':
                mass_Chirp = []
                if ev_type == 'mode0':
                    a_z = a[0] + a[1] * zs
                elif ev_type == 'mode1':
                    a_z = a[0] + a[1] * zs/(1+zs)
#                    print a_z[0]
                for i in range(len(a_z)):
#                    print a_z[i], zs[i], a[0] + a[1] * zs[i]/(1+zs[i])
                    bhmass_class = BH_mass_function(vol = 1, a=a_z[i], mbh_max=mbh_max, mbh_min=mbh_min)
                    mass_Chirp.append(bhmass_class.chirp_mass())
                mass_Chirp = np.asarray(mass_Chirp).T[0] #.reshape(3,len(a_z))
#                print mass_Chirp.shape
            #==============================================================================
            #   Calculate the observed events based on this vol of events 
            #   Rho = 8 Theta * r0/(dl) * (M_chirp_redshifted/1.2) **(5/6)
            #==============================================================================
            dlzs = (1+zs)*dH*dist_zs
            rhos = 8.*thetas * self.r0/dlzs * ((1+zs)*mass_Chirp[0]/1.2)**(5/6.)
            n_over_8 = np.sum([rhos>self.rho0])
            over_rate[j] = n_over_8/float(seed_vol)
            zs_detected = np.concatenate((zs_detected,zs[rhos>self.rho0]))
            if j == 0:
                masses = mass_Chirp.T[rhos>self.rho0]
            elif j >0:
                masses = np.concatenate((masses, mass_Chirp.T[rhos>self.rho0]))
            rhos_detected = np.concatenate((rhos_detected,rhos[rhos>self.rho0]))
            if j/5 > (j-1)/5:
                print "Total itera:", itera, "; Finish itera:", j
        av_over_rate = np.average(over_rate)
        event_rate = av_over_rate*norm_n
        '''
        One can also MC at each redshfit bins
        over_rate = np.zeros(len(z)-1)
        for i in range(len(z)-1):
            vol = 2001
            theta_class = Theta(vol = vol)
            thetas = theta_class.gene_theta_func()
            mass_Chirp = 6.7
            dist_zs = Dist[i]
            dlzs = (1+z[i])*dH*dist_zs
            rhos = 8.*thetas * self.r0/dlzs * ((1+z[i])*mass_Chirp/1.2)**(5/6.)
            n_over_8 = np.sum(rhos>self.rho0)
            over_rate[i] = float(n_over_8)/float(vol)
        n_detected = (p_dotN* over_rate).sum()
        return n_detected, over_rate, rhos
        '''
        return event_rate, zs_detected, masses, rhos_detected
#    def sim_BHBH(self, seed_vol= 200000, a0=2.35, a1=2.3, mbh_max=80., mbh_min=5., etype='linear'):
#        if etype='linear'
#        
#        self.a = a0 + a1 * z
        
        
        

#test = gene_BHBH()
##dis_dotN, year_rate, Ctheta2 = test.num_year_rate()
##event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(a=1.5, mbh_max=45.)
#event_rate, zs_detected, masses, rhos_detected = test.mc_year_rate(seed_vol= 1000, itera = 4,a=2.35, mbh_max=80.)
#print event_rate
#
#plt.hist(zs_detected)
#plt.show()
#plt.plot(test.z[1:], test.num_year_rate()[0])
#plt.show()
#
##plot the m1 intrinsic mass curve
#m = np.linspace(5,80,1000)
#a=2.35
#mbh_max=80.
#mbh_min=5.
#Norm = (mbh_max)**(-a+1)/(1-a) - (mbh_min)**(-a+1)/(1-a)
#dn_dm = m **(-a) / Norm
#print len(masses), dn_dm.sum() * (80-5)/1001.
#plt.hist(masses[:,1],bins=30, normed =1)
#plt.plot(m,dn_dm, 'r')
##plt.xlim(3,30)
##plt.yscale('log', nonposy='clip')
#plt.show()

#    '''
#    The Class to calculate the BHBH events yearly detection rate and randomly obtain the data format
#    ''' 
