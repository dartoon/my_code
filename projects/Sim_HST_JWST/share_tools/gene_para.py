#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 16:40:47 2017

@author: dxh

To generate the random parameters for simulation
"""
#==============================================================================
# Mag=17.7-5*(np.log10(cosmo.luminosity_distance(0.4546).value*10**6)-1)
#==============================================================================

from astropy.cosmology import FlatLambdaCDM
import numpy as np
from lenstronomy.Util import constants as cont
from lenstronomy.Cosmo.lens_cosmo import LensCosmo

class gene_para:
    """
    Generate the random lens parameter, including lens mass & light, source light, 
    para: fixh0 is the SEED for H0, if not set H0 SEED is the self.seed.
          fixz is to set the z, if None, generate randomly, if len()=2, [0]=lens_z,[1]=source_z
    """
    def __init__(self,seed,fixh0=None,fixz=None,rung3=False):
        if isinstance(fixh0,int):
            self.seed=fixh0
            np.random.seed(self.seed)
            self.H0=np.random.uniform(50,90)
        elif fixh0==None:
            self.seed=seed
            np.random.seed(self.seed)
            self.H0=np.random.uniform(50,90)
        else:
            raise ValueError("the fixh0 is not give correctly" % fixh0)
        self.seed=seed + 1000
        self.cosmo=FlatLambdaCDM(H0=self.H0, Om0=0.27)
        if not isinstance(fixz, np.ndarray):
            np.random.seed(self.seed)
            self.z_lens= np.random.normal(0.5,0.2)
            self.z_source=np.random.normal(2.0,0.4)
        elif isinstance(fixz, np.ndarray) and len(fixz)==2:
            self.z_lens=fixz[0]
            self.z_source=fixz[1]
        else:
            raise ValueError("the fixz is not give correctly" % fixz)
        lensunits=LensCosmo(z_lens=self.z_lens, z_source=self.z_source,cosmo=self.cosmo)
        self.D_l=lensunits.D_d
        self.D_s=lensunits.D_s
        self.D_ls=lensunits.D_ds
        
    def spemd(self):
        """
        Generate the random lens mass model as spemd type
        """
        np.random.seed(self.seed)
        self.velo_d=250+np.random.normal(0,25)   #velocity dispersion in Ein_r km/s
        theta_E=4.*np.pi*(self.velo_d/(cont.c/1000.))**2.*self.D_ls/self.D_s/cont.arcsec  #rough Einstein radius in arcsec
        gamma= 2+np.random.normal(0,0.1)
        q=np.random.uniform(0.7,1)
        phi_G=np.random.uniform(0,np.pi)
#        return theta_E,gamma,q,phi_G
        return {'theta_E': theta_E, 'gamma': gamma, 'center_x': 0, 'center_y': 0, 'q': q, 'phi_G': phi_G}
         
    def shear(self):
        """
        Generate the random lens shear
        """
        np.random.seed(self.seed)
        b = np.random.uniform(0,0.05)        #scaled shear
        pa = np.random.uniform(0,np.pi)      #orientation
        e1 = -b*np.cos(2*pa)
        e2 = -b*np.sin(2*pa)
        return {'e1': e1, 'e2': e2},{'b': b,'pa':pa/(np.pi)*180}  #use para.shear()[0] or para.shear()[1]
        
    def lens_light(self):
        np.random.seed(self.seed)
        mag=np.random.uniform(17,19)
        ratio=np.random.uniform(0.5,1)
        R_sersic=self.spemd()['theta_E']*ratio
        n_sersic=np.random.uniform(2,4)
        q=np.random.uniform(0.7,1)
        l_phi_G=self.spemd()['phi_G']*np.random.uniform(0.9,1.1)
        return {'mag_sersic': mag, 'R_sersic': R_sersic, 'n_sersic': n_sersic, 'phi_G': l_phi_G, 'q': q}
    
    def source_light(self):
        np.random.seed(self.seed)
        mag=np.random.uniform(20,22.5)
        if self.z_source>=1.0 and self.z_source<1.5:
            R_sersic=np.random.uniform(0.37,0.45)
        elif self.z_source>=1.5 and self.z_source<2.0:
            R_sersic=np.random.uniform(0.34,0.42)
        elif self.z_source>=2.0 and self.z_source<2.5:
            R_sersic=np.random.uniform(0.31,0.35)
        elif self.z_source>=2.5 and self.z_source<=3.0:
            R_sersic=np.random.uniform(0.23,0.42)
        else:
            raise ValueError("the self.z_source is out of range (1.0 ~ 3.0)" % self.z_source)
        n_sersic=np.random.uniform(2,4)
        q=np.random.uniform(0.7,1)
        phi_G=np.random.uniform(0,np.pi)
        return {'mag_sersic': mag, 'R_sersic': R_sersic, 'n_sersic': n_sersic, 'phi_G': phi_G, 'q': q}
    
#lens_p=gene_para(seed=4)
#print lens_p.spemd(),lens_p.velo_d

def para_hist(cata,name,num=500):
    '''
    To see the distributions of the random data
    input: para_hist(cata,name):
        cata= 'spemd'; 'shear'; 'lens_light'; 'source_light'
        'spemd'={'theta_E', 'gamma', 'center_x', 'center_y', 'q', 'phi_G'}
        'shear'={'e1', 'e2'}
        'light'={'mag_sersic', 'R_sersic', 'n_sersic', 'phi_G', 'q'}
    '''
#    num=500
    hist=np.zeros(num)
#    cata=catalog #input("The catalog 'spemd'; 'shear'; 'lens_light'; 'source_light':\n")
    for i in range(num):
        para=gene_para(seed=i)
        if cata=='spemd':
            hist[i]=para.spemd()[name]
        if cata=='shear':
            hist[i]=para.shear()[name]
        if cata=='lens_light':
            hist[i]=para.lens_light()[name]
        if cata=='source_light':
            hist[i]=para.source_light()[name]
    import matplotlib.pyplot as plt
    plt.hist(hist, bins='auto')
    plt.show()

    
        
        
