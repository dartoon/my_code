from scipy.integrate import quad
import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt
def tau(zs):       
      c=299790.      #speed of light [km/s]
      H=70
      om=0.3
      dH=c/H
      H0=H*1000./(3.0856e22)
      phi0=8.0e3-3*(H/100)**3
      sig0=161
      alpha=2.32
      bet=2.67
      TFactor = (sig0/c)**4*phi0*gamma((4.+alpha)/bet)/gamma(alpha/bet); 
      #print TFactor
      def EZ(z,om):
            return 1/(om*(1+z)**3+(1-om))**0.5
      def EE(z):
            return quad(EZ, 0, z , args=(om))[0]
      vec_EE=np.vectorize(EE)
      #print vec_Dist(np.arange(1.0, 4.0, 0.5))
      #zs=2              #input
      zl=np.linspace(0,zs,501)
      Dist=vec_EE(zl)
      #print np.shape(Dist),Dist[100]
      dtaudz=16*np.pi**3*dH**3*(Dist[500]-Dist)**2*Dist**2/Dist[500]**2*TFactor/EZ(zl,om) #*ymax^2
      dtaudz/= dtaudz.sum()
      #plt.plot(zl,dtaudz,'k')
      #plt.show()
      return (zl,dtaudz)
