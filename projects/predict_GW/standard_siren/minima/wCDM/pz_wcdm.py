import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from scipy import interpolate
   
def pz(om,w,h0):
   f0=open('NSNSrates.dat')
   f1=open('BHNSrates.dat')
   f2=open('BHBHrates.dat')
   nsns,nsbh,bhbh=np.loadtxt(f0),np.loadtxt(f1),np.loadtxt(f2)
   f0.close,f1.close,f2.close
   nsns=nsns[nsns[:,0].argsort(),]
   nsbh=nsbh[nsbh[:,0].argsort(),]
   bhbh=bhbh[bhbh[:,0].argsort(),]
   #print nsns
   
   c=299790.; # speed of light [km/s]
   dH=c/h0;
   rho0 = 8;
   
   # detector horizon 
   r0 = 1527.;	# in [Mpc] value for ET - polynomial Noise curve
   #r0 = 1591.;	# ET-D sensitivity
   #r0 = 1918.;    # ET - xylophone
   Mchirp0 = 1.2
   Mchirp1 = 3.2
   Mchirp2 = 6.7
   
   
   def Ez(z,om,w):
     return   1/np.sqrt(om*(1+z)**3+(1-om)*(1+z)**(3*(1+w)))
   def r(z,om,w):
      return integrate.quad(Ez, 0, z, args=(om,w))[0]    #use the cos distance r
   #print r(17,om,w)
   
   vec_r=np.vectorize(r)
   z=nsns[:,0]
   N=len(z)
   Dist = vec_r(z,om,w)
   
   x=(rho0/8.)*(1+z)**(1/6.)*dH*(Dist/r0)*(1.2)**(5/6.)
   x0=x/(Mchirp0**(5/6.))
   x1=x/(Mchirp1**(5/6.))
   x2=x/(Mchirp2**(5/6.))
   Ctheta0,Ctheta1,Ctheta2=np.zeros(N),np.zeros(N),np.zeros(N)
   for i in range (0,N):
     if x0[i]>=0 and x0[i]<=4:
        Ctheta0[i]=(1+x0[i])*(4-x0[i])**4/256.
     else:
        Ctheta0[i]=0
     if x1[i]>=0 and x1[i]<=4:
        Ctheta1[i]=(1+x1[i])*(4-x1[i])**4/256.
     else:
        Ctheta1[i]=0
     if x2[i]>=0 and x2[i]<=4:
        Ctheta2[i]=(1+x2[i])*(4-x2[i])**4/256.
     else:
        Ctheta2[i]=0
   #print Ctheta
   
   E=1/np.sqrt(om*(1+z)**3 + (1-om)*(1+z)**(3*(1+w)))
   
   n0 = nsns[:,2]*10**(-9) # standard low
   s0dz=4*np.pi*(dH)**3.*(n0/(1+z))*Dist**2*E*Ctheta0
   #plt.plot(z,sDNdz,'k')
   n1 = nsbh[:,2]*10**(-9) # standard low
   s1dz=4*np.pi*(dH)**3.*(n1/(1+z))*Dist**2*E*Ctheta1
   n2 = bhbh[:,2]*10**(-9) # standard low
   s2dz=4*np.pi*(dH)**3.*(n2/(1+z))*Dist**2*E*Ctheta2
  
   dp=np.zeros([N,4])
   dp[:,0]=z   
   dp[:,1]=s0dz   #for NS-NS note the real possibility should mutply by p(i)-p(i-1)
   dp[:,2]=s1dz   #for NS-BH
   dp[:,3]=s2dz   #for BH-BH
   ps=np.zeros([N,2])
   ps[:,0]=dp[:,0]
   ps[:,1]=dp[:,1]+dp[:,2]+dp[:,3]
   #p[:,1]=p[:,1]

   ps0=ps[ps[:,0]<2.98]
   ps1=ps[ps[:,0]>2.98]
   newx=np.linspace(ps1[:,0].min(),ps1[:,0].max(),144)
   tck=interpolate.splrep(ps1[:,0],ps1[:,1])
   newy=interpolate.splev(newx,tck)
   ps2=np.vstack((newx,newy)).T
   newps=np.vstack((ps0,ps2))
   #plt.plot(ps[:,0],ps[:,1],'k')
   #plt.plot(newps[:,0],newps[:,1])
   #plt.show()
   N=len(newps[:,0])
   R=np.zeros([N+1,3])
   R[:,0]=np.append(0,newps[:,0])
   R[:,1]=np.append(0,newps[:,1])
   R[0,2]=0
   for i in range (1,N):
	   R[i+1,2]=R[i,2]+R[i,1]*(R[i+1,0]-R[i,0])
   #print R
   R[:,1]=R[:,1]/R[N,2]
   newps[:,1]=newps[:,1]/R[N,2]   #normaized the possibility density
   R[:,2]=R[:,2]/R[N,2]
   return newps
#print pz(0.3,-1,70)
