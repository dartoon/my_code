from scipy.integrate import quad
import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)
H=70
om=0.3
c=299790.        #speed of light [km/s]
f=open('zs_zl')
data=np.loadtxt(f)
f.close
#print vec_EE(data[:,0])
dt=c*vec_EE(data[:,1])*(vec_EE(data[:,0])-vec_EE(data[:,1]))/vec_EE(data[:,0])/H
error= input('Which error level to set? in % unit:')
erbar=dt*error/100.                   #error bar level
bias=np.random.normal(0,dt*error/100.)   ##include system error to the data
dt_nois=dt+bias
#plt.errorbar(data[:,0],dt_nois,yerr=erbar,fmt='r.')
#plt.show()
DT=np.zeros([len(data),5])
DT[:,0]=data[:,0]     #zs
DT[:,1]=data[:,1]     #zs
DT[:,2]=dt            #real DT distance 
DT[:,3]=dt_nois       #after noise
DT[:,4]=erbar         #error
#print DT
g=open('TD_result','w')
np.savetxt(g,DT,fmt='%6.5f')
g.close

