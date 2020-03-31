import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import pickle

error= input('Which error level to set? in % unit:')
'''
c=299790.
def Ez(z,om):
  return   1/np.sqrt(om*(1+z)**3+(1-om))
def r(z,om):
  return integrate.quad(Ez, 0, z, args=(om))[0]    #use the cos distance r
H=70
om=0.3
c=299790.        #speed of light [km/s]

vec_r=np.vectorize(r)
f=open('zl_zs','r')
data=np.loadtxt(f)
f.close()
zl=data[:,0]
zs=data[:,1]
#y=data[:,3]
#err=data[:,4]

#####simulate the data######
y=c*vec_r(zl,om)*vec_r(zs,om)/(vec_r(zs,om)-vec_r(zl,om))/H
erbar=y*error/100.
bias=np.random.normal(0,erbar)   ##include system error to the data
dt=y+bias
#print zl, zs,y,erbar,dt

def lnlike(theta, zl, zs, dt, erbar):
    om, h0 = theta
    #print theta
    #model = om*z +h0
    model = c*vec_r(zl,om)*vec_r(zs,om)/(vec_r(zs,om)-vec_r(zl,om))/h0
    return -0.5*(np.sum((dt-model)**2./erbar**2.))
import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
#result = op.minimize(nll, [0.25,70], args=(zl, zs, y, erbar))
#m_om, m_h0 = result["x"]
bnds = ((-0.5, 1), (0, 100))
result = op.minimize(nll, (0.30, 70), method='SLSQP', bounds=bnds,args=(zl, zs, dt, erbar))
print result["x"].shape

group=10000
para=np.zeros([group,2])
for i in range(group):
  dt=y+np.random.normal(0,erbar)
  result = op.minimize(nll, (0.30, 70), method='SLSQP', bounds=bnds,args=(zl, zs, dt, erbar))
  para[i,:]=result["x"]
  if i%100==0:
    print i/100.
#print para[:,0]

###############save the data###############

pickle.dump(para, open("error_%s"%(error),'wb'))
'''

###############load the data##############
para=pickle.load(open("error_%s"%(error),'rb'))

plt.hist(para[:,0],25)
#plt.show()
#para[:,0]=(para[:,0]-0.3)/para[:,0]
#para[:,1]=(para[:,1]-70)/para[:,0]
import corner
fig = corner.corner(para, labels=["$\Omega_{m}$", "$H_{0}$"],
                    quantiles=[0.16, 0.5, 0.84], title_kwargs={"fontsize": 12},truths=[0.3,70],
                    plot_datapoints=False,smooth=1.0,smooth1d=1.0,plot_density=True, levels=(0.6826, 0.9544),
                    title_fmt='.3f',show_titles=True )
#fig.savefig("lcdm_%s.pdf"%(value))
#####  
plt.show()
