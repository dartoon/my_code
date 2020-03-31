import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
c=299790.
def Ez(z,om):
  return   1/np.sqrt(om*(1+z)**3+(1-om))
def r(z,om):
  return integrate.quad(Ez, 0, z, args=(om))[0]    #use the cos distance r

vec_r=np.vectorize(r)
f=open('TD_result_test','r')
data=np.loadtxt(f)
f.close()
zl=data[:,1]
zs=data[:,0]
y=data[:,3]
err=data[:,4]
'''
om=0.3
h0=70
model = c*vec_EE(data[:,1])*(vec_EE(data[:,0])-vec_EE(data[:,1]))/vec_EE(data[:,0])/h0
print model
'''

def lnlike(theta, zl, zs, y, err):
    om, h0 = theta
    #print theta
    #model = om*z +h0
    model = c*vec_r(zl,om)*(vec_r(zs,om)-vec_r(zl,om))/vec_r(zs,om)/h0
    return -0.5*(np.sum((y-model)**2/err**2))
#import scipy.optimize as op
#nll = lambda *args: -lnlike(*args)
#result = op.minimize(nll, [0.25,70], args=(zl, zs, y, err))
#m_om, m_h0 = result["x"]
#print m_om, m_h0
print lnlike(theta=[0.2,70], zl=data[:,1], zs=data[:,0], y=data[:,3], err=data[:,4])
print lnlike(theta=[0.3,70], zl=data[:,1], zs=data[:,0], y=data[:,3], err=data[:,4])
print lnlike(theta=[0.4,70], zl=data[:,1], zs=data[:,0], y=data[:,3], err=data[:,4])
