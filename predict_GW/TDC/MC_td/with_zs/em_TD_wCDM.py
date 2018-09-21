import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
c=299790.
def Ez(z,om,w):
  return   1/np.sqrt(om*(1+z)**3+(1-om)*(1+z)**(3*(1+w)))
def r(z,om,w):
  return integrate.quad(Ez, 0, z, args=(om,w))[0]    #use the cos distance r
vec_r=np.vectorize(r)
f=open('TD_result500','r')
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
    om, w , h0 = theta
    #print theta
    #model = om*z +h0
    model = c*vec_r(zl,om,w)*(vec_r(zs,om,w)-vec_r(zl,om,w))/vec_r(zs,om,w)/h0
    return -0.5*(np.sum((y-model)**2/err**2))

'''
import scipy.optimize as op
nll = lambda *args: -lnlike(*args)

bnds = ((0, None),(-2.5,-0.001), (0, None))
#result = op.minimize(nll, (0.30, -1, 70), method='SLSQP', bounds=bnds,args=(zl, y, err))

result = op.minimize(nll, [0.3,-1.0,70], method='SLSQP', bounds=bnds,args=(zl,zs, y, err))
m_om, m_w, m_h0 = result["x"]


print m_om, m_w, m_h0
'''

###MCMC###
def lnprior(theta):
    om, w, h0 = theta
    if 0.2 < om < 0.45 and -1.5<w<0.5 and 55 < h0 < 85:
        return 0.0
    return -np.inf

def lnprob(theta, zl, zs, y, yerr):
    lp = lnprior(theta)
    step=np.sum(sampler.flatchain!=0)/(ndim*nwalkers)
    print "step:",step,"percent",round(step/5.,2),"%","\r","theta:",theta,
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, zl, zs, y, err)
ndim, nwalkers = 3, 50
pos = [[0.3,-1,70] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(zl, zs, y, err))  ###the input args() should change
sampler.run_mcmc(pos, 500)

###############save the data###############
samples = sampler.chain.reshape((-1, ndim))   #sampler.chain[:, 50:, :] to cut the data
burn=sampler.chain[:,:,:].T
import pickle
pickle.dump([sampler.chain,burn], open('eas_500_wCDM.txt','wb'))


samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
import corner
fig = corner.corner(samples, labels=["$om$", "$w$", "$h0$"],
                      truths=[0.3, -1, 70])
#fig.savefig("triangle.png")

plt.show()


