import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from tau import *

##########to generate the lnPossible funtction, one need the Ka_square; P(zl|zs); P(zs).

#########to get the P(zs):#################
f=open('int_BHBH_sl.txt')
gw=np.loadtxt(f)
f.close
num=len(gw)
N=num
ps=np.zeros([num,2])
ps[:,0]=gw[:,0]
ps[:,1]=gw[:,1]   
def Ps(zl,zs):     #the zl is for define zl<zs
    p_bef=np.trim_zeros(ps[:,1]*(zl>=ps[:,0]))[-1]
    p_af=np.trim_zeros(ps[:,1]*(zs>=ps[:,0]))[-1]
    p=p_af/(p_af+p_bef)
    return p
#print Ps(0.0001)

vec_Ps=np.vectorize(Ps)
#######to get the P(zl|zs):##################

men=np.zeros([len(ps),2,501])
for i in range(len(ps)):    #creat a mem for fast cal
    men[i]=tau(ps[i,0])
def Pl(zl,zs):
    (a,b)=men[np.sum(zs>ps[:,0])-1]
    #(a,b)=tau(zs)
    p_va=np.trim_zeros(b*(zl>=a))[-1]/(len(a)-1)
    return p_va
print Pl(0.8,0.4)
vec_Pl=np.vectorize(Pl)

#######to get the ka^2############################
c=299790.
def Ez(z,om):
  return   1/np.sqrt(om*(1+z)**3+(1-om))
def r(z,om):
  return integrate.quad(Ez, 0, z, args=(om))[0]    #use the cos distance r
vec_r=np.vectorize(r)
f=open('TD_result50','r')
data=np.loadtxt(f)
f.close()
N_p=len(data)

zl=data[:,1]
zs=data[:,0]
y=data[:,3]
err=data[:,4]
##########def lnlike######################

def lnlike(theta, zl, y, err):    #allow zs to be 2D to get 2D sub_int, from (len(data)) to (len(data),133)
    om, h0 =theta[:2]
    zs=theta[2:]
    rzl=vec_r(zl,om)
    rzs=vec_r(zs,om)
    model = c*rzl*(rzs-rzl)/rzs/h0
    #print np.shape(rzl),np.shape(model)
    likh = (-0.5*(pow((y-model),2)/pow(err,2)))+np.log(vec_Pl(zl,zs))+np.log(vec_Ps(zl,zs))
    return np.sum(likh)

L1= [0.3, 70]
zs=data[:,0]
theta0=L1+zs.tolist()
#import time
#start=time.clock()
#print lnlike(theta=theta0, zl=data[:,1], y=data[:,3], err=data[:,4])
##print sub_like(theta=theta0, zl=data[1,1], y=data[1,3], err=data[1,4])
#elapsed=(time.clock()-start)
#print ('time1:',elapsed)

'''
###############minmize like#######################
import scipy.optimize as op
nll = lambda *args: -lnprob(*args)
L1= [0.27, 68]
zs=data[:,0]
result = op.minimize(nll, L1+zs.tolist(), args=(zl, y, err))
#m_om, m_h0, m_zs0 ,m_zs1 ,m_zs2 = result["x"]
#print result["x"]
'''

##############set prior and post post poss###############
def lnprior(theta,zl):
    [om, h0] =theta[:2]
    if 0.2 < om < 0.5 and 50 < h0 < 100 and all(theta[2:]>zl) and all(theta[2:]<20) :
        return 0.0
    return -np.inf

def lnprob(theta, zl, y, err):
    om, h0 =theta[:2]
    zs=theta[2:]
    step=np.sum(sampler.flatchain!=0)/(ndim*nwalkers)
    print "step:",step,"percent",round(step/3.,2),"%","\r",#"theta:",theta
    lp = lnprior(theta,zl)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, zl, y, err)

##################MCMC#########################################
ndim, nwalkers = N_p+2, 200
pos = [L1+zs.tolist() + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(zl, y, err))  ###the input args() should change
sampler.run_mcmc(pos, 300)

###############save the data###############
samples = sampler.chain.reshape((-1, ndim))   #sampler.chain[:, 50:, :] to cut the data
burn=sampler.chain[:,:,:].T
import pickle
pickle.dump([samples,burn], open('samp_nozs_50.txt','wb'))


samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
import corner
fig = corner.corner(samples[:,:2], labels=["$om$", "$h0$"],
                      truths=[0.3, 70])
#fig.savefig("triangle.png")

plt.show()

