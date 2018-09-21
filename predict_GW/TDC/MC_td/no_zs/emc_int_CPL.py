import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from tau import *

##########to generate the lnPossible funtction, one need the Ka_square; P(zl|zs); P(zs).

#########to get the P(zs):#################
f=open('../int_BHBH_sl.txt')
gw=np.loadtxt(f)
f.close
num=len(gw)
N=num
ps=np.zeros([num,2])
ps[:,0]=gw[:,0]
ps[:,1]=gw[:,1]   #not the real possibility should mutply by p(i)-p(i-1)

def Ps(zs):
    p_va=np.trim_zeros(ps[:,1]*(zs>=ps[:,0]))[-1]
    return p_va*(zs<30)
#print Ps(30)
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
vec_Pl=np.vectorize(Pl)
#print Pl(1,15)
############to inter the likehood with zs############
num=len(gw)
N=num
z_sca=np.zeros(num)
z_sca=gw[:,0]

f=open('TD_result500','r')
data=np.loadtxt(f)
f.close()
N_p=len(data)
zl=data[:,1]
zs=data[:,0]
y=data[:,3]
err=data[:,4]

def ct_zs(zl):                  #creat the zs array corrsponding to the zl, get the zs after zl in z_sca scale
    z_s=np.ones([len(zl),len(z_sca)])
    for i in range(len(zl)):
       z_s[i]=(z_sca*[zl[i]<z_sca])[0]
    z_s=z_s+30*np.ones(np.shape(z_s))*[z_s==0][0]
    return z_s

z_s=ct_zs(zl=data[:,1])
#print np.shape(z_s),np.shape(data[:,0])
#print ct_zs(zl=data[:,1])
#print sub_like(theta=[0.3,70], zs=ct_zs(zl=data[:,1]), zl=data[:,1], y=data[:,3], err=data[:,4])

################to get the ka^2 and combine together to get likehood for one data###################
c=299790.
def Ez(z,om,w0,wa):
  w=w0+wa*z/(1+z)
  return   1/np.sqrt(om*(1+z)**3+(1-om)*(1+z)**(3*(1+w)))
def r(z,om,w0,wa):
  if z < 20:
    return integrate.quad(Ez, 0, z, args=(om,w0,wa))[0]    #use the cos distance r
  else:
    return 3.33
vec_r=np.vectorize(r)

'''
def sub_like(theta, zs, zl, y, err):                 ###get the sub_lihod for one input data, !!!don't sum the ka in this step.
    om, h0 =theta
    rzl=vec_r(zl,om)
    rzs=vec_r(zs,om)
    model = c*rzl*(rzs-rzl)/rzs/h0
    #print np.shape((y-model)**2/err**2)
    likh = np.exp(-0.5*((y-model)**2/err**2))*vec_Ps(zs)*vec_Pl(zl,zs)
    #print likh
    return likh
#import time
#start=time.clock()
#print sub_like(theta=[0.3,70], zs=data[:,0], zl=data[:,1], y=data[:,3], err=data[:,4])
#print sub_like(theta=[0.3,70], zs=data[1,0], zl=data[1,1], y=data[1,3], err=data[1,4])
#elapsed=(time.clock()-start)
#print ('time1:',elapsed)
'''

z_s=ct_zs(zl)

def twod_like(theta, zl, y, err):    #allow zs to be 2D to get 2D sub_int, from (len(data)) to (len(data),133)
    om, w0 ,wa , h0 =theta
    rzl=vec_r(zl,om,w0,wa)[:,None]  #transfer rzl from (len(data),) to (len(data),1)
    #z_s=ct_zs(zl)
    rzs=vec_r(z_s,om,w0,wa)
    model = c*rzl*(rzs-rzl)/rzs/h0
    #print np.shape(rzl),np.shape(model)
    likh = np.exp(-0.5*(pow((y[:,None]-model),2)/pow(err[:,None],2)))*vec_Pl(zl[:,None],z_s)*vec_Ps(z_s)
    return likh
ddz=ps[:,0][1:]-ps[:,0][:-1]


import time
time0=time.clock()
print np.shape(np.sum(twod_like(theta=[0.3, -1, 0,70], zl=data[:,1], y=data[:,3], err=data[:,4])[:,:-1]*ddz.T,axis=1))
#print twod_like(theta=[0.3,70], zl=data[:,1], y=data[:,3], err=data[:,4])
time1=time.clock()
print ('time on generate the twod_like:',time1-time0)


def lnlike(theta, zl, y, err):
    om, w0, wa, h0 =theta
    #print theta
    likh=twod_like(theta, zl, y, err)[:,:-1]*ddz.T     #the sub of the intergral
    int_likh=np.sum(likh,axis=1)
    #print np.shape(int_likh)
    return np.sum(np.log(int_likh[int_likh!=0]))
#print lnlike(theta=[0.3,70], zl=data[:,1], y=data[:,3], err=data[:,4])

'''
###############minmize like#######################
import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
L1= [0.30, -1, 70]
#zs=data[:,0]
bnds = ((0, None), (0, None))
result = op.minimize(nll, (0.30, -1, 70), method='SLSQP', bounds=bnds,args=(zl, y, err))
print result["x"]
'''
######################MCMC#########################
def lnprior(theta):
    om, w0, wa, h0 = theta
    if 0.1 < om < 0.55 and -2 <w0< -0.5 and -2<wa<2 and 45 < h0 < 95:
        return 0.0
    return -np.inf

def lnprob(theta, zl, y, yerr):
    lp = lnprior(theta)
    step=np.sum(sampler.flatchain!=0)/(ndim*nwalkers)
    print "step:",step,"percent",round(step/5.,2),"%","\r",#"theta:",theta,
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, zl, y, err)

ndim, nwalkers = 4, 50
pos = [[0.30,-1, 0, 70] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(zl, y, err))  ###the input args() should change
sampler.run_mcmc(pos, 500)

###############save the data###############
samples = sampler.chain.reshape((-1, ndim))   #sampler.chain[:, 50:, :] to cut the data
burn=sampler.chain[:,:,:].T
import pickle
pickle.dump([sampler.chain,burn], open('int_CPL_500.txt','wb'))


samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
import corner
fig = corner.corner(samples, labels=["$om$", "$w0$", "wa", "$h0$"],
                      truths=[0.3, -1, 0, 70])
#fig.savefig("triangle.png")

plt.show()

