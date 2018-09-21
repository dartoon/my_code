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
    p_va=np.trim_zeros(ps[:,1]*(zs>=ps[:,0]))[-1]   #zs can be any number now, the Ps will be respesented the nearby zs
    return p_va*(zs<30)
#print Ps(13)
vec_Ps=np.vectorize(Ps)

#######to get the P(zl|zs):##################

men=np.zeros([len(ps),2,501])  #501 is the tau scale
for i in range(len(ps)):    #creat a mem for fast cal
    men[i]=tau(ps[i,0])
def Pl(zl,zs):
    (a,b)=men[np.sum(zs>ps[:,0])-1]  #the possibility distributino tau(zs)
    #(a,b)=tau(zs)
    p_va=np.trim_zeros(b*(zl>=a))[-1]/(len(a)-1)  #the possibility for zl when zs is given
    return p_va
vec_Pl=np.vectorize(Pl)
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

def ct_zs(zl):                 #creat a zs array above to the zl, get the zs after zl in z_sca scale, the other z_s set to 30
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
def Ez(z,om):
  return   1/np.sqrt(om*(1+z)**3+(1-om))
def r(z,om):
  if z < 20:
    return integrate.quad(Ez, 0, z, args=(om))[0]    #use the cos distance r
  else:
    return 3.33
vec_r=np.vectorize(r)

z_s=ct_zs(zl)

################to modify the P(zs)##################
R=np.zeros([N+1,3])
R[:,0]=np.append(0,ps[:,0])
R[:,1]=np.append(0,ps[:,1])
R[0,2]=0
for i in range (1,N):
	R[i+1,2]=R[i,2]+R[i,1]*(R[i+1,0]-R[i,0])
#print R
modf=np.zeros(N_p)
for i in range (N_p):
   modf[i]=1/(1-R[:,2][R[:,0]>zl[i]][0])
#print modf


def twod_like(theta, zl, y, err):    #set zs to be 2D to get 2D sub_int, from (len(data)) to (len(data),133)
    om, h0 =theta
    rzl=vec_r(zl,om)[:,None]  #transfer rzl from (len(data),) to (len(data),1)
    #z_s=ct_zs(zl)
    rzs=vec_r(z_s,om)
    model = c*rzl*(rzs-rzl)/rzs/h0
    #print np.shape(rzl),np.shape(model)
    likh = np.exp(-0.5*(pow((y[:,None]-model),2)/pow(err[:,None],2)))*vec_Pl(zl[:,None],z_s)*vec_Ps(z_s)
    return likh
#print twod_like(theta=[0.3,70], zl=data[:,1], y=data[:,3], err=data[:,4]).shape

ddz=ps[:,0][1:]-ps[:,0][:-1] #get the difference between each redshift grid

#import time
#time0=time.clock()
#int_likh=np.sum(twod_like(theta=[0.3,70], zl=data[:,1], y=data[:,3], err=data[:,4])[:,:-1]*ddz.T,axis=1)
#print int_likh.shape,np.where(int_likh==0)[0],int_likh[int_likh!=0].shape
##print twod_like(theta=[0.3,70], zl=data[:,1], y=data[:,3], err=data[:,4])
#time1=time.clock()
#print ('time on generate the twod_like:',time1-time0)

def lnlike(theta, zl, y, err):
    om, h0 =theta
    #print theta
    likh=twod_like(theta, zl, y, err)[:,1:]*ddz.T     #the sub of the intergral
    int_likh=np.sum(likh[:,:],axis=1)*modf
    #print np.shape(int_likh)
    return np.sum(np.log(int_likh[int_likh!=0]))
#print lnlike(theta=[0.3,70], zl=data[:,1], y=data[:,3], err=data[:,4])

print lnlike(theta=[0.2,70], zl=data[:,1], y=data[:,3], err=data[:,4])
print lnlike(theta=[0.3,70], zl=data[:,1], y=data[:,3], err=data[:,4])
print lnlike(theta=[0.4,70], zl=data[:,1], y=data[:,3], err=data[:,4])
print lnlike(theta=[0.5,70], zl=data[:,1], y=data[:,3], err=data[:,4])

