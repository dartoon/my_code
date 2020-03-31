import numpy as np
import matplotlib.pyplot as plt
from tau import *
'''
# The import of z_s Possibity dis
num=1000
a=sorted(np.random.uniform(0,3,num))
p=np.zeros([num,2])
p[:,0]=np.asarray(a)
p[:,1]=p[:,0]**2
p[:,1]/=p[:,1].sum()
#print p[:,1].sum()
#plt.plot(p[:,0],p[:,1],'k')
#plt.show()
########start to generate the data
'''
f=open('int_BHBH_sl.txt')
data=np.loadtxt(f)
f.close
#print len(data)
num=len(data)
N=num
p=np.zeros([num,2])
p[:,0]=data[:,0]
p[:,1]=data[:,1]   #not the real possibility should mutply by p(i)-p(i-1)
#plt.plot(p[:,0],p[:,1],'k')
#plt.show()
R=np.zeros([N+1,3])
R[:,0]=np.append(0,p[:,0])
R[:,1]=np.append(0,p[:,1])
R[0,2]=0
for i in range (1,N):
	R[i+1,2]=R[i,2]+R[i,1]*(R[i+1,0]-R[i,0])
#print R
d= input('How may lens time delay redshift data set to be generated?:')
da=np.zeros(d)
for i in range(d):
   k=np.random.random()
   #print k
   for j in range(1,N+1):
	if R[j-1,2]<k and R[j,2]>k:
		da[i]=np.random.uniform(R[j-1,0],R[j,0])

plt.hist(da,25)
plt.show()
#print da
######################################
##### for the zl with tau information
#####################################

####caculate the tau possiblity distribution and get a random zl from zs
daa=np.zeros(d)
N=len(np.asarray(tau(2)).T)
#print N
for i in range(len(da)):
  zs=da[i]         #input zs 
  Rl=np.zeros([N+1,3])
  (a,b)=tau(zs)
  Rl[:,0]=np.append(0,a)
  Rl[:,1]=np.append(0,b)
  Rl[0,2]=0
  for j in range(1,N):
	Rl[j+1,2]=Rl[j,2]+Rl[j,1]    
  k=np.random.random()       
  for j in range(1,N+1):
	if Rl[j-1,2]<k and Rl[j,2]>k:
	  daa[i]=np.random.uniform(Rl[j-1,0],Rl[j,0])
z=np.zeros([d,2])
z[:,0]=da
z[:,1]=daa
g=open('zs_zl','w')
np.savetxt(g,z,fmt='%6.5f')
g.close
#print z
#plt.hist(da,25)
#plt.show()

