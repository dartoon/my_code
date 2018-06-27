import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from pz import pz

########start to generate the data

#f=open('int_dis.txt')
#data=np.loadtxt(f)
#f.close
'''
num=len(data)
N=num
#print N
dp=np.zeros([num,4])
dp[:,0]=data[:,0]   # redshift assumption
dp[:,1]=data[:,1]   #for NS-NS note the real possibility should mutply by p(i)-p(i-1)
dp[:,2]=data[:,9]   #for NS-BH
dp[:,3]=data[:,17]   #for BH-BH
p=np.zeros([num,2])
p[:,0]=dp[:,0]
p[:,1]=dp[:,1]+dp[:,2]+dp[:,3]
'''
p=pz(0.3,70)
N=len(p[:,0])
R=np.zeros([N+1,3])
R[:,0]=np.append(0,p[:,0])
R[:,1]=np.append(0,p[:,1])
R[0,2]=0
for i in range (1,N):
	R[i+1,2]=R[i,2]+R[i,1]*(R[i+1,0]-R[i,0])
R[:,1]=R[:,1]/R[N,2]
R[:,2]=R[:,2]/R[N,2]
#print p
d= input('How may dataset to be generated?:')
#d=2000
zs=np.zeros(d)
seed=input('seed?:')
np.random.seed(seed)        #keep the random number to be the same every time
for i in range(d):
   k=np.random.random()
   #print k
   for j in range(1,N):
	if R[j,2]<k and R[j+1,2]>k:
		#zs[i]=np.random.uniform(R[j-1,0],R[j,0])     # should better not in this way
                zs[i]=R[j,0]                                  # the p_i(s) is for p_i to p_(i+1)


#plt.hist(zs,50)
#plt.show()
########################generate the luminorsity distance##############
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)
H=70             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]
dl=(1+zs)*c*vec_EE(zs)/H            # unit in Mpc 


'''
error= input('Which error level to set? in % unit:')
erbar=dl*error/100. 
'''

erbar=dl*1/100.   #1% level as the based
#erbar=((2*dl/Rho)**2+(0.05*zs*dl)**2)**0.5
#print np.mean(erbar/dl)
                  #error bar level
bias=np.random.normal(0,erbar)   ##include system error to the data, still with seed(1), normal and random is different type, it's linear to dl*error/100
#dl_nois=dl+bias
DL=np.zeros([len(zs),5])
DL[:,0]=zs            #zs
DL[:,1]=dl            #true dl distance 
#DL[:,2]=dl_nois       #after noise
#DL[:,3]=erbar         #error
#print DT
#plt.errorbar(DL[:,0],DL[:,2],yerr=DL[:,3],fmt='k.',mew=0,markersize=13)
#plt.show()
'''
DL[:,3]=erbar*1         #error
dl_nois=dl+bias*1
DL[:,2]=dl_nois       #after noise
g=open('TL_sim_1','w')
np.savetxt(g,DL,fmt='%6.5f')
g.close

DL[:,3]=erbar*3         #error
dl_nois=dl+bias*3
DL[:,2]=dl_nois       #after noise
g=open('TL_sim_3','w')
np.savetxt(g,DL,fmt='%6.5f')
g.close

DL[:,3]=erbar*5         #error
dl_nois=dl+bias*5
DL[:,2]=dl_nois       #after noise
g=open('TL_sim_5','w')
np.savetxt(g,DL,fmt='%6.5f')
g.close

DL[:,3]=erbar*10         #error
dl_nois=dl+bias*10
DL[:,2]=dl_nois       #after noise
g=open('TL_sim_10','w')
np.savetxt(g,DL,fmt='%6.5f')
g.close

DL[:,3]=erbar*15         #error
dl_nois=dl+bias*15
DL[:,2]=dl_nois       #after noise
g=open('TL_sim_15','w')
np.savetxt(g,DL,fmt='%6.5f')
g.close

DL[:,3]=erbar*25         #error
dl_nois=dl+bias*25
DL[:,2]=dl_nois       #after noise
g=open('TL_sim_25','w')
np.savetxt(g,DL,fmt='%6.5f')
g.close
'''
#DL[:,3]=erbar*10         #error
#dl_nois=dl+bias*10
#DL[:,2]=dl_nois       #after noise
#g=open('TL_sim_NO_100000','w')
#np.savetxt(g,DL,fmt='%6.5f')
#g.close
