import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'
plt.figure(figsize=(13.5,10))
#NAME=input('Which file to load?:')
f=open('TD_result')
data=np.loadtxt(f)
f.close
#plt.hist(data[:,0],25)
#plt.show()
m=np.amax(data[:,1])
c=data[:,1]/data[:,0]
#cc=[0,0,c]
plt.scatter(data[:,0],data[:,3],c=data[:,1]/data[:,0],s=95,zorder=100)
cl=plt.colorbar()
cl.set_label('redshift ratio between deflector and source',rotation=270,size=20)
cl.ax.get_yaxis().labelpad=35     #the distance of the colorbar titel from bar
cl.ax.tick_params(labelsize=30)   #the labe size
plt.errorbar(data[:,0],data[:,3],yerr=data[:,4],fmt='y.',mew=0,markersize=13)
plt.tick_params(labelsize=35)
plt.title("Simulated TD ",fontsize=35)
plt.xlabel("Source redshift",fontsize=35)
plt.ylabel("TD distance",fontsize=35)
plt.show()
