from scipy.integrate import quad
import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import interpolate
###load data first
c=299790. # speed of light [km/s]
H = 70.
om = 0.3
w = -1.
dH=c/H
rho0 = 8.
r0 = 6709 # in [Mpc] value for DECIGO - polynomial Noise curve
s_type_list = ['NSNS', 'BHNS', 'BHBH']
sce_list = [None, 'standard high', 'standard low', 'optimistic low', 'optimistic high', 'delayed SN low', 'delayed SN high', 'high BH kicks low', 'high BH kicks high']
def cal_GW_event(DCO_type, scenario, pho_upper=80, pho_grid=1001, I = 1):
    s_type = s_type_list[DCO_type]
    file = open(s_type + 'rates.dat')
    gwdata = np.loadtxt(file)
    gwdata=gwdata[gwdata[:,0].argsort(),]
    scenario = scenario #0:z 1:snh snl onl onh dnl dnh knl knh
    Mchirp0 = [1.22,3.2, 6.7][DCO_type] #type of the binaries 1.2 3.2 6.7
    z=gwdata[:,0]
    def Ez(z,om,w):
        return   1/np.sqrt(om*(1+z)**3+(1-om)*(1+z)**(3*(1+w)))
    def r(z,om,w):
        return integrate.quad(Ez, 0, z, args=(om,w))[0]    #use the cos distance r
    vec_r=np.vectorize(r)
    Dist = vec_r(z,om,w)
    n0 = gwdata[:,scenario]*10**(-9) # 2 is standard low
    rev_E=1/np.sqrt(om*(1+z)**3 + (1-om)*(1+z)**(3*(1+w)))   
    #The Marek 1409 way to calculate the total number:
    #N=len(z)
    #Ctheta0 =np.zeros(N) 
    #x0=(rho0/8)*(1+z)**(1/6.)*dH*(Dist/r0)*(1.2/Mchirp0)**(5/6.)
    #for i in range (0,N):
    #    if x0[i]>=0 and x0[i]<=4:
    #        Ctheta0[i]=(1+x0[i])*(4-x0[i])**4/256.
    #    else:
    #        Ctheta0[i]=0        
    #s0dz=4*np.pi*(dH)**3.*(n0/(1+z))*Dist**2*rev_E*Ctheta0     
    ##Calcualte the non lensed events:
    #plt.plot(z, s0dz)
    #plt.close()
    #print(s_type, ": total:", sce_list[scenario], np.sum(s0dz[:-1]*(z[1:]-z[:-1]))) 
    
    #Decompose into 2D as introduced by Ding2015
    pho_grid = pho_grid
    pho_upper = pho_upper
    Ptheta = np.zeros(len(z))
    sDrhoDz = np.zeros([len(z), pho_grid])
    pho_break_idx = int(8 /(80/pho_grid))
    for j in range(pho_grid):
        rho = 80*(j+1)/pho_grid
        x=(rho/8)*(1+z)**(1/6.)*dH*(Dist/r0)*(1.2/Mchirp0)**(5/6.)
        for i in range(len(z)):
            if x[i]>=0 and x[i]<=4:
                Ptheta[i]=5*x[i]*(4-x[i])**3/256.
            else:
                Ptheta[i]=0
        sDrhoDz[:,j]=4*np.pi*(dH)**3*(n0/(1+z))*Dist**2*rev_E*Ptheta*x/rho
#    #Calcualte the non lensed events:
#    plt.plot(z, np.sum(sDrhoDz[:,pho_break_idx:]*pho_upper/pho_grid, axis=1) )
#    plt.close()
#    print(s_type, ": total:", sce_list[scenario], np.sum(np.sum(sDrhoDz[:-1,pho_break_idx:], axis=1) *(z[1:]-z[:-1]))*pho_upper/pho_grid  )  
    #Calculate lening statistics:
    phi0=8.0e-3*(H/100)**3
    sig0=161
    alpha=2.32
    bet=2.67
    SNR = 8
    TFactor = (sig0/c)**4*phi0*gamma((4.+alpha)/bet)/gamma(alpha/bet)
    tau = np.zeros([len(z), pho_grid])
    for i in range(pho_grid): 
        rho0 =pho_upper*(i+1)/pho_grid
        ymax = 1./((SNR/rho0)**2 + I * 1)   #two images to decide
        if ymax>1:  #for the u_+ case ymax can be large than 1 
            ymax=1
        tau[:,i]=(16/30)*np.pi**3*dH**3*TFactor*ymax**2.*Dist**3
    #############################int z then rho######################
    Nd_rh=np.zeros_like(tau)
    for i in range(len(z)-1):
        Nd_rh[i+1,:] = Nd_rh[i,:]+sDrhoDz[i+1]*tau[i+1,:]*(z[i+1]-z[i])   #for continuous, calculated form the right end.
#    print(s_type, ": total:", sce_list[scenario],np.sum(Nd_rh[-1,:])*pho_upper/pho_grid)
#    plt.plot(np.linspace(0,80,pho_grid), Nd_rh[-1,:]/(np.sum(Nd_rh[-1,:])), 'k--')
#    plt.close()
    return sDrhoDz, Nd_rh ,z
#%%
#sDrhoDz, Nd_rh, z = cal_GW_event(DCO_type = 0, scenario = 1)
sce_list = [None, 'standard high', 'standard low', 'optimistic low', 'optimistic high', 'delayed SN low', 'delayed SN high', 'high BH kicks low', 'high BH kicks high']
read_list = [[2, 3, 5, 7], [1, 4, 6, 8]]
table_value = []
pho_upper = 80
pho_grid=1001
pho_break_idx = int(8 /(80/1000))
#%%
print("Calculated values in Table 2:")
for i in range(3):
    print('\n',s_type_list[i],':')
    for j in range(2):
        print("============================")
        for k in range(len(read_list[j])):
            sDrhoDz, Nd_rh, z = cal_GW_event(DCO_type = i, scenario = read_list[j][k])
#            table_value.append(np.sum(sDrhoDz[:,pho_break_idx:]*pho_upper/pho_grid, axis=1))
            print(sce_list[read_list[j][k]], round(np.sum(np.sum(sDrhoDz[:-1,pho_break_idx:], axis=1) *(z[1:]-z[:-1]))*pho_upper/pho_grid, 1))  
#%%
print("Calculated values in Table 4 (consider continuous observation time):")
for i in range(3):
    print('\n',s_type_list[i],':')
    for j in range(2):
        print("============================")
        for k in range(len(read_list[j])):
            sDrhoDz, Nd_rh, z = cal_GW_event(DCO_type = i, scenario = read_list[j][k])
#            table_value.append(np.sum(sDrhoDz[:,pho_break_idx:]*pho_upper/pho_grid, axis=1))
            print(sce_list[read_list[j][k]], round(np.sum(Nd_rh[-1,:])*pho_upper/pho_grid,1 ) )

#%%
print("Calculated values in Table 6 (consider continuous observation time):")
for i in range(3):
    print('\n',s_type_list[i],':')
    for j in range(2):
        print("============================")
        for k in range(len(read_list[j])):
            sDrhoDz, Nd_rh, z = cal_GW_event(DCO_type = i, scenario = read_list[j][k])
#            table_value.append(np.sum(sDrhoDz[:,pho_break_idx:]*pho_upper/pho_grid, axis=1))
            print(sce_list[read_list[j][k]], round(np.sum(Nd_rh[-1,:pho_break_idx])*pho_upper/pho_grid,1 ))

#%%
print("Calculated values in Table 7 (consider continuous observation time):")
for i in range(3):
    print('\n',s_type_list[i],':')
    for j in range(2):
        print("============================")
        for k in range(len(read_list[j])):
            sDrhoDz, Nd_rh, z = cal_GW_event(DCO_type = i, scenario = read_list[j][k], I = -1)
            print(sce_list[read_list[j][k]], round(np.sum(Nd_rh[-1,:pho_break_idx])*pho_upper/pho_grid,1 ))
