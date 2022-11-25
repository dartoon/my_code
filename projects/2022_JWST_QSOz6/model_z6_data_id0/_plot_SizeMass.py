#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 16:58:52 2022

@author: Dartoon

cp from my CEERS work 4_plot_SizeMass_relation.py
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'
cmap = matplotlib.cm.get_cmap('viridis')
import matplotlib as mpl
mpl.rc('image', cmap='jet')
import glob

from scipy.integrate import quad
h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)

folder = '/Users/Dartoon/Astro/Projects/QSO_decomposition/' #'/Comparsion/'

import pickle
# pickle.dump(results, open('size_mass_CANDELS_result.pkl', 'wb'))
results = pickle.load(open('../../2022_CEERS_checkup/real_analysis/model_JWST_stage3_Masafusa/size_mass_CANDELS_result.pkl','rb'))
da_result = 1/(1+results[:,0])*c*vec_EE(results[:,0])/h0  #in Mpc

#%% Plot and input my sample

relation = 0  # 0 M*- Reff ; 1 M*-color; 2 M*-n; 3 Reff-n relation
# z_range = [1.2,1.7]
z_range = [1.6,3.5]
z_cut = ((results[:,0]>z_range[0]) * (results[:,0]<z_range[1]))

def logR_mstar(mstar, logA, alpha):
    """
    The R_Mstar relation in
    http://adsabs.harvard.edu/abs/2014ApJ...788...28V
    mstar in unit of log(M*/Msun)
    """
#    mstar = np.log10(10**mstar/7.*10**10.)
#    logr = logA + alpha*(mstar) + 10*alpha + alpha*np.log10(1/7.)
    r = 10**logA*(10**mstar/(5.*10**10))**alpha
    logr = np.log10(r)
    return logr


#ssfr_break =-10.5
#blue_galaxy = ([results[:,5]>ssfr_break])[0]
#red_galaxy = ([results[:,5]<ssfr_break])[0]
all_galaxy = ([results[:,5]<100])[0]  #As all galaxy

#Define based on Color:
red_galaxy = []
blue_galaxy = []

all_color = []
for i in range(len(results)):
    all_color.append([results[:,9][i][0], results[:,9][i][1], results[:,9][i][2]])
all_color = np.asarray(all_color)
ColorUV = -(2.5* np.log10(all_color[:,0])-2.5* np.log10(all_color[:,1]))
ColorVJ = -(2.5* np.log10(all_color[:,1])-2.5* np.log10(all_color[:,2]))
for i in range(len(all_color)):
    if ColorVJ[i] < 0.8425:
        blue_galaxy.append((ColorUV[i] < 1.286))
        red_galaxy.append((ColorUV[i] > 1.286))
    else:
        k = 1.17
        b = 0.3
        line_p = k*ColorVJ[i]+b
        blue_galaxy.append((ColorUV[i] -line_p < 0))        
        red_galaxy.append((ColorUV[i] -line_p > 0))
blue_galaxy = np.asarray(blue_galaxy)
red_galaxy = np.asarray(red_galaxy)
#%%
#import scipy.optimize as opt
#def lfit(x,m,c):
#    return m*x+c
cmap_r = matplotlib.cm.get_cmap('RdBu_r')

fig, ax = plt.subplots(figsize=(14, 11))
Reff_kpc = da_result * 10 **3 * (results[:,1]/3600./180.*np.pi)
Reff_kpc = Reff_kpc.astype(float)

Mstar_candels = results[:,3]
Mstar_candels = Mstar_candels.astype(float)

mstar_cut = [(results[:,3]>9.5) * (results[:,3]<11.5)][0]
if relation == 0:
#    plt.scatter(results[:,3][z_cut* all_galaxy],np.log10(Reff_kpc[z_cut * all_galaxy]),
#                c=(results[:,7][z_cut* all_galaxy]/results[:,8][z_cut * all_galaxy]),s=280,marker=".",zorder=90, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)
#    mstar_line = np.linspace(10.5,11.5,20)
#    plt.plot(mstar_line, logR_mstar(mstar_line,logA=0.155 , alpha=0.76), 'r')
#    mstar_line = np.linspace(9,11.5,20)
#    plt.plot(mstar_line, logR_mstar(mstar_line,logA=0.675 , alpha=0.23), 'b')
    c1 = plt.scatter(Mstar_candels[z_cut* blue_galaxy],Reff_kpc[z_cut * blue_galaxy],   #!!!Dont need to care about the mstar_cut as would show in the figure.
                c='lightskyblue',s=280,marker=".",zorder=-90, alpha=0.6, edgecolors='white', cmap=cmap_r, label='CANDELS galaxy, star-forming')
    c2 = plt.scatter(Mstar_candels[z_cut* red_galaxy],Reff_kpc[z_cut * red_galaxy],
                c='darksalmon',s=280,marker=".",zorder=-90, alpha=0.6, edgecolors='white', cmap=cmap_r, label='CANDELS galaxy, quiescent')
    
    mstar_line = np.linspace(10.5,11.5,20)
#    m_cut = [(results[:,3]>10) * (results[:,3]<11.5)][0]
#    fit_red = opt.curve_fit(lfit, results[:,3][z_cut* red_galaxy * m_cut],np.log10(Reff_kpc[z_cut * red_galaxy * m_cut]))
#    plt.plot(mstar_line, 10**(lfit(mstar_line, fit_red[0][0],fit_red[0][1])), 'r',linewidth=3)
    # plt.plot(mstar_line, 10**(logR_mstar(mstar_line,logA=0.155 , alpha=0.76)), 'r',linewidth=3) #At z ~1.5
    plt.plot(mstar_line, 10**(logR_mstar(mstar_line,logA=-0.05 , alpha=0.76)), 'r',linewidth=3) #At z~2.5
    
    mstar_line = np.linspace(9,11.5,20)
    # plt.plot(mstar_line, 10**(logR_mstar(mstar_line,logA=0.675 , alpha=0.23)), 'b',linewidth=3)   #At z ~1.5
    plt.plot(mstar_line, 10**(logR_mstar(mstar_line,logA=0.55 , alpha=0.23)), 'b',linewidth=3)     #At z~2.5
    
    mstar_line = np.linspace(9,11.5,20)
#    plt.plot(mstar_line, 10**(0.54+ 0.57*(mstar_line-11.)), 'r--',linewidth=2,alpha=0.8)
#    plt.text(9.56, 10**(-0.34), 'z = 0.06', color='red', fontsize=35)    
elif relation == 1:
    plt.scatter(results[:,3][z_cut* all_galaxy],(results[:,7][z_cut* all_galaxy]/results[:,8][z_cut * all_galaxy]),
                c=(results[:,7][z_cut* all_galaxy]/results[:,8][z_cut * all_galaxy]),s=280,marker=".",zorder=90, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)
elif relation == 2:
    plt.scatter(results[:,3][z_cut* all_galaxy],results[:,2][z_cut* all_galaxy],
                c=(results[:,7][z_cut* all_galaxy]/results[:,8][z_cut * all_galaxy]),s=280,marker=".",zorder=90, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)
elif relation == 3:
    plt.scatter(np.log10(Reff_kpc[z_cut * all_galaxy * mstar_cut]),results[:,2][z_cut* all_galaxy* mstar_cut],
                c=(results[:,7][z_cut* all_galaxy* mstar_cut]/results[:,8][z_cut * all_galaxy* mstar_cut]),s=280,marker=".",zorder=90, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)
elif relation == 4:
    plt.scatter(results[:,0][z_cut* all_galaxy* mstar_cut], np.log10(Reff_kpc[z_cut * all_galaxy * mstar_cut]),
                c=(results[:,7][z_cut* all_galaxy* mstar_cut]/results[:,8][z_cut * all_galaxy* mstar_cut]),s=280,marker=".",zorder=90, vmin=0, vmax=7, alpha=0.6, edgecolors='white', cmap=cmap_r)
    z_line = np.linspace(1,2,21)
    plt.plot(z_line, np.log10(8.9*(1+z_line)**(-0.75)), 'b')
    plt.plot(z_line, np.log10(5.6*(1+z_line)**(-1.48)), 'r')
    
import sys
sys.path.insert(0,folder+'/py_tools')
from load_result import load_zs, load_mag, load_re, load_n, load_flux
ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281', 'CID255']
ALMA_ID = ['COSMOS-CID70','COSMOS-CID206','COSMOS-CID454','COSMOS-CID237','COSMOS-LID1273','COSMOS-XID2396','COSMOS-CID1174','COSMOS-XID2138','COSMOS-LID1538','COSMOS-XID2202','COSMOS-CID3242','COSMOS-CID543 ','COSMOS-CID607 ','COSMOS-CID597 ','COSMOS-CID3570','COSMOS-CID50  ','COSMOS-CID255 ','COSMOS-LID360 ','COSMOS-CID216 ','COSMOS-CID452 ']
ALMA_ID_str = 'COSMOS-CID70 COSMOS-CID206 COSMOS-CID454 COSMOS-CID237 COSMOS-LID1273 COSMOS-XID2396 COSMOS-CID1174 COSMOS-XID2138 COSMOS-LID1538 COSMOS-XID2202 COSMOS-CID3242 COSMOS-CID543  COSMOS-CID607  COSMOS-CID597  COSMOS-CID3570 COSMOS-CID50   COSMOS-CID255  COSMOS-LID360  COSMOS-CID216  COSMOS-CID452 '
in_list_ = [i for i in range(len(ID)) if ID[i] in ALMA_ID_str]
in_list = []
out_list = []
for i in range(len(ID)):
    if i in in_list_:
        in_list.append(True)
        out_list.append(False)
    else:
        in_list.append(False)
        out_list.append(True)
in_list = np.array(in_list)
out_list = np.array(out_list)
    
zs = np.asarray(load_zs(ID))
mags = np.array(load_mag(ID, folder = folder)[0])
Reffs = np.array(load_re(ID, folder = folder))[:,0]
Reffs_e = np.array(load_re(ID, folder = folder))[:,1]
indexs = np.array(load_n(ID, folder = folder))[:,0]
from load_result import load_host_p
Mstar = load_host_p(ID, folder = folder)[1]

for i in range(len(zs)):
    zs_ = zs[i]
    if zs_> 1.44:
        filt_wave = 14000
        key = 140
    else:
        filt_wave = 12500
        key = 125
    wave_rest = 5000
    wave = filt_wave / (1+zs_)
    z_p = int(key)/50  -1
    Reff_filt = Reffs[i]
    vdep = -0.35 + 0.12 * zs_ - 0.25 * np.log10( 10 ** Mstar[i] / 10**10)
    Reff_5000 = Reff_filt * ((1+zs_)/(1+z_p)) ** vdep
    Reffs[i] = Reff_5000


dl=(1+zs)*c*vec_EE(zs)/h0 *10**6   #in pc
da=1/(1+zs)*c*vec_EE(zs)/h0   #in Mpc
ID_Reff_kpc = da * 10 **3 * (Reffs/3600./180.*np.pi)
ID_Reff_kpc_e = da * 10 **3 * ((Reffs_e)/3600./180.*np.pi)

host_flux_WFC3 = np.array(load_flux(ID, folder = folder, flt = 'WFC3'))[:,0]
host_flux_ACS = []
for i in range(len(ID)):
    ifexit = glob.glob(folder + '/analysis_ACS/{0}'.format(ID[i]))
    if ifexit!= []:
        host_flux_ACS.append(load_flux([ID[i]], folder = folder, flt = 'ACS')[0][0])
    else:
        host_flux_ACS.append(-99)
host_flux_ACS = np.asarray(host_flux_ACS)

#cl=plt.colorbar()          #cl take the inforamtion from the lastest plt
#cl.set_label('filter flux ratio, WFC3 / ACS',rotation=270,size=30)
#cl.ax.get_yaxis().labelpad=35     #the distance of the colorbar titel from bar
#cl.ax.tick_params(labelsize=30)
f1 =folder + '/M_BH_relation/data/Bennert+2011_local.txt'
b11_l = np.loadtxt(f1)[:,1:]  #0 redshift; 1 M*; 2 BH mass;
b11_local_Reff = b11_l[:,8]
#    b11_local_mstar = b11_l[:,4]
b11_local_mstar = b11_l[:,9]  #!!! Change to total mass
#    plt.scatter(b11_local_mstar,b11_local_Reff,s=180, c ='black',
#                marker="o",zorder=100, vmin=0.5, vmax=5, edgecolors='white', label='local AGN (VB2011)')     
for i in range(len(Mstar)):
    if Reffs[i]-0.1 < 0.009 and zs[i]>1.2:
        plt.arrow(Mstar[i], ID_Reff_kpc[i], 0, -0.3, length_includes_head=True,
              head_width=0.08, head_length=0.05, zorder=102, color='black', linewidth=1.2)
   
log_Rerr = (np.log10(ID_Reff_kpc)-np.log10(ID_Reff_kpc-ID_Reff_kpc_e))
low_err = ID_Reff_kpc - 10**(np.log10(ID_Reff_kpc)-log_Rerr)
up_err = 10**(np.log10(ID_Reff_kpc)+log_Rerr) - ID_Reff_kpc

p2 = plt.scatter(Mstar[zs>1.2],ID_Reff_kpc[zs>1.2], s=170, linewidth=1.2, c =indexs[zs>1.2],
            marker="D",zorder=101, vmin=0.5, vmax=5, edgecolors='black', cmap=cmap_r)    
#%%
#Load M for 
JWST_z = np.array([1.646, 2.588, 3.442, 2.317, 3.465])
JWST_Reff = np.array([0.58063765, 0.0675265 , 0.16146611, 0.1249958 , 0.163813  ])
JWST_Reff_err = np.array([0.05988698, 0.04920861, 0.07036986, 0.02963798, 0.07025222])
da=1/(1+JWST_z)*c*vec_EE(JWST_z)/h0   #in Mpc
JWST_Reff_kpc = da * 10 **3 * (JWST_Reff/3600./180.*np.pi)
JWST_Reff_kpc_err = da * 10 **3 * (JWST_Reff_err/3600./180.*np.pi)
JWST_smass = np.array([11.255, 10.532, 10.07 , 10.191,  9.186])
JWST_smass_low = np.array([0.2  , 0.264, 0.449, 0.151, 0.125])
JWST_smass_high = np.array([0.13 , 0.154, 0.339, 0.153, 0.216])
JWST_n = np.array([2.07807855, 5.60323111, 0.96118981, 1.56714238, 1.41940508])

show = np.array([True, True, True, True, True])
jwst_p = p_jwst = plt.scatter(JWST_smass[show], JWST_Reff_kpc[show], s=480, linewidth=2.2, c =JWST_n[show],
            marker="H",zorder=201, vmin=0.5, vmax=5, edgecolors='black', cmap=cmap_r)   

show = np.array([True, True, True, True, True])
JWST_Reff_kpc_err[1] = 0.001
plt.errorbar(JWST_smass[show], JWST_Reff_kpc[show], 
             xerr=[JWST_smass_low[show],JWST_smass_high[show]], yerr=JWST_Reff_kpc_err[show], color='black', linewidth=3.2, 
             ecolor='black', fmt='.',zorder=200,markersize=1)


plt.arrow(JWST_smass[1], JWST_Reff_kpc[1], 0, -0.2, length_includes_head=True,
      head_width=0.08, head_length=0.05, zorder=102, color='black', linewidth=1.2)

# name_list = {1: 'SDSS1420A',
#  2: 'SDSS1420B',
#  0: 'SDSS1419',
#  51: 'AEGIS 477',
#  35: 'AEGIS 482'}
# idx_list = [1, 2, 0, 51, 35]
# for i,idx in enumerate(idx_list): 
#     idx = idx_list[i]
#     if idx != 51:
#         plt.text(JWST_smass[i]-0.5, JWST_Reff_kpc[i]*1.2, name_list[idx], fontsize = 21, zorder =202,
#                  bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 3})   
#     elif idx == 51:
#         plt.text(JWST_smass[i]+0.1, JWST_Reff_kpc[i]/1.2, name_list[idx], fontsize = 21, zorder =202,
#                  bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 3})    
#%%

#Target 0
smass = 10.63
R_filt_0 = np.array([1.87, 1.06 ])
q = 0.630
R_filt_0 = R_filt_0#*np.sqrt(q)
# wave = F356W
z = 6.34
key = 'F356W'
filt_wave_dic = {'F115W': 11623.89, 'F150W': 15104.23, 'F200W': 20028.15, 
           'F277W': 27844.64, 'F356W': 35934.49, 'F444W': 44393.50, 'F410M': 40886.55}
filt_wave = filt_wave_dic[key]
wave_rest = 5000
wave = filt_wave / (1+z)
z_p = int(key[1:4])/(wave_rest/100)  -1
Reff_filt = R_filt_0[0]
vdep = -0.35 + 0.12 * z - 0.25 * np.log10( 10 ** smass / 10**10)
R_0 = [Reff_filt * ((1+z)/(1+z_p)) ** vdep]
R_0.append(R_0[0]* R_filt_0[1]/R_filt_0[0] )
plt.scatter(smass, R_0[0],c=1.61,s=1280,marker="*",
            zorder=101, alpha=1, label = 'J2255+0251',  vmin=0.5, vmax=5, edgecolors='black', cmap=cmap_r)
plt.errorbar(smass, R_0[0], 
             yerr = [ R_0[1] ],
              # xerr = [ [smass - Muv_0[1]], [Muv_0[2] - smass] ],
              xerr = [ [0.6], [0.6] ],zorder = 100,
              c='black', elinewidth=3)

plt.text(smass+0., R_0[0]*1.2, 'J2255+0251', fontsize = 21, zorder =202,
                  bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 3})   

#Target 1
smass = 11.14
key = 'F150W'
# key = 'F356W'
if key == 'F356W':
    R_filt_0 = np.array([0.75, 0.11])
    q = 0.360
if key == 'F150W':
    R_filt_0 = np.array([0.59, 0.14])
    q = 0.269
R_filt_0 = R_filt_0 #*np.sqrt(q)
z = 6.4
filt_wave = filt_wave_dic[key]
wave = filt_wave / (1+z)
z_p = int(key[1:4])/(wave_rest/100)  -1
Reff_filt = R_filt_0[0]
vdep = -0.35 + 0.12 * z - 0.25 * np.log10( 10 ** smass / 10**10)
R_0 = [Reff_filt * ((1+z)/(1+z_p)) ** vdep]
R_0.append(R_0[0]* R_filt_0[1]/R_filt_0[0] )
plt.scatter(smass, R_0[0],c=1,s=1280,marker="*",
            zorder=101, alpha=1, label = 'J2236+0032',  vmin=0.5, vmax=5, edgecolors='black', cmap=cmap_r)
plt.errorbar(smass, R_0[0], 
             yerr = [ R_0[1] ],
              # xerr = [ [smass - Muv_0[1]], [Muv_0[2] - smass] ],
              xerr = [ [0.6], [0.6] ],zorder = 100,
              c='black', elinewidth=3)

plt.text(smass+0., R_0[0]*1.2, 'J2236+0032', fontsize = 21, zorder =202,
                  bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 3})  

#%%    
plt.xlim([8.5, 11.7])
plt.xlabel("log (M$_*$; units of M$_{\odot}$)",fontsize=35)
plt.tick_params(labelsize=25)
#plt.legend(loc='upper right',fontsize=21,numpoints=1)
from matplotlib.legend_handler import HandlerTuple
plt.legend([c1, c2, p2, jwst_p], ['CANDELS galaxy, star-forming', 'CANDELS galaxy, quiescent', 
                                  'Ding 2020 quasar 1.2<z<1.7', 'Ding 2022 quasar 1.6<z<3.5'],
               # handler_map={tuple: HandlerTuple(ndivide=None)},
               loc='upper left',fontsize=21,numpoints=1)
plt.ylabel(r"R$_{\rm eff, maj}$ (kpc)",fontsize=35)
plt.ylim([0.3, 20.5])
plt.yscale('log')
ax.tick_params(axis='both', which='major', length=12 , width = 2)
ax.tick_params(axis='y', which='minor', length=7, width=2)
ax.tick_params(axis='y', which='major', length=12, width=2)
ax.tick_params(axis='both', which='both', direction='in')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_ylabel('Sersic index', rotation=270, fontsize = 25, labelpad=25)
cbar.ax.tick_params(axis='both', which='both', direction='in')
cbar.ax.tick_params(axis='y', which='major', length=7, width=2)
plt.show()

