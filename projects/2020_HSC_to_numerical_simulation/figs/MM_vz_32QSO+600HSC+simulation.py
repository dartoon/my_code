#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 12:54:19 2020

@author: Dartoon
"""
import numpy as np
np.set_printoptions(precision=4)
import matplotlib.pyplot as plt
import matplotlib as mat
import matplotlib.lines as mlines
from matplotlib import colors
mat.rcParams['font.family'] = 'STIXGeneral'
host=plt.figure(figsize=(14.5,12))
ax=host.add_subplot(111)   #to get the log(1+z) and z label

import matplotlib as mpl
mpl.rc('image', cmap='jet')
# import sys
########## input local data ####
#==============================================================================
# The seleting for dm and host_total and dmag are in this local
#==============================================================================
########input 25 local by Bennert++2011 ############
bloc = np.array([[ 0.054 , 10.12  ,  0.24  ,  7.436 ,  0.4   ],
       [ 0.043 , 10.95  ,  0.23  ,  8.006 ,  0.4   ],
       [ 0.076 , 10.33  ,  0.22  ,  7.436 ,  0.4   ],
       [ 0.041 , 10.38  ,  0.23  ,  7.166 ,  0.4   ],
       [ 0.051 , 10.5   ,  0.23  ,  8.516 ,  0.4   ],
       [ 0.0524, 10.32  ,  0.23  ,  6.976 ,  0.4   ],
       [ 0.0475,  9.83  ,  0.24  ,  7.666 ,  0.4   ],
       [ 0.055 , 10.41  ,  0.23  ,  7.836 ,  0.4   ],
       [ 0.0355, 10.33  ,  0.22  ,  7.766 ,  0.4   ],
       [ 0.021 , 10.2   ,  0.22  ,  7.366 ,  0.4   ],
       [ 0.038 , 10.26  ,  0.22  ,  7.356 ,  0.4   ],
       [ 0.0229,  9.94  ,  0.24  ,  7.496 ,  0.4   ],
       [ 0.047 , 10.14  ,  0.22  ,  7.646 ,  0.4   ],
       [ 0.0559,  9.65  ,  0.22  ,  6.976 ,  0.4   ],
       [ 0.0501, 10.11  ,  0.23  ,  7.956 ,  0.4   ],
       [ 0.0541, 10.04  ,  0.23  ,  7.076 ,  0.4   ],
       [ 0.0558, 10.73  ,  0.21  ,  7.506 ,  0.4   ],
       [ 0.0365, 10.3   ,  0.24  ,  7.476 ,  0.4   ],
       [ 0.0304, 10.24  ,  0.24  ,  7.826 ,  0.4   ],
       [ 0.0481,  9.92  ,  0.22  ,  7.276 ,  0.4   ],
       [ 0.0483, 10.    ,  0.23  ,  7.626 ,  0.4   ],
       [ 0.0465,  9.82  ,  0.23  ,  7.546 ,  0.4   ],
       [ 0.0532,  9.95  ,  0.23  ,  7.776 ,  0.4   ],
       [ 0.0585, 10.33  ,  0.22  ,  7.216 ,  0.4   ],
       [ 0.0409, 10.33  ,  0.22  ,  7.34  ,  0.4   ]])
########input 30 local by Haring 04 ############
hloc = np.array([[3.7484e-03, 1.1778e+01, 1.8000e-01, 9.4771e+00, 2.1298e-01],
       [3.4930e-03, 1.0362e+01, 1.8000e-01, 7.1461e+00, 1.4262e-01],
       [2.4703e-03, 1.0833e+01, 1.8000e-01, 8.0000e+00, 1.0206e-01],
       [4.2822e-03, 1.1556e+01, 1.8000e-01, 8.6335e+00, 2.1306e-01],
       [7.3369e-03, 1.1556e+01, 1.8000e-01, 8.7160e+00, 3.2389e-01],
       [2.4297e-02, 1.1748e+01, 1.8000e-01, 8.7243e+00, 3.1358e-01],
       [1.3564e-02, 1.1462e+01, 1.8000e-01, 8.5185e+00, 2.0838e-01],
       [3.6091e-03, 9.7924e+00, 1.8000e-01, 7.1461e+00, 1.8221e-01],
       [5.6030e-03, 1.1114e+01, 1.8000e-01, 7.5682e+00, 1.6526e-01],
       [6.7825e-03, 1.1462e+01, 1.8000e-01, 9.3979e+00, 3.9967e-01],
       [1.7743e-04, 1.0568e+01, 1.8000e-01, 7.6532e+00, 1.8656e-01],
       [1.8910e-04, 8.9031e+00, 1.8000e-01, 6.3979e+00, 3.4062e-01],
       [2.6564e-03, 1.0839e+01, 1.8000e-01, 7.6435e+00, 4.6942e-01],
       [5.3251e-03, 1.0881e+01, 1.8000e-01, 7.1461e+00, 2.2797e-01],
       [2.2610e-03, 1.1079e+01, 1.8000e-01, 9.0000e+00, 2.3856e-01],
       [4.8618e-03, 1.0833e+01, 1.8000e-01, 8.3222e+00, 2.8354e-01],
       [2.6099e-03, 1.0204e+01, 1.8000e-01, 8.0000e+00, 6.1650e-01],
       [2.7029e-03, 1.0301e+01, 1.8000e-01, 7.2041e+00, 4.3571e-01],
       [5.3483e-03, 1.0987e+01, 1.8000e-01, 8.2788e+00, 2.5972e-01],
       [6.0889e-03, 1.1114e+01, 1.8000e-01, 8.4914e+00, 2.0586e-01],
       [3.5627e-03, 1.0079e+01, 1.8000e-01, 8.4771e+00, 2.4800e-01],
       [3.6556e-03, 1.0964e+01, 1.8000e-01, 8.0414e+00, 5.6735e-01],
       [3.4930e-03, 1.0643e+01, 1.8000e-01, 7.7482e+00, 4.3388e-01],
       [2.2842e-03, 1.1431e+01, 1.8000e-01, 9.0000e+00, 3.3450e-01],
       [3.9109e-03, 1.1690e+01, 1.8000e-01, 9.3010e+00, 4.8455e-02],
       [2.7262e-03, 1.1041e+01, 1.8000e-01, 8.2304e+00, 6.2621e-01],
       [6.0195e-03, 1.0568e+01, 1.8000e-01, 8.3802e+00, 1.0654e-01],
       [5.3483e-03, 1.0176e+01, 1.8000e-01, 7.1139e+00, 1.8447e-01],
       [3.0748e-03, 9.8451e+00, 1.8000e-01, 6.5441e+00, 1.4739e-01],
       [2.3349e-05, 1.0041e+01, 1.8000e-01, 6.5682e+00, 1.5707e-01]])

zs =  np.array([3.7484e-03, 3.4930e-03, 2.4703e-03, 4.2822e-03, 7.3369e-03,
       2.4297e-02, 1.3564e-02, 3.6091e-03, 5.6030e-03, 6.7825e-03,
       1.7743e-04, 1.8910e-04, 2.6564e-03, 5.3251e-03, 2.2610e-03,
       4.8618e-03, 2.6099e-03, 2.7029e-03, 5.3483e-03, 6.0889e-03,
       3.5627e-03, 3.6556e-03, 3.4930e-03, 2.2842e-03, 3.9109e-03,
       2.7262e-03, 6.0195e-03, 5.3483e-03, 3.0748e-03, 2.3349e-05])
hloc[:,0] = zs

#############################################################
###################fitting with MCMC#########################
x=np.append(bloc[:,1], hloc[:,1])
y=np.append(bloc[:,3], hloc[:,3])
yerr=(np.append(bloc[:,2], hloc[:,2])**2+np.append(bloc[:,4], hloc[:,4])**2)**0.5  # 0.2 is the uncertainty level for the L_R
def _lnlike(theta, x, y, yerr):
    m, b, sint= theta
    model = m * x + b
    sigma2 = (yerr**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf

import scipy.optimize as op
nll = lambda *args: -_lnlike(*args)
result = op.minimize(nll, [1.036, -1.947, 0.3], args=(x, y, yerr))
m_ml, b_ml,sint_ml= result["x"]
def _lnprior(theta):
    m, b, sint	 = theta
    if -5.0 < m < 5 and -10 < b < 10.0 and 0 < sint < 10:
        return 0.0
    return -np.inf
def _lnprob(theta, x, y, yerr):
    lp = _lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + _lnlike(theta, x, y, yerr)
ndim, nwalkers = 3, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, _lnprob, args=(x, y, yerr))
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

m_mid, b_mid, sint_mid =np.percentile(samples, 50,axis=0)

######################
xl = np.linspace(-0.9, 13, 100)
plt.errorbar(np.log10(hloc[:,0]+1),
             hloc[:,3]-(m_ml*hloc[:,1]+b_ml),yerr=(hloc[:,2]**2 + hloc[:,4]**2)**0.5 ,fmt='.',color='gray',markersize=10)
plt.errorbar(np.log10(bloc[:,0]+1),
             bloc[:,3]-(m_ml*bloc[:,1]+b_ml),yerr=(bloc[:,2]**2 + bloc[:,4]**2)**0.5 ,fmt='.',color='gray',markersize=10)
ty1=xl*0+np.std(y-(m_ml*x+b_ml))
ty2=xl*0-np.std(y-(m_ml*x+b_ml))
plt.fill_between(xl,ty1,ty2,color='lightgray',zorder=-50, alpha = 0.4)

Bkc=mlines.Line2D([], [], color='gray', ls='', marker='.', markersize=15)
# Hkc=mlines.Line2D([], [], color='gray', ls='', marker='.', markersize=15)
######################
#%%
#==============================================================================
# My new inference
#==============================================================================
# from load_result import load_host_p, load_MBH, load_err
# from load_result import load_zs, load_n
ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']
MB_ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'ECDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','LID1820','LID1622',\
'LID1878', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']
zs = np.array([1.63 , 1.301, 1.667, 1.447, 1.326, 1.57 , 1.552, 1.567, 1.618,
       1.532, 1.244, 1.407, 1.478, 1.239, 1.294, 1.617, 1.527, 1.579,
       1.325, 1.411, 1.276, 1.412, 1.585, 1.551, 1.516, 1.6  , 1.483,
       1.626, 1.337, 1.272, 1.445, 1.664])
host_n = np.array([4.8064, 0.4933, 3.616 , 2.0273, 0.4632, 2.2601, 1.8625, 6.1593,
       4.7361, 6.1424, 0.7278, 1.4017, 0.6239, 3.2167, 3.418 , 1.2204,
       2.8195, 0.7517, 1.9862, 1.6483, 5.58  , 2.3938, 2.0529, 1.2396,
       3.9565, 0.7716, 3.0717, 1.6708, 1.6177, 1.7523, 3.1462, 4.2099])
Mstar = np.array([10.2864, 10.4345, 10.5518, 10.7849, 10.6283, 11.1009, 10.6314,
       10.6296, 10.7542, 10.747 , 10.7134, 10.8578, 10.6966, 10.801 ,
       10.7502, 10.8897, 10.7078, 10.6594, 10.4861, 10.5363, 10.5055,
        9.6791, 10.6088, 10.4758, 10.7369, 10.6927, 10.4496, 10.7348,
        9.7876, 10.4582,  9.9995, 10.6496])
MBs = np.array([7.8378, 8.1868, 8.2689, 8.4403, 8.1997, 8.457 , 7.9915, 7.8539,
       8.2891, 8.4494, 7.8876, 8.2953, 8.3114, 8.4229, 8.5263, 8.5606,
       8.4694, 8.4497, 8.3326, 7.9417, 8.2039, 8.4646, 8.1956, 8.5519,
       8.4632, 8.455 , 8.5278, 8.1218, 8.0246, 7.8128, 7.7488, 8.2727])
Mstar_err = np.array([[-0.15,  0.19], [-0.15,  0.19], [-0.14,  0.16], [-0.14,  0.17], 
                      [-0.11,  0.12], [-0.2 ,  0.3 ], [-0.15,  0.18], [-0.1 ,  0.1 ], 
                      [-0.13,  0.14], [-0.15,  0.17], [-0.1 ,  0.1 ], [-0.1 ,  0.1 ], 
                      [-0.1 ,  0.1 ], [-0.21,  0.36], [-0.18,  0.25], [-0.12,  0.13], 
                      [-0.12,  0.13], [-0.11,  0.11], [-0.13,  0.14], [-0.13,  0.15], 
                      [-0.12,  0.12], [-0.24,  0.48], [-0.17,  0.23], [-0.12,  0.12], 
                      [-0.14,  0.16], [-0.19,  0.28], [-0.25,  0.53], [-0.14,  0.16], 
                      [-0.18,  0.25], [-0.18,  0.24], [-0.15,  0.18], [-0.15,  0.18]])
yerr_highz = [((m_ml*Mstar_err[:,0])**2+0.4**2)**0.5, ((m_ml*Mstar_err[:,1])**2+0.4**2)**0.5]
    
ding20_sample = np.log10(1+zs),MBs-(m_ml*Mstar+b_ml)
plt.scatter(np.log10(1+zs),MBs-(m_ml*Mstar+b_ml),c='lightsalmon',
            s=420,marker=".",zorder=10, vmin=0.3, vmax=5, edgecolors='k', alpha = 0.8)
plt.errorbar(np.log10(1+zs),MBs-(m_ml*Mstar+b_ml),
             yerr= yerr_highz,
             color='lightsalmon', fmt='.',markersize=1)    

#%%
HSC = {}
line_means = ['id','z','mbh','mbh_err','stellar_mass','lbol','spectra','bit','ps_gmag','ps_rmag','ps_imag','ps_rmag','ps_zmag','ps_ymag','host_gmag','host_rmag','host_imag','host_zmag','host_ymag']
infers  = np.loadtxt('../HSC_fitting/sdss_quasar_mbh.txt', dtype=str)
IDs_ = infers[:, 0]
HSC_z_overall = infers[:,1].astype(float)
HSC_Mstar_overall = infers[:,4].astype(float)
HSC_MBHs_overall = infers[:,2].astype(float)
HSC_ps_mag_overall = infers[:,10].astype(float)  #'ps_imag'
HSC_MBHs_err_overall  = infers[:,3].astype(float)
HSC_Lbol_overall = infers[:,5].astype(float)
HSC_i_mag_galaxy_overall = infers[:,16].astype(float)
core = True
if core == True:
    HSC_label_ = infers[:,7]
    HSC_z, HSC_Mstar, HSC_MBHs, HSC_ps_mag, HSC_MBHs_err, HSC_Lbol, HSC_i_mag_galaxy, HSC_label= [], [], [], [], [], [], [], []
    for i in range(len(IDs_)):
        if HSC_label_[i] in ['eboss_core', 'boss_core', 'ugri']:
            HSC_z.append(HSC_z_overall[i])
            HSC_Mstar.append(HSC_Mstar_overall[i])
            HSC_MBHs.append(HSC_MBHs_overall[i])
            HSC_ps_mag.append(HSC_ps_mag_overall[i])
            HSC_MBHs_err.append(HSC_MBHs_err_overall[i])
            HSC_Lbol.append(HSC_Lbol_overall[i])
            HSC_i_mag_galaxy.append(HSC_i_mag_galaxy_overall[i])
            HSC_label.append(HSC_label_[i])
    HSC_z_overall = np.array(HSC_z)
    HSC_Mstar_overall = np.array(HSC_Mstar)
    HSC_MBHs_overall = np.array(HSC_MBHs)   
    HSC_ps_mag_overall = np.array(HSC_ps_mag)   
    HSC_MBHs_err_overall  = np.array(HSC_MBHs_err)   
    HSC_Lbol_overall = np.array(HSC_Lbol)   
    HSC_i_mag_galaxy_overall = HSC_i_mag_galaxy
HSC['label'] = np.array(HSC_label)
HSC['HSC_z_overall'] = HSC_z_overall
HSC['HSC_Mstar_overall'] = HSC_Mstar_overall
HSC['HSC_MBHs_overall'] = HSC_MBHs_overall
HSC['HSC_ps_mag_overall'] = HSC_ps_mag_overall
HSC['HSC_MBHs_err_overall']  = HSC_MBHs_err_overall
HSC['HSC_Lbol_overall'] = HSC_Lbol_overall


z_range = np.arange(0.2, 1.0, 0.05)
mstar_cut_range = np.array([8.9, 9.1, 9.3, 9.4, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.3, 10.5, 10.5, 10.6, 10.7, 10.8])
mstar_cut = np.zeros_like(HSC['HSC_z_overall'])
i_mag_cut = np.zeros_like(HSC['HSC_z_overall'])
for i in range(len(mstar_cut)):
    mstar_cut[i] = mstar_cut_range[HSC['HSC_z_overall'][i]  > z_range][-1]
    if HSC['HSC_z_overall'][i]<0.5:
        i_mag_cut[i] = 20.5 
    else:
        i_mag_cut[i] = 22.0
    if HSC['label'][i] == 'ugri':
        i_mag_cut[i] = 19.1
    if HSC['label'][i] == 'eboss_core' or HSC['label'][i] == 'boss_core':
        i_mag_cut[i] = 22.0
    
select_bool = (HSC['HSC_Mstar_overall'] > mstar_cut) * (HSC['HSC_ps_mag_overall'] <i_mag_cut)
#%%
redshift_bool = (HSC_z_overall>0.2)*(HSC_z_overall<0.8)

HSC['HSC_Mstar'] = HSC_Mstar_overall[redshift_bool * select_bool]
HSC['HSC_MBHs'] = HSC_MBHs_overall[redshift_bool * select_bool]
HSC['HSC_ps_mag'] = HSC_ps_mag_overall[redshift_bool * select_bool]  #'ps_imag'
HSC['HSC_Lbol'] = HSC_Lbol_overall[redshift_bool * select_bool]
HSC['HSC_MBHs_err'] = HSC_MBHs_err_overall[redshift_bool * select_bool]
HSC['HSC_z'] = HSC_z_overall[redshift_bool * select_bool]

plt.errorbar(np.log10(1+HSC['HSC_z']),HSC['HSC_MBHs']-(m_ml*HSC['HSC_Mstar']+b_ml),
             yerr= 0.01,
             color='pink', fmt='.',markersize=10)    

obs_scatter = []

cal_z_range = ([0.2,0.4], [0.4,0.6], [0.6,0.8])
for i in range(len(cal_z_range)):
    s_bool = (HSC_z_overall>cal_z_range[i][0])*(HSC_z_overall<cal_z_range[i][1])
    cal_HSC_Mstar = HSC_Mstar_overall[s_bool * select_bool]
    cal_HSC_MBHs = HSC_MBHs_overall[s_bool * select_bool]
    obs_res = cal_HSC_MBHs-(m_ml*cal_HSC_Mstar+b_ml)
    obs_scatter.append( [np.mean(obs_res), np.std(obs_res)] )


#%%

# #####fit the evolution##########
# ################################
z_cosmos, y_cosmos = zs, MBs-(m_ml*Mstar+b_ml)
yerr_hz = (yerr_highz[0]+ yerr_highz[1])/2

obs_scatter.append( [np.mean(MBs-(m_ml*Mstar+b_ml)), np.std(MBs-(m_ml*Mstar+b_ml))] )
                              
# #if consider 32 AGN only:
# z=z_cosmos
# y=y_cosmos
# yerr = yerr_hz    
# yerr = np.sqrt(yerr**2 + sint_ml**2)

z=np.concatenate((z_cosmos, HSC['HSC_z']))
y=np.concatenate((y_cosmos, HSC['HSC_MBHs']-(m_ml*HSC['HSC_Mstar']+b_ml) ))
yerr = np.concatenate((yerr_hz, HSC['HSC_z']*0+0.438 ))
yerr = np.sqrt(yerr**2 + sint_ml**2)
#### fit with emcee ###############
x=np.log10(1+z)
y=y

def lnlike(theta, x, y, yerr):
    b, sint= theta
    model = b*x
    sigma2 = (yerr**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf

import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [1.8, 0.3], args=(x, y, yerr))
b_ml_offset,_= result["x"]

xp = np.array([5, 13])
def lnprior(theta):
    b, sint	 = theta
    if -10 < b < 10.0 and 0 < sint < 10:
        return 0.0
    return -np.inf
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)
ndim, nwalkers = 2, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

b_ml_offset, _ =np.percentile(samples, 50,axis=0)
#print "lnlike=",lnlike(theta=[b_ml_offset, sint_mid],x=x, y=y, yerr=yerr)
xl = np.linspace(0, 5, 100)
# plt.plot(xl, xl*0+xl*1, color="red", linewidth=4.0,zorder=0)
plt.plot(xl, xl*0+xl*b_ml_offset, color="k", linewidth=4.0,zorder=3)
# def find_n(array,value):           #get the corresponding b for a given m 
#     idx= (np.abs(array-value)).argmin()
#     return array[idx]
# b=np.percentile(samples,50,axis=0)[0]
# for i in range(100):
#     posi=np.random.uniform(16,84)
#     b=np.percentile(samples,posi,axis=0)[0]    
#     plt.plot(xl, xl*0+xl*b, color="pink", alpha=0.1,linewidth=7.0,zorder=-1+np.random.normal(0,0.02))
    
plt.plot(xl, xl*0, color="black", linewidth=2.0,zorder=0)
# value=round(b_ml_offset,2)
#####################
value,sig=round(b_ml_offset,2),round((np.percentile(samples,84,axis=0)[0]-np.percentile(samples,16,axis=0)[0])/2,2)
print(value,sig)



#%%Fill in the simulation data
if_int = 0
ls = '-'
if if_int == 0:
    SAM = np.array([(0.73, 0.49) , (0.65, 0.46)  , (0.51, 0.45)  , (0.51, 0.36) ])
    MBII =  np.array([(-0.15, 0.48) , (-0.16, 0.48)  , (0.14, 0.31)])
    Illustris =  np.array([(0.01, 0.52) , (0.08, 0.53)  , (0.06, 0.54)  , (0.07, 0.32) ])
    TNG100 =  np.array([(0.27, 0.48) , (0.24, 0.46)  , (0.24, 0.45)  , (0.38, 0.33)])
    TNG300 =  np.array([(0.26, 0.48) , (0.20, 0.48)  , (0.17, 0.48)  , (0.41, 0.34)])
    Horizon =  np.array([(0.16, 0.49) , (0.14, 0.47)  , (0.23, 0.47)  , (0.47, 0.35)])
elif if_int == 1:
    # SAM = np.array([(0.72, 0.20), (0.64, 0.18), (0.56, 0.16) , (0.08, 0.18) ]) #add noise but not select
    # MBII = np.array([(-0.15, 0.21), (-0.15, 0.22) , (0.08, 0.19)])
    # Illustris = np.array([(0.03, 0.32), (0.10, 0.36), (0.08, 0.36) , (0.04, 0.19) ])
    # TNG100 = np.array([ (0.27, 0.20), (0.27, 0.15), (0.26, 0.16) , (0.36, 0.15) ])
    # TNG300 = np.array([(0.25, 0.21), (0.19, 0.23), (0.20, 0.22) , (0.32, 0.16) ])
    # Horizon = np.array([(0.24, 0.21), (0.23, 0.22), (0.29, 0.19) , (0.37, 0.13) ])
    SAM = np.array([(0.75, 0.26),(0.74, 0.26),(0.94, 0.26) ,(0.59, 0.24) ])
    MBII = np.array([(-0.24, 0.22), (-0.26, 0.22), (-0.29, 0.22)])
    Illustris = np.array([(-0.78, 0.40),(-0.73, 0.40),(-0.69, 0.40) ,(-0.55, 0.36) ])
    TNG100 = np.array([ (0.33, 0.31),(0.31, 0.32),(0.29, 0.33) ,(0.22, 0.37) ])
    TNG300 = np.array([(0.42, 0.58),(0.40, 0.57),(0.39, 0.56) ,(0.35, 0.54) ])
    Horizon = np.array([(0.12, 0.30),(0.13, 0.30),(0.14, 0.28) ,(-0.00, 0.48) ])

sims = [SAM, MBII, Illustris, TNG100, TNG300, Horizon]
sim_label = [' SAM', ' MBII', ' Illustris', ' TNG100', ' TNG300', ' Horizon-AGN']
c = ['b', 'g', 'r', 'c', 'm', 'y']
obs_scatter = np.array(obs_scatter)
zs = np.array([0.3, 0.5, 0.7, 1.5])
if if_int == 0:
    plt.errorbar(np.log10(1+zs)-0.01/6*4, obs_scatter[:,0], obs_scatter[:,1], color = 'k', 
          zorder = 50, linewidth = 3, linestyle= '',fmt='o')
if if_int == 1:
    ls = '--'
for i in range(len(sims)):
    if i !=1 :
        zs = np.array([0.3, 0.5, 0.7, 1.5])
    elif i == 1:
        zs = np.array([0.3, 0.6, 1.5])
    # line = plt.plot(np.log10(1+zs), sims[i][:,0], linestyle = '--', color =c[i],  linewidth=3,zorder=10)
    plt.errorbar(np.log10(1+zs)+0.015/6*i, sims[i][:,0], sims[i][:,1], linestyle = ls, color =c[i],
                  zorder = 50, linewidth = 3)
    ty1=sims[i][:,0] + sims[i][:,1]
    ty2=sims[i][:,0] - sims[i][:,1]
    # zs[0] = zs[0]-0.1
    # zs[-1] = zs[-1]+0.3
    # ty1[0] = ty1[0] - (ty1[1] - ty1[0])*np.log10(0.5-0.3) / np.log10(0.3-0.2)
    # ty2[0] = ty2[0] - (ty2[1] - ty2[0])*np.log10(0.5-0.3) / np.log10(0.3-0.2)
    # ty1[-1] = ty1[-1] - (ty1[-2] - ty1[-1])*np.log10(1.5-0.7) / np.log10(1.8-1.5)
    # ty2[-1] = ty2[-1] - (ty2[-2] - ty2[-1])*np.log10(1.5-0.7) / np.log10(1.8-1.5)   
    # plt.fill_between(np.log10(1+zs),ty1,ty2, color=c[i], zorder=10, alpha = 0.2)
lines_legend = []
for color in c:
    lines_legend.append(mlines.Line2D([], [], color=color) )    

#%%
plt.xlabel(r"log(1+z)",fontsize=45)
plt.xticks(np.arange(-0.1,1,0.1))
xl=-0.01
xh=np.log10(1+2.5)
plt.yticks(np.arange(-5.5,6,0.5))
plt.axis([xl,xh,-2.0,3.5])
plt.ylim([-1.5,1.5])
plt.xlim([-0.01,0.5])
plt.ylabel(r"$\Delta$logM$_{\rm BH}$ (vs M$_*$)",fontsize=45)
plt.grid()
plt.tick_params(labelsize=35)
plt.tick_params(which='major', width=2, length=10, direction='in')
ax2=ax.twiny()
tticks=np.array([10**xl-1,0.5,1,1.5,2,10**xh-1])
ax2.set_xticks([np.log(t+1) for t in tticks])  # for the entire scale
ax2.set_xticklabels([0,0.5,1,1.5,2,2.5])  # 0 actuall is corresponds to 10**-0.01-1
ax2.set_xlabel('z',fontsize=45)
plt.tick_params(labelsize=35)
# SS13 = mlines.Line2D([], [], color='darkseagreen', ls='', marker='^', markersize=13)
ding_sample = mlines.Line2D([], [], color='lightsalmon', ls='', marker='.', markersize=20,markeredgecolor='k')
li_sample = mlines.Line2D([], [], color='pink', ls='', marker='.', markersize=10)

plt.legend([Bkc, 
            # SS13, 
            li_sample, ding_sample]+lines_legend,
[
'Local sample',\
# "Local by H&R",
# "Intermediate redshift AGNs",
"$0.2<$z$<0.8$ AGNs by HSC",
"$1.2<$z$<1.7$ AGNs by HST",
]+sim_label,scatterpoints=1,numpoints=1,loc=3,prop={'size':22,'family': 'Arial'},ncol=3,handletextpad=0)
if if_int == 0:
    plt.savefig("offset_summary_vz.pdf")
if if_int == 1:
    plt.savefig("offset_int_summary_vz.pdf")
plt.show()
