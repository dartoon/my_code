#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 20:32:34 2021

@author: Xuheng Ding
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import pickle
from matplotlib.colors import LogNorm

Bla_type = 'BL_Lac'
# Bla_type = 'FSRQ'

filename = 'model_result/{0}_host_result.txt'.format(Bla_type)
#Read files:
f = open(filename,"r")
string = f.read()
results = string.split('\n')   # Split in to \n
ids_file = glob.glob('../../data/{1}/*-{0}.fits'.format('I', Bla_type))
ids = [ids_file[i].split('/')[-1].split('_')[0] for i in range(len(ids_file))]

# for k in range(len(ids)):
#FSRQ
# for k in [2, 8, 13, 16, 23]:
#BL_Lac
for k in [25, 27, 28, 31, 34, 35, 39, 40, 42, 44, 45, 46, 50, 52, 54, 55, 57]:
    bands = ['G', 'R', 'I' ,'Z' , 'Y']
    mags = []
    chisqs = []
    residuals = []
    for band in bands:
        line = [i for i in range(len(results)) if ids[k]+'_HSC'+'-'+band in results[i]]
        if line != []:
            fit_result = results[line[0]].split("'magnitude': ")[1]
            if len(fit_result)>5:
                mags.append(float(fit_result[:-3]))
            else:
                mags.append(30)
        else:
            mags.append(30)
        
        pklfile_name = 'model_result/'+'{0}_HSC-{1}.pkl'.format(ids[k], band)
        if glob.glob(pklfile_name) != []:
            fitting_run_class = pickle.load(open(pklfile_name,'rb'))
            chisqs.append(fitting_run_class.reduced_Chisq)
            residuals.append(fitting_run_class.fitting_specify_class.kwargs_data['image_data'] - fitting_run_class.image_ps_list[0])#/fitting_run_class.fitting_specify_class.kwargs_data
        else:
            chisqs.append(99)
            residuals.append(np.zeros((60,60)) -10)
    mags = [round(mags[i],2) for i in range(len(mags))]
    chisqs = [round(chisqs[i],2) for i in range(len(chisqs))]
    # if np.max(chisqs) < 10 and np.min(mags)>0:
    if np.max(chisqs) > 0:
        print(ids[k], k) 
        print(mags)
        mag_err = [0.2, 0.2, 0.2, 0.2, 0.2]
        for i in range(len(bands)):
            band = bands[i]
            print("Chisq", chisqs[i])
            print("min residual", np.sum(residuals[i]<-0.5))
            pklfile_name = 'model_result/'+'{0}_HSC-{1}.pkl'.format(ids[k], band) 
            if glob.glob(pklfile_name) != []:
                fitting_run_class = pickle.load(open(pklfile_name,'rb'))
                fitting_run_class.plot_final_qso_fit(target_ID = ids[k]+'-'+band)
            else:
                print(pklfile_name.split('/')[1], " DOES NOT EXIST!!!")
            
        # In[ ]:
        #To estimate the AB mag for a given stellar template
        ID_org = ids[k]
        ID = str("".join(filter(str.isdigit, ID_org)))
        mag = mags
        fnu = [10 ** ((mag[i]-25)/(-2.5)) for i in range(len(mag))]
        fnu_up = [10 ** ((mag[i]-mag_err[i]-25)/(-2.5)) for i in range(len(mag))]
        fnu_dw = [10 ** ((mag[i]+mag_err[i]-25)/(-2.5)) for i in range(len(mag))]
        fnu_err = [(fnu_up[i]-fnu_dw[i])/2 for i in range(len(mag))]
        for i in range(len(bands)):
            if chisqs[i]>2:
                if np.sum(residuals[i]<-0.5)>25:
                    fnu_err[i] = 150
            if mags[i] > 23:
                fnu_err[i] = 150  
                
        import matplotlib as mat
        mat.rcParams['font.family'] = 'STIXGeneral'
        # from specutils import Spectrum1D
        if glob.glob("./SED_file/gsf_spec_{0}.fits".format(ID))!= []:
            filename_p  = './SED_file/SFH_{0}_PA00_param.fits'.format(ID)
            hdul = pyfits.open(filename_p)
            table = hdul[1].data
            name = hdul[1].columns
            age_idx = [i for i in range(len(name)) if 'T_MW' in str(name[i])][0]
            age = str(round(10**table[1][age_idx],3)) # 'AGE:' Gyr
            z_idx = [i for i in range(len(name)) if 'zmc' in str(name[i])][0]
            z = str(round(table[1][z_idx],2))
            mel_idx = [i for i in range(len(name)) if 'Z_MW' in str(name[i])][0]
            mel = str(round(10**table[1][mel_idx],3)) # 'Mel:' Z*/Z_sun
            
            spec1d = pyfits.open("./SED_file/gsf_spec_{0}.fits".format(ID))  
            name_spec = spec1d[1].columns
            table_spec = spec1d[1].data
            plt.figure(figsize=(10, 6))
            f = table_spec['f_model_50']            #This is the best-fit spectra (i.e. 50% by MCMC).
            wave = table_spec['wave_model']/10000.
            f_max = f.max()
            plt.plot(wave, f)
            
            f = table_spec['f_model_16']
            wave = table_spec['wave_model']/10000.
            plt.plot(wave, f)
            
            f = table_spec['f_model_84']
            wave = table_spec['wave_model']/10000.
            plt.plot(wave, f)
            
            HSC_g_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/314.fil')
            HSC_g_fil[:,2] = HSC_g_fil[:,2]/HSC_g_fil[:,2].max()* f_max/6
            plt.plot(HSC_g_fil[:,1]/10000., HSC_g_fil[:,2], label='HSC_G response')
            
            HSC_r_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/315.fil')
            HSC_r_fil[:,2] = HSC_r_fil[:,2]/HSC_r_fil[:,2].max()* f_max/6
            plt.plot(HSC_r_fil[:,1]/10000., HSC_r_fil[:,2], label='HSC_R response')
            
            HSC_i_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/316.fil')
            HSC_i_fil[:,2] = HSC_i_fil[:,2]/HSC_i_fil[:,2].max()* f_max/6
            plt.plot(HSC_i_fil[:,1]/10000., HSC_i_fil[:,2], label='HSC_I response')
            
            HSC_z_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/317.fil')
            HSC_z_fil[:,2] = HSC_z_fil[:,2]/HSC_z_fil[:,2].max()* f_max/6
            plt.plot(HSC_z_fil[:,1]/10000., HSC_z_fil[:,2], label='HSC_Z response')
            
            HSC_y_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf_v1.3/gsf/example/filter/318.fil')
            HSC_y_fil[:,2] = HSC_y_fil[:,2]/HSC_y_fil[:,2].max()* f_max/6
            plt.plot(HSC_y_fil[:,1]/10000., HSC_y_fil[:,2], label='HSC_Y response')
            
            yerr_l = np.array(fnu_err)/np.array(fnu)
            
            lam = np.array([4798.2, 6218.4, 7727.0, 8908.2, 9775.1])
            mag = -2.5*np.log10(fnu) + 25
            flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**19
            lam = lam/10000.
            
            plt.scatter(lam, flambda, c='r', zorder =100)
            plt.errorbar(lam, flambda, yerr=yerr_l*flambda,fmt='.',color='gray',markersize=1, zorder =90)
            #plt.plot(sov_jwst_f144w_fil[:,0]/10000., sov_jwst_f144w_fil[:,1], label='NIRCam F444W response curve')
            plt.xlim(0.25, 3)
            plt.ylim(-50,np.max(flambda) * 1.5)
            xmin, xmax, ymin, ymax = plt.axis()
            plt.title("photo-z of "+ID_org, fontsize=25)
            plt.text( (xmax-xmin)*0.7, (ymax-ymin)*0.6, 'stellar population with:', fontsize=17)
            plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.53, 'z={0}'.format(z), fontsize=17)
            plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.45, 'age={0:.3f} Gyr'.format(float(age)), fontsize=17)
            plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.37, 'Z$_*$/Z$_\odot$={0}'.format(mel), fontsize=17)
            plt.legend(prop={'size':15}, ncol=3, loc = 1)
            plt.tick_params(labelsize=15)
            plt.xlabel("um",fontsize=27)
            plt.ylabel(r"f$_\lambda$ 10$^{\rm -19}$ erg/s/cm$^2$/$\AA$",fontsize=27)
            # plt.savefig('./SED_file/photo-z_' + ID_org + '.pdf')
            plt.show()
        print(Bla_type, ids[k], k)
        hold = input("Next:")
# Bla_type, ids[k], k
# FSRQ J1520.5+4209 2
# FSRQ J1027.9+0252 8
# FSRQ J0217.8+0144 13  #Model image looks weird.
# FSRQ J0839.8+0105 16 #Image looks interesting, but SED Fitting Looks not good.
# FSRQ J0914.4+0249 23

# BL_Lac J0842.5+0251 2 #Host image is clear in the reder band.
# BL_Lac J1253.8+0327 4 #Bad SED, but very extended source
# BL_Lac J1117.2+0008 12 # See host clearly
# BL_Lac J0227.3+0201 16 # Bad fitting, but very dense field.
# BL_Lac J1311.0+0034 18 #Mergering system?
# BL_Lac J1056.0+0253 20 #Very clear host
# BL_Lac J0857.7+0137 21 # exciting image in Y band, clear host and merger
# BL_Lac J0022.0+0006 25 # Clear host
# BL_Lac J2204.3+0438 27 # Very extended iamge
# BL_Lac J2211.0-0003 28 # Detecte host in higher band
# BL_Lac J1154.0-0010 31 # Detecte host in higher band
# BL_Lac J1637.2+4327 34 # Clear host
# BL_Lac J1051.9+0103 35 # Mergering system?
# BL_Lac J0113.7+0225 39 # Very extended source,
# BL_Lac J1104.0+0020 40 # Close to a spiral galaxy, amazing merger system.
# BL_Lac J0152.6+0147 42 # Very extended source,
# BL_Lac J1506.4+4331 44 # Detecte host in higher band
# BL_Lac J0946.2+0104 45 # Mergering systems?
# BL_Lac J0201.1+0036 46 # Host image clear and mergering system?
# BL_Lac J1246.3+0112 50 # Detection the host see one D  profile
# BL_Lac J1428.5+4240 52 # Very extended source,
# BL_Lac J0831.8+0429 54 # Very dense field
# BL_Lac J1539.9+4220 55 # Detect host in Y band?
# BL_Lac J0849.5+0456 57 # Detection the host? see one D I band profile