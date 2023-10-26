#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 17:55:30 2019

@author: Dartoon
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner
import glob
import astropy.io.fits as pyfits
from gen_fit_id import gen_fit_id
import sys

band_list = ['G', 'R', 'I', 'Z', 'Y']
zp = 27.0

#QSO_id = sys.argv[1]
#QSO_RA = float(sys.argv[2])
#QSO_DEC = float(sys.argv[3])

QSO_id = "000017.88+002612.6"
QSO_RA = 0.07452999800443649
QSO_DEC = 0.4368380010128021

folder = 'result/'
QSO_id_address = glob.glob(folder + '*_HSC-R.pkl')
QSO_id_list = [QSO_id_address[i].split('_')[-2] for i in range(len(QSO_id_address))]

filename_ascii = folder + '/ascii_mag_error.txt'
if_file = glob.glob(filename_ascii)   
if if_file == []:
    ascii_err =  open(filename_ascii,'w') 
    ascii_err.write("#ID, mag_g, mag_g_err, mag_r, mag_r_err, mag_i, mag_i_err, mag_z, mag_z_err, mag_y, mag_y_err\
\n")
    ascii_err.close()
ascii_err = open(filename_ascii,"r+")
ascii_err.read()
for QSO_id in QSO_id_list:
    for band in band_list:
        picklename = folder +'fit_image_'+QSO_id+ '_HSC-{0}.pkl'.format(band)
        if_file = glob.glob(picklename)
        m_m, m_err = -99, -99
        if if_file == []:
            print "MCMC for {0}-{1} not performed!".format(QSO_id,band)
        else:
            result = pickle.load(open(picklename,'rb'))
            [best_fit,pso_fit,mcmc_fit, trans_paras] = result            
            source_result, image_host, ps_result, image_ps, _ =best_fit
            chain_list, param_list, _ = pso_fit
            samples_mcmc, param_mcmc, dist_mcmc, _ = mcmc_fit
            _, _, mcmc_new_list, labels_new, _ = trans_paras            
            mcmc_new_list = np.asarray(mcmc_new_list)
            idx = 1     #The translated flux for the host
            v_l=np.percentile(mcmc_new_list[:,idx],16,axis=0)
            v_m=np.percentile(mcmc_new_list[:,idx],50,axis=0)
            v_h=np.percentile(mcmc_new_list[:,idx],84,axis=0)
            #print labels_new[idx], ":", v_l, v_m, v_h
            m_l = -2.5 * np.log10(v_h) + zp
            m_m = -2.5 * np.log10(v_m) + zp
            m_h = -2.5 * np.log10(v_l) + zp
            m_err = ((m_h-m_m)**2 + (m_m-m_l)**2)**0.5
        if band == 'G':
            ascii_err.write("{0} {1} {2}".format(QSO_id, round(m_m,3), round(m_err,3)))
        else:
            ascii_err.write(" {0} {1}".format(round(m_m,3), round(m_err,3)))
    ascii_err.write("\n")
ascii_err.close()