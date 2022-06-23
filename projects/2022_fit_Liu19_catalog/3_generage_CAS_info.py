#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 14:22:19 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob, pickle
from galight.tools.astro_tools import plt_fits
# ID = 10004
# band = 'I'
folder = 'fit_result/'
# file_ = glob.glob(folder+"{0}-{1}.pkl".format(ID, band))

f = open("table_mag_frameflux.txt","r")
string = f.read()
lines = string.split('\n')   # Split in to \n

write_file = open('table_asy_concentration.txt','w') 
write_file.write("#ID, asymmetry, concentration by (GRIZY)\n")

for line in lines[1:-1]:
    results = line.split()
    ID, magG, magR, magI, magZ, magY = results
    write_file.write(ID + ' ')
    for band in 'GRIZY':
        file = glob.glob(folder+"{0}-{1}.pkl".format(ID, band))
        if file!=[]:
            file = file[0]
            fit_run = pickle.load(open(file,'rb'))
            # print(fit_run.final_result_galaxy)
            # host_image = fit_run.flux_2d_out['data'] - fit_run.image_ps_list[0]
            # plt_fits(host_image)
            #For the host image:
            from galight.tools.asymmetry_tools import CAS
            CAS_class = CAS(fit_run, seg_cal_reg = 'or', obj_id=0, rm_ps=True)
            cas = CAS_class.cal_CAS(mask_type='segm',extend=1, if_plot=False)
            # print('asymmetry:', cas[0])
            # print('smoothness (by abs diff and pos diff):', cas[1])
            # print('concentration:', cas[2])
            # print('Gini:', cas[3])
            # print('{0:.3f} {1:.3f}'.format(cas[0], cas[2] ) )
            write_file.write('{0:.3f} {1:.3f} '.format(cas[0], cas[2] ) )
        else:
            write_file.write('{0} {1} '.format(-99, -99 ) )
    write_file.write('\n')

write_file.close()