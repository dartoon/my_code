#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 13:49:21 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import os, glob
import shutil 
import pickle

hst_filt_id = {'F606W': '4', 'F814W':'6', 'F105W':'202', 'F125W':'203', 'F140W':'204', 'F160W':'205'}

jwst_filt_id = {'F115W': '352', 'F150W': '353', 'F200W': '354', 
           'F277W': '355', 'F356W': '356', 'F444W': '357', 'F410M': '362'}

def esti_smass(ID, mags_dict, z, folder = 'esti_smass/', flag = 0, if_run_gsf=True, band_as_upper = [],
               mag_err = [], just_run = False, metallicity = 0.0):
    from gsf import gsf
    if just_run == False:
        ID = ID
        z = z
        #%reate a cat file
        folder_path = folder + ID + '/'
        if glob.glob(folder_path) != []:   
            shutil.rmtree(folder_path)
        os.mkdir(path = folder_path)
        # text_temp = "# id F352 E352 F353 E353 F354 E354 F355 E355 F356 E356 F362 E362 F357 E357\n"
        text_temp = "# id "
        mags = []
        filterIDs=''
        filt_id  = jwst_filt_id | hst_filt_id
        # print(filt_id)
        if_hst = []
        for key in mags_dict.keys():
            text_temp = text_temp + ' F{0} E{0}'.format(filt_id[key])
            filterIDs = filterIDs + ','+filt_id[key]
            mags.append(mags_dict[key])
            if key in band_as_upper:
                if_hst.append(None)
            elif key in hst_filt_id.keys():
                if_hst.append(True)
            else:
                if_hst.append(False)
        filterIDs = filterIDs[1:]
        text_temp = text_temp + "\n"
        if mag_err == []:
            mag_err = [0.2] * len(mags)
        fnu = [10 ** ((mags[i]-25)/(-2.5)) for i in range(len(mags))]
        fnu_up = [10 ** ((mags[i]-mag_err[i]-25)/(-2.5)) for i in range(len(mags))]
        fnu_dw = [10 ** ((mags[i]+mag_err[i]-25)/(-2.5)) for i in range(len(mags))]
        fnu_err = [(fnu_up[i]-fnu_dw[i])/2 for i in range(len(mags))]
            
        write_file = open(folder_path+'sample.cat','w') 
        write_file.write(text_temp)
        _string = str(int(ID))
        for i in range(len(fnu)):
            if if_hst[i] == False:
                _string = _string + " {0:.8f} {1:.8f}".format(fnu[i], fnu_err[i])
            elif if_hst[i] == True:
                _string = _string + " {0:.8f} {1:.8f}".format(fnu[i], fnu_err[i]*2)
            else:
                _string = _string + " {0:.8f} {1:.8f}".format(0, fnu[i])
        write_file.write(_string)
        write_file.close()
    
        #Create a input file
        f = open("./SED_temp/sample_template.input","r")
        string = f.read()
        string = string.replace("idname", str(int(ID)))
        string = string.replace("zinfo", str(z))
        string = string.replace("folder/", folder_path)
        string = string.replace("filterIDs", filterIDs)
        string = string.replace("metaltemp", str(metallicity))
        string = string.replace("../../../../", "../../../")
        write_file = open(folder_path+'sample.input','w') 
        write_file.write(string)
        write_file.close()
        
    if if_run_gsf == True:
        gsf.run_gsf_all(folder+'{0}/sample.input'.format(ID), flag, idman=None)
        #Move things in position.
        mv_files = glob.glob('*_{0}_*'.format(ID)) + glob.glob('*_{0}.*'.format(ID))
        for mv_file in mv_files:
            shutil.move(mv_file, folder+'{0}/'.format(ID)+mv_file)
            
            
            
            