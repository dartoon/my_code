#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 16:11:45 2020

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
from astropy.cosmology import FlatLambdaCDM
import glob

from lenstronomy.ImSim.image_model import ImageModel
import lenstronomy.Util.param_util as param_util
import lenstronomy.Util.simulation_util as sim_util
from lenstronomy.Data.psf import PSF
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LightModel.light_model import LightModel

sys.path.insert(0, '../share_tools')
#%%Set up data basic information
from quasar_info import qso_info
import os
import rebin #From my share_tools    
oversample = 4
#seed = 0
zp_dic = {'F444W':27.3012, 'F356W':27.1841, 'F200W':27.0383, 'F150W':26.8627} #Using mirage
for seed in range(101,102):
    for ID in range(1,7): #int(input("input simulation ID:\n"))  #The ID for this simulation
        host_flux_rescale = False
        count = 0 
        for filt_i in [3,2,1,0]: #int(input("which filter 0: 'F444W', 1: 'F356W', 2: 'F200W', 3: 'F150W':\n"))
            filt  = ['F444W', 'F356W', 'F200W', 'F150W'][filt_i]
            numPix = [341, 341, 645, 645][filt_i]  # total frame pixel size #!!!Need to be changed for different filter
#            zp = [28., 27.9, 26.7, 27.75][filt_i]   #Using ETC, (616-556) total flux for 23.5 ab mag objects.  #Need to check
            zp = zp_dic[filt]   #Using mirage
            #properties:
            z_s = qso_info['ID'+repr(ID)]['z']       #AGN redshift
            point_mag = qso_info['ID'+repr(ID)]['AGN_{0}_mag'.format(filt)]
            host_mag = qso_info['ID'+repr(ID)]['galaxy_{0}_mag'.format(filt)]
            host_n = np.random.uniform(2,4)   #Host effective radius, unit: Kpc
            host_Reff_kpc = np.random.uniform(2,3)   #Host effective radius, unit: Kpc
            np.random.seed(seed = seed)
            #host_ratio = np.random.uniform(0.4, 0.7) #Set the random host flux ratio [40% - 70%].
            host_flux = 10**(0.4*(zp - host_mag))
            point_flux = 10**(0.4*(zp - point_mag))
            total_flux =  point_flux + host_flux
            host_ratio = host_flux/total_flux
            if count == 0:
                if host_ratio< 0.1:
                    host_ratio =np.random.uniform(0.1,0.2)
                    total_flux = point_flux/(1-host_ratio)
                    host_flux = total_flux - point_flux
                    host_flux_rescale = True
                elif host_ratio > 0.95:
                    host_ratio =np.random.uniform(0.8,0.95)
                    total_flux = point_flux/(1-host_ratio)
                    host_flux = total_flux - point_flux
                    host_flux_rescale = True
                host_flux_F150W = host_flux   
                ref_filt = filt
            if host_flux_rescale == True:
                temp_flux_F150w = 10**(0.4*(zp_dic[ref_filt] - qso_info['ID'+repr(ID)]['galaxy_{0}_mag'.format(ref_filt)]))
                temp_flux = 10**(0.4*(zp - qso_info['ID'+repr(ID)]['galaxy_{0}_mag'.format(filt)]))
                ratio = temp_flux/temp_flux_F150w
                host_flux = host_flux_F150W * ratio
                total_flux =  point_flux + host_flux
                host_ratio = host_flux/total_flux                
            count = count+1    
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
            scale_relation = cosmo.angular_diameter_distance(z_s).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
            #Take the PSF: Use PSF 1 - 8 to simulate, but use PSF0 to model the data
            psf_take_id = np.random.randint(1,8)
            psf_take_name = 'webPSF/highres_PSF_'+filt+'/PSF_id{0}.fits'.format(psf_take_id)
            psf = pyfits.open(psf_take_name)
            #Build up the simulation:
            pix_s = psf[0].header['PIXELSCL'] #* oversample
            print("Pixel scale for", filt, pix_s)
            host_Reff = host_Reff_kpc/scale_relation   #In arcsec
            sim_folder_name = 'sim'+'_ID'+repr(ID)+'_'+filt+'_seed'+repr(seed)
            print("Simulate: "+sim_folder_name+" host ratio:", host_ratio)
            if os.path.exists(sim_folder_name)==True:
                ifdel = input("Simulation of "+ sim_folder_name + " exist, delete the old folder?(Y/N):\n")
                if ifdel == 'Y':
                    import shutil 
                    shutil.rmtree(sim_folder_name)
                elif ifdel == 'N':
                    print("Simulation stop!")
                    from sys import exit
                    exit
                else:
                    raise ValueError("Input is not right, stop simulation.")        
            os.mkdir(sim_folder_name)
            psf_data = psf[0].data
            psf_data = psf_data[1:,1:]
            psf_data /= psf_data.sum()
            kwargs_psf_high_res = {'psf_type': 'PIXEL', 'kernel_point_source': psf_data, 'pixel_size': pix_s}
            kwargs_data_high_res = sim_util.data_configure_simple(numPix, pix_s)
            data_class = ImageData(**kwargs_data_high_res)
            psf_class = PSF(**kwargs_psf_high_res)
            center_x, center_y = np.random.uniform(-1.5, 1.5) * pix_s*oversample, np.random.uniform(-1.5,1.5)* pix_s*oversample
            point_amp = point_flux
            point_source_list = ['UNLENSED']
            pointSource = PointSource(point_source_type_list=point_source_list)
            kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]
            light_model_list = ['SERSIC_ELLIPSE']
            lightModel = LightModel(light_model_list=light_model_list)
            q = np.random.uniform(0.5,0.9)
            phi = np.random.uniform(0.,2*np.pi)
            e1, e2 = param_util.phi_q2_ellipticity(phi=phi, q=q)
            kwargs_numerics = {'supersampling_factor': 2, 'supersampling_convolution': False}
            kwargs_sersic_medi = {'amp': 1. , 'n_sersic': host_n, 'R_sersic': host_Reff/np.sqrt(q), 'e1': e1, 'e2': e2,
                             'center_x': center_x, 'center_y': center_y}
            kwargs_host_medi = [kwargs_sersic_medi]
            imageModel = ImageModel(data_class, psf_class, lens_light_model_class=lightModel,
                                            point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
            medi_host_flux = np.sum(imageModel.image(kwargs_lens_light=kwargs_host_medi, unconvolved=True))
            amp = 1. / medi_host_flux * host_flux
            kwargs_sersic = {'amp': amp, 'n_sersic': host_n, 'R_sersic': host_Reff/np.sqrt(q), 'e1': e1, 'e2': e2,
                             'center_x': center_x, 'center_y': center_y}
            kwargs_host = [kwargs_sersic]
            
            ## simulate image with the parameters we have defined above #
            total_highres = imageModel.image(kwargs_lens_light=kwargs_host, kwargs_ps=kwargs_ps, unconvolved=False)
            host_highres = imageModel.image(kwargs_lens_light=kwargs_host, kwargs_ps=[{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [0]}], unconvolved=False)
            point_highres = total_highres - host_highres
    #        print("AGN image:")
    #        plt.imshow(total_highres, origin='lower',cmap='gist_heat', norm=LogNorm())
    #        plt.colorbar()
    #        plt.show()
    #        print("HOST image:")    
    #        plt.imshow(host_highres, origin='lower',cmap='gist_heat', norm=LogNorm())
    #        plt.colorbar()
    #        plt.show()
    #        print("Point image:")    
    #        plt.imshow(point_highres, origin='lower',cmap='gist_heat', norm=LogNorm())
    #        plt.colorbar()
    #        plt.show()
            #%%
            ##==============================================================================
            ## #Bin the image res. from high to low. 
            ##==============================================================================
            factor=oversample
            pattern_x=[0,2,0,2,1,3,1,3]
            pattern_y=[0,0,2,2,3,3,1,1]      #from the info. given by observation
            ################Bin the lensed image################
            total_exp_grid=rebin.expend_grid(total_highres)
            total_cut_out=np.zeros([len(pattern_x),total_highres.shape[0]-5,total_highres.shape[1]-5])
            total_image_bin =np.zeros([len(pattern_x),int(total_highres.shape[0]/factor)-1,int(total_highres.shape[1]/factor)-1])
            
            host_exp_grid=rebin.expend_grid(host_highres)
            host_cut_out=np.zeros([len(pattern_x),host_highres.shape[0]-5,host_highres.shape[1]-5])
            host_image_bin =np.zeros([len(pattern_x),int(host_highres.shape[0]/factor)-1,int(host_highres.shape[1]/factor)-1])
            
            point_exp_grid=rebin.expend_grid(point_highres)
            point_cut_out=np.zeros([len(pattern_x),point_highres.shape[0]-5,point_highres.shape[1]-5])
            point_image_bin =np.zeros([len(pattern_x),int(point_highres.shape[0]/factor)-1,int(point_highres.shape[1]/factor)-1])
            
            
            for i in range(len(pattern_x)):
                total_cut_out[i]=total_exp_grid[pattern_x[i]:(numPix-5)+pattern_x[i],pattern_y[i]:(numPix-5)+pattern_y[i]]   #the size before bin
                total_image_bin[i]=rebin.block(total_cut_out[i],(int(numPix/factor)-1,int(numPix/factor)-1),factor=factor)
                pyfits.PrimaryHDU(total_image_bin[i]).writeto(sim_folder_name+'/non_drizzled-AGNclean-{0}.fits'.format(i+1),overwrite=False)
            
                host_cut_out[i]=host_exp_grid[pattern_x[i]:(numPix-5)+pattern_x[i],pattern_y[i]:(numPix-5)+pattern_y[i]]   #the size before bin
                host_image_bin[i]=rebin.block(host_cut_out[i],(int(numPix/factor)-1,int(numPix/factor)-1),factor=factor)
                pyfits.PrimaryHDU(host_image_bin[i]).writeto(sim_folder_name+'/non_drizzled-HOSTclean-{0}.fits'.format(i+1),overwrite=False)
            
                point_cut_out[i]=point_exp_grid[pattern_x[i]:(numPix-5)+pattern_x[i],pattern_y[i]:(numPix-5)+pattern_y[i]]   #the size before bin
                point_image_bin[i]=rebin.block(point_cut_out[i],(int(numPix/factor)-1,int(numPix/factor)-1),factor=factor)
                pyfits.PrimaryHDU(point_image_bin[i]).writeto(sim_folder_name+'/non_drizzled-POINTclean-{0}.fits'.format(i+1),overwrite=False)
    #        print("Rebin image:")    
    #        plt.imshow(total_image_bin[0], origin='lower',cmap='gist_heat', norm=LogNorm())
    #        plt.colorbar()
    #        plt.show()
            #==============================================================================
            # Add the noise same as Ding et al. 2017a 
            ######Since two long pics only ###########
            #==============================================================================
            bf_noz = total_image_bin#input simulate data to bf_noz
            rms = np.zeros_like(total_image_bin) #input rms
            noiz = np.zeros_like(total_image_bin) #input noiz
            image_data_noise=np.zeros_like(total_image_bin) #image after noiz
            exptim= 625.   #units of seconds. 
            stdd = [1.6, 1.06, 0.77, 0.75][filt_i] /np.sqrt(exptim)   #An empirical formula from ETC
            for i in range(len(pattern_x)):
                rms[i]=(bf_noz[i]/(exptim)+stdd**2)**0.5 #RMS not saved, thus its value not used here
                bkg_noise= stdd
                noiz[i]=np.random.normal(0, bkg_noise, size=rms[i].shape)
                image_data_noise[i]=noiz[i]+np.random.poisson(lam=bf_noz[i]*exptim)/(exptim)  #Non-drizzled imaged
                pyfits.PrimaryHDU(image_data_noise[i]).writeto(sim_folder_name+'/non_drizzled-image-{0}.fits'.format(i+1),overwrite=False)
                pyfits.PrimaryHDU(rms[i]).writeto(sim_folder_name + '/non_drizzled-noise_map-{0}.fits'.format(i+1),overwrite=False)
                pyfits.PrimaryHDU(rms[i]**2).writeto(sim_folder_name+'/rmsSQ-{0}.fits'.format(i+1),overwrite=False)
            filename_ascii = sim_folder_name+'/sim_info.txt'
            data_info =  open(filename_ascii,'w') 
            data_info.write("Filter:\t{0}\n".format(filt))
            data_info.write("Assumed exposure per frame:\t{0} (s)\n".format(exptim))
            data_info.write("oversample:\t{0}\n".format(oversample))     
            data_info.write("PSF used ID:\t{0}\n".format(psf_take_id))     
            data_info.write("z_s:\t{0}\n".format(z_s))
            data_info.write("host position (x, y):\t({0}, {1}) arcsec\n".format(round(center_x,5),round(center_y,5))) 
            data_info.write("host_flux (within frame):\t{0}\n".format(round(host_flux,3)))   
            data_info.write("host_flux_ratio:\t{0}%\n".format(round(host_ratio*100,3)))  
            data_info.write("host (phi, q):\t({0}, {1}) arcsec\n".format(round(phi,3),round(q,3))) 
            data_info.write("host_n:\t{0}\n".format(host_n))
            data_info.write("host_Reff_kpc:\t{0}\n".format(host_Reff_kpc)) 
            data_info.write("host_Reff:\t{0} arcsec\n".format(round(host_Reff,5))) 
            data_info.write("AGN_flux:\t{0}\n".format(point_flux))
            data_info.write("AGN position (x, y):\t({0}, {1}) arcsec\n".format(round(center_x,5),round(center_y,5))) 
            data_info.close()

