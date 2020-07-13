#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 22:34:48 2019

@author: Dartoon

Fit the Seyfert as Point source + two Sersic (Bluge + Disk).
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'./fitting_tools/')
from fit_qso import fit_qso
from transfer_to_result import transfer_to_result, transfer_obj_to_result
import glob
matplotlib.rcParams['font.family'] = 'STIXGeneral'
#import colormap
#cmap = matplotlib.cm.get_cmap('viridis')
from flux_profile import total_compare
import lenstronomy.Util.param_util as param_util
from fit_qso import condition_bulgediskbar

import numpy as np
import corner
from lenstronomy.Util import kernel_util
from lenstronomy.Util import image_util
from lenstronomy.Plots import model_plot

#namelist = ['2','5','10','11','70','71','74','76','79','99','102','103','106','109','126']#Gemini
namelist = ['10'] #individual
for i in range(len(namelist)):
    ID = namelist[i]

    rf = open("./fitting_smallpsf.dat",'r')
    for line in rf:
        ob = line.split()[0]
        if ob == ID:
            snratio = int(line.split()[1])
            numberpixels = int(line.split()[2])
            pixsz = float(line.split()[3])
            zero = float(line.split()[4])
            exptime = int(line.split()[5])
            fr = line.split()[6]
            kernel_size = int(line.split()[7])
            psfname = line.split()[8]
            if len(line.split())>9:
                mask = np.array(line.split()[9].split(',')).astype('int')
            else:
                mask=[]
            print (snratio, numberpixels, pixsz, zero, exptime, fr, kernel_size, psfname, mask)        

    ID = 'l{0}'.format(ID)             # Object ID
    fr = fr                           # frame size
    psf = pyfits.getdata('./psf/{0}.fits'.format(psfname))  # Which PSF to use
    block_id = mask
    deep_seed = True                  # False = trial run
    no_MCMC = True                     # True = without MCMC 
    kernel_size = kernel_size                   # psf kernel size (odd number)

    name_save = '{0}_sersicdiskbar'.format(ID)
    pix_sz = pixsz  # The pixel scale of the image.
    exp_time = exptime
    zp = zero
    pltshow = 1  #Change to 1 to show the plot while fitting.
    mpi_para = False

    frsize = int(fr)
    lowerR = pix_sz #lower limit radius
    upperR = pix_sz*frsize #upper limit radius
    lowerpos = -2*pix_sz #lower limit offset in x
    upperpos = 2*pix_sz #lower limit offset in y
    centerx = 2*pix_sz #step size in x
    centery = 2*pix_sz #step size in y
    
    agn_image = pyfits.getdata('./allscience/{0}_{1}_cutout.fits'.format(ID,fr))
    # agn_stdd = pyfits.getdata('./allscience/{0}_{1}_stdd.fits'.format(ID,fr))
    stdd =  4  #Measurement from empty retion, 0.016*0.08**2/0.13**2/np.sqrt(8)
    agn_stdd = (abs(agn_image/exp_time)+stdd**2)**0.5

    print ("The fitting image:")
    plt.imshow(agn_image, norm = LogNorm(), origin='lower')
    plt.savefig('./zoutput/{0}_agn.png'.format(ID))
    if pltshow == 0:
        plt.close()
    else:
        plt.show()

    print ("The adopted PSF:")

    psf = kernel_util.cut_psf(psf, psf_size=kernel_size)

    plt.imshow(psf, norm = LogNorm(), origin='lower')
    plt.savefig('./zoutput/{0}_psf.png'.format(ID))
    if pltshow == 0:
        plt.close()
    else:
        plt.show()

    from mask_objects import mask_obj
    _ , _, deblend_sources = mask_obj(agn_image, snr=snratio, npixels=numberpixels, return_deblend = True)

    print ("deblend image to find the ID for the Objects for the mask:")
    plt.imshow(deblend_sources, origin='lower',cmap=deblend_sources.cmap(random_state=12345))
    plt.colorbar()
    plt.savefig('./zoutput/{0}_mask.png'.format(ID))
    if pltshow == 0:
        plt.close()
    else:
        plt.show()

    if block_id == []:
        QSO_msk = np.ones_like(agn_image)
    else:
        for i in range(len(block_id)):
            if i ==0:
                mask = (np.array(deblend_sources)==block_id[i])
            else:
                mask += (np.array(deblend_sources)==block_id[i])
        
            QSO_msk = 1- mask
    print ("The QSO mask for the fitting:")
    plt.imshow(QSO_msk, origin='lower')
    plt.savefig('./zoutput/{0}_qso_mask.png'.format(ID))
    if pltshow == 0:
        plt.close()
    else:
        plt.show()

    fixed_source = []
    kwargs_source_init = []
    kwargs_source_sigma = []
    kwargs_lower_source = []
    kwargs_upper_source = []      
    fixed_source.append({})  
    kwargs_source_init.append({'R_sersic': 5.5, 'n_sersic': 1.4, 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
    kwargs_source_sigma.append({'n_sersic': 2, 'R_sersic': 2, 'e1': 0.2, 'e2': 0.2, 'center_x': centerx, 'center_y': centery})
    kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 4, 'n_sersic': 1, 'center_x': lowerpos, 'center_y': lowerpos})
    kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10, 'n_sersic': 5,'center_x': upperpos, 'center_y': upperpos})
    
    # kwargs_source_init.append({'R_sersic': pix_sz*4, 'n_sersic': 1.4, 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
    # kwargs_source_sigma.append({'n_sersic': 2, 'R_sersic': 2, 'e1': 0.2, 'e2': 0.2, 'center_x': centerx, 'center_y': centery})
    # kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': pix_sz*3, 'n_sersic': 1, 'center_x': lowerpos, 'center_y': lowerpos})
    # kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 1, 'n_sersic': 5,'center_x': upperpos, 'center_y': upperpos})
    # fixed_source.append({'n_sersic': 1.})  
    # kwargs_source_init.append({'R_sersic': 5.5, 'n_sersic': 1., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
    # kwargs_source_sigma.append({'n_sersic': 0.1, 'R_sersic': 3, 'e1': 0.2, 'e2': 0.2, 'center_x': centerx, 'center_y': centery})
    # kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': lowerR*3, 'center_x': lowerpos, 'center_y': lowerpos})
    # kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': upperR, 'center_x': upperpos, 'center_y': upperpos})
    # fixed_source.append({'n_sersic': 0.5})  
    # kwargs_source_init.append({'R_sersic': 2.78, 'n_sersic': 0.5, 'e1': 0.457, 'e2': -0.452, 'center_x': 0., 'center_y': 0.})
    # kwargs_source_sigma.append({'n_sersic': 0.1, 'R_sersic': 3, 'e1': 0.2, 'e2': 0.2, 'center_x': centerx, 'center_y': centery})
    # kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 2.78 * 0.5, 'center_x': lowerpos, 'center_y': lowerpos})
    # kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 2.78 * 2, 'center_x': upperpos, 'center_y': upperpos})


    source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]

    fixed_ps = []
    kwargs_ps_init = []
    kwargs_ps_sigma = []
    kwargs_lower_ps = []
    kwargs_upper_ps = []
    fixed_ps.append({})
    kwargs_ps_init.append({'ra_image': [0.], 'dec_image': [0.]})
    kwargs_ps_sigma.append({'ra_image': [centerx], 'dec_image': [centery]})
    kwargs_lower_ps.append({'ra_image': [lowerpos], 'dec_image': [lowerpos]})
    kwargs_upper_ps.append({'ra_image': [upperpos], 'dec_image': [upperpos]})
    ps_param = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]         

    source_result, ps_result, image_ps, image_host, error_map=fit_qso(agn_image, psf_ave=psf, ps_param = ps_param, 
                                                                  pix_sz = pix_sz, exp_time = exp_time, QSO_std = agn_stdd,  QSO_msk = QSO_msk,
                                                                  source_params=source_params, fixcenter=False, no_MCMC = no_MCMC,
                                                                  tag=name_save, deep_seed= deep_seed, pltshow = pltshow, flux_ratio_plot=True,
                                                                  dump_result = True)#, condition=condition_bulgediskbar)


    if len(image_host) == 1:
        host = image_host[0]
        label = ['data', 'QSO', 'host', 'model', 'normalized residual']
    elif len(image_host) >1:
        host = np.zeros_like(image_host[0])
        for i in range(len(image_host)):
            host += image_host[i]
        label = ['data', 'QSO', 'host as {0} components'.format(i+1), 'model', 'normalized residual']  #Print the numbers of objects
    flux_list = [agn_image, image_ps[0], host, error_map]
    fig = total_compare(label_list = label, flux_list = flux_list, target_ID = ID, pix_sz=pix_sz, zp = zp,
                    plot_compare = False, msk_image = QSO_msk)
    fig.savefig("./zoutput/{0}_SB_profile.pdf".format(name_save), bbox_inches = 'tight')
    if pltshow == 0:
        plt.close()
    else:
        fig.show()
    


    agn_image_hdu = pyfits.PrimaryHDU(agn_image)
    agn_image_hdu.header['targname'] = "The fitting image"
    image_ps_hdu = pyfits.ImageHDU(image_ps[0])
    image_ps_hdu.header['targname'] = "Best fit point source image"
    host_hdu = pyfits.ImageHDU(host)
    host_hdu.header['targname'] = "All objects image together"
    error_map_hdu = pyfits.ImageHDU(error_map)
    error_map_hdu.header['targname'] = "The error map of the fitting"
    QSO_msk_hdu = pyfits.ImageHDU(QSO_msk)
    QSO_msk_hdu.header['targname'] = "The QSO mask"
    bulge_hdu = pyfits.ImageHDU(image_host[0])
    bulge_hdu.header['targname'] = "Bulge"
    disk_hdu = pyfits.ImageHDU(image_host[1])
    disk_hdu.header['targname'] = "Disk"
    bar_hdu = pyfits.ImageHDU(image_host[2])
    bar_hdu.header['targname'] = "Bar"
    thdu_fluxlist = pyfits.HDUList([agn_image_hdu,image_ps_hdu,host_hdu, error_map_hdu, QSO_msk_hdu, bulge_hdu, disk_hdu, bar_hdu])
    thdu_fluxlist.writeto('./zoutput/'+ID+'_flux_list_sersicdiskbar.fits', overwrite=True)

    filename = './zoutput/{0}_fit_result.txt'.format(name_save)
    if_file = glob.glob(filename)   
    if if_file == []:
        fit_result =  open(filename,'w') 
    elif if_file is not []:
        fit_result = open(filename,"r+")
        fit_result.read()
    fit_result.write("Result for target " + ID + ":\n")
    ps_result_0 = ps_result[0]
    ps_result_0['PSF_mag'] = -2.5*np.log10(ps_result_0['point_amp']) + zp
    fit_result.write("point source result:\n")
    fit_result.write("    "+ repr(ps_result) + "\n")
    for i in range(len(source_result)):
        fit_result.write("obj {0} result:\n".format(i))
        result = transfer_obj_to_result(source_result =source_result[i] ,image_host=image_host[i], zp=zp)
        fit_result.write("    "+ repr(result) + "\n")
    fit_result.write('======================\n')
    fit_result.close()    

    print ("Total Magnitude Result:")
    print ("AGN magnitude:", -2.5*np.log10(image_ps[0].sum()) + zp)
    print ("Bulge magnitude:", -2.5*np.log10(image_host[0].sum()) + zp)
    print ("Disk magnitude:", -2.5*np.log10(image_host[1].sum()) + zp)
    print ("Bar magnitude:", -2.5*np.log10(image_host[2].sum()) + zp)
