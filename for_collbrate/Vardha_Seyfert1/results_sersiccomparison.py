#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 21:48:10 2019

@author: Dartoon
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
import corner
import astropy.io.fits as pyfits
from matplotlib.colors import LogNorm
import copy
import lenstronomy.Util.param_util as param_util
import corner
sys.path.insert(0,'./fitting_tools/')
from flux_profile import total_compare
ID = '11' #sys.argv[1]
namelist = [ID]

#namelist = ['2','10','11','70','71','74','76','79','99','102','103','106','109','126']
#namelist = ['99']
#namelist = ['2','10','11','106','109']
#namelist =  ['10','11','71','74','76','79','103','106','126']
print ("# ID, Run, Component, Reff, n, q")
for i in range(len(namelist)):
    ID = namelist[i]
    rf = open("./fitting_combinedpsf.dat",'r')
    for line in rf:
        ob = line.split( )[0]
        if ob == ID:
            snratio = int(line.split( )[1])
            numberpixels = int(line.split( )[2])
            pixsz = float(line.split( )[3])
            zero = float(line.split( )[4])
            exptime = int(line.split( )[5])
            fr = line.split( )[6]
            kernel_size = int(line.split( )[7])
            psfname = line.split( )[8]
            phibar = float(line.split( )[9])
            qbar = float(line.split( )[10])
            rbar = float(line.split( )[11])
            if len(line.split( ))>12:
                mask = np.array(line.split( )[12].split(',')).astype('int')
            else:
                mask=[]
            #print (snratio, numberpixels, pixsz, zero, exptime, fr, kernel_size, psfname, mask)

    pix_sz = pixsz
    zp = zero
    #runlist = ['1','2','3','4','5','6','7','8','9']
    #runlist = ['1','2','3']
    runlist = ['2']
    for i in range(len(runlist)):
        run = runlist[i]
        picklename = 'zoutput/l'+ID+'_sersic_'+run+'.pkl'
        result = pickle.load(open(picklename,'rb'))
        # best_fit, chain_list_result, trans_paras, material = result
        # #source_result, image_host, ps_result, image_ps, _ =best_fit
        # source_result, ps_result, image_ps, image_host, _=best_fit
        
        best_fit, chain_list_result, trans_paras, material = result
        source_result, image_host, ps_result, image_ps, _ = best_fit        
        chain_list, _ = chain_list_result
        multi_band_list, kwargs_model, kwargs_result, QSO_msk, kwargs_fixed_source, kwargs_fixed_ps, kwargs_constraints, kwargs_numerics, classes = material
        error_map = multi_band_list[0][0]['noise_map']
        phi0, q0 = param_util.ellipticity2phi_q(source_result[0]['e1'], source_result[0]['e2'])
        for i in range(len(kwargs_result['kwargs_source'])):
            print (ID, run, format(i), round(kwargs_result['kwargs_source'][i]['R_sersic'],2), round(kwargs_result['kwargs_source'][i]['n_sersic'],2), round(q0,2))
        print (error_map.shape)    

from lenstronomy.Plots import chain_plot
for i in range(len(chain_list)):
    f, axes = chain_plot.plot_chain_list(chain_list,i)
    plt.show()

agn_image = pyfits.getdata('./allscience/l{0}_{1}_cutout.fits'.format(ID,fr))
if len(image_host) == 1:
    host = image_host[0]
    label = ['data', 'QSO', 'host', 'model', 'normalized residual']
elif len(image_host) >1:
    host = np.zeros_like(image_host[0])
    for i in range(len(image_host)):
        host += image_host[i]
    label = ['data', 'QSO', 'host as {0} components'.format(i+1), 'model', 'normalized residual']  #Print the numbers
        
flux_list = [agn_image, image_ps[0], host, error_map]
fig = total_compare(label_list = label, flux_list = flux_list, target_ID = ID, pix_sz=pix_sz, zp = zp,
                    plot_compare = False, msk_image = QSO_msk)
from IPython.display import display
display(fig)


#%% Recover the plot
from lenstronomy.Plots.model_plot import ModelPlot
modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result,
                          arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[QSO_msk])
f, axes = plt.subplots(3, 3, figsize=(16, 16), sharex=False, sharey=False)
modelPlot.data_plot(ax=axes[0,0], text="Data")
modelPlot.model_plot(ax=axes[0,1])
modelPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
modelPlot.decomposition_plot(ax=axes[1,0], text='Host galaxy', source_add=True, unconvolved=True)
modelPlot.decomposition_plot(ax=axes[1,1], text='Host galaxy convolved', source_add=True)
modelPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
modelPlot.subtract_from_data_plot(ax=axes[2,0], text='Data - Point Source', point_source_add=True)
modelPlot.subtract_from_data_plot(ax=axes[2,1], text='Data - host galaxy', source_add=True)
modelPlot.subtract_from_data_plot(ax=axes[2,2], text='Data - host galaxy - Point Source', source_add=True, point_source_add=True)

f.tight_layout()
plt.show()