#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 08:01:23 2020

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'/Users/Dartoon/Astro/my_code/py_tools/')

from flux_profile import profiles_compare, SB_profile
import matplotlib

import lenstronomy.Util.util as util


def cal_gamma(convergence_map):
    r_array = np.linspace(len(convergence_map)/2*0.1, len(convergence_map)/2*0.7, len(convergence_map)/2)
    x_grid, y_grid = util.make_grid(numPix=len(convergence_map), deltapix=1)  #Creat the grid in 1 D (Copy from Lenstronomy)
    for r in r_array:
        R = np.sqrt(x_grid*x_grid + y_grid*y_grid)
        mask = np.empty_like(R)
        n = int(len(mask)**0.5)
        mask[R > r] = 0
        mask[R <= r] = 1
        mask_2d = mask.reshape(n,n)
        kappa_mean = np.sum(convergence_map*mask_2d)/np.sum(mask_2d)
        if kappa_mean < 1:
            Rein_r = r
            break
#        print "Rein = ", Rein_r * pixsize
    theta1 = (0.8*Rein_r)
    theta2 = (1.2*Rein_r)
    r_mid = len(convergence_map)/2    
    r_mid = len(convergence_map)/2
    k_profile, r_list = SB_profile(convergence_map, [r_mid,r_mid],radius=theta2+5,
                                   grids=150, if_annuli = True,gridspace='log',fits_plot=False,ifplot=False)
    p1 = np.where(abs(theta1-r_list) == abs(theta1-r_list).min())
    p2 = np.where(abs(theta2-r_list) == abs(theta2-r_list).min())
    k1 = k_profile[p1]
    k2 = k_profile[p2]
    s = np.log(k2/k1)/np.log(float(theta1)/theta2)  #Eq 2 in Dandan's paper
    return s+1, Rein_r

def compare_profile(prf_list,seed_id="?", gridspace = 'log', scal_list = [1,1], grids = 70, if_annuli = False, norm_pix = None, 
                    prf_name_list = ['Truth','Model'], y_log = True, radius = 45, start_p=0.5):
#    if gridspace == 'log':
#        radius = len(prf_list[1])/2
    fig, ax = plt.subplots(figsize=(10,7))
    prf_NO = len(prf_list)
    for i in range(prf_NO):
        b_c = len(prf_list[i])/2
        b_r = len(prf_list[i])/6
        center = np.reshape(np.asarray(np.where(prf_list[i]== prf_list[i][b_c-b_r:b_c+b_r,b_c-b_r:b_c+b_r].max())),(2))[::-1]
        scale = scal_list[i]
        r_SB, r_grids = SB_profile(prf_list[i], center, radius=radius*scale, start_p=start_p,
                                   grids=grids, gridspace=gridspace,if_annuli=if_annuli)
        r_grids /= scale
        if prf_name_list == None:
            plt.plot(r_grids, r_SB, 'x-', label="prf_list{0}".format(i))
        elif len(prf_name_list)==len(prf_list) and y_log == False:
            plt.plot(r_grids, r_SB, 'x-', label=prf_name_list[i])
        elif len(prf_name_list)==len(prf_list) and y_log == True:
            plt.plot(r_grids, np.log10(r_SB), 'x-', label=prf_name_list[i])
            plt.ylim([-1.5,1.25])
        else:
            raise ValueError("The profile name is not in right length")
    plt.title('Comparing for Seed {0}'.format(seed_id), fontsize = 25)
    plt.xlim([0.4,45])
    plt.tick_params(which='both', width=2)
    plt.tick_params(which='major', length=7)
    plt.tick_params(which='minor', length=4)#, color='r')
#    plt.grid()
    ax.set_ylabel("profiles of Convergence map (log10 space)", fontsize=20)
    ax.set_xlabel("Pixels",fontsize=20)
    if gridspace == 'log':
        ax.set_xscale('log')
    ax.set_xticks([1, 5, 10,20,30,40, 45])    
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())    
    plt.tick_params(labelsize=25)
#    plt.grid(which="minor")
#    s, Rein_r = cal_gamma(convergence_true_low_res)
#    print "True gamma", s
#    print "model gamma", cal_gamma(convergence_map_model)[0]
#    plt.plot(np.linspace(-1,2)*0+Rein_r, np.linspace(-1,2),'blue')
#    plt.text(Rein_r-8, -0.55, 'Einstein Radius',fontsize=20)
#    return ax, Rein_r

for i in range(0,6):
    map_folder = "/Users/Dartoon/Astro/Lens_Model_challenge/TDLMC_material/mock_data/rung2.5/lenstronomy/mock_with_kappa/f160w-seed{0}/".format(211+i)
    kappa = pyfits.open(map_folder+'kappa.fits')[0].data.copy()
    print "seed", 211+i , cal_gamma(kappa)