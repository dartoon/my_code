#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 21:20:09 2017

@author: dxh
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'


potential_mismatch_list = []
geometry_mismatch_list = []

    #folder_type = 'sim_lens_ID_'
    #file_type = 'model_result.pkl'

# folder_type = 'simulations_700_subg30/sim_lens_ID_subg30_'
# file_type = 'result_PSFerr001_PSFinter_subg3.pkl'

# folder_type = 'simulations_700_subg30/sim_lens_noqso_ID_subg30_'
# file_type = 'result_PSFerr001_subg3.pkl'

folder_type = 'AGN_result_folder/idx*_ID*_PSFnoerr_PSFinter.pkl'
folder_type = 'AGN_result_folder/idx*_ID*_PSFerr001_PSFinter.pkl'
folder_type = 'AGN_result_folder/idx*_ID*_PSFerr025_PSFinter.pkl'
# folder_type = 'AGN_result_folder_run3times/idx*_ID*noPSFerr_noPSFinter_morePSO*.pkl'

# folder_type = 'AGN_result_folder/idx*_ID*_PSFerr025_PSFinter.pkl'

import glob
folder_list = glob.glob(folder_type+'*')
folder_list.sort()
# test_numer = 50 #len(folder_list)
# folder_list = folder_list[:test_numer]

# outlier = [736, 728, 748, 710, 715]
# folder_list = [folder_list[i] for i in range(len(folder_list)) if int(folder_list[i][-3:]) not in outlier]

folder_list_list = ['AGN_result_folder_run3times/idx0_ID702_PSFerr025_notPSFinter_morePSO_0.pkl', #Best of 3 for PSF 50%
 'AGN_result_folder_run3times/idx10_ID713_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx11_ID714_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx12_ID715_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx13_ID717_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx14_ID718_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx15_ID719_PSFerr025_notPSFinter_morePSO_2.pkl',
 'AGN_result_folder_run3times/idx16_ID720_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx17_ID721_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx18_ID724_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx19_ID727_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx1_ID703_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx20_ID728_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx21_ID729_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx22_ID730_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx23_ID731_PSFerr025_notPSFinter_morePSO_2.pkl',
 'AGN_result_folder_run3times/idx24_ID733_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx25_ID734_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx26_ID735_PSFerr025_notPSFinter_morePSO_2.pkl',
 'AGN_result_folder_run3times/idx27_ID736_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx28_ID737_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx29_ID738_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx2_ID705_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx30_ID739_PSFerr025_notPSFinter_morePSO_2.pkl',
 'AGN_result_folder_run3times/idx31_ID740_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx32_ID741_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx33_ID742_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx34_ID743_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx35_ID744_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx36_ID745_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx37_ID747_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx38_ID748_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx39_ID749_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx3_ID706_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx40_ID751_PSFerr025_notPSFinter_morePSO_2.pkl',
 'AGN_result_folder_run3times/idx41_ID752_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx42_ID753_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx43_ID754_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx44_ID755_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx45_ID757_PSFerr025_notPSFinter_morePSO_2.pkl',
 'AGN_result_folder_run3times/idx46_ID759_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx47_ID760_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx4_ID707_PSFerr025_notPSFinter_morePSO_1.pkl',
 'AGN_result_folder_run3times/idx5_ID708_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx6_ID709_PSFerr025_notPSFinter_morePSO_0.pkl',
 'AGN_result_folder_run3times/idx7_ID710_PSFerr025_notPSFinter_morePSO_2.pkl',
 'AGN_result_folder_run3times/idx8_ID711_PSFerr025_notPSFinter_morePSO_2.pkl',
 'AGN_result_folder_run3times/idx9_ID712_PSFerr025_notPSFinter_morePSO_1.pkl']

for folder in folder_list_list:    
    #==============================================================================
    # ##### lens mass model 
    #==============================================================================
    from lenstronomy.LensModel.lens_model import LensModel
    lens_model_list = ['PEMD','SHEAR']
    lens_model_class = LensModel(lens_model_list)
#    print(folder)
    #Load the true parameter:
        
    ID = folder.split('ID')[1][:3]
    sim_folder = 'simulations_700_subg30/' + 'sim_lens_ID_subg30_' + ID +'/'
        
    model_lists, para_s, lens_info= pickle.load(open(sim_folder+'/sim_kwargs.pkl','rb'))  
    kwargs_lens_list = para_s[0]
    if len(para_s[-1]['ra_image']) < 4:
#        print(para_s[-1]['ra_image'])
        para_s[-1]['ra_image'], para_s[-1]['dec_image'] = para_s[-1]['ra_image'][:2], para_s[-1]['dec_image'][:2]
    x_image, y_image = para_s[-1]['ra_image'], para_s[-1]['dec_image']
    source_pos = [para_s[2][0]['center_x'], para_s[2][0]['center_y']]
    potential = lens_model_class.potential(x_image, y_image, kwargs=kwargs_lens_list)
    fer_pot = lens_model_class.fermat_potential(x_image, y_image, kwargs_lens_list)
    geometry = potential + fer_pot
    # geometry = ((x_image - source_pos[0])**2 + (y_image - source_pos[1])**2) / 2.
#    print("potential:", potential, "geometry:", geometry)
#    
    #Load the inferred parameter:  
    read_file = folder #[:-1]
    multi_band_list, kwargs_model, kwargs_result, chain_list, fix_setting, mcmc_new_list = pickle.load(open(read_file,'rb'))
    # print(kwargs_result['kwargs_lens'][0])
    # x_image_infe, y_image_infe = x_image, y_image
    # kwargs_lens_list_infe = kwargs_lens_list
    
    kwargs_lens_list_infe = kwargs_result['kwargs_lens']
    if kwargs_result['kwargs_ps'] != []:
        x_image_infe, y_image_infe = kwargs_result['kwargs_ps'][0]['ra_image'], kwargs_result['kwargs_ps'][0]['dec_image']
    # else:    
    #     ID = folder[-3:]
    #     qso_folder = folder_type[:-16] + 'ID_subg30_{0}/'.format(ID)
    #     _, _, kwargs_result_withQSO, _, _, _ = pickle.load(open(qso_folder+'model_result_calNoiseMap_modNoisemap_boostPossionx3_subg3.pkl','rb'))
    #     x_image_infe, y_image_infe = kwargs_result_withQSO['kwargs_ps'][0]['ra_image'], kwargs_result_withQSO['kwargs_ps'][0]['dec_image']
    
    source_pos_infe = [kwargs_result['kwargs_source'][0]['center_x'], kwargs_result['kwargs_source'][0]['center_y']]
    potential_infe = lens_model_class.potential(x_image_infe, y_image_infe, kwargs=kwargs_lens_list_infe)
    # geometry_infe = ((x_image_infe - source_pos_infe[0])**2 + (y_image_infe - source_pos_infe[1])**2) / 2.
    fer_pot_infe = lens_model_class.fermat_potential(x_image_infe, y_image_infe, kwargs_lens_list_infe)
    geometry_infe = potential_infe + fer_pot_infe
#    print(round(source_pos_infe[0],4), round(source_pos_infe[1],4))
#    print("potential:", potential_infe, "geometry:", geometry_inf)
    true_potential_diff = potential[0]- potential[1:]
    true_geometry_diff = geometry[0]- geometry[1:]
    infe_potential_diff = potential_infe[0]- potential_infe[1:]
    infe_geometry_diff = geometry_infe[0]- geometry_infe[1:]
#    print("A-BCD true: potential, geometry:", true_potential_diff, true_geometry_diff)
#    print("A-BCD infer: potential, geometry:",infe_potential_diff, infe_geometry_diff)
#    print("Mismatch:", true_potential_diff- infe_potential_diff, true_geometry_diff- infe_geometry_diff) 
    for i in range(len(true_potential_diff)):
        potential_mismatch_list.append(infe_potential_diff[i] - true_potential_diff[i])
        geometry_mismatch_list.append(infe_geometry_diff[i]- true_geometry_diff[i])
    # potential_mismatch_list.append(np.average( true_potential_diff - infe_potential_diff ))
    # geometry_mismatch_list.append(np.average( true_geometry_diff - infe_geometry_diff ))

#%%
if 'noqso' in folder_type:
    text = 'non-AGN'
else:
    text = 'AGN'
plt.figure(figsize=(8,8))
high0, x0, _ = plt.hist(potential_mismatch_list,density=False, histtype=u'step',
         label=('Shapiro delay mismatch'), linewidth = 2, color='orange')
high1, x1, _ = plt.hist(geometry_mismatch_list,density=False, histtype=u'step',
         label=('geometric delay mismatch'), linewidth = 2, color='green')

plt.title('Lensed {0} case'.format(text), fontsize=27)    
plt.xlabel("mismatch (inferred - truth)",fontsize=27)
plt.ylabel("Distribution density",fontsize=27)
plt.xlim(-0.025, 0.025)
plt.ylim(0,400)
plt.tick_params(which='both', width=2, length = 7, labelsize=20)
plt.legend(prop={'size':20})
plt.savefig('Fermat_bias_dis_{0}.pdf'.format(text), bbox_inches = "tight")
plt.show()

print("np.mean(geometry_mismatch_list):", np.mean(geometry_mismatch_list))
print("np.std(geometry_mismatch_list):", np.std(geometry_mismatch_list))

print("np.mean(potential_mismatch_list):", np.mean(potential_mismatch_list))
print("np.std(potential_mismatch_list):", np.std(potential_mismatch_list))
