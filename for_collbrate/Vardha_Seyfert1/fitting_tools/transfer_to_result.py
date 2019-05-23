#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 15:59:47 2018

@author: dartoon

Transfer the fit to a result dict
"""
import numpy as np
def transfer_to_result(data, source_result, ps_result, image_ps, image_host,
                       error_map, zp, fixcenter,ID, cut=0, plot_compare=True,
                       QSO_msk_list = "",pix_sz = 0.168, QSO_msk = None,
                       tag=None):
    #==============================================================================
    # Translate the e1, e2 to phi_G and q
    #==============================================================================
    import lenstronomy.Util.param_util as param_util
    source_result[0]['phi_G'], source_result[0]['q'] = param_util.ellipticity2phi_q(source_result[0]['e1'], source_result[0]['e2'])
    
    #==============================================================================
    # Save the result
    #==============================================================================
    from roundme import roundme
    import copy
    result = copy.deepcopy(source_result[0])
    del result['e1']
    del result['e2']
    result['QSO_amp'] = ps_result[0]['point_amp'][0]
    result['host_amp'] = image_host[0].sum()
    result['host_flux_ratio_percent']= result['host_amp']/(result['QSO_amp'] + result['host_amp'])*100
    zp = zp
    result['host_mag'] = - 2.5 * np.log10(result['host_amp']) + zp 
    #print "The host flux is ~:", image_host.sum()/(image_ps.sum() + image_host.sum())
    # =============================================================================
    # Save QSO position to result if not fix center
    # =============================================================================
    if fixcenter == False:
        result['qso_x'] = ps_result[0]['ra_image'][0]
        result['qso_y'] = ps_result[0]['dec_image'][0]
    
    ##==============================================================================
    ##Plot the images for adopting in the paper
    ##==============================================================================
    from flux_profile import total_compare
    data = data
    QSO = image_ps
    if QSO_msk is None:
        QSO_mask = data*0 + 1
    else:
        QSO_mask = QSO_msk
    if len(image_host) == 1:
        host = image_host[0]
        label = ['data', 'QSO', 'host', 'model', 'normalized residual']
    elif len(image_host) >1:
        host = np.zeros_like(image_host[0])
        for i in range(len(image_host)):
            host += image_host[i]
        label = ['data', 'QSO', 'host+{0}objs'.format(i), 'model', 'normalized residual']  #Print the numbers of objects
    flux_list = [data, QSO, host, error_map]
    import glob
    mask_list = glob.glob(QSO_msk_list)   # Read *.reg files in a list.
#    print "mask_list,muhahah", mask_list
    fig = total_compare(label_list = label, flux_list = flux_list, target_ID = ID, pix_sz=pix_sz,
                  data_mask_list = mask_list, data_cut = cut, zp = zp,
                  plot_compare = plot_compare, msk_image = QSO_mask)
    if tag is not None:
        fig.savefig("{0}_SB_profile.pdf".format(tag), bbox_inches = 'tight')
    fig1 = total_compare(label_list = label, flux_list = flux_list, target_ID = ID, pix_sz=pix_sz,
                  data_mask_list = mask_list, data_cut = cut, zp = zp,
                  plot_compare = plot_compare, msk_image = QSO_mask, if_annuli=True)
    if tag is not None:
        fig1.savefig("{0}_SB_profile_annuli.pdf".format(tag), bbox_inches = 'tight')
    
    # =============================================================================
    # Calculate reduced Chisq and save to result
    # =============================================================================
    chiq_map = ((data-image_ps-host)/(error_map))**2 * QSO_mask
    pixels=len(error_map)**2 - (1-QSO_mask).sum()
    reduced_Chisq = chiq_map.sum()/pixels
    result=roundme(result)
    result['redu_Chisq'] = round(reduced_Chisq,6)
    return result

def transfer_obj_to_result(source_result,image_host, zp):
    #==============================================================================
    # Translate the e1, e2 to phi_G and q
    #==============================================================================
    import lenstronomy.Util.param_util as param_util
    source_result['phi_G'], source_result['q'] = param_util.ellipticity2phi_q(source_result['e1'], source_result['e2'])
    from roundme import roundme
    import copy
    result = copy.deepcopy(source_result)
    del result['e1']
    del result['e2']
    result['host_amp'] = image_host.sum()
    result['host_mag'] = - 2.5 * np.log10(result['host_amp']) + zp 
    result=roundme(result)
    return result


def transfer_to_result_ps(data, ps_result, image_ps, image_host,
                       error_map, filt, fixcenter,ID, cut=0, plot_compare=True,
                       QSO_msk_list = "",pix_sz = 'drz06', QSO_msk = None,
                       tag=None):
    #==============================================================================
    # Translate the e1, e2 to phi_G and q
    #==============================================================================
    #==============================================================================
    # Save the result
    #==============================================================================
    from roundme import roundme
    result = {}
    ##==============================================================================
    ##Plot the images for adopting in the paper
    ##==============================================================================
    from flux_profile import total_compare
    data = data
    QSO = image_ps
    if QSO_msk is None:
        QSO_mask = data*0 + 1
    else:
        QSO_mask = QSO_msk
    if len(image_host) == 0:
        host = data*0
        label = ['data', 'QSO', 'host=0', 'model', 'normalized residual']
    elif len(image_host) >0:
        host = np.zeros_like(image_host[0])
        for i in range(len(image_host)):
            host += image_host[i]
        label = ['data', 'QSO', '{0}objs'.format(i), 'model', 'normalized residual']  #Print the numbers of objects
    flux_list = [data, QSO, host, error_map]
    import glob
    mask_list = glob.glob(QSO_msk_list)   # Read *.reg files in a list.
#    print "mask_list,muhahah", mask_list
    fig = total_compare(label_list = label, flux_list = flux_list, target_ID = ID, pix_sz=pix_sz,
                  data_mask_list = mask_list, data_cut = cut, zp = zp,
                  plot_compare = plot_compare, msk_image = QSO_mask)
    if tag is not None:
        fig.savefig("{0}_SB_profile.pdf".format(tag), bbox_inches = 'tight')
    # =============================================================================
    # Calculate reduced Chisq and save to result
    # =============================================================================
    chiq_map = ((data-image_ps-host)/(error_map))**2 * QSO_mask
    pixels=len(error_map)**2 - (1-QSO_mask).sum()
    reduced_Chisq = chiq_map.sum()/pixels
    result=roundme(result)
    result['redu_Chisq'] = round(reduced_Chisq,6)
    return result
