#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 14:51:15 2022

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from scipy import ndimage
import scipy.optimize as op
from galight.tools.astro_tools import plt_fits


def shift_img(img, shift_pix, order=1):
    shift_pix = shift_pix[::-1]  #uniform the yx to xy
    from scipy.ndimage.interpolation import shift
    shifted_digit_image=shift(img, shift_pix, order = order)
    return shifted_digit_image
def rotate_image(img, rotate_pix, order =1):
    shift_pix = [-rotate_pix[0]*2, -rotate_pix[1]*2]
    shift_ = shift_img(img, shift_pix, order=order)
    rotate = np.flip(shift_)
    return rotate

# %%
class Test_asy:
    """
    Measure the asymmetry of a image.
    """
    def __init__(self,img=None, segm = None, segm_id = 1, order=3, seg_cal_reg = 'or'):  #order relies on the background rms
        self.img = img
        if isinstance(segm, (np.ndarray)) or segm is None:
            self.segm = segm
        else:
            self.segm = segm.data
        self.segm_id = segm_id
        self.order = order
        self.seg_cal_reg = seg_cal_reg
    def residual(self, rotate_pix):
        rotate_ = rotate_image(self.img, rotate_pix, order = self.order)
        residual = self.img - rotate_  #Consider resdiual as data-model, where model is the rotation.
        return residual
    
    def abs_res(self, rotate_pix, if_plot=False):
        cal_areas, _, punish = self.segm_to_mask(rotate_pix)
        res_ = self.residual(rotate_pix)
        if if_plot == True:
            plt_fits(abs(res_*cal_areas),norm='log')
        if punish == False:
            return np.sum(abs(res_*cal_areas))
        else:
            return 10**6
        
    def find_pos(self, ini_pix = [0.0,0.0]):
        result = op.minimize(self.abs_res, ini_pix, method='nelder-mead',
                options={'xatol': 1e-8, 'disp': True})
        return result
    
    def segm_to_mask(self, rotate_pix, segm_id = None):
        if segm_id is None:
            segm_id = self.segm_id
        cal_area = self.segm == segm_id
        mask = (self.segm != segm_id) * (self.segm != 0)
        # plt.imshow(cal_area, origin='lower')
        # plt.show()
        rotate_pix = np.around(rotate_pix)
        cal_area_ = rotate_image(cal_area, rotate_pix,order =1)
        mask_ = rotate_image(mask, rotate_pix,order =1)
        punish = False
        if self.seg_cal_reg == 'and':
            cal_areas = cal_area * cal_area_
            mask_areas = mask * mask_
            if np.sum(cal_areas) < np.sum(cal_area)/3:
                punish = True
        elif self.seg_cal_reg == 'or':
            cal_areas = cal_area + cal_area_
            mask_areas = mask + mask_
        return cal_areas, mask_areas, punish
    
    def cal_asymmetry(self, rotate_pix):
        asy = self.abs_res(rotate_pix)
        cal_areas, masks, _ = self.segm_to_mask(rotate_pix)
        obj_flux = np.sum( self.img * cal_areas)
        obj_masks = cal_areas + masks
        obj_masks = obj_masks == False
        img = self.img * obj_masks
        print(obj_masks.sum())
        print(cal_areas.sum())
        bkg_asy = np.sum(abs(img - np.flip(img)) * obj_masks)  #!!! The mask for bkg also change with pos.
        plt.imshow(abs(img - np.flip(img)) * obj_masks, origin='lower')
        plt.show()
        # return asy, obj_flux, bkg_asy
        print(np.sum(obj_masks), np.sum(cal_areas))
        return asy/obj_flux - bkg_asy/np.sum(obj_masks) * np.sum(cal_areas)/obj_flux
        return asy, obj_flux, bkg_asy, np.sum(obj_masks), np.sum(cal_areas)


# Q:
# Which mask to use as total flux of object? The obj_mask + obj_mask_?
# When calculate residual image, which object total flux should be used?
# How to norm the bkg asy?

# %%
image = pyfits.open('./test_img.fits')[0].data.copy()
plt_fits(image)
pos = [0, 5.2]
image_shift = shift_img(image, pos,order = 2)
# plt_fits(image_shift)
image_rotate = rotate_image(image, pos,order = 2)
# plt_fits(image_rotate)
from galight.tools.measure_tools import detect_obj
apertures, segm_deblend = detect_obj(image, if_plot= True, segm_map= True, nsigma=1, auto_sort_center=True)
# plt.imshow(segm_mask_rot, origin='lower')
# plt.show()
img_class=Test_asy(img=image,segm=segm_deblend, order=1, segm_id=4, seg_cal_reg= 'or')
# cal_areas, mask_areas,_ = img_class.segm_to_mask([0,0])
result = img_class.find_pos(ini_pix = [0,0])
res = img_class.abs_res(rotate_pix = result["x"], if_plot=True)
asymmetry = img_class.cal_asymmetry(rotate_pix = result["x"])
print(result["x"])
print(res)
print(asymmetry)

# %%
res = img_class.abs_res(rotate_pix = [0,0], if_plot=True)
print(result["x"])
print(res)

residual = img_class.residual(result['x'])
plt_fits(abs(residual))
plt.imshow(abs(residual), origin='lower')
plt.show()

