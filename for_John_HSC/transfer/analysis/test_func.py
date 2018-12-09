#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:17:42 2018

@author: Dartoon
"""

from gen_fit_ingredient import gen_fit_ingredient
#==============================================================================
# Test the gen_fit_ingredient
#==============================================================================
zp_list = []
pix_scale_list = []
qso_center_list = []
for i in range(1,291):
    print i
    QSO_im, err_map, PSF, pix_scale, zp, qso_ID, qso_fr_center,_ = gen_fit_ingredient(i)
    zp_list.append(zp)
    pix_scale_list.append(pix_scale)
    qso_center_list.append(qso_fr_center)
#    print "qso_ID", qso_ID, 'pix_scale', pix_scale, '\n'

for i in range(290):
    print i, pix_scale_list[i]