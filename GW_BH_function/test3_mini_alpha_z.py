#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 23:15:38 2019

@author: Dartoon

Test 3, considering the alpha is evolving with redshift.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
from BH_mass_function import gene_BHBH, dl, solve_z

a, mbh_max, mbh_min = 2.35, 80., 5.
#filename = 'test3_sim_a_{0}_max_{1}_min_{2}'.format(round(a,2), round(mbh_max,1), round(mbh_min,1))
#if_file = glob.glob(filename)  

test = gene_BHBH(h0=70)
event_rate0, zs_detected0, masses0, rhos_detected0 = test.mc_year_rate(a=a, mbh_max=mbh_max, mbh_min=mbh_min, ev_type = 'const')
event_rate1, zs_detected1, masses1, rhos_detected1 = test.mc_year_rate(a=[2.35,0.1], mbh_max=mbh_max, mbh_min=mbh_min, ev_type = 'mode0')
event_rate2, zs_detected2, masses2, rhos_detected2 = test.mc_year_rate(a=[2.35,0.7], mbh_max=mbh_max, mbh_min=mbh_min, ev_type = 'mode1')

