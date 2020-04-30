#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 23:05:15 2020

@author: Dartoon
"""

file_list= ['0_fitting_084710.40.py','1_fitting_121405.12.py', '2_fitting_141637.44.py', '3_fitting_220906.91.py',  '4_fitting_233713.66.py', '5_fitting_021930.51.py', '6_fitting_022105.64.py']

# for i in [0]:
# for i in [1,3]:
# for i in [2]:
# for i in [4, 5]:    
for i in [6]:
    runfile('/Users/Dartoon/Astro/Projects/my_code/for_collbrate/Shenli_host_masses/MCMC_fits/{0}'.format(file_list[i]), wdir='/Users/Dartoon/Astro/Projects/my_code/for_collbrate/Shenli_host_masses/MCMC_fits')