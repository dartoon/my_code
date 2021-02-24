#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 09:36:10 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob

files = glob.glob('../dual_AGN/*/proof-2close.pdf')
from shutil import copyfile
for i in range(len(files)):
    file = files[i]
    ID = file.split('/')[2]
    newfile = ID + '_proof-2close.pdf'
    copyfile(file, newfile)
    
# /Users/Dartoon/Downloads/2close/011110.89+001709.2_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/113803.73+031457.8_HSC-I_proof-2close.pdf  #Lensing
# /Users/Dartoon/Downloads/2close/220713.35-000805.5_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/014630.53+001955.3_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/010748.49-000740.4_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/022258.94-034459.4_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/132441.58-015401.8_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/133222.62+034739.9_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/221101.45+001449.0_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/012110.93+010703.3_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/011935.29-002033.5_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/020721.92-005826.1_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/023600.28-010432.3_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/014536.01-000611.1_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/105458.01+043310.6_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/090654.53+021315.2_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/092604.90+012558.7_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/153008.91+425634.8_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/150216.66+025719.8_HSC-I_proof-2close.pdf #Like merging
# /Users/Dartoon/Downloads/2close/022906.04-051428.9_HSC-I_proof-2close.pdf #Like merging
# /Users/Dartoon/Downloads/2close/131618.38+045401.0_HSC-I_proof-2close.pdf #Like merging
# /Users/Dartoon/Downloads/2close/003659.43-001850.2_HSC-I_proof-2close.pdf #Like m
# /Users/Dartoon/Downloads/2close/124618.51-001750.2_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/093159.42+033622.0_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/084710.40-001302.6_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/101138.50+012344.6_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/090559.70+022408.2_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/122634.31+051024.0_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/000050.56-013055.2_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/235914.39-003428.4_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/122144.31-004144.1_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/110245.36+044806.3_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/095218.04-000459.1_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/230402.77-003855.4_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/022404.85+014941.9_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/162501.98+430931.6_HSC-I_proof-2close.pdf #WOW merging
# /Users/Dartoon/Downloads/2close/131512.46+015021.6_HSC-I_proof-2close.pdf #Wow WOw
# /Users/Dartoon/Downloads/2close/012727.51+003145.6_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/233718.07+002550.6_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/152445.62+440949.5_HSC-I_proof-2close.pdf #Lensing?
# /Users/Dartoon/Downloads/2close/090550.30+003948.1_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/224236.31-013121.1_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/012716.17-003557.6_HSC-I_proof-2close.pdf
# /Users/Dartoon/Downloads/2close/141637.44+003352.2_HSC-I_proof-2close.pdf  #Found before! 
# /Users/Dartoon/Downloads/2close/100043.13+020637.2_HSC-I_proof-2close.pdf #!!! Something new
# /Users/Dartoon/Downloads/2close/132951.99-020331.4_HSC-I_proof-2close.pdf #!!!
# /Users/Dartoon/Downloads/2close/092919.16+041414.2_HSC-I_proof-2close.pdf #Interesting
# /Users/Dartoon/Downloads/2close/022105.64-044101.5_HSC-I_proof-2close.pdf