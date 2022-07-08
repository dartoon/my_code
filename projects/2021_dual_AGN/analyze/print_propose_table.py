#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 17:58:21 2021

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from ID_list import ID_list
import glob

from tools import read_string_list, read_info, cal_oreination

f = open("material/ID_RA_DEC_z.txt","r")
string = f.read()
zlines = string.split('\n')   # Split in to \n
# def read_info(ID):
    # line = [zlines[i] for i in range(len(zlines)) if ID in zlines[i]]
    # if line != []:
    #     z = float(line[0].split(' ')[-1])
    # else:
    #     z = -99
    # return z
# Halpha = 6562.8

CIV, MgII, Hb, OIII = 1549, 2798, 4861, 5007
G102range = [8000, 11500]
G141range = [10750, 17000]
def av_filter(z):
    lines = np.array([CIV, MgII, Hb, OIII])
    redshift_lines = (1+z) * lines
    G102_bool =  (redshift_lines>G102range[0]+100) * (redshift_lines<G102range[1]-100)
    G141_bool =  (redshift_lines>G141range[0]+100) * (redshift_lines<G141range[1]-100)
    # return G102_bool, G141_bool
    s1 = np.array(['CIV', 'MgII', 'H$beta$', '[OIII]'])[G102_bool] 
    s2 = np.array(['CIV', 'MgII', 'H$beta$', '[OIII]'])[G141_bool] 
    s1 = [s1[i] for i in range(len(s1))]
    s2 = [s2[i] for i in range(len(s2))]
    # str1 = "G102: " + repr(s1)
    # str2 = " G141: " + repr(s2)    
    # s = str1 + str2
    if s2 != []:
        try:
            s = "G141 & " + s2[0] + '+' + s2[1]
        except:
            s = "G141 & " + s2[0]
    elif s2 == [] and s1 != []:
        try: 
            s = "G102 & " + s1[0] + '+' +  s1[1]
        except:
            s = "G102 & " + s1[0]
    else:
        s = "No fileter!!! & "
    return s

from astropy.cosmology import FlatLambdaCDM
# cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def deg2HMS(ra='', dec='', round=False):
  RA, DEC, rs, ds = '', '', '', ''
  if dec:
    if str(dec)[0] == '-':
      ds, dec = '-', abs(dec)
    deg = int(dec)
    decM = abs(int((dec-deg)*60))
    if round:
      decS = int((abs((dec-deg)*60)-decM)*60)
    else:
      decS = (abs((dec-deg)*60)-decM)*60
    DEC = '{0}{1} {2} {3}'.format(ds, deg, decM, decS)
  if ra:
    if str(ra)[0] == '-':
      rs, ra = '-', abs(ra)
    raH = int(ra/15)
    raM = int(((ra/15)-raH)*60)
    if round:
      raS = int(((((ra/15)-raH)*60)-raM)*60)
    else:
      raS = ((((ra/15)-raH)*60)-raM)*60
    RA = '{0}{1} {2} {3}'.format(rs, raH, raM, raS)
  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC

print("ID & RA & DEC & z & Separ.& Mag.& Grism & Lines & PA\\\\")
print("&&&&($''$, kpc)&(pair)&&&(deg)\\\\")
offset_kpc_h0_list, z_list = [], []
for i, ID in enumerate(ID_list):
    show_ID = ID[:4] + ID[9] +  ID[10:14]
    line = [zlines[i] for i in range(len(zlines)) if ID in zlines[i]]
    line[0] = line[0].replace('  ', ' ')
    z = float(line[0].split(' ')[-1])
    RA, Dec = line[0].split(' ')[1], line[0].split(' ')[2]
    # RA, Dec, z = read_info(ID)
    
    # files_1 = glob.glob('../proof2close_HSC_images_5band/*/' + ID + '/fit_result/')
    # files_2 = glob.glob('../extra/*/fit_result*/' + ID + '/')
    # files = files_1 + files_2
    
    add = '' #The first fit is better than this 'deep_'
    # add = 'deep_'
    folder_1 = glob.glob('../proof2close_HSC_images_5band/*/'+ ID+ '/')
    if folder_1 != []: # and 'z_below1' not in folder_1[0]:
        folder = folder_1[0] + add+'fit_result/'
        file = folder + add+'fit_result_I-band.txt'  #!!!
        # folder = folder_1[0] + 'fit_result/'
        # file = folder + 'fit_result_I-band.txt'
    # elif folder_1 != [] and 'z_below1' in folder_1[0]:
    #     folder = '../_John_fitted/'+ID+'_HSC-I/'   #For these z below 1(or z unkonwn), not fitted and use John's fit.
    #     file = folder + 'fit_result.txt'
    # else:
    #     folder_2 = glob.glob('../proofBHBH/model_Iband_zover1/'+ ID + '*/') 
    #     folder = folder_2[0]
    #     file = folder + 'fit_result_I-band.txt'
    
    # file = glob.glob(files[-1]+'fit_result_{0}-band.txt'.format('I'))
    if file != []:
        f = open(file,"r")    
    string = f.read()    
    lines = string.split('\n')   # Split in to \n
    trust = 2
    l1 = [i for i in range(len(lines)) if 'model_PS_result:' in lines[i]]
    AGN_dic = read_string_list(string = lines[l1[trust]].split('model_PS_result: ')[1])    
    AGN_pos = np.array([[-1*AGN_dic[i]['ra_image'], AGN_dic[i]['dec_image']] for i in range(len(AGN_dic))])    
    offset = np.sum( (AGN_pos[0] -  AGN_pos[1])**2)**0.5 
    scale_relation = cosmo.angular_diameter_distance(z).value * 10**3 * (1/3600./180.*np.pi)  #Kpc/arc
    offset_kpc = offset * scale_relation   #In kpc
    Mags = [AGN_dic[0]['magnitude'], AGN_dic[1]['magnitude']]
    APT_orie_1 = cal_oreination(ID,add=add,trust=trust)+135
    if APT_orie_1 >360:
        APT_orie_1 = APT_orie_1-360
    if APT_orie_1< 180:
        APT_orie_2 = APT_orie_1 + 180
    elif APT_orie_1> 180:
        APT_orie_2 = APT_orie_1 - 180
        
    mag1, mag2 = np.min(Mags), np.max(Mags) 
    zp = 27
    flux = 10**( 0.4*(zp-mag1)) + 10**( 0.4*(zp-mag2)) 
    tmag = -2.5*np.log10(flux) + zp
    print(i+1, show_ID, round(tmag,1), z, '& {0:.1f}'.format(cal_oreination(ID,add=add, trust = trust)), '%{0:.0f} {1:.0f}degree'.format(APT_orie_1, APT_orie_2) ) #'{0:.1f},{1:.1f}'.format(np.min(Mags), np.max(Mags)), '&', 
          # av_filter(z), '& {0:.1f} \\\\'.format(cal_oreination(ID)), '%{0:.0f} {1:.0f}degree'.format(APT_orie_1, APT_orie_2) )
    # print(deg2HMS(ra=float(RA), dec = float(Dec)) )
    
    offset_kpc_h0_list.append(offset_kpc * 70 / 100)
    z_list.append(z)
    
#%% Check how image looks
ID = [
    # '022404.85+014941.9',
    # '022906.04-051428.9',
    # '092532.13-020806.1',
  #  '095218.04-000459.1',
    # '105458.01+043310.6',
  #  '122144.31-004144.1',
    # '124618.51-001750.2',
    # '150216.66+025719.8',
    # '162501.98+430931.6',
    # '220642.82+003016.2',
    '230402.77-003855.4'
  ][0]

folder_1 = glob.glob('../proof2close_HSC_images_5band/*/'+ ID+ '/')
if folder_1 != []: # and 'z_below1' not in folder_1[0]:
    folder = folder_1[0] + 'deep_'+'fit_result/'
    file_glob = folder + 'fit_I-band_fit*pkl'  #!!!
    folder_ = folder_1[0] +'fit_result/'
    
f = open(folder_+'fit_result_I-band.txt',"r")    
string = f.read()    
lines = string.split('\n')   # Split in to \n
trust = 2
l1 = [i for i in range(len(lines)) if 'model_PS_result:' in lines[i]]
AGN_dic = read_string_list(string = lines[l1[trust]].split('model_PS_result: ')[1])    
AGN_pos = np.array([[-1*AGN_dic[i]['ra_image'], AGN_dic[i]['dec_image']] for i in range(len(AGN_dic))])    

file = glob.glob(file_glob)[0]
import pickle
fit_run = pickle.load(open(file,'rb'))
from galight.tools.astro_tools import plt_fits
plt_fits(fit_run.flux_2d_out['data'], hold=True)
c = len(fit_run.flux_2d_out['data'])/2
pixscale = fit_run.fitting_specify_class.deltaPix
for i in [0,1]:
    # x, y = fit_run.final_result_ps[i]['ra_image'], fit_run.final_result_ps[i]['dec_image']
    x, y = AGN_pos[i][0], AGN_pos[i][1]
    x, y = x/pixscale, y/pixscale
    plt.scatter( c+x, c+y )
plt.show()

AGN_dic = read_string_list(string = lines[l1[trust]].split('model_PS_result: ')[1])    
if AGN_dic[0]['flux_within_frame'] < AGN_dic[1]['flux_within_frame']:
    AGN_dic[0], AGN_dic[1] = AGN_dic[1], AGN_dic[0]
AGN_pos = np.array([[1*AGN_dic[i]['ra_image'], AGN_dic[i]['dec_image']] for i in range(len(AGN_dic))])    
dif = AGN_pos[1]-AGN_pos[0]
PA = np.arctan( dif[0]/dif[1] ) * 180 / np.pi
if dif[1]<0 and dif[0]>0:
    PA = 180 + PA
if dif[1]<0 and dif[0]<0:
    PA = 180 + PA
if dif[1]>0 and dif[0]<0:
    PA = 360 + PA
print(PA)

print(cal_oreination(ID, trust=trust))

#!!! I checked and confirm that the table in the proposal give perfect informaiton. Using trust = 2

#%%
# #%%Check V band mag (inbetween G and R)
# IDs = [
#    '022404.85+014941.9',
#     '022906.04-051428.9',
#     '092532.13-020806.1',
#     '095218.04-000459.1',
#     '105458.01+043310.6',
#     '122144.31-004144.1',
#     '124618.51-001750.2',
#     '150216.66+025719.8',
#     '162501.98+430931.6',
#     '220642.82+003016.2',
#     '230402.77-003855.4'
#   ]

# for ID in IDs:
#     fluxs = []
#     for band in ['G', 'R']:
#         folder_1 = glob.glob('../proof2close_HSC_images_5band/*/'+ ID+ '/')
#         # if folder_1 != []: # and 'z_below1' not in folder_1[0]:
#             # folder = folder_1[0] + add+'fit_result/'
#             # file_glob = folder + 'fit_-band_fit*pkl'  #!!!
#         folder_ = folder_1[0] +'fit_result/'
#         try:
#             f = open(folder_+'fit_result_{0}-band.txt'.format(band),"r")    
#             string = f.read()    
#             lines = string.split('\n')   # Split in to \n
#             trust = 2  
#             l1 = [i for i in range(len(lines)) if 'model_PS_result:' in lines[i]]
#             AGN_dic = read_string_list(string = lines[l1[trust]].split('model_PS_result: ')[1])    
#             l1 = [i for i in range(len(lines)) if 'model_Sersic_result:' in lines[i]]
#             sersic_dic = read_string_list(string = lines[l1[trust]].split('model_Sersic_result: ')[1])    
#             flux = AGN_dic[0]['flux_within_frame']+AGN_dic[1]['flux_within_frame']+sersic_dic[0]['flux_within_frame']
#             fluxs.append(flux)
#         except:
#             # print(ID, band, 'not exist')
#             None
#     # print(ID, fluxs)
#     print(ID, round(-2.5*np.log10(np.average(fluxs)) + 27.0,1 ) )
        
