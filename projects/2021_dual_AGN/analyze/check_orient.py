#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 15:49:36 2022

@author: Dartoon
"""

from ID_list import ID_list
import glob
import astropy.io.fits as pyfits
from tools import cal_oreination

f = open("material/ID_RA_DEC_z.txt","r")
string = f.read()
zlines = string.split('\n')   # Split in to \n
show_plot = True

def draw_line(x, y, angle):
    angle = angle-90 + 180
    angle = angle/180 * np.pi
    r = 40  # or whatever fits you
    plt.arrow(x, y,-r*np.cos(angle), -r*np.sin(angle),width=2)

# ID = [
#  #      '022404.85+014941.9',
#  # '022906.04-051428.9',
#   # '092532.13-020806.1',
#  # '095218.04-000459.1',
#  # '105458.01+043310.6',
#   # '122144.31-004144.1',
#  # '124618.51-001750.2',
#  # '150216.66+025719.8',
#  # '162501.98+430931.6',
#   '220642.82+003016.2',
#  # '230402.77-003855.4'
#  ][0]
for ID in ID_list:
# for ID in [ID]:
    adding_ori = 0
    if ID == '022906.04-051428.9':
        adding_ori = 2
    if ID == '092532.13-020806.1':
        adding_ori = 4 #+5
    if ID == '124618.51-001750.2':
        adding_ori = 7
    if ID == '162501.98+430931.6':
        adding_ori = 5
    if ID == '230402.77-003855.4':
        adding_ori = 5 +180
    if ID == '122144.31-004144.1':
        adding_ori = 0 #- 30
    if ID == '162501.98+430931.6':
        adding_ori = 0 +8 # To avoid the left source
    pair_ori = cal_oreination(ID)
    APT_orie_1 = pair_ori+135
    APT_orie_1 = round(APT_orie_1) + adding_ori  #The value to input to APT to get the aperture ori as shown
    if APT_orie_1 >360:
        APT_orie_1 = APT_orie_1-360
    if APT_orie_1< 180:
        APT_orie_2 = APT_orie_1 + 180
    elif APT_orie_1> 180:
        APT_orie_2 = APT_orie_1 - 180
    search_ori = pair_ori + 90  + adding_ori  
    folder_1 = glob.glob('../proof2close_HSC_images_5band/*/'+ ID+ '/')
    folder = folder_1[0]
    
    files = glob.glob(folder+'*_HSC-I.fits')
    fitsFile = pyfits.open(files[0])
    fov_image= fitsFile[1].data
    header = fitsFile[1].header # if target position is add in WCS, the header should have the wcs information, i.e. header['EXPTIME']
    
    
    from galight.tools.measure_tools import measure_bkg
    bkglight = measure_bkg(fov_image, if_plot=False)
    fov_image = fov_image-bkglight
    #Check out the flux distribution at the particular angle.
    
    #!!! Need to check the define of ori
    # Search ori should be the +90 to the pair_ori. Thus, when search_ori = 0, it should be to North, so that
    # so when 
    from galight.tools.astro_tools import plt_fits
    from astropy.coordinates import Angle
    from regions import PixCoord
    from regions import RectanglePixelRegion, CirclePixelRegion
    
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    norm = LogNorm(vmin=None, vmax=None)#np.max(img[~np.isnan(img)]))
    show_im = fov_image[40:-40,40:-40]
    reg = RectanglePixelRegion(center=PixCoord(x=len(show_im)/2, y=len(show_im)/2),
                                      width=7, height=155*0.13/0.168*2,   #155*0.13
                                      angle=Angle(search_ori, 'deg'))
    _aperture = CirclePixelRegion(PixCoord(len(show_im)/2, len(show_im)/2), 10.)
    _mask = _aperture.to_mask()
    aperture = _mask.to_image(show_im.shape)
    _mask = reg.to_mask()
    mask = _mask.to_image(show_im.shape)
    
    import numpy as np
    fluxs_within = show_im * mask * (1-aperture)
    print(ID, APT_orie_1)
    print('exit flux rate:', round(np.sum(fluxs_within)/np.sum(show_im* aperture ),2), "\nflux max pixel:", round(np.max(fluxs_within) ,2))
    if show_plot == True:
        fig, ax = plt.subplots(figsize=None)
        plt.imshow(show_im, norm=norm, origin='lower') 
        patch = reg.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                         label='Rectangle')
        ax.legend(handles=(patch,), loc='upper center')
        patch = _aperture.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                         label='aperture')
        ax.set_aspect('equal')
        draw_line(len(show_im)/2+10, y=len(show_im)/2, angle = search_ori)
        plt.show()
            # plt.imshow(fluxs_within, origin='lower')
            # plt.show()
    # print((show_im * mask).max(),fluxs_within.max() )
    # print( np.sum(show_im* aperture ), np.sum(fluxs_within))
    # plt.hist(np.reshape(fluxs_within, (1,-1))[0])
    # plt.show()
    
