#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 10:52:29 2023

@author: Dartoon
"""

import astropy.io.fits as pyfits
import numpy as np
import os, shutil
import matplotlib.pyplot as plt
import glob
# %matplotlib inline/''

# from hstaxe import axetasks

cwd = os.getcwd()
print("We are in %s" % (cwd))

#%%
name = "SDSSJ1246-0017"  #G102
# name = "SDSSJ1502+0257"    #G141
# name = "SDSSJ1625+4309"    #G141
# name = "SDSSJ2304-0038"    #G102
# name = "SDSSJ2206+0030"    #G141
# name = "SDSSJ0224+0149"    #G102


z_dic = {"SDSSJ1246-0017": 2.559, "SDSSJ1502+0257": 1.484, "SDSSJ1625+4309":1.647, 
            "SDSSJ2304-0038": 2.775, "SDSSJ2206+0030": 1.444, "SDSSJ0224+0149":2.721}

z = z_dic[name]

CIV, MgII, Hb, OIII, Halpha = 1549, 2798, 4861, 5007, 6563

#%%
run_reduce = False
if run_reduce == True:
    from drizzlepac import astrodrizzle
    os.chdir(cwd)
    if os.path.isdir(name):
        shutil.rmtree(name)
    os.mkdir(name)

    os.system("cp ../MAST_download/*_{0}/iexo*/*flt.fits ./{0}/".format(name))
    os.chdir(name)
        
    spec_filenames = []
    image_filenames = []
    for filename in glob.glob("*flt.fits"):
        fitsfile = pyfits.open(filename)
        print(fitsfile[0].header['FILTER'])
        if fitsfile[0].header['FILTER'][0] == 'F':   #F105W or F140W
            image_filenames.append(filename)
            filt_image = fitsfile[0].header['FILTER']
        elif fitsfile[0].header['FILTER'][0] == 'G': #G102 or G104
            spec_filenames.append(filename)
            filt_spec = fitsfile[0].header['FILTER']
            
    if filt_spec == 'G141':
        x1, x2 = 11500, 16500
    else:
        x1, x2 = 8000, 11500

    os.mkdir(filt_spec)
    os.mkdir(filt_image)
    os.chdir('./'+filt_spec)
    import shutil
    write_file = open("{0}.lis".format(filt_spec),'w') 
    for filename in spec_filenames:
        write_file.write("{0}\n".format(filename))
        shutil.move('../'+filename, './')
    write_file.close()
    astrodrizzle.AstroDrizzle("@{0}.lis".format(filt_spec),output=name+'_'+filt_spec,build=True)
    
    os.chdir('../'+filt_image)
    write_file = open("{0}.lis".format(filt_image),'w') 
    for filename in image_filenames:
        write_file.write("{0}\n".format(filename))
        shutil.move('../'+filename, './')
    write_file.close()
    ref = '../'+filt_spec+'/'+name+'_'+filt_spec + '_drz.fits[1]'
    astrodrizzle.AstroDrizzle("@{0}.lis".format(filt_image),output=name+'_'+filt_image,
                              in_memory=False,skysub="yes",build=True,driz_cr_corr=True,driz_cr=True,
                              final_wcs=True,driz_separate=True,driz_sep_wcs=True,driz_sep_refimage=ref,
                              final_refimage=ref)
    os.chdir(cwd+'/'+name)
    from hstaxe import axetasks
    
    if os.path.isdir("CONF"):
        shutil.rmtree("CONF")
    os.mkdir("CONF")
    
    os.system("cp ../../config_file/WFC3.IR.{0}.cal.V4.32/* CONF/".format(filt_spec))
    #%%
    os.chdir(cwd+'/'+name)
    
    if os.path.isdir("DATA"):
        shutil.rmtree("DATA")
    os.mkdir("DATA")
    os.environ['AXE_IMAGE_PATH'] = './DATA/' 
    print ('--> variable AXE_IMAGE_PATH   set to "./DATA/"')
    
    os.environ['AXE_CONFIG_PATH'] = './CONF/'
    print ('--> variable AXE_CONFIG_PATH  set to "./CONF/"')
    
    if os.path.isdir("OUTPUT"):
        shutil.rmtree("OUTPUT")
    os.mkdir("OUTPUT")
    os.environ['AXE_OUTPUT_PATH'] = './OUTPUT/'
    print ('--> variable AXE_OUTPUT_PATH  set to "./OUTPUT/"')
    
    print ("Length of AXE_IMAGE_PATH is",len(os.environ['AXE_IMAGE_PATH']),"characters")

    # dimension_info = "0,0,0,0"
    dimension_info = "183,85,50,50"
    os.system("cp {0}/*flt.fits DATA/".format(filt_spec))
    os.system("cp {0}/*flt.fits DATA/".format(filt_image))
    
    os.chdir(cwd+'/'+name)
    os.chdir(filt_image)
    os.system("cp ../../../cookbook_{0}.cat ./cookbook.cat".format(name))
    
    axetasks.iolprep(drizzle_image='{0}_{1}_drz.fits'.format(name, filt_image),
                         input_cat='cookbook.cat',
                         dimension_in=dimension_info)
    os.system("cp *_1.cat ../DATA/")
    
    os.chdir(cwd+'/'+name + '/DATA')
    # write_file = open("aXe.lis",'w') 
    img_times = []
    for filename in image_filenames:
        fitsfile = pyfits.open(filename)
        img_times.append(fitsfile[0].header['EXPEND'])
        # print(fitsfile[1].header['CRVAL1'], fitsfile[1].header['CRVAL2'])  # To confirm if the position is in right order.
    image_filenames = [image_filenames[np.argsort(img_times)[i]] for i in range(4)]
    
    spec_times = []
    for filename in spec_filenames:
        fitsfile = pyfits.open(filename)
        spec_times.append(fitsfile[0].header['EXPEND'])
        # print(fitsfile[1].header['CRVAL1'], fitsfile[1].header['CRVAL2'])  #
    spec_filenames = [spec_filenames[np.argsort(spec_times)[i]] for i in range(4)]
    # write_file.close()
    
    os.chdir(cwd+'/'+name)
    write_file = open("aXe.lis",'w') 
    for i in range(4):
        write_file.write("{0} {1} {2} \n".format(spec_filenames[i], image_filenames[i].replace('.fits', '_1.cat'), image_filenames[i]))
    write_file.close()
    
    axetasks.axeprep(inlist="aXe.lis",
                         configs="{0}.{1}.V4.32.conf".format(filt_spec, filt_image),
                         backgr=True,
                         backims="WFC3.IR.{0}.sky.V1.0.fits".format(filt_spec),
                         norm=False,
                         mfwhm=3.0)
    
    axetasks.axecore('aXe.lis',
                     "{0}.{1}.V4.32.conf".format(filt_spec, filt_image),
                     extrfwhm=4.,
                     drzfwhm=3.,
                     backfwhm=0.,
                     orient=False,
                     weights=True,
                     slitless_geom=False,
                     cont_model='gauss',  #fluxcube #gauss
                     sampling='drizzle',
                     exclude=True)

    
    ID = 1
    
    filenames = glob.glob('./OUTPUT/*STP.fits')
    
    plt.rcParams["figure.figsize"] = (15,3)
    plt.subplot(2,2,1)
    try:
        d1 = pyfits.open(filenames[0])["BEAM_%dA" % (ID)].data
        im1 = plt.imshow(d1,origin="lower")
        im1.set_clim(0,.1)
        print(np.sum(d1))
    except:
        pass
    
    plt.subplot(2,2,2)
    try:
        d1 = pyfits.open(filenames[1])["BEAM_%dA" % (ID)].data
        im1 = plt.imshow(d1,origin="lower")
        im1.set_clim(0,.1)
        # print(np.sum(d1))
    except:
        pass
    
    plt.subplot(2,2,3)
    try:
        d1 = pyfits.open(filenames[2])["BEAM_%dA" % (ID)].data
        im1 = plt.imshow(d1,origin="lower")
        im1.set_clim(0,.1)
        # print(np.sum(d1))
    except:
        pass
    
    plt.subplot(2,2,4)
    try:
        d1 = pyfits.open(filenames[3])["BEAM_%dA" % (ID)].data
        im1 = plt.imshow(d1,origin="lower")
        im1.set_clim(0,.1)
        # print(np.sum(d1))
    except:
        pass
    
    plt.show()
    
    import glob
    for s in glob.glob("OUTPUT/*2.SPC.fits"):
        print( s)
        d1 = pyfits.open(s)["BEAM_%dA" % (ID)].data
        w = d1["LAMBDA"]
        f = d1["FLUX"]
        e = d1["FERROR"]
        # plt.errorbar(w,f,e)
        # vg = (w>11500) & (w<16500)
        vg = (w>x1) & (w<x2)
        plt.errorbar(w[vg],f[vg],e[vg])
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel(r'Flux ($erg/s/cm^2/\AA/s$)');
    plt.show()
    
    # os.mkdir('./DRIZZLE/tmp')
    
    os.chdir(cwd+'/'+name)
    opt_extr=False
    
    axetasks.drzprep(inlist = "aXe.lis",
                configs = "{0}.{1}.V4.32.conf".format(filt_spec, filt_image),
                back = False,opt_extr=opt_extr)
    axetasks.axecrr(inlist="aXe.lis",
        configs="{0}.{1}.V4.32.conf".format(filt_spec, filt_image),
        infwhm = 4.0,
        outfwhm = 3.0,
        back = False,
        driz_separate = 'yes',
        opt_extr=opt_extr
        )
#%%
if run_reduce == False:
    os.chdir(cwd+'/'+name)
    try:
        filt_spec = glob.glob('G102')[0]
        filt_image = glob.glob('F105W')[0]
    except:
        filt_spec, filt_image = 'G141', 'F140W'

from galight.tools.astro_tools import plt_fits
from matplotlib.colors import LogNorm
import copy, matplotlib
ID = 2
my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')
d = pyfits.open("./DRIZZLE/aXeWFC3_{0}_2.STP.fits".format(filt_spec))["BEAM_%dA" % (ID)].data
plt.imshow(d, norm=LogNorm(), origin='lower',cmap = my_cmap)
# im.set_clim(0,0.1)
# ID = 2
# d = pyfits.open("./DRIZZLE/aXeWFC3_G102_2.STP.fits")["BEAM_%dA" % (ID)].data
# im = plt.imshow(d)
# im.set_clim(0,0.1)
plt.show()

fin = pyfits.open("./DRIZZLE/aXeWFC3_{0}_2.SPC.fits".format(filt_spec))
tdata = fin["BEAM_%dA" % (ID)].data
x = tdata["LAMBDA"]
f = tdata["FLUX"]
e = tdata["FERROR"]
c = tdata["CONTAM"]

if filt_spec == 'G141':
    x1, x2 = 11500, 16500
else:
    x1, x2 = 8000, 11500

vg = (x>x1) & (x<x2)
plt.plot(x[vg],f[vg])
plt.axvline(x = MgII * (1+z), color = 'b', label = 'MgII')
plt.axvline(x = Hb * (1+z), color = 'b', label = 'Hb')
plt.axvline(x = OIII * (1+z), color = 'b', label = 'OIII')
plt.errorbar(x[vg],f[vg],e[vg])
plt.plot(x[vg],c[vg])
plt.xlim([x1-100,x2+100])
plt.show()


import glob

# for s in glob.glob("OUTPUT/*2.SPC.fits"):
#     print (s)
#     d1 = pyfits.open(s)["BEAM_%dA" % (ID)].data
#     w = d1["LAMBDA"]
#     f = d1["FLUX"]
#     e = d1["FERROR"]
#     c = d1["CONTAM"]
#     wg = (w>x1) & (w<x2)
#     plt.errorbar(w[wg],f[wg],e[wg])
    # plt.plot(w[vg],c[vg])
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel(r'Flux ($erg/s/cm^2/\AA/s$)');


fin = pyfits.open("./DRIZZLE/aXeWFC3_{0}_2.SPC.fits".format(filt_spec))
tdata = fin["BEAM_%dA" % (ID)].data
x = tdata["LAMBDA"]
f = tdata["FLUX"]
e = tdata["FERROR"]
c = tdata["CONTAM"]
#plt.errorbar(x[vg],y[vg],e[vg])
plt.plot(x[vg],f[vg],color='k',lw=2)
plt.errorbar(x[vg],f[vg],e[vg],color='k',lw=2)
plt.xlim([x1-100,x2+100])
plt.show()
#%%
image_file = glob.glob(filt_image+'/*_drz.fits')[0]
image = pyfits.open(image_file)['SCI'].data
pos = np.loadtxt(filt_image+'/cookbook.cat')[0][:2]
plt.imshow(image[int(pos[1])-20:int(pos[1])+20,int(pos[0])-20:int(pos[0])+20], 
          norm=LogNorm(), origin='lower',cmap = my_cmap)
plt.show()


fits = pyfits.open("./DRIZZLE/aXeWFC3_{0}_mef_ID{1}.fits".format(filt_spec, ID))
d = fits['SCI'].data
header = fits['SCI'].header

# fits = pyfits.open("./DRIZZLE/aXeWFC3_{0}_2.STP.fits".format(filt_spec))["BEAM_%dA" % (ID)]
# d = fits.data
# wcs.all_pix2world([[0,0]])

from astropy.wcs import WCS
wcs = WCS(header, naxis=1, relax=False, fix=False)
lam = wcs.wcs_pix2world(np.arange(len(d.T)), 0)[0]

# plt.imshow(d[10:-2,:], norm=LogNorm(), origin='lower',cmap = my_cmap)
y = 10
y = 7
d = d[y:y+4,:]
fig, ax = plt.subplots(1,1,figsize=(12.5, 2))
ax.imshow(d, norm=LogNorm(), origin='lower',cmap = my_cmap)
ticks = [25,50,75,100,125,150,175]
ax.set_xticks(ticks)
ax.set_xticklabels([round(lam[i]*10**8)*10**2 for i in ticks])
plt.show()

# plt.plot(np.linspace(0,len(np.sum(d,axis=0))-1,len(np.sum(d,axis=0))),np.sum(d,axis=0))
# plt.plot((lam*10**10)[(lam*10**10>x1) & (lam*10**10<x2)],(np.sum(d[10:-7,:],axis=0)/(lam*10**7))[(lam*10**10>x1) & (lam*10**10<x2)])
plt.plot((lam*10**10)[(lam*10**10>x1) & (lam*10**10<x2)],(np.sum(d,axis=0)/(lam*10**7))[(lam*10**10>x1) & (lam*10**10<x2)])
lines = CIV, MgII, Hb, OIII, Halpha
for line in lines:
    if line * (1+z) > x1 and line * (1+z) < x2:
        plt.axvline(x = line * (1+z), color = 'b', label = str(line))
        plt.text(line * (1+z), np.max((np.sum(d,axis=0)/(lam*10**7))[(lam*10**10>x1) & (lam*10**10<x2)]), str(line))
# plt.text(MgII * (1+z), np.max((np.sum(d,axis=0)/(lam*10**7))[(lam*10**10>x1) & (lam*10**10<x2)]), 'MgII')
# plt.axvline(x = Hb * (1+z), color = 'b', label = 'Hb')
# plt.text(Hb * (1+z), np.max((np.sum(d,axis=0)/(lam*10**7))[(lam*10**10>x1) & (lam*10**10<x2)]), 'Hb')
# plt.axvline(x = OIII * (1+z), color = 'b', label = 'OIII')
# plt.text(OIII * (1+z), np.max((np.sum(d,axis=0)/(lam*10**7))[(lam*10**10>x1) & (lam*10**10<x2)]), 'OIII')
plt.xlim([x1-100,x2+100])
plt.tight_layout()
plt.show()

