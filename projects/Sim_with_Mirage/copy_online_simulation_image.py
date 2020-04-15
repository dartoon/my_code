#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 22:49:33 2020

Copy from: https://github.com/spacetelescope/mirage/blob/master/examples/Imaging_simulator_use_examples.ipynb
"""
import os

os.environ["MIRAGE_DATA"] = '/Volumes/Seagate_Expansion_Drive/Mirage_data/'
#os.environ["CRDS_DATA"] = "/user/myself/crds_cache"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

# For examining outputs
from glob import glob
from scipy.stats import sigmaclip
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
#%matplotlib inline

# mirage imports
from mirage import imaging_simulator
from mirage.seed_image import catalog_seed_image
from mirage.dark import dark_prep
from mirage.ramp_generator import obs_generator
from mirage.yaml import yaml_generator

#%%
# Specify the xmnl and pointing files exported from APT
xml_file = 'example/example_imaging_program.xml'
pointing_file = 'example/example_imaging_program.pointing'

# Source catalogs to be used.
cat_dict = {'GOODS-S-FIELD': {'point_source':'imaging_example_data/ptsrc_catalog.cat'}}

# Set reference file values. 
# Setting to 'crds_full_name' will search for and download needed
# calibration reference files (commonly referred to as CRDS reference files) when
# the yaml_generator is run. 
# 
# Setting to 'crds' will put placeholders in the yaml files and save the downloading
# for when the simulated images are created.
reffile_defaults = 'crds'

# Optionally set the cosmic ray library and rate
cosmic_rays = {'library': 'SUNMAX', 'scale': 1.0}

# Optionally set the background signal rates to be used
background = 'medium'

# Optionally set the telescope roll angle (PAV3) for the observations
pav3 = 12.5
# Optionally set the observation date to use for the data. Note that this information
# is placed in the headers of the output files, but not used by Mirage in any way.
dates = '2022-10-31'

#%%
# Set the directory into which the yaml files will be written
output_dir = './imaging_example_data/'

# You can also set a separate directory where the simulated data
# will eventually be saved to
simulation_dir = './imaging_example_data/'

#%%
datatype = 'linear, raw'
# Run the yaml generator
yam = yaml_generator.SimInput(input_xml=xml_file, pointing_file=pointing_file,
                              catalogs=cat_dict, cosmic_rays=cosmic_rays,
                              background=background, roll_angle=pav3,
                              dates=dates, reffile_defaults=reffile_defaults,
                              verbose=True, output_dir=output_dir,
                              simdata_output_dir=simulation_dir,
                              datatype=datatype)
yam.create_inputs()
yfiles = glob(os.path.join(output_dir,'jw*.yaml'))

#%%
# Choose one of the yaml files just created
yamlfile = 'imaging_example_data/jw42424001001_01101_00001_nrcb1.yaml'

# Run all steps of the imaging simulator for yaml file #1
img_sim = imaging_simulator.ImgSim()
img_sim.paramfile = yamlfile
img_sim.create()

#%%
def show(array,title,min=0,max=1000):
    plt.figure(figsize=(12,12))
    plt.imshow(array,clim=(min,max))
    plt.title(title)
    plt.colorbar().set_label('DN$^{-}$/s')
    
# First, look at the noiseless seed image
show(img_sim.seedimage,'Seed Image',max=20)    

#%%
lin_file = 'imaging_example_data/jw42424001001_01101_00001_nrcb1_linear.fits'
with fits.open(lin_file) as hdulist:
    linear_data = hdulist['SCI'].data
print(linear_data.shape)

show(linear_data[0, 3, :, :], "Final Group", max=250)

raw_file = 'imaging_example_data/jw42424001001_01101_00001_nrcb1_uncal.fits'
with fits.open(raw_file) as hdulist:
    raw_data = hdulist['SCI'].data
print(raw_data.shape)

show(raw_data[0, 3, :, :], "Final Group", max=15000)

show(1. * raw_data[0, 3, :, :] - 1. * raw_data[0, 0, :, :], "Last Minus First Group", max=200)