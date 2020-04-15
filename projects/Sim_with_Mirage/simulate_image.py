#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 22:49:33 2020

Copy from: https://github.com/spacetelescope/mirage/blob/master/examples/Imaging_simulator_use_examples.ipynb
"""
import os

os.environ["MIRAGE_DATA"] = '/Volumes/Seagate_Expansion_Drive/Mirage_data/mirage_data/nircam'
#os.environ["CRDS_DATA"] = "/user/myself/crds_cache"
#os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

# For examining outputs
import glob
from scipy.stats import sigmaclip
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
#%matplotlib inline

# mirage imports
#from mirage import imaging_simulator
#from mirage.seed_image import catalog_seed_image
#from mirage.dark import dark_prep
#from mirage.ramp_generator import obs_generator
#from mirage.yaml import yaml_generator

#%%
#from mirage.catalogs import catalog_generator
#if glob.glob('example/point_sources.cat') == []:
#    print("Generate point source cat")
#    ra_list = np.random.random(10) + 53.5
#    dec_list = np.random.random(10) - 67.2
#    
#    nrc_f200w_mag = np.random.random(10) + 16.
#    nrc_f212n_mag = np.random.random(10) + 19.
#    nis_f090w_mag = np.random.random(10) + 15.5
#    
#    ptsrc = catalog_generator.PointSourceCatalog(ra=ra_list, dec=dec_list)
#    ptsrc.add_magnitude_column(nrc_f200w_mag, instrument='nircam', filter_name='F200W', magnitude_system='abmag')
#    ptsrc.add_magnitude_column(nrc_f212n_mag, instrument='nircam', filter_name='F212N', magnitude_system='abmag')
#    ptsrc.add_magnitude_column(nis_f090w_mag, instrument='niriss', filter_name='F090W', magnitude_system='abmag')
#    ptsrc.save('example/point_sources.cat')
#    
#    x_list = np.random.random(10) * 2048
#    y_list = np.random.random(10) * 2048
#    
#    ptsrc = catalog_generator.PointSourceCatalog(x=x_list, y=y_list)
#    ptsrc.add_magnitude_column(nrc_f200w_mag, instrument='nircam', filter_name='F200W', magnitude_system='abmag')
#    ptsrc.add_magnitude_column(nrc_f212n_mag, instrument='nircam', filter_name='F212N', magnitude_system='abmag')
#    ptsrc.add_magnitude_column(nis_f090w_mag, instrument='niriss', filter_name='F090W', magnitude_system='abmag')
#    ptsrc.save('example/point_sources_xy.cat')

#if glob.glob('example/galaxies.cat') == []:
#    print("Generate galaxy cat")    
#    ra_list = np.random.random(10) + 53.5
#    dec_list = np.random.random(10) - 67.2
#    ellipticity = np.random.random(10) * 0.75
#    sersic_index = np.random.random(10) * 4.
#    position_angle = np.random.random(10) * 359.
#    
#    nrc_f200w_mag = np.random.random(10) + 16.
#    nrc_f212n_mag = np.random.random(10) + 19.
#    nis_f090w_mag = np.random.random(10) + 15.5
#    
#    gal = catalog_generator.GalaxyCatalog(ra=ra_list, dec=dec_list, ellipticity=ellipticity,
#                                          sersic_index=sersic_index, position_angle=position_angle)
#    gal.add_magnitude_column(nrc_f200w_mag, instrument='nircam', filter_name='F200W', magnitude_system='abmag')
#    gal.add_magnitude_column(nrc_f212n_mag, instrument='nircam', filter_name='F212N', magnitude_system='abmag')
#    gal.add_magnitude_column(nis_f090w_mag, instrument='niriss', filter_name='F090W', magnitude_system='abmag')
#    gal.save('example/galaxies.cat')

#%%
from mirage.imaging_simulator import ImgSim
sim = ImgSim(paramfile='example/my_yaml_file.yaml')
print("Simulation setting up DONE!")
sim.create()