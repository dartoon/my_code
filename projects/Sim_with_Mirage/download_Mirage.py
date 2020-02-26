#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 22:49:33 2020

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from mirage.reference_files import downloader
download_path = '/Volumes/Seagate_Expansion_Drive/Mirage_data/'
downloader.download_reffiles(download_path, instrument='NIRCam', dark_type='linearized', skip_darks=False, skip_cosmic_rays=False, skip_psfs=False, skip_grism=False)
#
#from mirage.imaging_simulator import ImgSim
