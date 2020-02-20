#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 09:48:36 2018

@author: Dartoon
"""

import numpy as np
import matplotlib.pyplot as plt
#If you have trouble accessing the image you can download it straight away using Python:

from astropy.extern.six.moves.urllib import request
url = "http://python4astronomers.github.com/_downloads/image2.fits"
open("image2.fits", "wb").write(request.urlopen(url).read())
#ls
from sherpa.astro.ui import *
image_data()