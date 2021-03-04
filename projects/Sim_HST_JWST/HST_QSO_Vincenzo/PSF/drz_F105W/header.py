
from numpy import *
import astropy.io.fits as pyfits
import glob
num = 8 #8 dither pattern
for i in range(num):
	psf=pyfits.open('non_drizzled_psf-{0}.fits'.format(i+1),mode = 'update')
	psf[0].header["EXPTIME"]=1
	psf.flush()