
from numpy import *
import astropy.io.fits as pyfits
import glob
num = 8 #8 dither pattern
id_num = 9
for i in range(id_num):
	for j in range(num):
		psf=pyfits.open('non_drizzled_psf_id{0}-{1}.fits'.format(i, j+1),mode = 'update')
		psf[0].header["EXPTIME"]=1
		psf.flush()