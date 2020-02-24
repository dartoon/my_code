
from numpy import *
import astropy.io.fits as pyfits
for i in range(8):
	hdu=pyfits.open('non_drizzled-lens-image-{0}.fits'.format(i+1),mode = 'update')
	hdu[0].header["EXPTIME"]=1
	hdu.flush()
	psf=pyfits.open('non_drizzled_psf-{0}.fits'.format(i+1),mode = 'update')
	psf[0].header["EXPTIME"]=1
	psf.flush()