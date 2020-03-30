
from numpy import *
import astropy.io.fits as pyfits
import glob
num = len(glob.glob('non_drizzled-image-*.fits'))
for i in range(num):
	hdu=pyfits.open('non_drizzled-image-{0}.fits'.format(i+1),mode = 'update')
	hdu[0].header["EXPTIME"]=1
	hdu.flush()
#	psf=pyfits.open('non_drizzled_psf-{0}.fits'.format(i+1),mode = 'update')
#	psf[0].header["EXPTIME"]=1
#	psf.flush()
	hdu=pyfits.open('non_drizzled-AGNclean-{0}.fits'.format(i+1),mode = 'update')
	hdu[0].header["EXPTIME"]=1
	hdu.flush()
	hdu=pyfits.open('non_drizzled-HOSTclean-{0}.fits'.format(i+1),mode = 'update')
	hdu[0].header["EXPTIME"]=1
	hdu.flush()
	hdu=pyfits.open('non_drizzled-POINTclean-{0}.fits'.format(i+1),mode = 'update')
	hdu[0].header["EXPTIME"]=1
	hdu.flush()	