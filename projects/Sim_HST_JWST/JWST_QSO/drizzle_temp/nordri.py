import numpy as np
import astropy.io.fits as pyfits

import glob
nimage = len(glob.glob('non_drizzled-image-*.fits'))

#scaling back to cps
filename='drimage.fits'
im = pyfits.open(filename)[0].data.copy()
im /=nimage
pyfits.PrimaryHDU(im).writeto('Drz_QSO_image.fits')#,overwrite=True)

filename='dripsf.fits'
psf = pyfits.open(filename)[0].data.copy()
psf /= np.sum(psf)  
pyfits.PrimaryHDU(psf).writeto('Drz_PSF.fits')#,overwrite=True)

filename='drAGNclean.fits'
AGN = pyfits.open(filename)[0].data.copy()
AGN /= nimage
pyfits.PrimaryHDU(AGN).writeto('Drz_AGNclean_image.fits')#,overwrite=True)

filename='drHOSTclean.fits'
HOST = pyfits.open(filename)[0].data.copy()
HOST /=nimage
pyfits.PrimaryHDU(HOST).writeto('Drz_HOSTclean_image.fits')#,overwrite=True)

filename='drPOINTclean.fits'
POINT = pyfits.open(filename)[0].data.copy()
POINT /=nimage
pyfits.PrimaryHDU(POINT).writeto('Drz_POINTclean_image.fits')#,overwrite=True)
