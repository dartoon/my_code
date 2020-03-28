import numpy as np
import astropy.io.fits as pyfits

filename='drimage.fits'
nimage=8.
im = pyfits.open(filename)[0].data.copy()
#scaling back to cps
im /=nimage
pyfits.PrimaryHDU(im).writeto('Drz_QSO_image.fits')#,overwrite=True)

filename='dripsf.fits'
psf = pyfits.open(filename)[0].data.copy()
psf /= np.sum(psf)  
pyfits.PrimaryHDU(psf).writeto('Drz_PSF.fits')#,overwrite=True)