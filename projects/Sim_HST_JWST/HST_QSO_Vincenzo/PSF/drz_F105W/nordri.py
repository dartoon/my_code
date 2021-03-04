import numpy as np
import astropy.io.fits as pyfits

import glob
nimage = len(glob.glob('non_drizzled-image-*.fits'))

#scaling back to cps

filename='dripsf.fits'
psf = pyfits.open(filename)[0].data.copy()
psf /= np.sum(psf)  
pyfits.PrimaryHDU(psf).writeto('Drz_PSF.fits')#,overwrite=True)
