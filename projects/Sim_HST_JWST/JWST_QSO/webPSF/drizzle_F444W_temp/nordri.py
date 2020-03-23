import numpy as np
import astropy.io.fits as pyfits

import glob
nimage = len(glob.glob('non_drizzled-image-*.fits'))

#scaling back to cps
id_num = 9
for i in range(id_num):
	filename='dripsf_id{0}.fits'.format(i)
	psf = pyfits.open(filename)[0].data.copy()
	psf /= np.sum(psf)  
	pyfits.PrimaryHDU(psf).writeto('Drz_PSF_id{0}.fits').format(i)#,overwrite=True)
