#!/usr/bin/env python
import drizzlepac
from drizzlepac import astrodrizzle
#unlearn astrodrizzle
astrodrizzle.AstroDrizzle('*flt.fits',output='final', build='Yes', static='Yes',
                          skysub='No',driz_separate='Yes', median='Yes', blot='Yes',
                          driz_cr='Yes',driz_combine='Yes',final_wcs='Yes',
                          final_bits=576,final_scale=0.0642,final_pixfrac=0.8,
                          final_kernel='gaussian',final_rot=None, overwrite=True)

#final_kernel='gaussian',


