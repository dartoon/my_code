#!/usr/bin/env python
# import drizzlepac
# from drizzlepac import astrodrizzle
# #unlearn astrodrizzle
# astrodrizzle.AstroDrizzle('*i2d.fits',output='final', build='Yes', static='Yes',
#                           skysub='No',driz_separate='Yes', median='Yes', blot='Yes',
#                           driz_cr='Yes',driz_combine='Yes',final_wcs='Yes',
#                           final_bits=576,final_scale=0.031,final_pixfrac=0.8,
#                           final_kernel='gaussian',final_rot=None, overwrite=True)

# #final_kernel='gaussian',

# from jwst.resample import ResampleStep
# ResampleStep(name = '*i2d.fits')


#%%
# import astropy.io.fits as pyfits
# from jwst.resample.resample import ResampleData
# from jwst import datamodels
# import glob
# input_models = []
# files= glob.glob('*i2d.fits')
# for file in files:
#     fitsFile = pyfits.open(file)
#     input_models.append(fitsFile)
    
# ResampleData(input_models, output='final',single=False,
#             blendheaders=True, pixfrac=1.0, kernel='square',
#             fillval='INDEF', weight_type='ivm', good_bits=0,
#             pscale_ratio=1.0, pscale=None)


#%%
import os
os.environ["CRDS_PATH"] = "/Users/Dartoon/Downloads/crds_cache"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

from jwst.pipeline import calwebb_image3

# jw02736-o001_20220712t161855_image3_003_asn

asn_file = 'jw02736-o001_20220712t161855_image3_003_asn.json'

# Create an instance of the pipeline class
image3 = calwebb_image3.Image3Pipeline()

# Set some parameters that pertain to the entire pipeline
image3.output_dir = 'output_dir'
image3.save_results = True

# image3.schema_url = None
# schema_url = 'http://stsci.edu/schemas/jwst_datamodel/image.schema'

# Set some parameters that pertain to some of the individual steps
# Turn off TweakRegStep
image3.tweakreg.skip = True  
# Turn off SkyMatchStep
image3.skymatch.skip = True
# Set the ratio of input to output pixels to create an output mosaic 
# on a 0.015"/pixel scale
image3.resample.pixel_scale_ratio = 0.48

# Call the run() method
image3.run(asn_file)