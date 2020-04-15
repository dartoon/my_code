#!/usr/bin/env python
# coding: utf-8

# # Example: Using MIRAGE to Generate Imaging Exposures

# This notebook shows the general workflow for creating simulated data with Mirage, beginning with an APT file. For users without an APT file, Mirage will work with manually-created instrument/exposure parameter files, generally referred to as [input yaml files](#yaml_example). 

# Under the hood, the `Mirage` simulator is broken up into three basic stages:
# 
# 1. **Creation of a "seed image".**<br>
#    This is generally a noiseless countrate image that contains signal
#    only from the astronomical sources to be simulated. Currently, the 
#    mirage package contains code to produce a seed image starting
#    from object catalogs.<br><br>
#    
# 2. **Dark current preparation.**<br>
#    The simualted data will be created by adding the simulated sources
#    in the seed image to a real dark current exposure. This step
#    converts the dark current exposure to the requested readout pattern
#    and subarray size requested by the user.<br><br>
#    
# 3. **Observation generation.**<br>
#    This step converts the seed image into an exposure of the requested
#    readout pattern and subarray size. It also adds cosmic rays and 
#    Poisson noise, as well as other detector effects (IPC, crosstalk, etc).
#    This exposure is then added to the dark current exposure from step 2.<br><br>
#    
# For imaging mode observations, these steps are wrapped by the `imaging_simulator.py` module, as shown below.

# *Table of Contents:*
# * [Generating `yaml` files](#make_yaml)
# * [Create Simulated Data](#run_steps_together)
# * [Simulating Multiple Exposures](#mult_sims)
# * [Running Simulation Steps Independently](#run_steps_independently)
# * [Example `yaml` file](#yaml_example)

# ---
# ## Getting Started
# 
# <div class="alert alert-block alert-warning">
# **Important:** 
# Before proceeding, ensure you have set the MIRAGE_DATA environment variable to point to the directory that contains the reference files associated with MIRAGE.  
# <br/><br/>
# If you want JWST pipeline calibration reference files to be downloaded in a specific directory, you should also set the CRDS_DATA environment variable to point to that directory. This directory will also be used by the JWST calibration pipeline during data reduction.
# <br/><br/>
# You may also want to set the CRDS_SERVER_URL environment variable set to https://jwst-crds.stsci.edu. This is not strictly necessary, and Mirage will do it for you if you do not set it, but if you import the crds package, or any package that imports the crds package, you should set this environment variable first, in order to avoid an error.
# </div>

# In[3]:


import os


# In[4]:


os.environ["MIRAGE_DATA"] = "/Volumes/Seagate_Expansion_Drive/Mirage_data/mirage_data"
os.environ["CRDS_DATA"] = "/user/Dartoon/crds_cache"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"


# In[5]:


# For examining outputs
from glob import glob
from scipy.stats import sigmaclip
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# In[6]:


# mirage imports
from mirage import imaging_simulator
from mirage.seed_image import catalog_seed_image
from mirage.dark import dark_prep
from mirage.ramp_generator import obs_generator
from mirage.yaml import yaml_generator


# ---
# <a id='make_yaml'></a>
# # Generating input yaml files

# For convenience, observing programs with multiple pointings 
# and detectors can be simulated starting with the program's 
# APT file. The xml and pointings files must be exported from 
# APT, and are then used as input to the *yaml_generator*, which will
# generate a series of yaml input files.

# In[7]:


# Specify the xmnl and pointing files exported from APT
xml_file = 'imaging_example_data/example_imaging_program.xml'
pointing_file = 'imaging_example_data/example_imaging_program.pointing'


# ## Optional user inputs

# See Mirage's [Mirage's yaml_generator documentation](https://mirage-data-simulator.readthedocs.io/en/latest/yaml_generator.html#additional-yaml-generator-inputs "Yaml Generator Inputs")
# for details on the formatting options for the inputs listed below. The formats will vary based on the complexity of your inputs.

# In[8]:


# Source catalogs to be used.
cat_dict = {'GOODS-S-FIELD': {'point_source':'imaging_example_data/ptsrc_catalog.cat'}}


# In[9]:


# Set reference file values. 
# Setting to 'crds_full_name' will search for and download needed
# calibration reference files (commonly referred to as CRDS reference files) when
# the yaml_generator is run. 
# 
# Setting to 'crds' will put placeholders in the yaml files and save the downloading
# for when the simulated images are created.
reffile_defaults = 'crds'


# In[10]:


# Optionally set the cosmic ray library and rate
cosmic_rays = {'library': 'SUNMAX', 'scale': 1.0}


# In[11]:


# Optionally set the background signal rates to be used
background = 'medium'


# In[12]:


# Optionally set the telescope roll angle (PAV3) for the observations
pav3 = 12.5


# In[13]:


# Optionally set the observation date to use for the data. Note that this information
# is placed in the headers of the output files, but not used by Mirage in any way.
dates = '2022-10-31'


# ## Run the yaml_generator

# This will create a collection of yaml files that will be used as inputs when creating the simulated data. There will be one yaml file for each detector and exposure, so there can be quite a few files created if your program has lots of exposures or dithers.

# In[14]:


# Set the directory into which the yaml files will be written
output_dir = './imaging_example_data/'


# In[15]:


# You can also set a separate directory where the simulated data
# will eventually be saved to
simulation_dir = './imaging_example_data/'


# You can specify the data reduction state of the Mirage outputs.
# Options are 'raw', 'linear', or 'linear, raw'. 
# 
# If 'raw' is specified, the output is a completely uncalibrated file, with a filename ending in "uncal.fits"
# 
# If 'linear' is specified, the output is a file with linearized signals, ending in "linear.fits". This is equivalent to having been run through the dq_init, saturation flagging, superbias subtraction, reference pixel subtraction, and non-linearity correction steps of the calibration pipeline. Note that this product does not include dark current subtraction.
# 
# If 'linear, raw', both outputs are saved.
# 
# In order to fully process the Mirage output with the default steps used by the pipeline, it would be best to use the 'raw' output and run the entire calibration pipeline.

# In[16]:


datatype = 'linear, raw'


# In[ ]:


# Run the yaml generator
yam = yaml_generator.SimInput(input_xml=xml_file, pointing_file=pointing_file,
                              catalogs=cat_dict, cosmic_rays=cosmic_rays,
                              background=background, roll_angle=pav3,
                              dates=dates, reffile_defaults=reffile_defaults,
                              verbose=True, output_dir=output_dir,
                              simdata_output_dir=simulation_dir,
                              datatype=datatype)
yam.create_inputs()


# In[ ]:


yfiles = glob(os.path.join(output_dir,'*.yaml'))
print(yfiles)


# ---
# <a id='run_steps_together'></a>
# # Create simulated data

# ### The imaging simulator class

# The imaging_simulator.ImgSim class is a wrapper around the three main steps of the simulator (detailed in the [Running simulator steps independently](#run_steps_independently) section below). This convenience function is useful when creating simulated imaging mode data. WFSS data will need to be run in a slightly different way. See the WFSS example notebook for details.

# In[55]:


# Choose one of the yaml files just created
yamlfile = 'imaging_example_data/jw42424001001_01101_00001_nrcb1.yaml'


# In[56]:


# Run all steps of the imaging simulator for yaml file #1
img_sim = imaging_simulator.ImgSim()
img_sim.paramfile = yamlfile
img_sim.create()


# ### Examine the Output

# In[57]:


def show(array,title,min=0,max=1000):
    plt.figure(figsize=(12,12))
    plt.imshow(array,clim=(min,max))
    plt.title(title)
    plt.colorbar().set_label('DN$^{-}$/s')


# #### Noiseless Seed Image

# This image is an intermediate product. It contains only the signal from the astronomical sources and background. There are no detector effects, nor cosmic rays added to this count rate image.

# In[58]:


# First, look at the noiseless seed image
show(img_sim.seedimage,'Seed Image',max=20)


# #### Final Output Product

# Next examine the final output product. The linear output will make a prettier picture.

# In[ ]:


lin_file = 'imaging_example_data/jw42424001001_01101_00001_nrcb1_linear.fits'
with fits.open(lin_file) as hdulist:
    linear_data = hdulist['SCI'].data
print(linear_data.shape)


# In[ ]:


show(linear_data[0, 3, :, :], "Final Group", max=250)


# Examine the raw output. First a single group, which is dominated by noise and detector artifacts. 

# In[ ]:


raw_file = 'imaging_example_data/jw42424001001_01101_00001_nrcb1_uncal.fits'
with fits.open(raw_file) as hdulist:
    raw_data = hdulist['SCI'].data
print(raw_data.shape)


# In[ ]:


show(raw_data[0, 3, :, :], "Final Group", max=15000)


# Many of the instrumental artifacts can be removed by looking at the difference between two groups. Raw data values are integers, so first make the data floats before doing the subtraction.

# In[ ]:


show(1. * raw_data[0, 3, :, :] - 1. * raw_data[0, 0, :, :], "Last Minus First Group", max=200)


# This raw data file is now ready to be run through the [JWST calibration pipeline](https://jwst-pipeline.readthedocs.io/en/stable/) from the beginning. If dark current subtraction is not important for you, you can use Mirage's linear output, skip some of the initial steps of the pipeline, and begin by running the [Jump detection](https://jwst-pipeline.readthedocs.io/en/stable/jwst/jump/index.html?highlight=jump) and [ramp fitting](https://jwst-pipeline.readthedocs.io/en/stable/jwst/ramp_fitting/index.html) steps.

# ---
# <a id='mult_sims'></a>
# ## Simulating Multiple Exposures

# Each yaml file will simulate an exposure for a single pointing using a single detector. To simulate multiple exposures, or a single exposure with multiple detectors, multiple calls to the *imaging_simulator* must be made.

# ### In Series
# ```python
# paramlist = [yaml_a1,yaml_a2,yaml_a3,yaml_a4,yaml_a5]
# 
# def many_sim(paramlist):
#     '''Function to run many simulations in series
#     '''
#     for file in paramlist:
#         m = imaging_simulator.ImgSim()
#         m.paramfile = file
#         m.create()
# ```
# 
# ### In Parallel
# 
# Since each `yaml` simulations does not depend on the others, we can parallelize the process to speed things up:
# ```python
# from multiprocessing import Pool
# 
# n_procs = 5 # number of cores available
# 
# with Pool(n_procs) as pool:
#     pool.map(make_sim, paramlist)
# ```

# ---
# <a id='run_steps_independently'></a>
# # Running simulation steps independently

# The steps detailed in this section are wrapped by the `imaging_simulator` mentioned above. General users will not need to worry about the details of these three steps.

# Mirage is composed of three main steps:
# <br></br>
# <br></br>
# Seed image creation
# <br></br>
# Dark current preparation
# <br></br>
# Observation creation (combining simulated sources and dark current)
# <br></br><br></br>
# This section shows how to call the three steps independently. The `imaging_simulator` function above is a wrapper around these three steps. Most users will want simply to call `imaging_simulator`.

# ## First generate the "seed image" 
# 
# This is generally a 2D noiseless countrate image that contains only simulated astronomical sources.
# 
# A seed image is generated based on a `.yaml` file that contains all the necessary parameters for simulating data. For this exercise, use the same yaml file that was used in the [Create Simulated Data](#run_steps_together) section as input.

# In[ ]:


cat = catalog_seed_image.Catalog_seed()
cat.paramfile = yamlfile
cat.make_seed()


# ### Look at the seed image

# In[ ]:


show(cat.seedimage,'Seed Image',max=20)


# ## Prepare the dark current exposure
# This will serve as the base of the simulated data.
# This step will linearize the dark current (if it 
# is not already), and reorganize it into the 
# requested readout pattern and number of groups.

# In[ ]:


d = dark_prep.DarkPrep()
d.paramfile = yamlfile
d.prepare()


# ### Look at the dark current 
# For this, we will look at an image of the final group
# minus the first group

# In[ ]:


exptime = d.linDark.header['NGROUPS'] * cat.frametime
diff = (d.linDark.data[0,-1,:,:] - d.linDark.data[0,0,:,:]) / exptime
show(diff,'Dark Current Countrate',max=0.1)


# ## Create the final exposure
# Turn the seed image into a exposure of the 
# proper readout pattern, and combine it with the
# dark current exposure. Cosmic rays and other detector
# effects are added. 
# 
# The output can be either this linearized exposure, or
# a 'raw' exposure where the linearized exposure is 
# "unlinearized" and the superbias and 
# reference pixel signals are added, or the user can 
# request both outputs. This is controlled from
# within the yaml parameter file.

# In[ ]:


obs = obs_generator.Observation()
obs.linDark = d.prepDark
obs.seed = cat.seedimage
obs.segmap = cat.seed_segmap
obs.seedheader = cat.seedinfo
obs.paramfile = yamlfile
obs.create()


# ### Examine the final output image
# Again, we will look at the last group minus the first group

# In[ ]:


with fits.open(obs.linear_output) as h:
    lindata = h[1].data
    header = h[0].header


# In[ ]:


exptime = header['EFFINTTM']
diffdata = (lindata[0,-1,:,:] - lindata[0,0,:,:]) / exptime
show(diffdata,'Simulated Data',min=0,max=20)


# In[ ]:


# Show on a log scale, to bring out the presence of the dark current
# Noise in the CDS image makes for a lot of pixels with values < 0,
# which makes this kind of an ugly image. Add an offset so that
# everything is positive and the noise is visible
offset = 2.
plt.figure(figsize=(12,12))
plt.imshow(np.log10(diffdata+offset),clim=(0.001,np.log10(80)))
plt.title('Simulated Data')
plt.colorbar().set_label('DN$^{-}$/s')


# --- 
# <a id='yaml_example'></a>
# ## Example yaml input file
# 
# For an example of a yaml file, see the [example yaml file](https://mirage-data-simulator.readthedocs.io/en/latest/example_yaml.html "Example Yaml File") page
# in the Mirage documentation.
# 
# Entries listed as 'config' have default files that are present in the 
# config directory of the repository. The scripts are set up to 
# automatically find and use these files. The user can replace 'config'
# with a filename if they wish to override the default.
# 
# In general, if 'None' is placed in a field, then the step that uses
# that particular file will be skipped.
# 
# Note that the linearized_darkfile entry overrides the dark entry, unless
# linearized_darkfile is set to None, in which case the dark entry will be
# used.
# 
# Use of a valid readout pattern in the readpatt entry will cause the 
# simulator to look up the values of nframe and nskip and ignore the 
# values given in the yaml file.

# In[ ]:




