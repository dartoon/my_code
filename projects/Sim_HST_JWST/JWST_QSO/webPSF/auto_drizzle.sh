#!/bin/bash
#Run using iraf27 environment
read -p 'Simulation ID:' filter
cp -r drizzle_${filter}_temp/* drizzle_PSF_${filter}
cd drizzle_PSF_${filter}
python header.py
#echo "cl<dri_psf.cl; logout" | cl
python nordri.py
rm -r uparm
rm *.cl *.py
rm *w.fits
rm dr*.fits