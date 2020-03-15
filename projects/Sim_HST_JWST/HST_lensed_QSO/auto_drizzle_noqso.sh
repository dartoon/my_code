#!/bin/bash
#Run using iraf27 environment
read -p 'Simulation ID (noqso):' id
cp -r ../drizzle_temp/* sim_lens_noqso_ID_${id}
cd sim_lens_noqso_ID_${id}
python header.py
echo "cl<dri_img.cl; logout" | cl
echo "cl<dri_psf.cl; logout" | cl
python nordri.py
rm -r uparm
rm *.cl *.py
rm *w.fits
rm dri*.fits