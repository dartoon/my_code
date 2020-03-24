#!/bin/bash
#Run using iraf27 environment
read -p 'Drizzling filter:' filt
read -p 'Simulation ID:' id
cp -r ${filt}_drizzle_temp/* sim_ID_${id}
cd sim_${filt}_ID_${id}
python header.py
echo "cl<dri_img.cl; logout" | cl
echo "cl<dri_AGNclean.cl; logout" | cl
echo "cl<dri_HOSTclean.cl; logout" | cl
echo "cl<dri_POINTclean.cl; logout" | cl
python nordri.py
rm -r uparm
rm *.cl *.py
rm *w.fits
rm dr*.fits