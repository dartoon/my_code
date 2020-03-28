#!/bin/bash
#Run using iraf27 environment
# read -p 'Simulation ID:' id
for id in 501 502 503 504 505 506 507 508 509 510 511 512 513 514 515 516 517 518 519 520 521
do
cp -r ../drizzle_temp/* sim_lens_ID_${id}
cd sim_lens_ID_${id}
python header.py
echo "cl<dri_img.cl; logout" | cl
echo "cl<dri_psf.cl; logout" | cl
python nordri.py
rm -r uparm
rm *.cl *.py
rm *w.fits
rm dri*.fits
cd ../
done