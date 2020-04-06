#!/bin/bash
#Run using iraf27 environment
#read -p 'Simulation ID (noqso):' id
for id in 605 606 607 608 609 610 611 612 613 614 615 616 617 618 619 620 621
do
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
cd ../
done