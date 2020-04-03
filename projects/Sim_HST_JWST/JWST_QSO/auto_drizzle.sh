#!/bin/bash
#Run using iraf27 environment
#read -p 'Drizzling filter:' filt
#read -p 'Simulation ID:' id
for filt in F444W F356W F200W F150W
do
for seed in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
do
for id in 1 2 3 4 5
do
# echo sim_ID${id}_${filt}_seed${seed}
cp -r drizzle_temp_${filt}/* sim_ID${id}_${filt}_seed${seed}
cd sim_ID${id}_${filt}_seed${seed}
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
cd ../
done
done
done