#!/bin/bash
#Run using iraf27 environment
#read -p 'Drizzling filter:' filt
#read -p 'Simulation ID:' id
for filt in F356W F200W
do
for seed in {301..350}
do
for id in 3
do
# echo sim_ID${id}_${filt}_seed${seed}
cp -r drizzle_temp_${filt}/* sim_ID${id}_${filt}_seed${seed}
cd sim_ID${id}_${filt}_seed${seed}
python header.py
echo "cl<dri_img.cl; logout" | cl
echo "cl<dri_AGNclean.cl; logout" | cl
echo "cl<dri_HOSTclean.cl; logout" | cl
echo "cl<dri_POINTclean.cl; logout" | cl
echo "cl<dri_rmsSQ.cl; logout" | cl
python nordri.py
rm -r uparm
rm *.cl *.py
rm *w.fits
rm dr*.fits
cd ../
done
done
done