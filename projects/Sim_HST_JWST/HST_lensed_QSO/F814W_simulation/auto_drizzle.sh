#!/bin/bash
#Run using iraf27 environment
# read -p 'Simulation ID:' id
#for id in 701 702 703 704 705 706 707 708 709 710 711 712 713 714 715 716 717 718 719 720 721 722 723 724 725 726 727 728 729 730 731 732 733 734 735 736 737 738 739 740 741 742 743 744 745 746 747 748 749 750
for id in {702..799}
do
cp -r ../../drizzle_temp/* sim_lens_ID_${id}
cd sim_lens_ID_${id}
python header.py
echo "cl<dri_img.cl; logout" | cl
echo "cl<dri_psf.cl; logout" | cl
echo "cl<dri_rmsSQ.cl; logout" | cl
python nordri.py
rm -r uparm
rm *.cl *.py
rm *w.fits
rm dri*.fits
cd ../
done