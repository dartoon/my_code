#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 10:13:22 2020

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

taka_fil = np.loadtxt('204.fil')
sov_fil = np.loadtxt('HST_WFC3_IR.F140W.dat')

plt.plot(taka_fil[:,1], taka_fil[:,2])
plt.plot(sov_fil[:,0], sov_fil[:,1])

sov_fil = np.loadtxt('HST_WFC3_IR.F140W.dat')

print(taka_fil[:,2].max())
print(sov_fil[:1].max())
#%%
#Produce a fil as Taka's template
sov_jwst_f444w_fil = np.loadtxt('../esti_filter_magnitude/JWST_NIRCam.F444W.dat')
sov_jwst_f444w_fil[:,1] /= sov_jwst_f444w_fil[:,1].max()
plt.plot(sov_jwst_f444w_fil[:,0], sov_jwst_f444w_fil[:,1])
file_text = []
for i in range(len(sov_jwst_f444w_fil)):
    file_text.append([i+1, sov_jwst_f444w_fil[i,0], sov_jwst_f444w_fil[i,1]])
file_text = np.asarray(file_text)
#Write information for good team
write_file = open('f444w.fil','w') 
write_file.write("#File normed by Xuheng from SOV http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=JWST/NIRCam.F444W&&mode=browse&gname=JWST&gname2=NIRCam#filter\n")
				 
count = 1				 
for i in range(len(file_text)):
#	if float(i/3) == int(i/3):
#		print(i)
	write_file.write("{0} {1} {2}\n".format(count, float(file_text[i][1]), float(file_text[i][2])))
	count = count+1
write_file.close()

#%%
#Produce a fil as Taka's template
sov_jwst_f356w_fil = np.loadtxt('../esti_filter_magnitude/JWST_NIRCam.F356W.dat')
sov_jwst_f356w_fil[:,1] /= sov_jwst_f356w_fil[:,1].max()
plt.plot(sov_jwst_f356w_fil[:,0], sov_jwst_f356w_fil[:,1])
file_text = []
for i in range(len(sov_jwst_f356w_fil)):
    file_text.append([i+1, sov_jwst_f356w_fil[i,0], sov_jwst_f356w_fil[i,1]])
file_text = np.asarray(file_text)
#Write information for good team
write_file = open('f356w.fil','w') 
write_file.write("#File normed by Xuheng from SOV http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=JWST/NIRCam.F356W&&mode=browse&gname=JWST&gname2=NIRCam#filter\n")
				 
count = 1				 
for i in range(len(file_text)):
#	if float(i/3) == int(i/3):
#		print(i)
	write_file.write("{0} {1} {2}\n".format(count, float(file_text[i][1]), float(file_text[i][2])))
	count = count+1
write_file.close()

#%%
#Produce a fil as Taka's template
sov_jwst_f200w_fil = np.loadtxt('../esti_filter_magnitude/JWST_NIRCam.F200W.dat')
sov_jwst_f200w_fil[:,1] /= sov_jwst_f200w_fil[:,1].max()
plt.plot(sov_jwst_f200w_fil[:,0], sov_jwst_f200w_fil[:,1])
file_text = []
for i in range(len(sov_jwst_f200w_fil)):
    file_text.append([i+1, sov_jwst_f200w_fil[i,0], sov_jwst_f200w_fil[i,1]])
file_text = np.asarray(file_text)
#Write information for good team
write_file = open('f200w.fil','w') 
write_file.write("#File normed by Xuheng from SOV http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=JWST/NIRCam.F200W&&mode=browse&gname=JWST&gname2=NIRCam#filter\n")
				 
count = 1				 
for i in range(len(file_text)):
#	if float(i/3) == int(i/3):
#		print(i)
	write_file.write("{0} {1} {2}\n".format(count, float(file_text[i][1]), float(file_text[i][2])))
	count = count+1
write_file.close()

#%%
#Produce a fil as Taka's template
sov_jwst_f150w_fil = np.loadtxt('../esti_filter_magnitude/JWST_NIRCam.F150W.dat')
sov_jwst_f150w_fil[:,1] /= sov_jwst_f150w_fil[:,1].max()
plt.plot(sov_jwst_f150w_fil[:,0], sov_jwst_f150w_fil[:,1])
file_text = []
for i in range(len(sov_jwst_f150w_fil)):
    file_text.append([i+1, sov_jwst_f150w_fil[i,0], sov_jwst_f150w_fil[i,1]])
file_text = np.asarray(file_text)
#Write information for good team
write_file = open('f150w.fil','w') 
write_file.write("#File normed by Xuheng from SOV http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=JWST/NIRCam.F150W&&mode=browse&gname=JWST&gname2=NIRCam#filter\n")
				 
count = 1				 
for i in range(len(file_text)):
#	if float(i/3) == int(i/3):
#		print(i)
	write_file.write("{0} {1} {2}\n".format(count, float(file_text[i][1]), float(file_text[i][2])))
	count = count+1
write_file.close()