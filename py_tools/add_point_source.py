#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 21:20:18 2017

@author: dxh
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from scipy import ndimage
def exp_grid(img,nums,drc):
    '''
    expand the image frame with zero, with nums, direct to expand
        img: 
            2d array_like
        num:
            number of pixels to input
        drc:
            direction. 1: -x 2: -y 3: x 4: -y (i.e.: from -x, anticlockwise)
    '''
    if drc==1:
        exp_img=np.concatenate((np.zeros([len(img), nums]),img), axis=1)
    if drc==2:
        exp_img=np.concatenate((np.zeros([nums, len(img.T)]), img), axis=0)
    if drc==3:
        exp_img=np.concatenate((img,np.zeros([len(img), nums])), axis=1)
    if drc==4:
        exp_img=np.concatenate((img,np.zeros([nums, len(img.T)])), axis=0)
    return exp_img

class Add_psf:
    """
    Interpolate, Shift and resize the PSF, to add as point source
    """
    def __init__(self,img=None,order=1,exp_numPix=601, numPix=241):  #Matt said the order is 1
        self.exp_numPix=exp_numPix
        self.numPix=numPix
        self.img = img
        self.dy,self.dx = img.shape
        self.x0 = self.dx/2
        self.y0 = self.dy/2
        self.order = order
        self.nx, self.ny = None, None   #For coordinate the pixel position
        if  self.order==1:
            self.model = self.img
        else:
            self.model = ndimage.spline_filter(self.img,output=np.float64,order=order)

    def getPix(self,sx,sy):
        x = (sx-self.nx)+self.x0
        y = (sy-self.ny)+self.y0
        return x,y

    def evaluateSource(self,sx,sy):
        ix,iy = self.getPix(sx,sy)
        return ndimage.map_coordinates(self.model,[[iy],[ix]],order=self.order)[0,:,:]
  
    def shift_psf(self, shift_pos):
        '''
        Shift the PSF 
        psf_info=Class_shifted_psf(img=psf,order=1)
        psf_shift=psf_info.shift_psf(shift_pos=[0.56,0.64])
        '''
        y,x = np.indices((int(len(self.img)),int(len(self.img)))).astype(np.float64)
        self.nx=len(self.img)/2+shift_pos[0]  # 0.x if the residul is *.x
        self.ny=len(self.img)/2+shift_pos[1]  # 0.y if the residul is *.y
#        print "shift_pos:", shift_pos
        return self.evaluateSource(x,y)
  
    def exp_psf(self, int_pos, shift_pos):
        '''
        Expand the grid for the Shiftted the PSF 
        '''
        img_size=len(self.img)
        exp_img=self.shift_psf(shift_pos)
        exp_img=exp_grid(exp_img,int_pos[0]+(self.exp_numPix-img_size)/2-self.numPix/2,1)   #to make start from zero
        exp_img=exp_grid(exp_img,int_pos[1]+(self.exp_numPix-img_size)/2-self.numPix/2,2)
        exp_img=exp_grid(exp_img,self.exp_numPix-len(exp_img),3)
        exp_img=exp_grid(exp_img,self.exp_numPix-len(exp_img),4)
        return exp_img

    def add_psf(self, pos):
        '''
        Directly produce a interpolated PSF to add as Point source
        '''
        int_pos=pos.astype(int)
        shift_pos=pos-int_pos
        if shift_pos.min()>=0:
            exp=self.exp_psf(int_pos=int_pos,shift_pos=shift_pos) 
            cnt1=(self.exp_numPix-self.numPix)/2
            cnt2=(self.exp_numPix-self.numPix)/2+self.numPix
            trim_exp=exp[cnt1:cnt2,cnt1:cnt2]
            trim_exp /= trim_exp.sum() 
            return trim_exp
        else:
            raise ValueError('Have negative position value when interpolated PSF as Point source')

##==============================================================================
## example of shift a position at (x,y)=(120.56, 120.64)
##==============================================================================
#psf = pyfits.open('../material/PSF/f160w/psf_sub4.fits')[0].data.copy()
#psf=psf[1:,1:]
#psf /= psf.sum()   
#psf_info=Add_psf(img=psf,order=1,exp_numPix= 600, numPix=241)
#trim_exp1=psf_info.add_psf(pos=np.array([120.56, 120]))
#trim_exp2=psf_info.add_psf(pos=np.array([120.56, 120.64]))
#plt.imshow(np.log10(trim_exp2-trim_exp1),origin='lower')
#plt.colorbar()
#plt.show()


