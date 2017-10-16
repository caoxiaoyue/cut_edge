#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 14:55:24 2017
take an aggressive cutting-method
all the object belong to edge-type are cutted
(even the non-edge-type object which have same id with edge-object
are cutted) 

@author: cao
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np



'''
#test row/col in fits image
fname = '/home/cao/Desktop/easy_test/test_B.fits'#'../outputs_mem_M/0_lensed_M_fits.gz'
fitsfile = fits.open(fname)
imdata = fitsfile[0].data  #posx/posy in imdata is switched
imdata[...] = 0
imdata[20:30,2:7] = 1
plt.imshow(imdata)  #origin -- top-left   left->right top->down
row,col = imdata.shape
fits.writeto('/home/cao/Desktop/easy_test/11111.fits',imdata)
#origin -- down-left   left->right  down->top
'''


fname = '../outputs_mem_M/0_lensed_M_fits.gz'
fitsfile = fits.open(fname)
imdata = fitsfile[0].data  #posx/posy in imdata is switched
id1 = imdata.reshape(imdata.size)
row,col = imdata.shape


top_edge = imdata[0,:]
down_edge = imdata[-1,:] 
left_edge = imdata[:,0]
right_edge = imdata[:,-1]

search_top = np.unique(top_edge)
search_down = np.unique(down_edge)
search_left = np.unique(left_edge)
search_right = np.unique(right_edge)

if search_top.size > 1 :
    ss= search_top[1:]
    for idtmp in ss:
        idx = (id1 == idtmp).nonzero()[0]
        if (idx.size != 0) :
            id1[idx] = 0
            
if search_down.size > 1 :
    ss= search_down[1:]
    for idtmp in ss:
        idx = (id1 == idtmp).nonzero()[0]
        if (idx.size != 0) :
            id1[idx] = 0
            
            
if search_left.size > 1 :
    ss= search_left[1:]
    for idtmp in ss:
        idx = (id1 == idtmp).nonzero()[0]
        if (idx.size != 0) :
            id1[idx] = 0
            
            
if search_right.size > 1 :
    ss= search_right[1:]
    for idtmp in ss:
        idx = (id1 == idtmp).nonzero()[0]
        if (idx.size != 0) :
            id1[idx] = 0
            
            
fitsfile.writeto('./111.fits',overwrite=True)
        
        
    