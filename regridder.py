#import packages
from scipy import interpolate
import pyfits as pyfits
import numpy as np
import matplotlib.pyplot as plt

#filament to be regridded
fil=9

#open fits file
image=pyfits.open('mosaicv3_ch4_336.fits')
imagemask=pyfits.open('imagemask_filament_'+str(fil)+'_with_header_new.fits')
data=image[0].data #image data
header=image[0].header #image header
im_data=imagemask[0].data #imagemask data
im_header=imagemask[0].header #imagemask header

#measuring the size of the file
nrows=len(data[:,0]) #number of rows of fits file
ncols=len(data[0]) #number of columns of fits file

n=5 #factor to be regridded by

#initialize maps for storing regridded images
reduced=np.zeros((nrows/n,ncols/n)) #map of data
reduced1=np.zeros((nrows/n,ncols/n)) #map of imagemask

#regridding
for i in range(nrows/n):
    for j in range(ncols/n):
        reduced[i,j] = data[n*i,n*j]
        reduced1[i,j] = im_data[n*i, n*j]
        
#saving files
pyfits.writeto('/Users/Cohn/Desktop/summer15research/filament_'+str(fil)+'_with_header_regridded_by5.fits', reduced,header) #regridded data
pyfits.writeto('/Users/Cohn/Desktop/summer15research/imagemask_filament_'+str(fil)+'_with_header_regridded_by5.fits', reduced1, im_header) #regridded imagemask
                                
