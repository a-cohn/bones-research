#amy cohn
#program for cropping
#created: june 11, 2015

#import packages
import pyfits as pyfits
import numpy as np
import matplotlib.pyplot as plt

#open fits file
f=pyfits.open('mosaicv3_ch4_333.fits')
fim=pyfits.open('imagemask_filament_10.fits')
dt=f[0].data
imdt=fim[0].data
data=np.asarray(dt)
imdata=np.asarray(imdt)

#define region of interest (ellipse)
h=6726.3267 
k=3156.819
a=1727.069
b=418.38837

#measuring the size of the file
nrows=len(data[:,0]) #number of rows of fits file
ncols=len(data[0]) #number of columns of fits file


#initializing cropped maps
cropbox=np.zeros((3*b,3*a)) #map for cropped data
imcropbox=np.zeros((3*b,3*a)) #map for cropped imagemask


#cropping
for row in range(nrows):
    print "row", row, "nrows", nrows
    for col in range(ncols):
        if k-1.5*b<row<k+1.5*b:
            if h-1.5*a<col<h+1.5*a:
                cropbox[row-k+1.5*b,col-h+1.5*a]=data[row,col]
                imcropbox[row-k+1.5*b,col-h+1.5*a]=imdata[row,col]

#saving files                
hdu = pyfits.PrimaryHDU(cropbox)
imhdu = pyfits.PrimaryHDU(imcropbox)
hdulist = pyfits.HDUList([hdu])
imhdulist = pyfits.HDUList([imhdu])
hdulist.writeto('/Users/Cohn/Desktop/summer15research/filament_10_cropped.fits')
imhdulist.writeto('/Users/Cohn/Desktop/summer15research/imagemask_filament_10_cropped.fits')

