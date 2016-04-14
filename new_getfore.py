#import packages
import pyfits as pyfits
import numpy as np
import matplotlib.pyplot as plt

#open data file
f=pyfits.open('/Users/Cohn/Desktop/summer15research/mosaicv3_ch4_018_cropped.fits')
dt=f[0].data
data=np.asarray(dt)

#open image mask
im=pyfits.open('imagemask_mosaicv3_ch4_018_cropped.fits')
imdt=im[0].data
imdata=np.asarray(imdt)

#measuring the size of the file
nrows=len(data[:,0]) #number of rows of fits file
ncols=len(data[0]) #number of columns of fits file

minval=1000
total_flux=0
total_pixels=0
#find minimum value and mean flux
for row in range(nrows):
    #print "row: ", row, "nrows: ", nrows
    for col in range(ncols):
        if imdata[row,col]==0:
            data[row,col]=np.nan_to_num(data[row,col])
            total_flux += data[row,col]
            total_pixels += 1
            if data[row,col]<minval:
                minval=data[row,col]
mean_flux=total_flux/total_pixels
sum=0 #sum of square differences
#find sigma
for row in range(nrows):
    #print "row: ", row, "nrows: ", nrows
    for col in range(ncols):
        if imdata[row,col]==0:
            sum += (data[row,col]-mean_flux)**2
n=total_pixels
std=(sum/(n-1))**0.5
sat=np.zeros((nrows,ncols))
for row in range(nrows):
    #print "row: ", row, "nrows: ", nrows
    for col in range(ncols):
        value=0
        if imdata[row,col]==0:
            if minval<data[row,col]<minval+2*std:
                for row_new in range(row-8, row+8):
                    #print range(row-8, row+8)
                    for col_new in range(col-8, col+8):
                        #print range(col-8, col+8)
                        if data[row_new,col_new]==minval:
                            value=1
                if value==0:
                    sat[row,col]=data[row,col]
total_flux_sat=0
total_pix_sat=0
for row in range(nrows):
    #print "row: ", row, "nrows: ", nrows
    for col in range(ncols):
        if imdata[row,col]==0:
            if sat[row,col]>0:
                total_flux_sat += sat[row,col]
                total_pix_sat += 1
mean_sat=total_flux_sat/total_pix_sat
fi=mean_sat-2*std
print "foreground intensity: ", fi
                       
            


#define the region of interest

#find minimum value of foreground radiation inside this region

#search for all pixels within the region between the minimum value and the minimum value + two sigma

#select all pixels outside of this range as the saturated pixels

#find the mean value of the saturated pixels

#foreground intensity raito = mean value of saturated pixels - two sigma
