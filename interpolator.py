#amy cohn
#program for interpolation
#created: june 11, 2015
#updated: june 18, 2015

import pyfits as pyfits
import numpy as np
import matplotlib.pyplot as plt

#open fits file
f=pyfits.open('/Users/Cohn/Desktop/summer15research/mosaicv3_ch4_018_cropped_50pxrmedian_regridded_by10.fits')
dt=f[0].data
data=np.asarray(dt)

#measuring the size of the file
nrows=len(data[:,0]) #number of rows of fits file
ncols=len(data[0]) #number of columns of fits file

#open image mask
im=pyfits.open('imagemask_cropped_regridded_by10.fits')
imdt=im[0].data
imdata=np.asarray(imdt)

#defining variables
doe=np.zeros((nrows,ncols)) #data outside ellipse
weightab=np.zeros((nrows,ncols)) #table of weights
interpo=np.zeros((nrows,ncols)) #interpolation
distances=np.zeros((nrows,ncols)) #table of distances
newdistances=np.zeros((nrows,ncols)) #shifted table of distances
rowarray=np.zeros((nrows,ncols)) #row array
colarray=np.zeros((nrows,ncols)) #col array

#collecting data outside ellipse
for row in range(nrows):
    for col in range(ncols):
        data[row,col]=np.nan_to_num(data[row,col])
        if imdata[row,col]==1:
            doe[row,col]=data[row,col]
            interpo[row,col]=doe[row,col]

#creating row and column arrays
for row in range(nrows):
    for col in range(ncols):
        rowarray[row,col]=row
        colarray[row,col]=col
        
#interpolating on ellipse
for row in range(nrows):
    loop=0
    distances=np.zeros((nrows,ncols))
    newdistances=np.zeros((nrows,ncols))
    weightab=np.zeros((nrows,ncols))
    print row, nrows
    for col in range(ncols):
        if imdata[row,col]>0:
          pass
        if imdata[row,col]==0:
            if loop==0:
                for i in range(nrows):
                    for j in range(ncols):
                        if imdata[i,j]==1:
                            distances[i,j]=((i-row)**2+(j-col)**2)**0.5
                            distances[i,j]=np.nan_to_num(distances[i,j])
                            weightab[i,j]=1/(distances[i,j]**2*10**5)
                            weightab[i,j]=np.nan_to_num(weightab[i,j])
                        if imdata[i,j]==0:
                            pass
                for i in range(nrows):
                  for j in range(ncols):
                    interpo[row,col] += doe[i,j]*weightab[i,j]/np.sum(weightab)
                interpo[row,col]=np.nan_to_num(interpo[row,col])
                loop += 1
            if loop>0:
                weightab=np.zeros((nrows,ncols))
                newdistaneces=np.zeros((nrows,ncols))
                for j in range(ncols):
                      if imdata[i,j]==1:
                         if j<col:
                                newdistances[:,j]=distances[:,j]+loop
                                newdistances[:,j]=np.nan_to_num(newdistances[:,j])
                                for i in range(nrows):
                                    weightab[i,j]=1/(newdistances[i,j]**2*10**5)
                                    weightab[i,j]=np.nan_to_num(weightab[i,j])

                         if j>col:
                               newdistances[:,j]=distances[:,j]-loop
                               newdistances[:,j]=np.nan_to_num(newdistances[:,j])
                               for i in range(nrows):
                                  weightab[i,j]=1/(newdistances[i,j]**2*10**5)
                                  weightab[i,j]=np.nan_to_num(weightab[i,j])
                      if imdata[i,j]==0:
                        pass
                
                for i in range(nrows):
                  for j in range(ncols):
                    interpo[row,col] += doe[i,j]*weightab[i,j]/np.sum(weightab)
                interpo[row,col]=np.nan_to_num(interpo[row,col])
                if interpo[row,col]>10**5:
                      interpo[row,col]=np.average(doe)
                loop += 1
hdu=pyfits.PrimaryHDU(interpo)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('/Users/Cohn/Desktop/summer15research/interpolation_regridded_by10_with_new_mean_routine_and_weightab_fix_and_right_edge_fix.fits')
                                
                           
