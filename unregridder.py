#import packages
import numpy as np
import scipy as sp
import scipy.interpolate as interp
import pyfits as pyfits

#filament
n=5

#open fits file
f=pyfits.open('/Users/Cohn/Desktop/summer15research/filament_'+str(n)+'_cropped_with_header_regridded_by5_50pxrmedian.fits')
data=f[0].data
header=f[0].header

#measuring the size of the file
nrows=len(data[:,0]) #number of rows of fits file
ncols=len(data[0]) #number of columns of fits file

#unregridding
factor=5
rows = np.linspace(0,nrows,nrows)
cols = np.linspace(0,ncols,ncols)
myspline = interp.RectBivariateSpline(rows,cols,data)
new_rows = np.linspace(0,nrows,nrows*factor)
new_cols = np.linspace(0,ncols,ncols*factor)
splineoutput = myspline(new_rows,new_cols)

#save file
pyfits.writeto('/Users/Cohn/Desktop/summer15research/filament_'+str(n)+'_cropped_with_header_50pxrmedian_unregridded.fits',splineoutput,header)
                      
