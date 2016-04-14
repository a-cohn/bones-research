#amy cohn
#program for masking an ellipse
#creates map of a FITS file where unmasked pixels are marked with 0's and masked pixels (inside of the ellipse of interest) are marked with 1's
#created: june 10, 2015

import pyfits as pyfits
import numpy as np
import matplotlib.pyplot as plt
f=pyfits.open('/Users/Cohn/Desktop/summer15research/gaussfit_iter_beta175_column_l026_itervar_conv36.fits')
data=f[0].data
header=f[0].header
header.remove('Equinox')

#define region of interest (ellipse)
h=1129.2329
k=461.27555
a=12.554057
b=a
tilt=0*np.pi/180


#measuring the size of the file
nrows=len(data[:,0]) #number of rows of fits file
ncols=len(data[0]) #number of columns of fits file

#create map of file
map=np.ones((nrows,ncols))

#picking pixels
for row in range(nrows):
    print row, nrows
    for col in range(ncols):
        if k-(a**2+b**2)**0.5<row<k+(a**2+b**2)**0.5:
            if h-(a**2+b**2)**0.5<col<h+(a**2+b**2)**0.5:
                if (((col-h)*np.cos(tilt)+(row-k)*np.sin(tilt))**2/a**2)+(((col-h)*np.sin(tilt)-(row-k)*np.cos(tilt))**2/b**2)<1:
                    map[row,col]=0

pyfits.writeto('/Users/Cohn/Desktop/summer15research/imagemask_filament_1_herschel_off_region.fits',map,header)


