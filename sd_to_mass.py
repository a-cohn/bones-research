#import packages
import pyfits as pyfits
import numpy as np
import matplotlib.pyplot as plt

#import distances to filaments
distances=np.array([4.6,3.82,3.24,4.26,3.60,3.34,3.1,3.0,3.2,3.3])

#for filaments
for n in range(1,11):
    if n>11: 
        pass
    else:
        #import image data
        f=pyfits.open('/Users/Cohn/Desktop/summer15research/masses with corrected foreground intensity ratios/filament_'+str(n)+'_mass_above_5_background_std_above_mean_map.fits')
        data=f[0].data
        header=f[0].header

        #measuring the size of the file
        nrows=len(data[:,0]) #number of rows of fits file
        ncols=len(data[0]) #number of columns of fits file

        #create map
        new_map=np.zeros((nrows,ncols))

        #initialize variables
        mass=0
        pix=0
        sum_sd=0

        #generate surface density map
        for row in range(nrows):
            for col in range(ncols):
                if data[row,col]>0:
                    new_map[row,col]=data[row,col]*(distances[n-1]*3.0857*10**21*7.8/206265)**2/(1.989*10**33)
                    mass += new_map[row,col]
                    sum_sd += data[row,col]
                    pix += 1
        pyfits.writeto('/Users/Cohn/Desktop/summer15research/masses with corrected foreground intensity ratios/filament_'+str(n)+'_mass_above_5_background_std_above_mean_mass_map.fits',new_map,header)
        if pix==0:
            avg_sd=0
            avg_mass=0
        else:
            avg_sd=sum_sd/pix
            avg_mass=mass/pix
        differences2=0
        differences2_mass=0
        for row in range(nrows):
            for col in range(ncols):
                if data[row,col]>0:
                    differences2 += (data[row,col]-avg_sd)**2
                    differences2_mass += (data[row,col]*(distances[n-1]*3.0857*10**21*7.8/206265)**2/(1.989*10**33)-avg_mass)**2
        std=(differences2/(pix-1))**0.5
        std_mass=(differences2_mass/(pix-1))**0.5
        sstd=pix**0.5*std
        sstd_mass=pix**0.5*std_mass
        print "filament: ", n, "mass above 5 background std above mean: ", mass , "standard deviation: ", sstd_mass

