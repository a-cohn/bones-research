#import packages
import pyfits as pyfits
import numpy as np
import matplotlib.pyplot as plt

#load distances
distances=np.array([4.6,3.82,3.24,4.26,3.60,3.34,3.1,3.0,3.2,3.3])

#for filaments
for n in range(1,11):
    #load data
    f=pyfits.open('/Users/Cohn/Desktop/summer15research/masses with corrected foreground intensity ratios/filament_'+str(n)+'_surface_density_map_with_header_new.fits')
    data=f[0].data
    sd=f[0].data
    header=f[0].header
    f1=pyfits.open('/Users/Cohn/Desktop/summer15research/imagemask_filament_'+str(n)+'_background.fits')
    imdata=f1[0].data
    f2=pyfits.open('/Users/Cohn/Desktop/summer15research/imagemask_filament_'+str(n)+'_cropped_with_header.fits')
    regdata=f2[0].data
    
    #measuring the size of the file
    nrows=len(data[:,0]) #number of rows of fits file
    ncols=len(data[0]) #number of columns of fits file

    pix=0
    sum_sd=0
    for row in range(nrows):
        for col in range(ncols):
            if imdata[row,col]==0:
                if sd[row,col]>0:
                    pix += 1
                    sum_sd += sd[row,col]
    avg_sd=sum_sd/pix
    differences2=0
    for row in range(nrows):
        for col in range(ncols):
            if imdata[row,col]==0:
                if sd[row,col]>0:
                    differences2 += (sd[row,col]-avg_sd)**2
    std=(differences2/(pix-1))**0.5
    mass0=0
    mass0_map=np.zeros((nrows,ncols))
    mass0_mask=np.zeros((nrows,ncols))
    mass3=0
    mass3_map=np.zeros((nrows,ncols))
    mass3_mask=np.zeros((nrows,ncols))
    mass5=0
    mass5_map=np.zeros((nrows,ncols))
    mass5_mask=np.zeros((nrows,ncols))
    pix=0
    for row in range(nrows):
        for col in range(ncols):
            if regdata[row,col]==0:
                if sd[row,col]>0:
                    mass0 += sd[row,col]*(distances[n-1]*3.0857*10**21*1.2/206265)**2/(1.989*10**33)
                    mass0_map[row,col]=sd[row,col]*(distances[n-1]*3.0857*10**21*1.2/206265)**2/(1.989*10**33)
                    mass0_mask[row,col]=1
                    pix += 1
                if sd[row,col]>avg_sd+3*std:
                    mass3 += sd[row,col]*(distances[n-1]*3.0857*10**21*1.2/206265)**2/(1.989*10**33)
                    mass3_map[row,col]=sd[row,col]*(distances[n-1]*3.0857*10**21*1.2/206265)**2/(1.989*10**33)
                    mass3_mask[row,col]=1
                if sd[row,col]>avg_sd+5*std:
                    mass5 += sd[row,col]*(distances[n-1]*3.0857*10**21*1.2/206265)**2/(1.989*10**33)
                    mass5_map[row,col]=sd[row,col]*(distances[n-1]*3.0857*10**21*1.2/206265)**2/(1.989*10**33)
                    mass5_mask[row,col]=1
    avg_mass=mass0/pix
    differences2_mass=0
    for row in range(nrows):
        for col in range(ncols):
            if regdata[row,col]==0:
                if sd[row,col]>0:
                    differences2_mass += (sd[row,col]*(distances[n-1]*3.0857*10**21*1.2/206265)**2/(1.989*10**33)-avg_mass)**2
    std_mass=(differences2_mass/(pix-1))**0.5
    sstd_mass=std_mass*pix**0.5
    '''pyfits.writeto('masses with corrected foreground intensity ratios/filament_'+str(n)+'_mass_above_0_mask.fits',mass0_mask,header)
    pyfits.writeto('masses with corrected foreground intensity ratios/filament_'+str(n)+'_mass_above_0_map.fits',mass0_map,header)
    pyfits.writeto('masses with corrected foreground intensity ratios/filament_'+str(n)+'_mass_above_3_background_std_above_mean_map.fits',mass3_map,header)
    pyfits.writeto('masses with corrected foreground intensity ratios/filament_'+str(n)+'_mass_above_3_background_std_above_mean_mask.fits',mass3_mask,header)
    pyfits.writeto('masses with corrected foreground intensity ratios/filament_'+str(n)+'_mass_above_5_background_std_above_mean_map.fits',mass5_map,header)
    pyfits.writeto('masses with corrected foreground intensity ratios/filament_'+str(n)+'_mass_above_5_background_std_above_mean_mask.fits',mass5_mask,header)'''
    print "filament n: ", n, "mass above 0: ", mass0, "mass above 3 bacgkround std above mean: ", mass3, "mass above 5 background std above mean: ", mass5, "standard deviation: ", sstd_mass

