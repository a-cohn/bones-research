#import packages
import pyfits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy

#for each filament
for n in range(1,11):
    if 0<n<5 or 8<n<11: #select which filaments 
        pass
    else:
        #open data file
        f=pyfits.open('/Users/Cohn/Desktop/summer15research/filament_'+str(n)+'_surface_density_map_with_header.fits')
        data=f[0].data
        header=f[0].header

        #open image mask
        im=pyfits.open('imagemask_filament_'+str(n)+'_cropped_with_header.fits')
        imdt=im[0].data
        imdata=np.asarray(imdt)

        #open median filtered file
        mf=pyfits.open('/Users/Cohn/Desktop/summer15research/filament_'+str(n)+'_surface_density_map_with_header_10pxrmedian.fits')
        mfdt=mf[0].data
        mfdata=np.asarray(mfdt)

        #measuring the size of the file
        nrows=len(data[:,0]) #number of rows of fits file
        ncols=len(data[0]) #number of columns of fits file

        #measuring average flux
        total_flux=0
        total_pixels=0
        for row in range(nrows):
            for col in range(ncols):
                data[row,col]=np.nan_to_num(data[row,col])
                if imdata[row,col]==0:
                    total_flux += data[row,col]
                    total_pixels += 1
        mean=total_flux/total_pixels

        #mesausring standard deviation
        sum2=0
        for row in range(nrows):
            for col in range(ncols):
                sum2 += (data[row,col]-mean)**2
        std=(sum2/(total_pixels-1))**0.5

        #generating new maps
        new_map=np.zeros((nrows,ncols))
        cloud_map=np.zeros((nrows,ncols))
        new_map2=np.zeros((nrows,ncols))
        cloud_map2=np.zeros((nrows,ncols))
        new_map0=np.zeros((nrows,ncols))
        cloud_map0=np.zeros((nrows,ncols))
        for row in range(nrows):
            for col in range(ncols):
                if imdata[row,col]==0:
                    if data[row,col]>2*std+mean:
                        new_map2[row,col]=data[row,col]
                        cloud_map2[row,col]=1
                    if data[row,col]>std+mean:
                        new_map[row,col]=data[row,col]
                        cloud_map[row,col]=1
                    if data[row,col]>0:
                        new_map0[row,col]=data[row,col]
                        cloud_map0[row,col]=1
        #saving new maps
        pyfits.writeto('/Users/Cohn/Desktop/summer15research/filament_'+str(n)+'_cloud_mask_above_1_sigma_above_mean_with_header.fits',cloud_map,header)
        pyfits.writeto('/Users/Cohn/Desktop/summer15research/filament_'+str(n)+'_cloud_mask_above_2_sigma_above_mean_with_header.fits',cloud_map2,header)
        pyfits.writeto('/Users/Cohn/Desktop/summer15research/filament_'+str(n)+'_surface_density_above_1_sigma_above_mean_with_header.fits',new_map,header)
        pyfits.writeto('/Users/Cohn/Desktop/summer15research/filament_'+str(n)+'_surface_density_above_2_sigma_above_mean_with_header.fits',new_map2,header)

        ''' #surface density path along spine (does not work well)
        path=np.zeros((nrows,ncols))
        for col in range(ncols):
            sum=0
            maxval=0
            maxval_row=0
            avg=0
            print col, ncols
            for row in range(nrows):
                if new_map[row,col]>maxval:
                    maxval=new_map[row,col]
                    maxval_row=row
            avg=mfdata[maxval_row,col]
            path[maxval_row,col]=avg
        #pyfits.writeto('/Users/Cohn/Desktop/summer15research/filament_'+str(n)+'_surface_density_path_along_spine_with_averaging_with_header.fits',path,header)

        x=np.zeros((ncols))
        y=np.zeros((ncols))
        weights=np.zeros((ncols))
        for col in range(ncols):
            x[col]=col
            for row in range(nrows):
                if path[row,col]>0:
                    y[col]=path[row,col]
                    weights[col]=1
        plt.scatter(x,y)
        #plt.savefig('filament_'+str(n)+'_mass_per_unit_length_scatter_plot_path_along_spine_with_averaging.png')'''
