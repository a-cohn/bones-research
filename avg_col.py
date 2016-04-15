#import packages
import pyfits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import astropy as astropy
from astropy.wcs import WCS
#data
col_name_array=np.array(['gaussfit_iter_beta175_column_l026_itervar_conv36.fits','gaussfit_iter_beta175_column_l026_itervar_conv36.fits','gaussfit_iter_beta175_column_l024_itervar_conv36.fits','gaussfit_iter_beta175_column_l022_itervar_conv36.fits','gaussfit_iter_beta175_column_inner40_conv_mosaics_itervar_conv36.fits','gaussfit_iter_beta175_column_inner40_conv_mosaics_itervar_conv36.fits','gaussfit_iter_beta175_column_inner40_conv_mosaics_itervar_conv36.fits','gaussfit_iter_beta175_column_inner40_conv_mosaics_itervar_conv36.fits','gaussfit_iter_beta175_column_l334_itervar_conv36.fits','gaussfit_iter_beta175_column_l332_itervar_conv36.fits'])
label_name_array=np.array(['gaussfit_iter_beta175_column_l026_itervar_conv36_label.fits','gaussfit_iter_beta175_column_l026_itervar_conv36_label.fits','gaussfit_iter_beta175_column_l024_itervar_conv36_label.fits','gaussfit_iter_beta175_column_l022_itervar_conv36_label.fits','gaussfit_iter_beta175_column_inner40_conv_mosaics_itervar_conv36_label.fits','gaussfit_iter_beta175_column_inner40_conv_mosaics_itervar_conv36_label.fits','gaussfit_iter_beta175_column_inner40_conv_mosaics_itervar_conv36_label.fits','gaussfit_iter_beta175_column_inner40_conv_mosaics_itervar_conv36_label.fits','gaussfit_iter_beta175_column_l334_itervar_conv36_label.fits','gaussfit_iter_beta175_column_l332_itervar_conv36_label.fits'])
for n in range(1,11):
    if n>11:
        pass
    else:
        #open data file
        f=pyfits.open('/Users/Cohn/Desktop/summer15research/'+col_name_array[n-1])
        data=f[0].data
        header=f[0].header
        if 4<n<9:
            pass
        else:
            header.remove('Equinox')

        #open cloud image mask
        im=pyfits.open('/Users/Cohn/Desktop/summer15research/masses with corrected foreground intensity ratios/filament_'+str(n)+'_mass_above_3_background_std_above_mean_mask.fits')
        imdata=im[0].data

        #open label mask
        lb=pyfits.open('/Users/Cohn/Desktop/summer15research/'+label_name_array[n-1])
        lbdata=lb[0].data

        #measuring the size of the file
        nrows_im=len(imdata[:,0]) #number of rows of fits file
        ncols_im=len(imdata[0]) #number of columns of fits file

        nrows_dat=len(data[:,0])
        ncols_dat=len(data[0])

        w = WCS('/Users/Cohn/Desktop/summer15research/'+col_name_array[n-1])
        wim = WCS('/Users/Cohn/Desktop/summer15research/masses with corrected foreground intensity ratios/filament_'+str(n)+'_mass_above_3_background_std_above_mean_mask.fits')
        
        lon_array=np.zeros((nrows_im,ncols_im))
        lat_array=np.zeros((nrows_im,ncols_im))
        for row in range(nrows_im):
            for col in range(ncols_im):
                lon, lat = wim.all_pix2world(col,row,0)
                lon_array[row,col]=lon
                lat_array[row,col]=lat

        x_array=np.zeros((nrows_im,ncols_im))
        y_array=np.zeros((nrows_im,ncols_im))
        for row in range(nrows_im):
            for col in range(ncols_im):
                if imdata[row,col]==1:
                    x, y = w.all_world2pix(lon_array[row,col],lat_array[row,col],0)
                    x_array[row,col] = x
                    y_array[row,col] = y

        data=lbdata*data
        new_map=np.zeros((nrows_dat,ncols_dat))
        total_col=0
        total_pix=0
        for row in range(nrows_im):
            for col in range(ncols_im):
                x=x_array[row,col]
                y=y_array[row,col]
                data[y,x]=np.nan_to_num(data[y,x])
                if data[y,x]>0:
                    new_map[y,x]=data[y,x]
                    total_col += data[y,x]
                    total_pix += 1
        avg_col=total_col/total_pix
        print "filament: ", n, "avg_col: ", avg_col
        difference2=0
        for row in range(nrows_im):
            for col in range(ncols_im):
                x=x_array[row,col]
                y=y_array[row,col]
                data[y,x]=np.nan_to_num(data[y,x])
                if data[y,x]>0:
                    difference2 += (data[y,x]-avg_col)**2
        std=(difference2/(total_pix-1))**0.5
        print "filament: ", n, "standard deviation: ", std
        #pyfits.writeto('/Users/Cohn/Desktop/summer15research/masses with corrected foreground intensity ratios/filament_'+str(n)+'_temp_map_from_5_background_std_masks.fits',new_map,header)
