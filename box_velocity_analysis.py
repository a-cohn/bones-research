#import packages
import pyfits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import math as math
import astropy.wcs
from astropy.wcs import WCS
from scipy.optimize import curve_fit
from statistics import *
from pylab import *


distances=np.array([4.6,3.82,3.24,4.26,3.60,3.34,3.1,3.0,3.2,3.3])
vel_array=(["/Users/Cohn/Desktop/summer15research/GRS files/filament 1/GRS_40to100kms.fits", "/Users/Cohn/Desktop/summer15research/GRS files/filament 2/GRS_30to90kms (1).fits","/Users/Cohn/Desktop/summer15research/GRS files/filament 3/GRS_20to80kms.fits", "/Users/Cohn/Desktop/summer15research/GRS files/filament 4/GRS_40to100kms.fits","/Users/Cohn/Desktop/summer15research/GRS files/filament 5/grs_30_60kms.fits",0,0,"/Users/Cohn/Desktop/summer15research/GRS files/filament 8/thrumms_neg30to30kms.fits","/Users/Cohn/Desktop/summer15research/GRS files/filament 9/thrumms_neg80toneg20kms.fits","/Users/Cohn/Desktop/summer15research/GRS files/filament 10/thrumms_neg80toneg20kms.fits"])

for n in range(10):
    if n>11:
        pass
    elif vel_array[n]==0:
        pass
    else:
        center_l=np.array([26.94,25.24,24.95,21.25,18.88,11.13,4.14,357.62,335.31,332.21])
        center_x=np.array([796,940,929,983,1746,934,1239,1113,2857,5484])
        center_b=np.array([-0.30,-0.45,-0.17,-0.15,-0.09,-0.12,-0.02,-0.33,-0.29,-0.04])
        center_y=np.array([255,396,543,374,290,347,355,437,1409,1586])
        length_deg=np.array([0.16,0.63,0.18,0.20,0.70,0.38,0.69,0.40,0.60,0.90])
        length_pc=np.array([13,47,14,14,45,22,37,21,34,52])
        pc_per_deg=length_pc/length_deg
        arcsec_per_pix=1.2
        pc_per_pix=length_pc[n]*arcsec_per_pix/(length_deg[n]*60*60)
        num_boxes=length_pc[n]/5.0
        f=pyfits.open('/Users/Cohn/Desktop/summer15research/masses with corrected foreground intensity ratios/filament_'+str(n+1)+'_surface_density_map_with_header_new.fits')
        data=f[0].data
        header=f[0].header

        fm=pyfits.open('/Users/Cohn/Desktop/summer15research/masses with corrected foreground intensity ratios/filament_'+str(n+1)+'_mass_above_5_background_std_above_mean_map.fits')
        mdata=fm[0].data

        nrows=len(data[:,0]) #number of rows of fits file
        ncols=len(data[0]) #number of columns of fits file

        im=pyfits.open('imagemask_filament_'+str(n+1)+'_cropped_with_header.fits')
        imdt=im[0].data
        imdata=np.asarray(imdt)
        if vel_array[n]==0:
            pass
        else:
            v=pyfits.open(vel_array[n])
            vdata=v[0].data
            vheader=v[0].header
            vcols=vheader["NAXIS1"]
            vrows=vheader["NAXIS2"]
            vdepth=vheader["NAXIS3"]
        #print(vcols,vrows,vdepth)
        #print(len(vdata[0,0,:]),len(vdata[0,:,0]),len(vdata[:,0,0]))
        vchannels=len(vdata[:,0,0])
        vrows=len(vdata[0,:,0])
        vcols=len(vdata[0,0,:])
        v0=vheader["CRVAl3"]
        vpix=vheader["CRPIX3"]
        vdelt=vheader["CDELT3"]
                
        w=WCS('/Users/Cohn/Desktop/summer15research/filament_'+str(n+1)+'_cropped_with_header.fits')

        new_map=np.zeros([nrows,ncols])
        for col in range(ncols):
            brightest_val=0
            brightest_val_row=0
            for row in range(nrows):
                if imdata[row,col]==0:
                    if mdata[row,col]>brightest_val:
                        brightest_val=mdata[row,col]
                        brightest_val_row=row
            new_map[brightest_val_row,col]=brightest_val
                    
        w=WCS('/Users/Cohn/Desktop/summer15research/filament_'+str(n+1)+'_cropped_with_header.fits')
        #pyfits.writeto('/Users/Cohn/Desktop/summer15research/filament_'+str(n+1)+'_brightest_val_point_per_column.fits',new_map,header)
        x=center_x[n]
        y=center_y[n]
        #generate boxes
        num_boxes=int(math.ceil(num_boxes))
        x_lims=np.zeros(num_boxes+1)
        l_lims=np.zeros(num_boxes+1)
        for i in range(num_boxes+1):
            x_lims[i]=center_x[n]+(i-num_boxes/2.0)*5.0/pc_per_pix
            l_lims[i]=center_l[n]+(i-num_boxes/2.0)*5.0/pc_per_deg[n]
        y_vals=np.zeros(num_boxes)
        b_vals=np.zeros(num_boxes)
        for j in range(num_boxes):
            y_sum=0
            y_n=1
            y_avg=0
            b_sum=0
            b_n=1
            sd_sum=0
            for col in range(ncols):
                if x_lims[j]<col<x_lims[j+1]:
                    for row in range(nrows):
                        lon,lat=w.all_pix2world(col,row,0)
                        if new_map[row,col]>0:
                            if row>0:
                                y_sum += row*new_map[row,col]
                                y_n += 1
                                b_sum += lat*new_map[row,col]
                                b_n += 1
                                if new_map[row,col]>0:
                                    sd_sum += new_map[row,col]
            if sd_sum>0:
                y_avg=y_sum/sd_sum
                b_avg=b_sum/sd_sum
            y_vals[j]=y_avg
            b_vals[j]=b_avg
        x_vals=np.zeros(num_boxes)
        l_vals=np.zeros(num_boxes)
        for k in range(num_boxes):
            x_vals[k]=(x_lims[k]+x_lims[k+1])/2
            l_vals[k]=(l_lims[k]+l_lims[k+1])/2
        print("filament: ", n+1)
        for l in range(num_boxes):
            print("box("+str(l_vals[l])+","+str(b_vals[num_boxes-1-l])+","+str(5.0/pc_per_deg[n])+","+str(5.0/pc_per_deg[n])+",0)") 
        print("----------")
        
        #bvlon, bvlat = w.all_pix2world(col, row, 0)
        #w=WCS(vel_array[n])
        #vy, vx, vd = w.all_world2pix(bvlon, bvlat, 0, 0)
        w=WCS(vel_array[n])
        x, y, d = w.all_world2pix(l_vals,b_vals,0,0)
        for i in range(len(x)):
            plt.clf()
            avg_intensity_array=np.zeros(vdepth)
            avg_intensity_array_2=np.zeros(vdepth)
            avg_intensity_array_3=np.zeros(vdepth)
            #print "filament: ", n, "box: ", len(x)-i, x[i]-(x[0]-x[1])/2, x[i]+(x[0]-x[1])/2
            col_counter=0
            row_counter=0
            col_counter_2=0
            row_counter_2=0
            col_counter_3=0
            row_counter_3=0
            print(vchannels,vrows,vcols)
            #box
            for col in range(vcols):
                if x[i]-(x[0]-x[1])/2<col<x[i]+(x[0]-x[1])/2:
                    col_counter += 1
                    for row in range(vrows):
                        y_height=1
                        if y[len(x)-1-i]-y_height*(x[0]-x[1])<row<y[len(x)-1-i]+y_height*(x[0]-x[1]):
                            row_counter += 1
                            for z in range(vchannels):
                                #print z, vdata[z,row,col]
                                avg_intensity_array[z] += vdata[z,row,col]
            tot_sum=0
            for z in range(len(avg_intensity_array)):
                tot_sum += avg_intensity_array[z]
            avg_intensity_array=avg_intensity_array/(col_counter*row_counter)
            v_array=np.zeros(vchannels)
            for z in range(vchannels):
                v_array[z]=(z-vpix)*vdelt+v0
            #sky subtraction
            '''for col in range(vcols):
                if x[i]-(x[0]-x[1])/2<col<x[i]+(x[0]-x[1])/2:
                    col_counter_2 += 1
                    for row in range(vrows):
                        y_height_2=0.5
                        if y[len(x)-1-i]-(y_height_2*3+y_height)*(x[0]-x[1])<row<y[len(x)-1-i]-(y_height_2*2+y_height)*(x[0]-x[1]):
                            row_counter_2 += 1
                            for z in range(vchannels):
                                #print z, vdata[z,row,col]
                                avg_intensity_array_2[z] += vdata[z,row,col]

            tot_sum_2=0
            for z in range(len(avg_intensity_array)):
                tot_sum_2 += avg_intensity_array_2[z]
            avg_intensity_array_2=avg_intensity_array_2/(col_counter_2*row_counter_2)
            v_array_2=np.zeros(vchannels)
            for z in range(vchannels):
                v_array_2[z]=(z-vpix)*vdelt+v0
            for col in range(vcols):
                if x[i]-(x[0]-x[1])/2<col<x[i]+(x[0]-x[1])/2:
                    col_counter_3 += 1
                    for row in range(vrows):
                        y_height_3=0.5
                        if y[len(x)-1-i]+(y_height+y_height_3*2)*(x[0]-x[1])<row<y[len(x)-1-i]+(y_height+y_height_3*3)*(x[0]-x[1]):
                            row_counter_3 += 1
                            for z in range(vchannels):
                                #print z, vdata[z,row,col]
                                avg_intensity_array_3[z] += vdata[z,row,col]
            tot_sum_3=0
            for z in range(len(avg_intensity_array)):
                tot_sum_3 += avg_intensity_array_2[z]
            avg_intensity_array_3=avg_intensity_array_3/(col_counter_2*row_counter_2)
            v_array_3=np.zeros(vchannels)
            for z in range(vchannels):
                v_array_3[z]=(z-vpix)*vdelt+v0
            plt.plot(v_array,avg_intensity_array_2)
            plt.savefig("/Users/Cohn/Desktop/summer15research/spectra/sky subtracted/bottom box filament"+str(n+1)+"box"+str(len(x)-i)+".png")
            plt.clf()
            plt.plot(v_array,avg_intensity_array_3)
            plt.savefig("/Users/Cohn/Desktop/summer15research/top box filament"+str(n+1)+"box"+str(len(x)-i)+".png")
            plt.clf()
            avg_intensity_array=avg_intensity_array-avg_intensity_array_2-avg_intensity_array_3'''
            plt.plot(v_array,avg_intensity_array)
            plt.savefig("/Users/Cohn/Desktop/summer15research/filament"+str(n+1)+"box"+str(len(x)-i)+".png")
       
        '''c=pyfits.open('/Users/Cohn/Desktop/summer15research/masses with corrected foreground intensity ratios/filament_'+str(n+1)+'_column_density_above_3_background_std_above_mean_to_mass_map.fits')
        coldensedata=c[0].data
        coldenseheader=c[0].header
        ccols=len(coldensedata[:,0])
        crows=len(coldensedata[0])
        lon, lat = w.all_world2pix(l_vals, b_vals, 0)
        for step in range(-int((x_lims[1]-x_lims[0])/2),int((x_lims[1]-x_lims[0])/2)+int((x_lims[1]-x_lims[0])/(2*2)),int((x_lims[1]-x_lims[0])/(2*2))):
            print "x-offset from box center: ", step
            total_mass=0
            total_boxes=0
            sstd=0
            std=0
            for l in range(num_boxes):
                mass_per_box=0
                for col in range(ncols):
                    for row in range(nrows):
                        if imdata[row,col]==0:
                            if x_lims[l]+step<col<x_lims[l+1]+step:
                                mass_per_box += mdata[row,col]
                                total_mass += mdata[row,col]
                if y_vals[l]>0:
                    print "filament: ", n+1, "box number: ", l, "mass per box: ", mass_per_box
                    total_boxes += 1
            avg_mass_per_box=total_mass/total_boxes
            for l in range(num_boxes):
                mass_per_box=0
                for col in range(ncols):
                    for row in range(nrows):
                        if imdata[row,col]==0:
                            if x_lims[l]+step<col<x_lims[l+1]+step:
                                mass_per_box += mdata[row,col]
                if y_vals[l]>0:
                    sstd += (mass_per_box-avg_mass_per_box)**2
            std=(sstd/(total_boxes-1))**0.5
            print"filament: ", n+1, "total mass: ", total_mass, "average mass per box: ", total_mass/total_boxes, "standard deviation: ", std
            print("-------")'''
    
            

