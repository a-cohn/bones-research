#amy cohn
#program for finding the surface density and mass of a dust cloud
#created: june 18, 2015

#import packages
import pyfits as pyfits
import numpy as np
import matplotlib.pyplot as plt

distances=np.array([4.6,3.82,3.24,4.26,3.60,3.34,3.1,3.0,3.2,3.3]) #load distances
fir=np.array([0.228110503699,0.166105172217,0.129696069108,0.18066353617,0.132235411746,0.100190511107,0.0794047216383,0.0740902435082,0.126654172247,0.142069364321]) #load foreground intensity ratios

#for filaments 
for n in range(1,11):
    if n>11:
        pass
    else:
        #open data file
        f=pyfits.open('/Users/Cohn/Desktop/summer15research/filament_'+str(n)+'_cropped_with_header.fits')
        data=f[0].data
        header=f[0].header

        #open image mask
        im=pyfits.open('/Users/Cohn/Desktop/summer15research/imagemask_filament_'+str(n)+'_cropped_with_header.fits')
        imdt=im[0].data
        imdata=np.asarray(imdt)

        #open interpolation file
        intp=pyfits.open('/Users/Cohn/Desktop/summer15research/interpolation_filament_'+str(n)+'_cropped_with_header_unregridded.fits')
        intpdt=intp[0].data
        intpdata=np.asarray(intpdt)

        #measuring the size of the file
        nrows=len(data[:,0]) #number of rows of fits file
        ncols=len(data[0]) #number of columns of fits file

        #import measurements
        I_f=np.zeros((nrows,ncols)) #radiation in front of the cloud
        I_b=np.zeros((nrows,ncols)) #radiation behind the cloud
        for row in range(nrows):
            #print row, nrows
            for col in range(ncols):
                if imdata[row,col]==0:
                    I_f[row,col]=data[row,col] #radiation in front of the cloud
                    I_b[row,col]=intpdata[row,col] #radiation behind the cloud

        #define constants
        do=11.7 #dust opacity
        s=0.3 #scattering
        f_f=fir[n-1] #foreground intensity ratio
    
        sd=np.zeros((nrows,ncols)) #surface density

        #creating surface density map
        for row in range(nrows):
            for col in range(ncols):
                if imdata[row,col]==0:
                    sd[row,col]=-1/do*np.log(((s+1)*I_f[row,col]-(s+f_f)*I_b[row,col])/((1-f_f)*I_b[row,col]))
            
        #save surface density map to fits file

        pyfits.writeto('/Users/Cohn/Desktop/summer15research/masses with corrected foreground intensity ratios/filament_'+str(n)+'_surface_density_map_with_header_new.fits',sd,header)

        '''#find mass
        mass0=0
        pix=0
        sum_sd=0
        for row in range(nrows):
            for col in range(ncols):
                #print "mass round", row, nrows
                if imdata[row,col]==0:
                    if sd[row,col]<0.0:
                        sd[row,col]=0
                    elif sd[row,col]>0.0:
                        mass0 += (sd[row,col]*(distances[n-1]*3.0857*10**21*1.2/206265)**2)/(1.989*10**33)
                        sum_sd += sd[row,col]
                        pix += 1
        avg_sd=sum_sd/pix
        avg_mass=mass0/pix
        differences2=0
        differences2_mass=0
        for row in range(nrows):
            for col in range(ncols):
                if imdata[row,col]==0:
                    if sd[row,col]>0:
                        differences2 += (sd[row,col]-avg_sd)**2
                        differences2_mass += ((sd[row,col]*(distances[n-1]*3.0857*10**21*1.2/206265)**2)/(1.989*10**33)-avg_mass)**2
        std=(differences2/(pix-1))**0.5
        std_mass=(differences2_mass/(pix-1))**0.5
        sstd=pix**0.5*std
        sstd_mass=pix**0.5*std_mass
        #print n, sstd_mass
        mass1=0
        mass2=0
        for row in range(nrows):
            for col in range(ncols):
                if imdata[row,col]==0:
                    if sd[row,col]>avg_sd+std:
                        mass1 += (sd[row,col]*(distances[n-1]*3.0857*10**21*1.2/206265)**2)/(1.989*10**33)
                    if sd[row,col]>avg_sd+2*std:
                        mass2+= (sd[row,col]*(distances[n-1]*3.0857*10**21*1.2/206265)**2)/(1.989*10**33)
        print "filament: ", n
        print "number of pix: ", pix
        print "avg_mass: ", avg_mass
        print "std_mass: ", std_mass
        print "mass above 0: ", mass0
        print "standard deviation of summed mass: ", sstd_mass
        print "mass above 1 sigma: ", mass1
        print "mass above 2 sigma: ", mass2'''
   
