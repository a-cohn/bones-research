#Find an estimate of the fraction of foreground emssision at 
# a distance 'dist' and longitude 'l' assuming an exponentially 
# decreasing galactic surface density
# Distance is in 'kpc' and longitude (l) is in degrees

#Adapted from Butler & Tan 2009
# Written by C. Battersby March 8, 2009

import numpy as np
import math

l_array=np.array([26.94, 25.25, 24.95, 21.25, 18.88, 11.13, 4.14, 357.62, 335.31, 332.21])
dist_array=np.array([4.6,3.82,3.24,4.26,3.6,3.34,3.1,3.0,3.2,3.3])

for n in range(1,11):
    dgc = 8.  #galactocentric distance to sun
    l=l_array[n-1]
    dist=dist_array[n-1]
    ll = l
    #Convert l to radians, and calculate it's cosine
    ll = ll*3.14159625/180. #will need to make statement 
    cosl = np.cos(ll)  #of l more general perhaps

    #increment for numerical summation (pseudo-column-integration)
    inc = 0.001 #distance increment in kpc
    hr = 3.5 #scale length in kpc
    dd = 0.001 # distance along the column
    fore = 0. #initialize foreground estimate
    back = 0. # and background estimate

    #Between the sun and the cloud (at distance, d, calculate the 
    # galactocentric radius (using the law of cosines) at each
    # distance increment, and sum these to get the total 
    # foreground emission estimate (not normalized, only
    # care about the ratio)
    while dd<dist:
        r = (dd**2 + dgc**2 - 2.*dd*dgc*cosl)**0.5
        fore = fore + math.exp(-r/hr)*inc
        dd = dd + inc


    dd = 0.001
    r = 0.
    #Between the sun and a galactocentric radius of 16 kpc, calculate
    # the galactocentric radius (using the law of cosines) at each
    # distance increment and sum these to get the total background
    # emission along that column.
    while r<16:
        r = (dd**2 + dgc**2 - 2.*dd*dgc*cosl)**0.5
        back = back + math.exp(-r/hr)*inc
        dd = dd + inc


    #Take ratio, print, and return
    ffore = fore/back
    print "filament", n, "fore, back, ffore: ", fore, back, ffore

