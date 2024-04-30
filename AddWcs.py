# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 11:09:49 2023

@author: 35385
"""
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

import numpy as np

from astropy import wcs
from astropy.io import fits

# get image info 

filename = get_pkg_data_filename('C:/Users/35385/OneDrive/Desktop/4thyearProject/WCS_test/AlignedCombinedH_test.FITS')
fits.info(filename)  
hdulist = fits.open(filename)  
hdu = hdulist[0]  




# reference: https://keflavich-astropy.readthedocs.io/en/latest/wcs/

# Set the WCS information manually by setting properties of the WCS
# object.


w = wcs.WCS(naxis=2)



w.wcs.crpix = [599.764,413.701]
#w.wcs.cdelt = np.array([-0.000053472,0.000053472]) #binning 2 used so twice this
w.wcs.cdelt = np.array([-0.0002138,0.0002138])
w.wcs.crval = [76.95606708,30.401327]   '05 07 49.4561' ,'30 24 04.7772'
w.wcs.ctype = ["RA---TAN", "DEC--TAN"] # ra--air, dec--air
w.wcs.set_pv([(2, 1,180 )]) # axis no. , Parameter number, value
#w.wcs.set_pv([(1,1,360)])

# reference pixels to test
pixcrd = np.array([[513.833, 533], [711.333, 374.667], [535.5, 690.5]], dtype=np.float64)


# for 0-based coordinates
world = w.wcs_pix2world(pixcrd, 0)
print(world)

# Convert the same coordinates back to pixel coordinates.
#pixcrd2 = w.wcs_world2pix(world, 0)
#print(pixcrd2)





#write to header
header = w.to_header()
#s_header = w.to_header_string()

hdu = fits.PrimaryHDU(header=header)
# Save to FITS file
#hdu.writeto('AlignedCombinedH.FITS',overwrite = True) # doesnt work

with fits.open('C:/Users/35385/OneDrive/Desktop/4thyearProject/AlignedCombinedH.FITS',mode = 'update') as hdul:
    hdr = hdul[0].header
    #del hdr[8]
    
    for i in range(len(header)):
       
        
        #hdul[0].header.append(str(header._keyword_from_index(i)[0]), str(header[i]),end=True)
        hdr[str(header._keyword_from_index(i)[0])] =  (str(header[i]))
        
        #print(i)
        #print(hdr[str(header._keyword_from_index(i)[0])])
        hdul.flush()
    
    
    

#print(len(header))

#print(header)
#print(header[0],[1])


#print(header._keyword_from_index(0)[0])
