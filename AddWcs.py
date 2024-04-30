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
#wcs = WCS(hdu.header)  
#print(wcs)  

#print(wcs.pixel_n_dim)  
#print(wcs.world_n_dim)  
#print(wcs.array_shape)  

#print(wcs.world_axis_physical_types)  

#coord = wcs.pixel_to_world([1, 2], [4, 3])  
#print(coord)  

# add the wcs

# reference: https://keflavich-astropy.readthedocs.io/en/latest/wcs/

# Set the WCS information manually by setting properties of the WCS
# object.

# Create a new WCS object.  The number of axes must be set
# from the start
w = wcs.WCS(naxis=2)

# Set up an "Airy's zenithal" projection
# Vector properties may be set with Python lists, or Numpy arrays
w.wcs.crpix = [599.764,413.701]
#w.wcs.cdelt = np.array([-0.000053472,0.000053472]) # 0.385 arc seconds to degrees  buuuuuuttt binning 2 was used so maybe each pixel is actually twice that???????????????????????????????
w.wcs.cdelt = np.array([-0.0002138,0.0002138])
w.wcs.crval = [76.95606708,30.401327]  #5.130404472,30.4 '05 07 49.4561' ,'30 24 04.7772'
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
w.wcs.set_pv([(2, 1,180 )]) # axis no. , Parameter number, value
#w.wcs.set_pv([(1,1,360)])

# Three pixel coordinates of interest.
# The pixel coordinates are pairs of [X, Y].
# The "origin" argument indicates whether the input coordinates
# are 0-based (as in Numpy arrays) or
# 1-based (as in the FITS convention, for example coordinates
# coming from DS9).
pixcrd = np.array([[513.833, 533], [711.333, 374.667], [535.5, 690.5]], dtype=np.float64)

# Convert pixel coordinates to world coordinates.
# The second argument is "origin" -- in this case we're declaring we
# have 0-based (Numpy-like) coordinates.
world = w.wcs_pix2world(pixcrd, 0)
print(world)

# Convert the same coordinates back to pixel coordinates.
#pixcrd2 = w.wcs_world2pix(world, 0)
#print(pixcrd2)

#These should be the same as the original pixel coordinates, modulo
# some floating-point error.
#assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6

"""
# The example below illustrates the use of "origin" to convert between
# 0- and 1- based coordinates when executing the forward and backward
# WCS transform.
x = 0
y = 0
origin = 0
assert (w.wcs_pix2world(x, y, origin) ==
        w.wcs_pix2world(x + 1, y + 1, origin + 1))
"""

# Now, write out the WCS object as a FITS header
header = w.to_header()
#s_header = w.to_header_string()

# header is an astropy.io.fits.Header object.  We can use it to create a new
# PrimaryHDU and write it to a file.
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
