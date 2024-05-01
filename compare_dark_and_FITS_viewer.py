# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 04:40:08 2023

@author: 35385
"""



import matplotlib.pyplot as plt
from astropy import units as u

import ccdproc
import specutils
import numpy as np
from astropy.io import fits
from astropy.nddata import CCDData

# from ccdproc and astropy documents


"""

# OPEN SPECTRA
hdu_list = fits.open('D:/Spectra/Uves_495-707_RwAur.2020-06-26T08 44 09.229.fits')


#hdu_list.info()
image_data = hdu_list[0].data
spectra = hdu_list[1].data 
wavelength = spectra.WAVE[0]
flux = spectra.FLUX_REDUCED[0]
#print(type(image_data))
#print(image_data.shape)
#hdu_list.close()

plt.figure()
plt.plot(wavelength,flux)
"""


hdu_list = fits.open('D:/4thyearProject/OHP-2023-10-18/Darksnight2/CCD Image 2414.fits')

print(hdu_list)
hdu_list.info()
Darkimage_data = hdu_list[0].data
dark = CCDData(Darkimage_data,unit = 'adu')
#image_data = CCDData(image_data,unit = 'adu')
#image = hdu_list[0].data 
#wavelength = spectra.WAVE[0]
#flux = spectra.FLUX_REDUCED[0]
#print(type(image_data))
#print(image_data.shape)
#headdata = hdu_list[0].header
#print(headdata)
hdu_list.close()
#plt.imshow(image_data, cmap='gray')

#plt.figure()
#plt.plot(wavelength,flux)


hdu_list = fits.open('D:/4thyearProject/OHP-2023-10-16/DBF night1/Bias/CCD Image 2379.fits')
#print(hdu_list)
#hdu_list.info()
biasimage_data = hdu_list[0].data
bias = CCDData(biasimage_data,unit = 'adu')
#image_data = CCDData(image_data,unit = 'adu')
#image = hdu_list[0].data 
#wavelength = spectra.WAVE[0]
#flux = spectra.FLUX_REDUCED[0]
#print(type(image_data))
#print(image_data.shape)
headdata = hdu_list[0].header
#print(headdata)
hdu_list.close()




d = []
b = []


diff_list = []
for i in range(1024):
    for j in range(1024):
        d .append(int(Darkimage_data[i][j]))
        b.append(int(biasimage_data[i][j]))
    
diff = []        
for i in range(len(d)):
    diff.append(b[i]-d[i])
    

      



mean_val = np.mean(diff)
print('Mean difference Value: ', mean_val)


#print(biasimage_data[0][3]-Darkimage_data[0][3])

"""
# OPEN IMAGE
image_data = fits.getdata('D:/4thyearProject/Isacc_newton_telescope/axy.fits')


print(type(image_data)) 
print(image_data.shape)
#plt.imshow(image_data, cmap='gray')
#plt.colorbar()

#print('Min:', np.min(image_data))
#print('Max:', np.max(image_data))
#print('Mean:', np.mean(image_data))
#print('Stdev:', np.std(image_data))
"""


"""
CCD Image 2339.fits Blue filter - 100s - ok
CCD Image 2342.fits Visual filter - 35s - ok
CCD Image 2343.fits Sulphur 2 - 120s - ok 
CCD Image 2344.fits Sulphur 2 - 220s - ok
CCD Image 2347.fits Sulphur 2 - 300s - ok
CCD Image 2349.fits Halpha - 200s - ok 
CCD Image 2350.fits halpha - 250s -ok

"""
