# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:18:00 2024

@author: 35385
"""

from scipy.fftpack import fft
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import julian
from datetime import datetime
import math
from astropy import units as u
import ccdproc
import specutils
from astropy.io import fits
from astropy.nddata import CCDData
from specutils import Spectrum1D


###############################################################################
#reference of jd converter https://gist.github.com/jiffyclub/1294443 
def jd_to_date(jd):
    """
    Convert Julian Day to date.
    
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
        4th ed., Duffet-Smith and Zwart, 2011.
    
    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
        
    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
    
    day : float
        Day, may contain fractional part.
        
    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.
    
    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)
    
    """
    jd = jd + 0.5
    
    F, I = math.modf(jd)
    I = int(I)
    
    A = math.trunc((I - 1867216.25)/36524.25)
    
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        
    C = B + 1524
    
    D = math.trunc((C - 122.1) / 365.25)
    
    E = math.trunc(365.25 * D)
    
    G = math.trunc((C - E) / 30.6001)
    
    day = C - E + F - math.trunc(30.6001 * G)
    
    day = math.floor(day) # adjusted to get rid of decimal place.
    
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
        
    return year, month, day
############################################################################### 

magnitude_data = pd.read_csv("aavso_light_data.txt",
sep = ",",header = None, names =[ "JD","Magnitude","Uncertainty",
"HQuncertainty","Band,""Observer Code","Comment Code(s)","Comp Star 1",
"Comp Star 2","Charts","Comments","Transfomed","Airmass","Validation Flag",
"Cmag","Kmag","HJD","Star Name","Observer Affiliation","Measurement Method",
"Grouping Method","ADS Reference","Digitizer","Credit","extra"])


Vmags = []
Vdates = []


Bmags = []
Bdates = []


for i in range(1,len(magnitude_data)):
    #print(str((magnitude_data.iloc[i,3])))
    if(str(magnitude_data.iloc[i,4]) == 'V' ):#or str(magnitude_data.iloc[i,4]) == 'Vis.'):
        Vmags.append(magnitude_data.iloc[i,1])
        Vdates.append(magnitude_data.iloc[i,0])

    if(str(magnitude_data.iloc[i,4]) == 'B'):
        Bmags.append(magnitude_data.iloc[i,1])
        Bdates.append(magnitude_data.iloc[i,0])  
 
   


for i in range(len(Bmags)):
    Bdates[i] = float(Bdates[i])
    Bmags[i] = float(Bmags[i])
    
    
    
    
for i in range(len(Vmags))    :    
    Vdates[i] = float(Vdates[i])
    Vmags[i] = float(Vmags[i])
    
###############################################################################

rdateB = []
rdateV = []

a = datetime(2020,9,14)

# convert to regular date. 
for i in range(len(Bmags)):
    temp = jd_to_date(Bdates[i])
    rdateB.append(datetime(temp[0],temp[1],temp[2]))
    
    
    
    
    
for i in range(len(Vmags))    :    
    temp = jd_to_date(Vdates[i])
    rdateV.append(datetime(temp[0],temp[1],temp[2]))
    
    
    
###############################################################################   
owndataV = 11.0857 # average V magnitude
uncertaintyV = 0.0499 # uncertainty
jdV = 2460234.58565

owndataB = 11.9638
uncertaintyB = 0.1502
jdB = 2460234.57955
    
###############################################################################
"""
plt.figure()
plt.plot(Bdates,Bmags,'bx',label = 'B filter')
plt.plot(Vdates,Vmags,'gx',label = 'V filter')
plt.plot(jdV,owndataV,'rx',label ='Own V data')
plt.plot(jdB,owndataB,'ro',label = 'Own B data')


plt.gca().invert_yaxis()
plt.legend()
plt.xlabel('Julian Date')
plt.ylabel('Magnitude')
plt.grid(True)
"""
###############################################################################
hdu_list = fits.open('C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/Uves_495-707_RwAur.2020-06-26T08 44 09.229.fits')
spec_date = datetime(2008,12,7)
specy = 10.5
#hdu_list.info()
image_data = hdu_list[0].data
spectra = hdu_list[1].data 
wavelength = spectra.WAVE[0]
flux = spectra.FLUX_REDUCED[0]
#headdata = hdu_list[1].header
#print(headdata)

plt.figure()
plt.plot(wavelength,flux,'k')
plt.title('uves 495-707nm')
plt.xlabel('wavelength')
hdu_list.close()
###############################################################################

hdu_list = fits.open('C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/X_Shooter_RwAur_533-1020nm_ADP.2017-05-12T16_14_35.495.fits')

# modified julian date mid obs = 57661.373
spec_date2 = datetime(2016,9,30)
specy2 = 10.8
#hdu_list.info()
image_data = hdu_list[0].data
spectra = hdu_list[1].data 
wavelength = spectra.WAVE[0]
flux = spectra.FLUX_REDUCED[0]
#headdata = hdu_list[1].header
#print(headdata)

plt.figure()
plt.plot(wavelength,flux,'r')
plt.title('Xshooter 533-1020nm 2016')
plt.xlabel('wavelength')
hdu_list.close()

###############################################################################
hdu_list = fits.open('C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/Xshooter_533-1020_Rw_AurADP.2015-04-17T12_28_15.697 (1).fits')

# modified julian date mid obs = 57101.0147619471
spec_date3 = datetime(2015,3,20)
specy3 = 12.5
#hdu_list.info()
image_data = hdu_list[0].data
spectra = hdu_list[1].data 
wavelength = spectra.WAVE[0]
flux = spectra.FLUX_REDUCED[0]
#headdata = hdu_list[1].header
#print(headdata)

plt.figure()
plt.plot(wavelength,flux,'m')
plt.title('Xshooter 533-1020nm 2015')
plt.xlabel('wavelength')
hdu_list.close()

###############################################################################
hdu_list = fits.open('C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/UVES_665-1043nm_Rw_Aur_2009-11-17_ADP.2020-06-19T11_26_47.264.fits')
# modified julian date mid obs = 55152.23429014246
spec_date4 = datetime(2009,11,17)
specy4 = 10.6
#hdu_list.info()
image_data = hdu_list[0].data
spectra = hdu_list[1].data 
wavelength = spectra.WAVE[0]
flux = spectra.FLUX_REDUCED[0]
headdata = hdu_list[1].header
print(headdata)

plt.figure()
plt.plot(wavelength,flux,'orange')
plt.title('UVES 2009 665-1043nm')
plt.xlabel('wavelength')
hdu_list.close()

#TTYPE1  = 'WAVE    '
#TUNIT1  = 'angstrom'   
#TTYPE4  = 'BGFLUX_REDUCED'                                                        
#TUNIT6  = '10**(-16)erg.cm**(-2).s**(-1).angstrom**(-1)'

###############################################################################

V = jd_to_date(jdV)
B  =jd_to_date(jdB)

spec_date5 = datetime(2016,9,30)
specy5 = 11

spec_date6 = datetime(2010,1,3)
specy6 = 10.5
datV = datetime(V[0],V[1],V[2])
datB  =datetime(B[0],B[1],B[2])

plt.figure()
plt.plot_date(rdateB,Bmags,'bx',label = 'B filter')
plt.plot_date(rdateV,Vmags,'gx',label = 'V filter')
#plt.plot_date(datV,owndataV,'rx',label ='Own V data')
#plt.plot_date(datB,owndataB,'ro',label = 'Own B data')
plt.scatter(spec_date,specy,c= 'k',marker = 'o',s = 100,zorder = 100, label = 'UVES 459-707 spectra')
plt.scatter(spec_date2,specy2,c = 'r',marker = 'o',s = 100,zorder = 100,label = 'Xshooter 533-1020 spectra')
plt.scatter(spec_date3,specy3,c = 'm',marker = 'o',s = 100,zorder = 100,label = 'Xshooter 533-1020 spectra')
plt.scatter(spec_date4,specy4,c = 'orange',marker = 'o',s = 100,zorder = 100,label = 'UVES 665-1043nm spectra')
plt.scatter(spec_date5,specy5,c = 'cyan',marker = 'o',s = 100,zorder = 100,label = 'Xshooter 298-556 spectra')
plt.scatter(spec_date6,specy6,c = 'yellow',marker = 'o',s = 100,zorder = 100,label = 'UVES 373-500nm spectra')
plt.xlabel('Date',fontsize = 20)
plt.ylabel('Magnitude',fontsize = 20)
plt.xticks(size = 20)
plt.yticks(size = 20)
plt.title('Spectra used plotted on their corresponding date and magnitude of the light curve.',fontsize = 25)
plt.gca().invert_yaxis()
plt.legend(fontsize = 20)
plt.grid(True)



# from ccdproc and astropy documents

