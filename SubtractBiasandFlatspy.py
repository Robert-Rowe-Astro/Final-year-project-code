# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 23:45:38 2023


'Hardcoded' image reduction

@author: Robert Rowe
"""
import astroalign as aa

import matplotlib.pyplot as plt
from astropy import units as u

import ccdproc

import numpy as np
from astropy.io import fits
from astropy.nddata import CCDData


# TAKE IN IMAGES
BIm1 = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-16/RR CTT/CCD Image 2339.fits') # Blue filter image 1 100s

VIm1 = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-16/RR CTT/CCD Image 2342.fits') # V filter image 1 35s

S2Im1 = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-16/RR CTT/CCD Image 2343.fits') # S2 Image 1 120s

S2Im2 = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-16/RR CTT/CCD Image 2344.fits') # S2 image 2 220s

S2Im3 = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-16/RR CTT/CCD Image 2347.fits') # S2 image 3 300s

HalphaIm1 = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-16/RR CTT/CCD Image 2349.fits') # Ha image 1 200s

HalphaIm2 = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-16/RR CTT/CCD Image 2350.fits') # Ha image 2 250s

#*****************************************************************************#


BIm1 = CCDData(BIm1,unit = 'adu')

VIm1 = CCDData(VIm1,unit = 'adu')

S2Im1 = CCDData(S2Im1,unit = 'adu')

S2Im2 = CCDData(S2Im2,unit = 'adu')

S2Im3 = CCDData(S2Im3,unit = 'adu')

HIm1 = CCDData(HalphaIm1,unit = 'adu')

HIm2 = CCDData(HalphaIm2,unit = 'adu')

#*****************************************************************************#

drk = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-18/Darksnight2MasterDark.FITS')

drk = CCDData(drk,unit = 'adu')

#*****************************************************************************#
# take in the masterbias and masterflat from 16/10/23
# that were created from python scripts projectdata.py and
# astropy-methods.py


# B FILTER
MBias = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-16/DBF night1/BiasMasterbias.FITS')

MBias = CCDData(MBias,unit = 'adu')


MFlatB = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-16/DBF night1/Flat/MasterFlats/B flatsMasterBflat.FITS')

MFlatB = CCDData( MFlatB,unit = 'adu' )


#*****************************************************************************#
# 

BiasSubtractedDark = ccdproc.subtract_bias(drk,MBias)




#*****************************************************************************#

BiasSubtractedFlatB = ccdproc.subtract_bias(MFlatB, MBias) 
darkSubtractedFlatB = ccdproc.subtract_dark(BiasSubtractedFlatB,drk, dark_exposure=(3600 * u.second), data_exposure = (0.5 * u.second), scale = True)



BiasSubtractedBIm1 = ccdproc.subtract_bias(BIm1,MBias)
darkSubtractedBIm1 = ccdproc.subtract_dark(BiasSubtractedBIm1,drk, dark_exposure=(3600 * u.second), data_exposure = (100 * u.second), scale = True)



FinalBIm1 = ccdproc.flat_correct(BiasSubtractedBIm1,BiasSubtractedFlatB )

#FinalBIm1.write('FinalBIm1incdark.fits') # already written

#*****************************************************************************#




# V FILTER


MFlatV = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-16/DBF night1/Flat/MasterFlats/V flatsMasterVflat.FITS')

MFlatV = CCDData(MFlatV,unit = 'adu')



#*****************************************************************************#

BiasSubtractedFlatV = ccdproc.subtract_bias(MFlatV, MBias) 
darkSubtractedFlatV = ccdproc.subtract_dark(BiasSubtractedFlatV,drk, dark_exposure=(3600 * u.second), data_exposure = (0.5 * u.second), scale = True)



BiasSubtractedVIm1 = ccdproc.subtract_bias(VIm1,MBias)
darkSubtractedVIm1 = ccdproc.subtract_dark(BiasSubtractedVIm1,drk, dark_exposure=(3600 * u.second), data_exposure = (35 * u.second), scale = True)



FinalVIm1 = ccdproc.flat_correct(BiasSubtractedVIm1,BiasSubtractedFlatV )

#FinalVIm1.write('FinalVImincdark.fits')

#*****************************************************************************#




# HALPHA FILTER


MFlatH = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-16/DBF night1/Flat/MasterFlats/Halpha flatsMasterHALPHAflat.FITS')

MFlatH = CCDData(MFlatH,unit = 'adu')



#*****************************************************************************#

BiasSubtractedFlatH = ccdproc.subtract_bias(MFlatH, MBias) 
darkSubtractedFlatH = ccdproc.subtract_dark(BiasSubtractedFlatH,drk, dark_exposure=(3600 * u.second), data_exposure = (40 * u.second), scale = True)



BiasSubtractedHIm1 = ccdproc.subtract_bias(HIm1,MBias)
darkSubtractedHIm1 = ccdproc.subtract_dark(BiasSubtractedHIm1,drk, dark_exposure=(3600 * u.second), data_exposure = (200 * u.second), scale = True)



FinalHIm1 = ccdproc.flat_correct(BiasSubtractedHIm1,BiasSubtractedFlatH )

#FinalHIm1.write('FinalHIm1Woutdark.fits')

#*****************************************************************************#




# HALPHA FILTER Image 2


#BiasSubtractedFlatH = ccdproc.subtract_bias(MFlatH, MBias) 
#darkSubtractedFlatH2 = ccdproc.subtract_dark(BiasSubtractedFlatH,drk, dark_exposure=(3600 * u.second), data_exposure = (40 * u.second), scale = True)



BiasSubtractedHIm2 = ccdproc.subtract_bias(HIm2,MBias)
darkSubtractedHIm2 = ccdproc.subtract_dark(BiasSubtractedHIm2,drk, dark_exposure=(3600 * u.second), data_exposure = (250 * u.second), scale = True)



FinalHIm2 = ccdproc.flat_correct(BiasSubtractedHIm2,BiasSubtractedFlatH )

#FinalHIm2.write('FinalHIm2Woutdark.fits')

#*****************************************************************************#
















# S2 FILTER


MFlatS2 = fits.getdata('C:/Users/35385/OneDrive/Desktop/4thyearProject/OHP-2023-10-16/DBF night1/Flat/MasterFlats/S2 flatsMasterS2flat.FITS')

MFlatS2 = CCDData(MFlatS2,unit = 'adu')



#*****************************************************************************#

BiasSubtractedFlatS2 = ccdproc.subtract_bias(MFlatS2, MBias) 
darkSubtractedFlatS2 = ccdproc.subtract_dark(BiasSubtractedFlatS2,drk, dark_exposure=(3600 * u.second), data_exposure = (40 * u.second), scale = True)



BiasSubtractedS2Im1 = ccdproc.subtract_bias(S2Im1,MBias)
darkSubtractedS2Im1 = ccdproc.subtract_dark(BiasSubtractedS2Im1,drk, dark_exposure=(3600 * u.second), data_exposure = (120 * u.second), scale = True)



FinalS2Im1 = ccdproc.flat_correct( BiasSubtractedS2Im1,BiasSubtractedFlatS2 )

#FinalS2Im1.write('FinalS2Im1Woutdark.fits')

#*****************************************************************************#


# S2 FILTER IMAGE 2


BiasSubtractedS2Im2 = ccdproc.subtract_bias(S2Im2,MBias)
darkSubtractedS2Im2 = ccdproc.subtract_dark(BiasSubtractedS2Im2,drk, dark_exposure=(3600 * u.second), data_exposure = (220 * u.second), scale = True)



FinalS2Im2 = ccdproc.flat_correct(BiasSubtractedS2Im2,BiasSubtractedFlatS2 )

#FinalS2Im2.write('FinalS2Im2Woutdark.fits')

#*****************************************************************************#

# S2 FILTER IMAGE 3


BiasSubtractedS2Im3 = ccdproc.subtract_bias(S2Im3,MBias)
darkSubtractedS2Im3 = ccdproc.subtract_dark(BiasSubtractedS2Im3,drk, dark_exposure=(3600 * u.second), data_exposure = (300 * u.second), scale = True)



FinalS2Im3 = ccdproc.flat_correct(BiasSubtractedS2Im3,BiasSubtractedFlatS2 )


#FinalS2Im3.write('FinalS2Im3Woutdark.fits')


#*****************************************************************************#

#S2incD = [FinalS2Im1,FinalS2Im2,FinalS2Im3]

#S2FinalIncDark = ccdproc.combine(S2incD,'C:/Users/35385/OneDrive/Desktop/4thyearProject/CombinedS2IncDark.FITS','median')

#HincD = [FinalHIm1,FinalHIm2]

#HFinalIncDark = ccdproc.combine(S2incD,'C:/Users/35385/OneDrive/Desktop/4thyearProject/CombinedHIncDark.FITS','median')

# align the images



# ref https://pypi.org/project/astroalign/

s2alignedim1, s2footprint1 = aa.register(FinalS2Im1,FinalS2Im3)
s2alignedim1 =  CCDData(s2alignedim1,unit = 'adu')



s2alignedim2, s2footprint2 = aa.register(FinalS2Im2,FinalS2Im3)
s2alignedim2 =  CCDData(s2alignedim2,unit = 'adu')

Haligned, Hfootprint = aa.register(FinalHIm1,FinalHIm2)
Haligned = CCDData(Haligned,unit = 'adu')




S2list = [s2alignedim1,s2alignedim2,FinalS2Im3]
alignedSulphurFinal = ccdproc.combine(S2list,'C:/Users/35385/OneDrive/Desktop/4thyearProject/AlignedCombinedS2.FITS','median')


hlist = [Haligned,FinalHIm2]
HFinalIncDark = ccdproc.combine(hlist,'C:/Users/35385/OneDrive/Desktop/4thyearProject/AlignedCombinedH.FITS','median')













