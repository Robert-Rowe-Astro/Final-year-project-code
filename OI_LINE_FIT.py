# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:46:17 2024

@author: 35385
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 00:37:58 2024

@author: 35385
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 11:21:38 2024

@author: 35385
"""




import numpy as np
import matplotlib.pyplot as plt


from astropy import units as u
#import ccdproc
import specutils
from astropy.io import fits


from specutils import Spectrum1D, SpectralRegion
from specutils.fitting import fit_generic_continuum
from specutils.analysis import equivalent_width 

from astropy.visualization import quantity_support

from specutils.analysis import line_flux

from specutils.manipulation import extract_region


from astropy.modeling import models ,fitting
#from astropy.modeling.models import Blackbody
from dust_extinction.parameter_averages import G23

from scipy.stats import norm

from scipy.integrate import simpson
from numpy import trapz


import math

###############################################################################

#Halpha_vacuum = 6564.614 * u.AA

OI_vacuum = 6365.536 * u.AA

CaII_k_vacuum = 3933.66 * u.AA



OI557 = 5577 * u.AA

OI630 = 6300 * u.AA # get ref

OI636 = 6360 * u.AA # NOT VIABLE

SII406 = 4060.9 * u.AA

SII671 = 67100 * u.AA

SII673 = 6731 * u.AA # get ref

NII658 = 6580 * u.AA



L_solar_erg = 3.839e33 # erg s**-1
L_solar_watt = 3.839e26 # watts
R_solar_m = 6.95508e8 # metres
M_solar_kg = 1.9891e30 #kg  

G = 6.67e-11
e = 2.7

###############################################################################
quantity_support()

ext = G23(Rv = 3.1)

#******************************************************************************

# gaussian fitting and line fluxes reference:
# ref https://lukeholden.com/blog/measuring_emission_line_fluxes_python.html

def gaussfit(a,b,c,x):
    
    gauss_fit = []
    for i in range(len(x)):
        
       gauss_fit.append( a * np.exp(-((x-b)**2)/(2*c**2)))
       
       
    return gauss_fit   


def flux_calc(peak,width): # when doing error ,peak_err,width_err
    flux = peak * width * np.sqrt(2*np.pi)
    return flux



def Spectra_Analyse(filename,instrument,plot_colour,plot_title,lam_vacuum,v_correct):
    hdu_list = fits.open(filename)
    #spec_date = datetime(2008,12,7)
    #specy = 10.5
    #hdu_list.info()
     
    # read the data 
    image_data = hdu_list[0].data
    spectra = hdu_list[1].data 
    wavelength1 = spectra.WAVE[0] 
    flux1 = spectra.FLUX[0] 
    
    # x shooter wavlength is in nanometres by default
    if(instrument != 'UVES'):
        wavelength1 = wavelength1 * 10
    else:
         flux1 = flux1 * 1e-16    
        
        
    wavelength1 = wavelength1 * u.AA
    flux1 = flux1 * u.Unit('erg cm**-2 s**-1 nm**-1')
    
    flux1 = flux1 / ext.extinguish(wavelength1,Av = 1.0)


    

    #create spectrum1d object 
    spec1 = Spectrum1D(spectral_axis = wavelength1,flux = flux1,radial_velocity = v_correct * u.Unit('km/s'),velocity_convention ='optical',rest_value = lam_vacuum)
    
    #continuum subtract the spectrum1d object
    fit = fit_generic_continuum(spec1,exclude_regions=[SpectralRegion((lam_vacuum - 10* u.AA) ,(lam_vacuum + 10* u.AA))])
    yfit = fit(spec1.wavelength)
    spec1sub = spec1 - yfit 
    #spec1sub = spec1 - yfit
    
    #take region around peak
    
    
    #--------------------------------------------------------------------------
    region = SpectralRegion((lam_vacuum -10 * u.AA),(lam_vacuum + 10 * u.AA))

    
    line_regionwav = extract_region(spec1,region)
    
    line_regionflux = extract_region(spec1sub,region)
    
    
    #extcorr = line_regionflux.flux / ext.extinguish(line_regionwav.wavelength,Av = 1.0)
    
    """
    sson = simpson(line_regionflux.flux,line_regionwav.wavelength)
    print('\n')
    print('\n')
    #print('The flux calculated by specutils for ', plot_title,' is:')
    #print(sson)
    
    line_lum = Macc(sson,183,1.13,1.74,1.5,1.8)
    print('\n')
    print('\n')
    print('The line luminosity of ', plot_title,' at ', lam_vacuum, ' is:')
    print(line_lum)
    print('\n')
    print('\n')
    """
    
    #line_lum = Macc(specutilsflux,183,1.13,1.74,1.5,1.8)
    
    
    return spec1,spec1sub, line_regionwav,line_regionflux, 7 #,extcorr


###############################################################################
def Macc(fline,dist,a,b,mass,radius):
    
    #convert dist from parsecs to m then cm
    
    
    dm = dist * 3.08567758e16 # dist in metres
    
    dcm = dm * 100  # dist in centimetres
    
    Lline = (4 * np.pi *  (dcm**2)) * (fline) # fline in erg/s cm
    
    divline = Lline/L_solar_erg
    
    #Log_Lacc = math.log10(divline)
    
    #Log_Lacc = (Log_Lacc * a) + b
    
    #print(Log_Lacc)
    
    #Lacc = 10**Log_Lacc
    
    #Macc = ((Lacc * L_solar_watt))*((radius * R_solar_m )) / ((G)*(mass * M_solar_kg))
    
    #Macc = Macc * 1.25 # result is in kg/s 
    
    #Macc = Macc * 31536000 # get in terms of years
    
    #Macc = Macc / 1.9819e30 # get in terms of solar mass 
    
    #log_Macc = math.log10(Macc)
    
    #return log_Macc
    return Lline
    
    #print(str(filename) + ' spectra complete')

    #--------------------------------------------------------------------------    
        
       

   


#******************************************************************************
###############################################################################

#spectra files to process
#******************************************************************************
file_list = ['C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/Uves_495-707_RwAur.2020-06-26T08 44 09.229.fits',
             'C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/X_Shooter_RwAur_533-1020nm_ADP.2017-05-12T16_14_35.495.fits',
             'C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/Xshooter_533-1020_Rw_AurADP.2015-04-17T12_28_15.697 (1).fits',
             'C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/UVES_res58640_snr185.8_2010-01-03_373.2-500nm_.fits',
             'C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/Xshooter_RwAur_res9861_snr300.7_298-556nm_2016-09-30.fits']

Instrument_list = ['UVES','XSH','Xsh','UVES','Xsh']

rad_vel = [1.416323309,27.850689942181,-29.461835387817,-12.126655856364,27.850832009443]

unit_scale = [0,1e-16,1e-15,1e-16,1e-16]
#******************************************************************************



# initialise extinction profile
#ext = G23(Rv = 3.1)



s4wav, s5flux, reg4wav,reg4flux, regext4 = Spectra_Analyse(file_list[1],Instrument_list[1], 'red' , 'Xshooter 2016 533-1020nm', OI630,rad_vel[1])


s5wav, sflux, reg5wav,reg5flux, regext5 = Spectra_Analyse(file_list[2],Instrument_list[2], 'm' , 'Xshooter 2015 533-1020', OI630,rad_vel[2])


# Gaussian 5

#gauss1init4 = models.Gaussian1D(amplitude = 1.21e-13, mean =  -105.41, stddev = 145)#150 #145
#gauss2init4 = models.Gaussian1D(amplitude = 1.66e-13, mean =  -30, stddev = 15) #20 #10 #10 #15
#gauss3init4 = models.Gaussian1D(amplitude = 1.38e-13, mean =  75.45, stddev = 10) #70 # 100 #120 #150 # 5 #10
#gauss4init4 = models.Gaussian1D(amplitude = 1.24e-13, mean =  104.04, stddev = 15)#40 # 60 #30   # 15 #15

gauss1init4 = models.Gaussian1D(amplitude = 1.655e-13, mean =  -38.73, stddev = 50)#150 #145
gauss2init4 = models.Gaussian1D(amplitude = 1.655e-13, mean =  -19.8, stddev = 15) #20 #10 #10 #15
gauss3init4 = models.Gaussian1D(amplitude = 1.316e-13, mean =  74.45, stddev = 10) #70 # 100 #120 #150 # 5 #10
gauss4init4 = models.Gaussian1D(amplitude = 1.213e-13, mean =  106.2, stddev = 30)#40 # 60 #30   # 15 #15



#gtot1init4 = gauss1init4  +  gauss2init4 + gauss3init4 + gauss4init4

gtot1init4 = gauss1init4  +  gauss2init4 + gauss3init4 + gauss4init4

fitgauss4 = fitting.LevMarLSQFitter()

gtot4 = fitgauss4(gtot1init4,reg4wav.velocity,reg4flux.flux)



a1 = gtot4.parameters[0]
b1 = gtot4.parameters[1]
c1 = gtot4.parameters[2]

a2 = gtot4.parameters[3]
b2 = gtot4.parameters[4]
c2 = gtot4.parameters[5]

a3 = gtot4.parameters[6]
b3 = gtot4.parameters[7]
c3 = gtot4.parameters[8]

a4 = gtot4.parameters[9]
b4 = gtot4.parameters[10]
c4 = gtot4.parameters[11]

#r4 = float(reg4wav.wavelength)


#g1 = gaussfit(a1,b1,c1,r4)

#g2 = gaussfit(a2,b2,c2,r4)

#g3 = gaussfit(a3,b3,c3,r4)

#g4 = gaussfit(a4,b4,c4,r4)








###############################################################################

plt.figure()
plt.plot(reg4wav.velocity,reg4flux.flux,'r',label='Xshooter 533-1020nm 2016')
plt.plot(reg4wav.velocity, gtot4(reg4wav.velocity),'blue',linestyle = '--',label = 'gaussian fit 4')

#plt.plot(reg4wav.wavelength, g1)
#plt.plot(reg4wav.wavelength, g2)
#plt.plot(reg4wav.wavelength, g3)
#plt.plot(reg4wav.wavelength, g4)

plt.xlabel('wavelength (angstroms)',fontsize = 20)
plt.ylabel('flux erg $cm^{-2} s^{-1}angstrom^{-1}$',fontsize = 20)
plt.title('Gaussians fitted to the X-shooter 2016 [OI] line. ', fontsize = 20)

plt.legend(fontsize = 20)
plt.xticks(size = 20)
plt.yticks(size = 20)
plt.grid()

diff = reg4flux.flux - gtot4(reg4wav.velocity) 
square_diff = diff**2
divsqdiff = square_diff/reg4flux.flux
tot = 0
for i in range(len(divsqdiff)):
    tot = tot + divsqdiff[i]
    
print(tot)
plt.figure()
plt.plot(reg4wav.velocity,diff)




plt.figure()
plt.xticks(size = 20)
plt.yticks(size = 20)
gtot_sep = norm.pdf(reg4wav.velocity,-119.5, 76.5)
plt.plot(reg4wav.velocity,gtot_sep )

gtot_sep2 = norm.pdf(reg4wav.velocity,-32, 32)
plt.plot(reg4wav.velocity,gtot_sep2 )

gtot_sep3 = norm.pdf(reg4wav.velocity,71.4,3.28)
plt.plot(reg4wav.velocity,gtot_sep3 )

gtot_sep4 = norm.pdf(reg4wav.velocity,89.7, 62)
plt.plot(reg4wav.velocity,gtot_sep4)

###############################################################################


###############################################################################


print('oi gaussian fit  paramters: ')
print(gtot4.parameters)
print('\n')
print('oi gaussian covariance matrix')
print(fitgauss4.fit_info['param_cov'])
print('\n')


flux = flux_calc(8.85058022e-14,3.20605845e+01)

print(Macc(flux,183,1.13,1.74,1.5,1.8))


chisq = 0
for i in range(len(reg4wav.velocity)):
     cons = ((reg4flux.flux[i] - gtot4(reg4wav.velocity[i]))**2)/reg4flux.flux[i]
     chisq = chisq + cons
    

print(chisq)
    
print(len(reg4wav.velocity))  

reducedchi = chisq/(len(reg4flux.flux))

print(reducedchi)

    


