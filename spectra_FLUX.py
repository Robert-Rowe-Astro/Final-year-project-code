# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 18:20:43 2024

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

Halpha_vacuum = 6564.614 * u.AA # ref https://classic.sdss.org/dr7/products/spectra/vacwavelength.php

OI_vacuum = 6365.536 * u.AA # REF https://classic.sdss.org/dr6/algorithms/linestable.php

CaII_k_vacuum = 3933.66 * u.AA # ref (for now) wikipedia frauenhofer liines


#ref carrol ostlie

L_solar_erg = 3.839e33 # erg s**-1
L_solar_watt = 3.839e26 # watts
R_solar_m = 6.95508e8 # metres
M_solar_kg = 1.9891e30 #kg  

G = 6.67e-11

###############################################################################
quantity_support()

def Spectra_Analyse(filename,instrument,plot_colour,plot_title,lam_vacuum,v_correct,unitscale):
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
    
    flux1 = flux1 * u.Unit('erg cm**-2 s**-1 AA**-1')

   
    

    #create spectrum1d object 
    spec1 = Spectrum1D(spectral_axis = wavelength1,flux = flux1,radial_velocity = v_correct * u.Unit('km/s'),velocity_convention ='optical',rest_value = lam_vacuum)
    
    #continuum subtract the spectrum1d object
    fit = fit_generic_continuum(spec1,exclude_regions=[SpectralRegion((lam_vacuum - 10* u.AA) ,(lam_vacuum + 10* u.AA))])
    yfit = fit(spec1.wavelength)
    spec1sub = spec1 - yfit  #-1
    #spec1sub = spec1 - yfit
    

    
    #take region around peak
    #--------------------------------------------------------------------------
    region = SpectralRegion((lam_vacuum-15 * u.AA),(lam_vacuum + 15 * u.AA))

    
    line_regionwav = extract_region(spec1,region)
    
    line_regionflux = extract_region(spec1sub,region)
    
    
    specutilsflux = line_flux(line_regionflux)
    print('\n')
    print('The flux calculated by specutils for ', plot_title,' is:')
    print(specutilsflux)
    print('\n')
    
 
    
    return spec1,spec1sub, line_regionwav,line_regionflux


#******************************************************************************
###############################################################################
def Macc(fline,dist,a,b,mass,radius):
    
    #convert dist from parsecs to m then cm
    
    
    dm = dist * 3.08567758e16 # dist in metres
    
    dcm = dm * 100  # dist in centimetres
    
    Lline = (4 * np.pi *  (dcm**2)) * (fline) # fline in erg/s cm
    
    divline = Lline/L_solar_erg
    
    Log_Lacc = math.log10(divline)
    
    Log_Lacc = (Log_Lacc * a) + b
    
    print(Log_Lacc)
    
    Lacc = 10**Log_Lacc
    
    Macc = ((Lacc * L_solar_watt))*((radius * R_solar_m )) / ((G)*(mass * M_solar_kg))
    
    Macc = Macc * 1.25 # result is in kg/s 
    
    Macc = Macc * 31536000 # get in terms of years
    
    Macc = Macc / 1.9819e30 # get in terms of solar mass 
    
    log_Macc = math.log10(Macc)
    
    return log_Macc
    
    
    
    





#spectra files to process
#******************************************************************************
file_list = ['C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/Uves_495-707_RwAur.2020-06-26T08 44 09.229.fits',
             'C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/X_Shooter_RwAur_533-1020nm_ADP.2017-05-12T16_14_35.495.fits',
             'C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/Xshooter_533-1020_Rw_AurADP.2015-04-17T12_28_15.697 (1).fits',
             'C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/UVES_res58640_snr185.8_2010-01-03_373.2-500nm_.fits',
             'C:/Users/35385/OneDrive/Desktop/4thyearProject/Spectra_Rwaur/Xshooter_RwAur_res9861_snr300.7_298-556nm_2016-09-30.fits']

Instrument_list = ['UVES','XSH','Xsh','UVES','Xsh']

rad_vel = [1.416323309,27.850689942181,-29.461835387817,-12.126655856364,27.850832009443]

unit_correct = [0,1e-16,1e-16,1e-16,1e-16]
#******************************************************************************
# initialise extinction profile
ext = G23(Rv = 3.1)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#s1wav, s1flux, reg1wav,reg1flux = Spectra_Analyse(file_list[0],Instrument_list[0],'k','uves 495-707nm 2008',Halpha_vacuum,rad_vel[0])
#reg1flux.new_flux_unit(u.Unit('erg cm**-2 s**-1 nm**-1'))
#reg1extflux = reg1flux.flux / ext.extinguish(reg1wav.wavelength,Av = 0.39)

s2wav, s2flux, reg2wav,reg2flux = Spectra_Analyse(file_list[1],Instrument_list[1],'r','Xshooter 533-1020nm 2016',Halpha_vacuum,rad_vel[1],unit_correct[1])
#reg2flux.new_flux_unit(u.Unit('erg cm**-2 s**-1 nm**-1'))
reg2extflux = reg2flux.flux / ext.extinguish(reg2wav.wavelength,Av = 1.0)


s3wav, s3flux, reg3wav,reg3flux = Spectra_Analyse(file_list[2],Instrument_list[2],'m','Xshooter 533-1020nm 2015',Halpha_vacuum,rad_vel[2],unit_correct[2])
#reg3flux.new_flux_unit(u.Unit('erg cm**-2 s**-1 nm**-1'))
reg3extflux = reg3flux.flux / ext.extinguish(reg3wav.wavelength,Av = 4.0)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

s4wav, s5flux, reg4wav,reg4flux = Spectra_Analyse(file_list[3],Instrument_list[3], 'yellow' , 'UVES 2010  373-500nm', CaII_k_vacuum,rad_vel[3],unit_correct[3])
reg4extflux = reg4flux.flux / ext.extinguish(reg4wav.wavelength,Av = 1.0)

s5wav, sflux, reg5wav,reg5flux = Spectra_Analyse(file_list[4],Instrument_list[4], 'cyan' , 'Xshooter 2016 298-556nm', CaII_k_vacuum,rad_vel[4],unit_correct[4])
reg5extflux = reg5flux.flux / ext.extinguish(reg5wav.wavelength,Av = 1.0)



plt.figure()
#plt.plot(reg1wav.velocity,reg1flux.flux,'k',label='Uves 495-707nm')
plt.plot(reg2wav.velocity,reg2flux.flux,'r',label='Xshooter 533-1020nm 2016')
plt.plot(reg3wav.velocity,reg3flux.flux,'m',label='Xshooter 533-1020nm 2015')
plt.xlabel('Velocity km$s^{-1}',fontsize = 20)
plt.xticks(size = 20)
plt.ylabel('flux erg $cm^{-2} s^{-1}angstrom^{-1}$',fontsize = 20)
plt.yticks(size = 20)
plt.title(' velocity centered at ' + str(Halpha_vacuum), fontsize = 25)
plt.legend()
plt.grid()


plt.figure()
plt.plot(reg4wav.velocity,reg4flux.flux,'yellow',label='UVES 2010  373-500nm')
plt.plot(reg5wav.velocity,reg5flux.flux,'cyan',label='Xshooter 2016 298-556nm')
plt.xlabel('Velocity km$s^{-1}$',fontsize = 20)
plt.ylabel('flux erg $cm^{-2} s^{-1}angstrom^{-1}$',fontsize = 20)
plt.title(' velocity centered at ' + str(CaII_k_vacuum),fontsize = 25)
plt.xticks(size = 20)
plt.yticks(size = 20)
plt.legend()
plt.grid()


# plot in wavelength 
plt.figure()
#plt.plot(reg1wav.wavelength,reg1flux.flux,'k',linestyle = '--',label='Uves 495-707nm')
plt.plot(reg2wav.wavelength,reg2flux.flux,'r',linestyle = '--',label='Xshooter 533-1020nm 2016')
plt.plot(reg3wav.wavelength,reg3flux.flux,'m',linestyle = '--',label='Xshooter 533-1020nm 2015')

#plt.plot(reg1wav.wavelength,reg1extflux,'k',label='Uves 495-707nm extinction corrected')
plt.plot(reg2wav.wavelength,reg2extflux,'r',label='Xshooter 533-1020nm 2016 extinction corrected')
plt.plot(reg3wav.wavelength,reg3extflux,'m',label='Xshooter 533-1020nm 2015 extinction corrected')
plt.xlabel('wavelength (angstroms)',fontsize = 20)
plt.ylabel('flux erg $cm^{-2} s^{-1}angstrom^{-1}$',fontsize = 20)
plt.title(' extinction corrected wavelength centered at ' + str(Halpha_vacuum),fontsize = 25)
plt.xticks(size = 20)
plt.yticks(size = 20)
plt.legend(fontsize = 20)
plt.grid()

plt.figure()
#plt.plot(reg1wav.wavelength,reg1flux.flux,'k',linestyle = '--',label='Uves 495-707nm')
plt.plot(reg4wav.wavelength,reg4flux.flux,'green',linestyle = '--',label= 'UVES 2010  373-500nm')
plt.plot(reg5wav.wavelength,reg5flux.flux,'blue',linestyle = '--',label= 'Xshooter 2016 298-556nm')

#plt.plot(reg1wav.wavelength,reg1extflux,'k',label='Uves 495-707nm extinction corrected')
plt.plot(reg4wav.wavelength,reg4extflux,'green',label=' UVES 2010  373-500nm extinction corrected')
plt.plot(reg5wav.wavelength,reg5extflux,'blue',label='Xshooter 2016 298-556nm extinction corrected')
plt.xlabel('wavelength (angstrom)',fontsize = 20)
plt.ylabel('flux erg $cm^{-2} s^{-1}angstrom^{-1}$',fontsize = 20.)
plt.title(' extinction corrected wavelength centered at ' + str(CaII_k_vacuum),fontsize = 25)
plt.xticks(size = 20)
plt.yticks(size = 20)
plt.legend(fontsize = 20)
plt.grid()

###############################################################################


###############################################################################
#find the flux with numerical integration of the gaussian

#gtot_simp = simpson(reg1flux.flux,reg1wav.wavelength)

print('The flux of the Uves 495-707nm H alpha line from the simpsons rule integration is: ')
#print(gtot_simp)
#gtot_simpext = simpson(reg1extflux,reg1wav.wavelength)
print('The extinction corrected flux of the Uves 495-707nm H alpha line from the simpsons rule integration is: ')
#print(gtot_simpext)
print('\n')

gtot_simp2 = simpson(reg2flux.flux,reg2wav.wavelength)
print('The flux of Xshooter 533-1020nm 2016 H alpha line from the simpsons rule integration is: ')
print(gtot_simp2)
gtot_simpext2 = simpson(reg2extflux,reg2wav.wavelength)
print('The extinction corrected flux of the Xshooter 533-1020nm 2016 H alpha line from the simpsons rule integration is: ')
print(gtot_simpext2)
print('\n')

gtot_simp3 = simpson(reg3flux.flux,reg3wav.wavelength)
print('The flux of the Xshooter 533-1020nm 2015 H alpha line from the simpsons rule integration is: ')
print(gtot_simp3)
gtot_simpext3 = simpson(reg3extflux,reg3wav.wavelength)
print('The extinction corrected flux of the Xshooter 533-1020nm 2015 H alpha line from the simpsons rule integration is: ')
print(gtot_simpext3)
print('\n')
#------------------------------------------------------------------------------
gtot_simp4 = simpson(reg4flux.flux,reg4wav.wavelength)
print('The flux of the Uves 373-500nm CAII line from the simpsons rule integration is: ')
print(gtot_simp4)
gtot_simpext4 = simpson(reg4extflux,reg4wav.wavelength)
print('The extinction corrected flux of the Uves 373-500nm CAII line from the simpsons rule integration is: ')
print(gtot_simpext4)
print('\n')

gtot_simp5 = simpson(reg5flux.flux,reg5wav.wavelength)
print('The flux of the Xshooter 2016 298-556nm CAII line from the simpsons rule integration is: ')
print(gtot_simp5)
gtot_simpext5 = simpson(reg5extflux,reg5wav.wavelength)
print('The extinction corrected flux of the Xshooter 2016 298-556nm line from the simpsons rule integration is: ')
print(gtot_simpext5)
print('\n')



gtot_simpext2 = gtot_simpext2 #* 1e-14 #2016 xshoo 533-1020 
gtot_simpext3 = gtot_simpext3 #* 1e-16 # 2015 533-1020 xshoo

gtot_simpext4 = gtot_simpext4  #* 1e-16
gtot_simpext5 = gtot_simpext5 #* 1e-14 # 298-5556 2015 xshoo

print('Ha',Macc(gtot_simpext2,183,1.13,1.74,1.5,1.8))
print('\n')
print('Ha',Macc(gtot_simpext3,183,1.13,1.74,1.5,1.8))
print('\n')

print('CAII UVES',Macc(gtot_simpext4,183,1.13,1.74,1.5,1.8))
print('\n')
print('CAII XSHOOTER',Macc(gtot_simpext5,183,1.13,1.74,1.5,1.8))
print('\n')

###############################################################################







    


