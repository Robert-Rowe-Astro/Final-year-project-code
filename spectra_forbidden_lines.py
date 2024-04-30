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

#Halpha_vacuum = 6564.614 * u.AA # ref https://classic.sdss.org/dr7/products/spectra/vacwavelength.php

OI_vacuum = 6365.536 * u.AA # REF https://classic.sdss.org/dr6/algorithms/linestable.php

CaII_k_vacuum = 3933.66 * u.AA # ref (for now) wikipedia frauenhofer liines



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

###############################################################################
quantity_support()

ext = G23(Rv = 3.1)

#******************************************************************************

# gaussian fitting and line fluxes reference:
# ref https://lukeholden.com/blog/measuring_emission_line_fluxes_python.html

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
        
       
        
    #plot the continuum subtraction
    #plt.figure()
    #plt.plot(spec1.wavelength,yfit)
    """
    plt.figure()
    plt.plot(spec1.velocity,spec1sub.flux,color = plot_colour)
    #plt.plot(spec1.velocity,yfit)
    plt.title(plot_title)
    plt.xlabel('Velocity km/s')
    plt.ylabel(r"$\frac{10^-16 erg}{cm^2 s angstrom}$")
    #plt.xlim(lam_vacuum - 10* u .AA,lam_vacuum + 10* u.AA)
    hdu_list.close()
    
    plt.figure()
    plt.plot(line_regionwav.velocity,line_regionflux.flux,color = plot_colour)
    plt.grid()
    plt.title(plot_title + ' centered at ' + str(lam_vacuum) ) 
    """
    
   


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


s1wav, s1flux, reg1wav,reg1flux, regext1 = Spectra_Analyse(file_list[3],Instrument_list[3],'green','uves 373-500nm ',SII406, rad_vel[3])


s2wav, s2flux, reg2wav,reg2flux, regext2 = Spectra_Analyse(file_list[4],Instrument_list[4],'cyan','Xshooter 298-556',SII406, rad_vel[4])





s4wav, s5flux, reg4wav,reg4flux, regext4 = Spectra_Analyse(file_list[1],Instrument_list[1], 'red' , 'Xshooter 2016 533-1020nm', OI630,rad_vel[1])


s5wav, sflux, reg5wav,reg5flux, regext5 = Spectra_Analyse(file_list[2],Instrument_list[2], 'm' , 'Xshooter 2015 533-1020', OI630,rad_vel[2])

###############################################################################

f1wav,f1flux, freg1wav, freg1flux, fregext1 = Spectra_Analyse(file_list[1],Instrument_list[1],'r', 'xshooter', OI557 , rad_vel[0])

##############################################################################
#f4wav,f4flux, freg4wav, freg4flux = Spectra_Analyse(file_list[0],Instrument_list[0],'blue', 'uves 495-707 2008', SII673 , rad_vel[0],unit_scale[1])

f5wav,f5flux, freg5wav, freg5flux, fregext5 = Spectra_Analyse(file_list[1],Instrument_list[1],'r', 'Xshooter 533-1020nm 2016', SII673 , rad_vel[1])

f6wav,f6flux, freg6wav, freg6flux, fregext6 = Spectra_Analyse(file_list[2],Instrument_list[2],'m', 'Xshooter 533-1020nm 2015', SII673 , rad_vel[2])
###############################################################################

"""

plt.figure()
plt.plot(reg1wav.velocity,reg1flux.flux,'green',linestyle = '--',label='Uves 373-500nm')
plt.plot(reg1wav.velocity,regext1,'green',label = 'Uves 373-500nm extinction corrected')

plt.plot(reg2wav.velocity,reg2flux.flux,'cyan',linestyle = '--',label='Xshooter 298-556 2016')
plt.plot(reg2wav.velocity,regext2,'cyan',label = 'Xshooter 298-556 2016 extinction corrected')
#plt.plot(reg3wav.velocity,reg3flux.flux,'m',label='Xshooter 533-1020nm 2015')
plt.xlabel('Velocity km$s^{-1}$',fontsize = 20)
plt.ylabel('flux erg $cm^{-2} s^{-1}angstrom^{-1}$',fontsize = 20)
plt.title('Spectra velocity centered at ' + str(SII406),fontsize = 25)
plt.xticks(size = 20)
plt.yticks(size = 20)
plt.legend(fontsize = 20)
plt.grid()


plt.figure()
plt.plot(reg4wav.velocity,reg4flux.flux,'r',linestyle = '--',label='Xshooter 2016 522-1020')
plt.plot(reg4wav.velocity,regext4,'r',label = 'Xshooter 533-10202nm 2016 extinction corrected')

plt.plot(reg5wav.velocity,reg5flux.flux,'m',linestyle = '--',label='Xshooter 2015 522-1020nm')
plt.plot(reg5wav.velocity,regext5,'m',label = 'Xshooter 533-10202nm 2015 extinction corrected')
plt.xlabel('Velocity  km$s^{-1}$',fontsize = 20)
plt.ylabel('flux erg $cm^{-2} s^{-1}angstrom^{-1}$',fontsize = 20)
plt.title('Spectra velocity centered at ' + str(OI630),fontsize = 25)
plt.xticks(size = 20)
plt.yticks(size = 20)
plt.legend(fontsize = 15)
plt.grid()


plt.figure()
plt.plot(freg1wav.velocity,freg1flux.flux,'r',linestyle = '--',label='Xshooter 2016 522-1020nm')
plt.plot(freg1wav.velocity,fregext1,'r',label = 'Xshooter 533-1020nm 2016 extinction corrected')
#plt.plot(freg2wav.velocity,freg2flux.flux,'pink',label='Xshooter 533-1020nm 2016')
#plt.plot(freg3wav.velocity,freg3flux.flux,'orange',label='Xshooter 533-1020nm 2015')
plt.xlabel('Velocity  km$s^{-1}$',fontsize = 20)
plt.ylabel('flux erg $cm^{-2} s^{-1}angstrom^{-1}$',fontsize = 20)
plt.title('Spectra velocity centered at ' + str( OI557),fontsize = 25)
plt.legend(fontsize = 20)
plt.xticks(size = 20)
plt.yticks(size = 20)
plt.grid()


plt.figure()
#plt.plot(freg4wav.velocity,freg4flux.flux,'blue',label='Uves 495-707nm')
plt.plot(freg5wav.velocity,freg5flux.flux,'r',linestyle = '--',label='Xshooter 533-1020nm 2016')
plt.plot(freg5wav.velocity,fregext5,'r',label = 'Xshooter 533-1020nm 2016 extinction corrected')

plt.plot(freg6wav.velocity,freg6flux.flux,'m',linestyle = '--',label='Xshooter 533-1020nm 2015')
plt.plot(freg6wav.velocity,fregext6,'m',label = 'Xshooter 533-1020nm 2015 extinction corrected')
plt.xlabel('Velocity km$s^{-1}$',fontsize = 20)
plt.ylabel('flux  erg $cm^{-2} s^{-1}angstrom^{-1}$',fontsize = 20)
plt.title('Spectra  velocity centered at ' + str( SII673),fontsize = 25)
plt.xticks(size = 20)
plt.yticks(size = 20)
plt.legend(fontsize = 20)
plt.grid()
"""

# plot in wavelength 


###############################################################################




#-----------------------------------------------------------------------------#
# Gaussian 1


gauss1init = models.Gaussian1D(amplitude = 6.02e-13, mean = 201.7, stddev = 100 )



gtotinit = gauss1init #  + gauss2init + gauss3init




fitgauss = fitting.LevMarLSQFitter()
gtot = fitgauss(gtotinit,reg1wav.velocity,reg1flux.flux)
#-----------------------------------------------------------------------------#

# Gaussian 2


gauss1init2 = models.Gaussian1D(amplitude = 1.47e-13, mean = 256.3, stddev = 130 )



gtot2init = gauss1init2  # + gauss2init2 # + gauss3init

fitgauss2 = fitting.LevMarLSQFitter()
gtot2 = fitgauss2(gtot2init,reg2wav.velocity,reg2flux.flux)
#-----------------------------------------------------------------------------#

# Gaussian 4



#gauss1init4 = models.Gaussian1D(amplitude = 2.84e-13, mean = 106.77 , stddev = 50)

#gauss2init3 = models.Gaussian1D(amplitude = 0.1182, mean = 106.77, stddev = 50 )
#gauss3init3 = models.Gaussian1D(amplitude = 2.1, mean = -264, stddev = 100 )

#gtot4init = gauss1init4 # +  gauss3init3


#fitgauss4 = fitting.LevMarLSQFitter()
#gtot4 = fitgauss(gtot4init,freg4wav.velocity,freg4flux.flux)
#-----------------------------------------------------------------------------#
# Gaussian 5



#gauss1init5 = models.Gaussian1D(amplitude = 1.15e-14, mean =  86.8, stddev = 25)



#gtot5init = gauss1init5 # +  gauss3init3


#fitgauss5 = fitting.LevMarLSQFitter()
#gtot5 = fitgauss(gtot5init,freg5wav.velocity,freg5flux.flux)

#-----------------------------------------------------------------------------#

# Gaussian 5



fgauss1init5 = models.Gaussian1D(amplitude = 6.22e-14, mean =  54.6, stddev = 80)



fgtot1init5 = fgauss1init5 # +  gauss3init3


ffitgauss5 = fitting.LevMarLSQFitter()
fgtot5 = ffitgauss5(fgtot1init5,freg5wav.velocity,freg5flux.flux)
#-----------------------------------------------------------------------------#

# Gaussian 5



fgauss1init6 = models.Gaussian1D(amplitude = 2.331e-14, mean =  86.8, stddev = 40)



fgtot1init6 = fgauss1init6 # +  gauss3init3


ffitgauss6 = fitting.LevMarLSQFitter()
fgtot6 = ffitgauss6(fgtot1init6,freg6wav.velocity,freg6flux.flux)
#-----------------------------------------------------------------------------#

# Gaussian 5

#gauss1init4 = models.Gaussian1D(amplitude = 1.21e-13, mean =  -105.41, stddev = 145)#150 #145
#gauss2init4 = models.Gaussian1D(amplitude = 1.66e-13, mean =  -30, stddev = 15) #20 #10 #10 #15
#gauss3init4 = models.Gaussian1D(amplitude = 1.38e-13, mean =  75.45, stddev = 10) #70 # 100 #120 #150 # 5 #10
#gauss4init4 = models.Gaussian1D(amplitude = 1.24e-13, mean =  104.04, stddev = 15)#40 # 60 #30   # 15 #15

gauss1init4 = models.Gaussian1D(amplitude = 1.655e-13, mean =  -38.73, stddev = 20)#150 #145
gauss2init4 = models.Gaussian1D(amplitude = 1.655e-13, mean =  -19.8, stddev = 15) #20 #10 #10 #15
gauss3init4 = models.Gaussian1D(amplitude = 1.316e-13, mean =  74.45, stddev = 7) #70 # 100 #120 #150 # 5 #10
gauss4init4 = models.Gaussian1D(amplitude = 1.213e-13, mean =  106.2, stddev = 10)#40 # 60 #30   # 15 #15



#gtot1init4 = gauss1init4  +  gauss2init4 + gauss3init4 + gauss4init4

gtot1init4 = gauss1init4  +  gauss2init4 + gauss3init4 + gauss4init4

fitgauss4 = fitting.LevMarLSQFitter()
gtot4 = fitgauss4(gtot1init4,reg4wav.velocity,reg4flux.flux)



#-----------------------------------------------------------------------------#
plt.figure()
plt.plot(reg1wav.velocity,reg1flux.flux,'green',label='Uves 376-500nm')
plt.plot(reg1wav.velocity, gtot(reg1wav.velocity),'grey',linestyle = '--',label = 'gaussian fit 1')

plt.plot(reg2wav.velocity,reg2flux.flux,'cyan',label='Xshooter 298-556nm 2015')
plt.plot(reg2wav.velocity, gtot2(reg2wav.velocity),'orange',linestyle = '--',label = 'gaussian fit 2')

plt.legend()
plt.grid()

#*****************************************************************************#

###############################################################################


###############################################################################





plt.figure()
plt.plot(freg5wav.velocity,freg5flux.flux,'r',label='Xshooter 533-1020nm 2016')
plt.plot(freg5wav.velocity, fgtot5(freg5wav.velocity),'blue',linestyle = '--',label = 'gaussian fit 3')


plt.plot(freg6wav.velocity,freg6flux.flux,'m',label='Xshooter 533-1020nm 2015')
plt.plot(freg6wav.velocity, fgtot6(freg6wav.velocity),'lightblue',linestyle = '--',label = 'gaussian fit 3')
plt.legend()
plt.grid()

#******************************************************************************#

plt.figure()
plt.plot(reg4wav.velocity,reg4flux.flux,'r',label='Xshooter 533-1020nm 2016')
plt.plot(reg4wav.velocity, gtot4(reg4wav.velocity),'blue',linestyle = '--',label = 'gaussian fit 4')
#gtot_sep = norm.pdf(reg1wav.velocity,-32.9, 30.65)
#plt.plot(reg1wav.velocity,gtot_sep,alpha = 0.6 )

#gtot_sep2 = norm.pdf(reg1wav.velocity,-11.86, 7.57)
#plt.plot(reg1wav.velocity,gtot_sep2,alpha = 0.6 )

#gtot_sep3 = norm.pdf(reg1wav.velocity,46.1, 1.17e-38)
#plt.plot(reg1wav.velocity,gtot_sep3,alpha = 0.6 )

#gtot_sep4 = norm.pdf(reg1wav.velocity,86.99, 62.1)
#plt.plot(reg1wav.velocity,gtot_sep4,alpha = 0.6 )
plt.legend()
plt.grid()





#*****************************************************************************#

print('gaussian  fit 1 paramters: ')
print(gtot.parameters)
print('\n')

print('gaussian  fit 2 paramters: ')
print(gtot2.parameters)
print('\n')


print('gaussian  fit 3 paramters: ')
print(fgtot5.parameters)
print('\n')

print('gaussian fit 4 paramters: ')
print(fgtot6.parameters)
print('\n')

print('oi gaussian fit  paramters: ')
print(gtot4.parameters)
print('\n')

###############################################################################
#plt.figure()
#plt.plot(freg4wav.velocity,freg4flux.flux,'blue',label='Uves 495-707nm')
#plt.plot(freg4wav.velocity, gtot4(freg4wav.velocity),'purple',linestyle = '--',label = 'gaussian fit 1')

#plt.plot(freg6wav.velocity,freg6flux.flux,'r',label='Xshooter 533-1020nm 2015')
#plt.plot(freg6wav.velocity, gtot6(freg6wav.velocity),'cyan',linestyle = '--',label = 'gaussian fit 2')

#plt.legend()
#plt.grid()


###############################################################################


###############################################################################










    


