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




###############################################################################

Halpha_vacuum = 6564.614 * u.AA # ref https://classic.sdss.org/dr7/products/spectra/vacwavelength.php

OI_vacuum = 6365.536 * u.AA # REF https://classic.sdss.org/dr6/algorithms/linestable.php

CaII_k_vacuum = 3933.66 * u.AA # ref (for now) wikipedia frauenhofer liines

OI = 6300 * u.AA # get ref

SII = 6731 * u.AA # get ref

###############################################################################
quantity_support()

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
    
    # x shooter wavlength is in nanometres by default
    if(instrument != 'UVES'):
        wavelength1 = wavelength1 * 10
        
        
    wavelength1 = wavelength1 * u.AA
    flux1 = spectra.FLUX_REDUCED[0] * u.Unit('erg cm**-2 s**-1 nm**-1')


    

    #create spectrum1d object 
    spec1 = Spectrum1D(spectral_axis = wavelength1,flux = flux1,radial_velocity = v_correct * u.Unit('km/s'),velocity_convention ='optical',rest_value = lam_vacuum)
    
    #continuum subtract the spectrum1d object
    fit = fit_generic_continuum(spec1,exclude_regions=[SpectralRegion((lam_vacuum - 10* u.AA) ,(lam_vacuum + 10* u.AA))])
    yfit = fit(spec1.wavelength)
    spec1sub = (spec1/yfit)-1
    
    
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
    
    return spec1,spec1sub, line_regionwav,line_regionflux


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
#******************************************************************************



# initialise extinction profile
ext = G23(Rv = 3.1)


s1wav, s1flux, reg1wav,reg1flux = Spectra_Analyse(file_list[0],Instrument_list[0],'k','uves 495-707nm 2008',Halpha_vacuum,rad_vel[0])
#reg1flux.new_flux_unit(u.Unit('erg cm**-2 s**-1 nm**-1'))
reg1extflux = reg1flux.flux / ext.extinguish(reg1wav.wavelength,Av = 1.0)

s2wav, s2flux, reg2wav,reg2flux = Spectra_Analyse(file_list[1],Instrument_list[1],'r','Xshooter 533-1020nm 2016',Halpha_vacuum,rad_vel[1])
#reg2flux.new_flux_unit(u.Unit('erg cm**-2 s**-1 nm**-1'))
reg2extflux = reg2flux.flux / ext.extinguish(reg2wav.wavelength,Av = 1.0)


s3wav, s3flux, reg3wav,reg3flux = Spectra_Analyse(file_list[2],Instrument_list[2],'m','Xshooter 533-1020nm 2015',Halpha_vacuum,rad_vel[2])
#reg3flux.new_flux_unit(u.Unit('erg cm**-2 s**-1 nm**-1'))
reg3extflux = reg3flux.flux / ext.extinguish(reg3wav.wavelength,Av = 1.0)

s4wav, s5flux, reg4wav,reg4flux = Spectra_Analyse(file_list[3],Instrument_list[3], 'yellow' , 'UVES 2010  373-500nm', CaII_k_vacuum,rad_vel[3])

s5wav, sflux, reg5wav,reg5flux = Spectra_Analyse(file_list[4],Instrument_list[4], 'cyan' , 'Xshooter 2016 298-556nm', CaII_k_vacuum,rad_vel[4])
###############################################################################

f1wav,f1flux, freg1wav, freg1flux = Spectra_Analyse(file_list[0],Instrument_list[0],'green', 'uves 495-707 2008', OI , rad_vel[0])

f2wav,f2flux, freg2wav, freg2flux = Spectra_Analyse(file_list[1],Instrument_list[1],'pink', 'Xshooter 533-1020nm 2016', OI , rad_vel[1])

f3wav,f3flux, freg3wav, freg3flux = Spectra_Analyse(file_list[2],Instrument_list[2],'orange', 'Xshooter 533-1020nm', OI , rad_vel[2])

##############################################################################
f4wav,f4flux, freg4wav, freg4flux = Spectra_Analyse(file_list[0],Instrument_list[0],'blue', 'uves 495-707 2008', SII , rad_vel[0])

f5wav,f5flux, freg5wav, freg5flux = Spectra_Analyse(file_list[1],Instrument_list[1],'purple', 'Xshooter 533-1020nm 2016', SII , rad_vel[1])

f6wav,f6flux, freg6wav, freg6flux = Spectra_Analyse(file_list[2],Instrument_list[2],'r', 'Xshooter 533-1020nm', SII , rad_vel[2])
###############################################################################



plt.figure()
plt.plot(reg1wav.velocity,reg1flux.flux,'k',label='Uves 495-707nm')
plt.plot(reg2wav.velocity,reg2flux.flux,'r',label='Xshooter 533-1020nm 2016')
plt.plot(reg3wav.velocity,reg3flux.flux,'m',label='Xshooter 533-1020nm 2015')
plt.xlabel('Velocity')
plt.ylabel('flux')
plt.title('spectra comparison velocity centered at ' + str(Halpha_vacuum))
plt.legend()
plt.grid()


plt.figure()
plt.plot(reg4wav.velocity,reg4flux.flux,'yellow',label='UVES 2010  373-500nm')
plt.plot(reg5wav.velocity,reg5flux.flux,'cyan',label='Xshooter 2016 298-556nm')
plt.xlabel('Velocity')
plt.ylabel('flux')
plt.title('spectra comparison velocity centered at ' + str(CaII_k_vacuum))
plt.legend()
plt.grid()


plt.figure()
plt.plot(freg1wav.velocity,freg1flux.flux,'green',label='Uves 495-707nm')
plt.plot(freg2wav.velocity,freg2flux.flux,'pink',label='Xshooter 533-1020nm 2016')
plt.plot(freg3wav.velocity,freg3flux.flux,'orange',label='Xshooter 533-1020nm 2015')
plt.xlabel('Velocity')
plt.ylabel('flux')
plt.title('spectra comparison velocity centered at ' + str( OI))
plt.legend()
plt.grid()


plt.figure()
plt.plot(freg4wav.velocity,freg4flux.flux,'blue',label='Uves 495-707nm')
plt.plot(freg5wav.velocity,freg5flux.flux,'purple',label='Xshooter 533-1020nm 2016')
plt.plot(freg6wav.velocity,freg6flux.flux,'red',label='Xshooter 533-1020nm 2015')
plt.xlabel('Velocity')
plt.ylabel('flux')
plt.title('spectra comparison velocity centered at ' + str( SII))
plt.legend()
plt.grid()


# plot in wavelength 
plt.figure()
plt.plot(reg1wav.wavelength,reg1flux.flux,'k',linestyle = '--',label='Uves 495-707nm')
plt.plot(reg2wav.wavelength,reg2flux.flux,'r',linestyle = '--',label='Xshooter 533-1020nm 2016')
plt.plot(reg3wav.wavelength,reg3flux.flux,'m',linestyle = '--',label='Xshooter 533-1020nm 2015')

plt.plot(reg1wav.wavelength,reg1extflux,'k',label='Uves 495-707nm extinction corrected')
plt.plot(reg2wav.wavelength,reg2extflux,'r',label='Xshooter 533-1020nm 2016 extinction corrected')
plt.plot(reg3wav.wavelength,reg3extflux,'m',label='Xshooter 533-1020nm 2015 extinction corrected')
plt.xlabel('wavelength')
plt.ylabel('flux')
plt.title('spectra comparison extinction corrected wavelength centered at ' + str(Halpha_vacuum))
plt.legend()
plt.grid()

###############################################################################




#-----------------------------------------------------------------------------#
# Gaussian 1


gauss1init = models.Gaussian1D(amplitude = 6.6, mean = 108, stddev = 100 )
gauss2init = models.Gaussian1D(amplitude = 3.482, mean = -13, stddev = 70 )
gauss3init = models.Gaussian1D(amplitude = 4.359, mean = -205, stddev = 80 )

gtotinit = gauss1init + gauss2init + gauss3init

fitgauss = fitting.LevMarLSQFitter()
gtot = fitgauss(gtotinit,reg1wav.velocity,reg1flux.flux)
#-----------------------------------------------------------------------------#

# Gaussian 2


gauss1init2 = models.Gaussian1D(amplitude = 5.088, mean = 36, stddev = 130 )
gauss2init2 = models.Gaussian1D(amplitude = 2.14, mean = -310, stddev = 80 )
#gauss3init2 = models.Gaussian1D(amplitude = 4.5, mean = -215, stddev = 1 )

gtot2init = gauss1init2 + gauss2init2 # + gauss3init

fitgauss2 = fitting.LevMarLSQFitter()
gtot2 = fitgauss2(gtot2init,reg2wav.velocity,reg2flux.flux)
#-----------------------------------------------------------------------------#

# Gaussian 3 



gauss1init3 = models.Gaussian1D(amplitude = 6.98, mean = 31, stddev = 60 )
#gauss2init3 = models.Gaussian1D(amplitude = 3.5, mean = -66, stddev = 1.2 )
gauss3init3 = models.Gaussian1D(amplitude = 2.1, mean = -264, stddev = 100 )

gtot3init = gauss1init3  +  gauss3init3

fitgauss3 = fitting.LevMarLSQFitter()
gtot3 = fitgauss(gtot3init,reg3wav.velocity,reg3flux.flux)
#-----------------------------------------------------------------------------#


plt.figure()
plt.plot(reg1wav.velocity,reg1flux.flux,'k',label='Uves 495-707nm')
plt.plot(reg1wav.velocity, gtot(reg1wav.velocity),'purple',linestyle = '--',label = 'gaussian fit 1')

plt.plot(reg2wav.velocity,reg2flux.flux,'r',label='Xshooter 533-1020nm 2016')
plt.plot(reg2wav.velocity, gtot2(reg2wav.velocity),'green',linestyle = '--',label = 'gaussian fit 2')

plt.plot(reg3wav.velocity,reg3flux.flux,'m',label='Xshooter 533-1020nm 2015')
plt.plot(reg3wav.velocity, gtot3(reg3wav.velocity),'blue',linestyle = '--',label = 'gaussian fit 3')

plt.legend()
plt.grid()

plt.figure()

plt.plot(reg1wav.velocity, gtot(reg1wav.velocity),'purple',linestyle = '--',label = 'Uves 495-707nm gaussian fit 1')


plt.plot(reg2wav.velocity, gtot2(reg2wav.velocity),'green',linestyle = '--',label = 'Xshooter 533-1020nm 2016 gaussian fit 2')


plt.plot(reg3wav.velocity, gtot3(reg3wav.velocity),'blue',linestyle = '--',label = 'Xshooter 533-1020nm 2015 gaussian fit 3')

plt.legend()
plt.grid()
###############################################################################

gtot_sep1 = norm.pdf(reg1wav.velocity,gtot.parameters[2],gtot.parameters[3])
plt.figure()
plt.plot(reg1wav.velocity,reg1flux.flux,'k',label='Uves 495-707nm')
plt.plot(reg1wav.velocity, gtot(reg1wav.velocity),'red',linestyle = '--', label = 'gaussian fit 1')
gtot_sep = norm.pdf(reg1wav.velocity,gtot.parameters[2],gtot.parameters[3])
plt.plot(reg1wav.velocity,gtot_sep1,alpha = 0.6 )
plt.legend()
plt.grid()

print('\n')
print('gaussian fit 1 paramters: ')
print(gtot.parameters)
print('\n')


plt.figure()
plt.plot(reg2wav.velocity,reg2flux.flux,'r',label='Xshooter 533-1020nm 2016')
plt.plot(reg2wav.velocity, gtot2(reg2wav.velocity),'green',linestyle = '--',label = 'gaussian fit 2')
plt.legend()
plt.grid()

print('gaussian fit 2 paramters: ')
print(gtot2.parameters)
print('\n')



plt.figure()
plt.plot(reg3wav.velocity,reg3flux.flux,'m',label='Xshooter 533-1020nm 2015')
plt.plot(reg3wav.velocity, gtot3(reg3wav.velocity),'blue',linestyle = '--',label = 'gaussian fit 3')
plt.legend()
plt.grid()

print('gaussian fit 3 paramters: ')
print(gtot3.parameters)
print('\n')


###############################################################################

"""
flux1_1 = flux_calc(gtot.parameters[0],gtot.parameters[2])
flux1_2 = flux_calc(gtot.parameters[3],gtot.parameters[5])
flux1_3 = flux_calc(gtot.parameters[6],gtot.parameters[8])

flux1tot = flux1_1 + flux1_2 + flux1_3


flux2_1 = flux_calc(gtot2.parameters[0],gtot2.parameters[2])
flux2_2 = flux_calc(gtot2.parameters[3],gtot2.parameters[5])

flux2tot = flux2_1 + flux2_2


flux3_1 = flux_calc(gtot3.parameters[0],gtot3.parameters[2])
flux3_2 = flux_calc(gtot3.parameters[3],gtot3.parameters[5])

flux3tot = flux3_1 + flux3_2

print('the flux of the Uves 495-707nm H alpha line from adding the gaussians is: ')
print(flux1tot)
print('\n')

print('the flux of the Xshooter 533-1020nm 2016 H alpha line from adding the gaussians is: ')
print(flux2tot)
print('\n')

print('the flux of the Xshooter 533-1020nm 2015 H alpha line from adding the gaussians is: ')
print(flux3tot)
print('\n')
"""

###############################################################################
#find the flux with numerical integration of the gaussian

gtot_simp = simpson(reg1flux.flux,reg1wav.wavelength)

print('The flux of the Uves 495-707nm H alpha line from the simpsons rule integration is: ')
print(gtot_simp)
gtot_simpext = simpson(reg1extflux,reg1wav.wavelength)
print('The extinction corrected flux of the Uves 495-707nm H alpha line from the simpsons rule integration is: ')
print(gtot_simpext)
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

###############################################################################
# do the same with numpy trapezoidal rule for comparison

gtot_trap = trapz(reg1flux.flux,reg1wav.wavelength)
print('The flux of the Uves 495-707nm H alpha line from the trapezoidal rule integration is: ')
print(gtot_trap)
gtot_trapext = trapz(reg1extflux,reg1wav.wavelength)
print('The extinction corrected flux of the Uves 495-707nm H alpha line from the trapezoidal rule integration is: ')
print(gtot_trapext)
print('\n')

gtot_trap2 = trapz(reg2flux.flux,reg2wav.wavelength)
print('The flux of Xshooter 533-1020nm 2016 H alpha line from the trapezoidal rule integration is: ')
print(gtot_trap2)
gtot_trapext2 = trapz(reg2extflux,reg2wav.wavelength)
print('The extinction corrected flux of theXshooter 533-1020nm 2016 H alpha line from the trapezoidal rule integration is: ')
print(gtot_trapext2)
print('\n')

gtot_trap3 = trapz(reg3flux.flux,reg3wav.wavelength)
print('The flux of the Xshooter 533-1020nm 2015 H alpha line from the trapezoidal rule integration is: ')
print(gtot_trap3)
gtot_trapext3 = trapz(reg3extflux,reg3wav.wavelength)
print('The extinction corrected flux of the Xshooter H alpha line from the trapezoidal rule integration is: ')
print(gtot_trapext3)
print('\n')

###############################################################################




# paramters returns amplitude of peak,mean,width

# gtot = norm.pdf(x values, mean value, sigma)
# plot(xvalues, gtot alpha 0.6)


# 







    


