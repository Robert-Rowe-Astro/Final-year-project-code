# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 16:00:43 2024

@author: 35385
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 16:13:30 2024

@author: 35385
"""

from scipy.fftpack import fft
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import julian
import datetime
from astropy.timeseries import LombScargle
from astropy import units as u 

    

magnitude_data = pd.read_csv("aavso_light_data.txt",sep = ",",header = None, names =[ "JD","Magnitude","Uncertainty","HQuncertainty","Band,""Observer Code","Comment Code(s)","Comp Star 1","Comp Star 2","Charts","Comments","Transfomed","Airmass","Validation Flag","Cmag","Kmag","HJD","Star Name","Observer Affiliation","Measurement Method","Grouping Method","ADS Reference","Digitizer","Credit","extra"])


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
 
"""   
test = float(Bdates[0])
dt = julian.from_jd(test, fmt='mjd')
print(dt)
"""    
 
"""        
Vcalender = [] 
for i in Vdates:
    Vcalender.append(Vdates[i])  
"""


for i in range(len(Bmags)):
    Bdates[i] = float(Bdates[i])
    Bmags[i] = float(Bmags[i])
    
    
    
for i in range(len(Vmags))    :    
    Vdates[i] = float(Vdates[i])
    Vmags[i] = float(Vmags[i])
    
    
    
"""

for i in range(len(Bmags)):
    Bdates[i] = (Bdates[i]) * u.day
    Bmags[i] = (Bmags[i]) * u.mag
    
    
    
for i in range(len(Vmags))    :    
    Vdates[i] = (Vdates[i]) * u.day
    Vmags[i] = (Vmags[i])     * u.mag
"""
    
owndataV = 11.0857 # average V magnitude
uncertaintyV = 0.0499 # uncertainty
jdV = 2460234.58565

owndataB = 11.9638
uncertaintyB = 0.1502
jdB = 2460234.57955
    


#frequencyV, powerV = LombScargle(Vdates, Vmags).autopower()

#plt.plot(frequencyV,powerV)
#frequencyB, powerB = LombScargle(Bdates, Bmags).autopower()


#bestfreq = 


# VISUAL
#freqVday = np.linspace(2451484.0,2460274.858211,8426)
#powerVday = LombScargle(Vdates, Vmags).power(freqVday)
"""
freqVday = np.arange(1,7,1/7)


powerVday = LombScargle(Vdates, Vmags).power(freqVday)

#freqVday.unit
#powerVday.unit
"""

"""
freqVday = np.arange(1,7,1/7)
t = Vdates * u.day
mags = Vmags * u.mag
ls = LombScargle(t,mags)
power = ls.power(freqVday/u.week)
plt.plot(freqVday,power)
plt.xlabel('Frequency (1/d)')
plt.ylabel('Power/Frequency')
plt.title('Periodogram')
plt.grid()
"""

freqBday = np.arange(1,11,1/11)
tB = Bdates * u.day
magsB = Bmags * u.mag
ls2 = LombScargle(tB,magsB)
powerB = ls2.power(freqBday/u.year)

plt.plot(freqBday,powerB)
plt.xlabel('Frequency (1/d)')
plt.ylabel('Power/Frequency')
plt.title('Periodogram')
plt.grid()


bestfreq = freqBday[np.argmax(powerB)]
tfit = np.linspace(2452000,2460275)
ls = LombScargle(Bdates,Bmags)
yfit = ls.model(tfit,bestfreq)

"""
freqVweek = np.linspace(2451484.0,2460274.858211,1204) 
powerVweek = LombScargle(Vdates, Vmags).power(freqVweek)

freqVmonth = np.linspace(2451484.0,2460274.858211,301)
powerVmonth = LombScargle(Vdates, Vmags).power(freqVmonth)

freqVyear = np.linspace(2451484.0,2460274.858211,23)
powerVyear = LombScargle(Vdates, Vmags).power(freqVyear)



#BLUE
freqBday = np.linspace(2456251.89071,2460271.91817,4020)
powerBday = LombScargle(Bdates, Bmags).power(freqBday)

freqBweek = np.linspace(2456251.89071,2460271.91817,574)
powerBweek = LombScargle(Bdates, Bmags).power(freqBweek)

freqBmonth = np.linspace(2456251.89071,2460271.91817,144)
powerBmonth = LombScargle(Bdates, Bmags).power(freqBmonth)

freqByear = np.linspace(2456251.89071,2460271.91817,11)
powerByear = LombScargle(Bdates, Bmags).power(freqByear)
"""

"""
bestfreq = freqVday[np.argmax(powerVday)]
tfit = np.linspace(2452000,2460275)
ls = LombScargle(Vdates,Vmags)
yfit = ls.model(tfit,bestfreq)
"""
#Vdates = Vdates/(np.max(Vdates))    
#Vmags = Vmags/((np.max(Vmags)/10))
"""
# VISUAL
fig, ax =  plt.subplots()

ax1 = plt.subplot(141)
ax2 = plt.subplot(142)
ax3 = plt.subplot(143)
ax4 = plt.subplot(144)

ax1.plot(freqVday, powerVday)
ax1.set_title('daily')
ax1.set_xlabel('frequency')
ax1.set_ylabel('power')

ax2.plot(freqVweek, powerVweek,'g')
ax2.set_title('weekly')
ax2.set_xlabel('frequency')
ax2.set_ylabel('power')

ax3.plot(freqVmonth, powerVmonth,'r')
ax3.set_title('monthly')
ax3.set_xlabel('frequency')
ax3.set_ylabel('power')



ax4.plot(freqVyear, powerVyear,'y')
ax4.set_title('yearly')
ax4.set_xlabel('frequency')
ax4.set_ylabel('power')
"""



"""
# BLUE
fig, ax =  plt.subplots()

ax1 = plt.subplot(141)
ax2 = plt.subplot(142)
ax3 = plt.subplot(143)
ax4 = plt.subplot(144)

ax1.plot(freqBday, powerBday)
ax1.set_title('daily')
ax1.set_xlabel('frequency')
ax1.set_ylabel('power')

ax2.plot(freqBweek, powerBweek,'g')
ax2.set_title('weekly')
ax2.set_xlabel('frequency')
ax2.set_ylabel('power')

ax3.plot(freqBmonth, powerBmonth,'r')
ax3.set_title('monthly')
ax3.set_xlabel('frequency')
ax3.set_ylabel('power')



ax4.plot(freqByear, powerByear,'y')
ax4.set_title('yearly')
ax4.set_xlabel('frequency')
ax4.set_ylabel('power')
"""




#VISUAL
#plt.plot(freqVday,powerVday,label='dayV')
#plt.plot(freqVweek,powerVweek,label = 'weekV')
#plt.plot(freqVmonth,powerVmonth,label = 'monthV')
#plt.plot(freqVyear,powerVyear,label='yearV')



#BLUE
#plt.plot(freqBday,powerBday,label='dayB')
#plt.plot(freqBweek,powerBweek,label = 'weekB')
#plt.plot(freqBmonth,powerBmonth,label = 'monthB')
#plt.plot(freqByear,powerByear,label='yearB')
#plt.legend()

"""
invBdate = fft(Bdates)
invBmag = fft(Bmags)

invVdate = fft(Vdates)

plt.plot(invBdate,invBmag,'yx',label='fouriertransformedB')
"""


"""
#plt.plot(Bdates,Bmags,'bx',label = 'B filter')
plt.plot(Vdates,Vmags,'gx',label = 'V filter')
#plt.plot(tfit,yfit)
#plt.plot(jdV,owndataV,'rx',label ='Own V data')
#plt.plot(jdB,owndataB,'ro',label = 'Own B data')

plt.gca().invert_yaxis()
plt.legend()
plt.xlabel('Julian Date')
plt.ylabel('Magnitude')
plt.grid(True)
"""







"""
for i in range(len(frequencyB)):
    print('V\n')
    print(frequencyV[i])
    
for i in range(len(frequencyV)):
    print('B\n')
    print(frequencyB[i])    
"""








