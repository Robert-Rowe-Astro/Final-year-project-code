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

    

magnitude_data = pd.read_csv("aavso_light_data.txt",sep = ",",header = None,
names =[ "JD","Magnitude","Uncertainty","HQuncertainty","Band,""Observer Code",
"Comment Code(s)","Comp Star 1","Comp Star 2","Charts","Comments","Transfomed",
"Airmass","Validation Flag","Cmag","Kmag","HJD","Star Name",
"Observer Affiliation","Measurement Method","Grouping Method","ADS Reference",
"Digitizer","Credit","extra"])


Vmags = []
Vdates = []


Bmags = []
Bdates = []


for i in range(1,len(magnitude_data)):
    #print(str((magnitude_data.iloc[i,3])))
    if(str(magnitude_data.iloc[i,4]) == 'V' ):#or str(magnitude_data.iloc[i,4]) == 'Vis.'):
        Vmags.append(float(magnitude_data.iloc[i,1]))
        Vdates.append(float(magnitude_data.iloc[i,0]))

    if(str(magnitude_data.iloc[i,4]) == 'B'):
        Bmags.append(float(magnitude_data.iloc[i,1]))
        Bdates.append(float(magnitude_data.iloc[i,0]))  
 

    
owndataV = 11.0857 # average V magnitude
uncertaintyV = 0.0499 # uncertainty
jdV = 2460234.58565

owndataB = 11.9638
uncertaintyB = 0.1502
jdB = 2460234.57955
    

freqVday = np.arange(0.01,2,0.01)
freqVday = freqVday/u.year
t = Vdates * u.day
mags = Vmags * u.mag
ls = LombScargle(t,mags)
power = ls.power(freqVday)
plt.figure()
plt.plot(freqVday,power,'g')
plt.xlabel('Frequency (1/year)')
plt.ylabel('Power/Frequency')
plt.title('Periodogram V filter')
plt.grid()




freqBday = np.arange(0.01,2,0.01)
freqBday = freqBday/u.year
tB = Bdates * u.day
magsB = Bmags * u.mag
ls2 = LombScargle(tB,magsB)
powerB = ls2.power(freqBday)
plt.figure()
plt.plot(freqBday,powerB)
plt.xlabel('Frequency (1/year)')
plt.ylabel('Power/Frequency')
plt.title('Periodogram B filter')
plt.grid()




#V
bestfreqV = freqVday[np.argmax(power)]
tfitV = np.linspace(t[0],t[-1])
ls4 = LombScargle(t,mags)
yfitV = ls4.model(tfitV,bestfreqV)
yfitV = yfitV * u.mag

"""
m = freqVday[np.argmax(power)]
for i in range(len(freqVday)):
    if( freqVday[i] == m):
        freqVday[i] = 0/u.year
        power[np.argmax(power)] = 0
"""      
#repeat
bestfreqV2 = 1.2/u.year
tfitV2 = np.linspace(t[0],t[-1])
ls42 = LombScargle(t,mags)
yfitV2 = ls42.model(tfitV2,bestfreqV2)
yfitV2 = yfitV2 * u.mag


#B
bestfreqB = freqBday[np.argmax(powerB)]
tfitB = np.linspace(tB[0],tB[-1])
ls5 = LombScargle(tB,magsB)
yfitB = ls5.model(tfitB,bestfreqB)
yfitB = yfitB * u.mag

"""
m = freqBday[np.argmax(powerB)]
for i in range(len(freqBday)):
    if( freqBday[i] == m):
        freqBday[i] = 0/u.year
        powerB[np.argmax(powerB)] = 0
 """
       
bestfreqB2 = 0.67/u.year
tfitB2 = np.linspace(tB[0],tB[-1])
ls52 = LombScargle(tB,magsB)
yfitB2 = ls52.model(tfitB2,bestfreqB2)
yfitB2 = yfitB2 * u.mag

###############################################################################
#PER DAY  to compare to paper



freqVcomp = np.arange(0.0,1.5,0.01)
freqVcomp = freqVcomp/u.day
tcomp = Vdates * u.day
magscomp = Vmags * u.mag
lscomp = LombScargle(tcomp,magscomp)
powercomp = ls.power(freqVcomp)
plt.figure()
plt.plot(freqVcomp,powercomp,'g')
plt.xlabel('Frequency (1/day)')
plt.ylabel('Power/Frequency')
plt.title('Periodogram V filter days')
plt.grid()
###############################################################################











plt.figure()
plt.plot(t,mags,'gx',label = 'V filter')
plt.plot(tfitV,yfitV,'r',label='highest power frequency')
#plt.plot(tfitV2,yfitV2,'k--',label=' second highest power frequency')
plt.gca().invert_yaxis()
plt.legend(fontsize = 20)
plt.xlabel('Julian Date',fontsize = 20)
plt.ylabel('Magnitude',fontsize = 20)
plt.title('Best fit frequency sinusoid V',fontsize = 25)
plt.xscale('log')
#plt.yscale('log')
plt.grid(True)
plt.show()


plt.figure()
plt.plot(tB,magsB,'bx',label = 'B filter')
plt.plot(tfitB,yfitB,'r' ,label='highest power frequency')
#plt.plot(tfitB2,yfitB2,'k--' ,label='second highest power frequency')
plt.gca().invert_yaxis()
plt.legend(fontsize = 20)
plt.title('Best fit frequency sinusoid B',fontsize = 25)
plt.xlabel('Julian Date',fontsize = 20)
plt.ylabel('Magnitude',fontsize = 20)
plt.xscale('log')
#plt.yscale('log')

plt.grid(True)
plt.show()




###############################################################################







