# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 17:11:39 2023

@author: 35385
"""
from scipy.fftpack import fft
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import julian
from datetime import datetime
import math

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
plt.figure()
plt.plot(Bdates,Bmags,'bx',label = 'B filter')
plt.plot(Vdates,Vmags,'gx',label = 'V filter')
plt.plot(jdV,owndataV,'rx',label ='Own V data')
plt.plot(jdB,owndataB,'ro',label = 'Own B data')


plt.gca().invert_yaxis()
plt.legend()
plt.title('RW Aurigae lightcurve AAVSO Data',fontsize = 25)
plt.xlabel('Julian Date',fontsize = 20)
plt.xticks(size = 20)
plt.ylabel('Magnitude',fontsize = 20)
plt.yticks(size = 16)
plt.grid(True)
#


#

V = jd_to_date(jdV)
B  =jd_to_date(jdB)

datV = datetime(V[0],V[1],V[2])
datB  =datetime(B[0],B[1],B[2])
plt.figure()
plt.plot_date(rdateB,Bmags,'bx',label = 'B filter')
plt.plot_date(rdateV,Vmags,'gx',label = 'V filter')
plt.scatter(datV,owndataV,c = 'orange',marker = 'o',s = 300,zorder = 1000,label = 'Calculated V')
plt.scatter(datB,owndataB,c = 'r',marker = 'o',s = 300,zorder = 1000,label = 'Calculated B')
#plt.xticks(size = 45,rotation = 85)


plt.gca().invert_yaxis()
plt.legend(fontsize = 20)
plt.title('RW Aurigae light curve AAVSO Data',fontsize = 25)
plt.xlabel('Year',fontsize = 20)
plt.xticks(size = 20)
plt.ylabel('Magnitude',fontsize = 20)
plt.yticks(size = 20)
plt.grid(True)












