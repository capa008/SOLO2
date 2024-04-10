#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 13:36:58 2024

@author: cperezal
"""

from sunpy.time import parse_time
import pandas as pd
from pylab import *
from scipy import *

# Define function for calculating a power law
powerlaw = lambda x, amp, index: amp * (x**index)


####################################################
urMAG1hr = '/Users/cperezal/Documents/SoLo/data/SOLO_COHO1HR_MERGED_MAG_PLASMA_1699652.txt'
icMAG1hr = pd.read_table(urMAG1hr, delim_whitespace=True, skiprows=62)
icMAG1hr.columns = ['year', 'time', 'R', 'B', 'Vp', 'Np', 'Tp']

dateMAG1hr = ["" for x in range(len(icMAG1hr))]
for x in range(len(icMAG1hr)):
    dateMAG1hr[x] = icMAG1hr['year'][x][6:10]+'-'+icMAG1hr['year'][x][3:5]+'-'+icMAG1hr['year'][x][0:2]+' '+icMAG1hr['time'][x][0:8]
  
srDateMAG1hr = pd.Series(dateMAG1hr)
icMAG1hr = icMAG1hr.assign(srDateMAG1hr=srDateMAG1hr.values)

ICME_ti = '2022-09-06 10:01'
ME_ti = '2022-09-07 06:47'
ME_tf = '2022-09-08 04:10'

punto_inicial = parse_time(ME_ti).plot_date
punto_final = parse_time(ME_tf).plot_date

DateMAGVal1hr = parse_time(srDateMAG1hr).plot_date
srDateMAGVal1hr = pd.Series(DateMAGVal1hr)
icMAG1hr = icMAG1hr.assign(srDateMAGVal1hr=srDateMAGVal1hr.values)
interval2 = icMAG1hr[icMAG1hr.srDateMAGVal1hr.ge(punto_inicial ) & icMAG1hr.srDateMAGVal1hr.le(punto_final)]
    

urMAG = '/Users/cperezal/Documents/SoLo/data/SOLO_L2_MAG-RTN-NORMAL-1-MINUTE_2613226.txt'
#urMAG = '/home/arturo/Documentos/SOLO/data/SOLO_L2_MAG-RTN-NORMAL-1-MINUTE_2613226.txt'
icMAG = pd.read_table(urMAG, delim_whitespace=True, skiprows=67)
icMAG.columns = ['year', 'time', 'Br', 'Bt', 'Bn']
Bar = np.arange(len(icMAG), dtype=float)
dateMAG = ["" for x in range(len(icMAG))]
for x in range(len(icMAG)):
    dateMAG[x] = icMAG['year'][x][6:10]+'-'+icMAG['year'][x][3:5]+'-'+icMAG['year'][x][0:2]+' '+icMAG['time'][x][0:8]
    Bar[x] =  np.sqrt(icMAG['Br'][x]**2 + icMAG['Bt'][x]**2 + icMAG['Bn'][x]**2)
    
srDateMAG = pd.Series(dateMAG)
icMAG = icMAG.assign(srDateMAG=srDateMAG.values)
B = pd.Series(Bar)
icMAG = icMAG.assign(B=B.values)    
DateMAGVal = parse_time(srDateMAG).plot_date
srDateMAGVal = pd.Series(DateMAGVal)
icMAG = icMAG.assign(srDateMAGVal=srDateMAGVal.values)
interval = icMAG[icMAG.srDateMAGVal.ge(punto_inicial ) & icMAG.srDateMAGVal.le(punto_final)]


###################################################



##########
# Generate data points with noise
##########
num_points = len(interval)

# Note: all positive, non-zero data
xdata = linspace(1.1, 10.1, num_points) 
ydata = interval['B'] #powerlaw(xdata, 10.0, -2.0)     # simulated perfect data
yerr = 0.2 * ydata                      # simulated errors (10%)

#ydata += randn(num_points) * yerr       # simulated noisy data

##########
# Fitting the data -- Least Squares Method
##########

# Power-law fitting is best done by first converting
# to a linear equation and then fitting to a straight line.
#
#  y = a * x^b
#  log(y) = log(a) + b*log(x)
#

logx = log10(xdata)
logy = log10(ydata)
logyerr = yerr / ydata

# define our (line) fitting function
fitfunc = lambda p, x: p[0] + p[1] * x   
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

pinit = [1.0, -1.0]
out = optimize.leastsq(errfunc, pinit,
                       args=(logx, logy, logyerr), full_output=1)

pfinal = out[0]
covar = out[1]
print(pfinal)
print(covar)

index = pfinal[1]
amp = 10.0**pfinal[0]

indexErr = sqrt( covar[0][0] ) 
ampErr = sqrt( covar[1][1] ) * amp

##########
# Plotting data
##########

clf()
subplot(2, 1, 1)
plot(xdata, powerlaw(xdata, amp, index))     # Fit
#errorbar(xdata, ydata, yerr=yerr, fmt='k.')  # Data
#text(5, 6.5, 'Ampli = %5.2f +/- %5.2f' % (amp, ampErr))
#text(5, 5.5, 'Index = %5.2f +/- %5.2f' % (index, indexErr))
title('Best Fit Power Law')
xlabel('X')
ylabel('Y')
#xlim(1, 11)

subplot(2, 1, 2)
loglog(xdata, powerlaw(xdata, amp, index))
#errorbar(xdata, ydata, yerr=yerr, fmt='k.')  # Data
xlabel('X (log scale)')
ylabel('Y (log scale)')
#xlim(1.0, 11)



#source: http://scipy.github.io/old-wiki/pages/Cookbook/FittingData#:~:text=simulated%20noisy%20data-,Fitting%20the%20data,positional%20uncertainties%20during%20the%20fit.

