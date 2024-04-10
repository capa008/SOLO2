#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 14:58:08 2024

@author: cperezal
"""
import numpy as np
from sunpy.time import parse_time
import pandas as pd
from pylab import *
from scipy import *
from scipy.optimize import curve_fit

ICME_ti = '2022-09-06 10:01'
ME_ti = '2022-09-07 06:47'
ME_tf = '2022-09-08 04:10'

punto_inicial = parse_time(ME_ti).plot_date
punto_final = parse_time(ME_tf).plot_date

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


num_points = len(interval)
xdata = linspace(1.1, 10.1, num_points) 
ydata = interval['B']

# logarithmic function
def func(x, p1,p2):
  return p1*np.log(x)+p2


#def func(x, p1, p2):
#    return p1 * x**(p2) 

popt, pcov = curve_fit(func, xdata, ydata)

# curve params
p1 = popt[0]
p2 = popt[1]

subplot(2, 1, 1)
plot(xdata, ydata)

subplot(2, 1, 2)
# plot curve
curvex=xdata#np.linspace(15,85,1000)
curvey=func(curvex,p1,p2)
plot(curvex,curvey,'r')

# Source: https://stats.stackexchange.com/questions/190107/curve-fit-with-logarithmic-regression-in-python

