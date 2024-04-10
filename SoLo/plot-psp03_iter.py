#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 15:14:31 2024

@author: cperezal
"""

import matplotlib.dates as dates
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np
from sunpy.time import parse_time
import pandas as pd
from pylab import *
from scipy import *
from scipy.optimize import curve_fit

# -------------------------------------------------------------------------------------------------------
#	
#	Plots Parker Solar Probe PSP
#   Fitt using power-law
#   Final version del script plots-psp03.py pero iterativa
# --------------------------------------------------------------------------------------------------------

url = '/Users/cperezal/Documents/SoLo/lists/PSP_events.txt'
dt = pd.read_table(url, header=0)

for i in range(len(dt)):
    name = '/Users/cperezal/Documents/SoLo/data/'+dt['name_file'][i]+'.txt'
    icMAG = pd.read_table(name, delim_whitespace=True, skiprows=62)
    icMAG.columns = ['year', 'time', 'Br', 'Bt', 'Bn']



#PSP_FLD_L2_MAG_RTN_1MIN_4036897.txt ----> PSP_FLD_L2_MAG_RTN_1MIN_20210926.txt
#PSP_FLD_L2_MAG_RTN_1MIN_4044937.txt ----> PSP_FLD_L2_MAG_RTN_1MIN_20220905.txt

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
    
    # ------ initial times
    # ---------------------------------------
    ICME_ti = dt['ICME_ti'][i]
    ME_ti = dt['ME_ti'][i]
    ME_tf = dt['ME_tf'][i]
    # ---------------------------------------
    
    # -------- PLOT ---------
    #
    fig = plt.figure(figsize=(8.0,9.5))
    top=0.965
    bottom=0.055
    left=0.120
    right=0.975
    hspace=0.070
    wspace=0.6

    fig.subplots_adjust(top=top,bottom=bottom,left=left,right=right,hspace=hspace, wspace=wspace)
    plt.rc('font', size=12)
    # ---------------------------------------------
    #B_avg = icMAG['B'].rolling(window=10).mean()
    ax1 = plt.subplot(7, 1,1)
    gs = gridspec.GridSpec(7, 3)
    icMAG['srDateMAG'] = pd.to_datetime(icMAG['srDateMAG'], format='%Y-%m-%d %H:%M:%S')
    plt.plot(icMAG['srDateMAG'], icMAG['B'], linewidth='0.5')
    #plt.plot(icMAG['srDateMAG'], B_avg, linewidth='0.5')
    formatter = dates.DateFormatter('%Y-%m-%d %H:%M:%S') 
    plt.title('PSP/MAG', fontsize=18)
    plt.ylabel("|B| [nT]")
    plt.axvline(pd.to_datetime(ICME_ti), color='black', linestyle='-', lw=2)
    plt.axvline(pd.to_datetime(ME_ti), color='black', linestyle='--', lw=2)
    plt.axvline(pd.to_datetime(ME_tf), color='black', linestyle='--', lw=2)
    plt.axvspan(ME_ti, ME_tf, alpha=0.2)
    ax1.axes.xaxis.set_ticklabels([])	
    plt.xlim(icMAG['srDateMAG'][0], icMAG['srDateMAG'][len(icMAG['srDateMAG'])-1])
    #-----
    ax2 = plt.subplot(7, 1,2)
    plt.plot(icMAG['srDateMAG'], icMAG['Br'], color='black', linewidth='0.5')
    plt.plot(icMAG['srDateMAG'], icMAG['Bt'], 'r-',linewidth='0.5')
    plt.plot(icMAG['srDateMAG'], icMAG['Bn'], 'b-',linewidth='0.5')
    plt.ylabel(r'$B_{r,t,n} [nT]$')
    plt.axvline(pd.to_datetime(ICME_ti), color='black', linestyle='-', lw=2)
    plt.axvline(pd.to_datetime(ME_ti), color='black', linestyle='--', lw=2)
    plt.axvline(pd.to_datetime(ME_tf), color='black', linestyle='--', lw=2)
    plt.axvspan(ME_ti, ME_tf, alpha=0.2)
    ax2.axes.xaxis.set_ticklabels([])
    plt.xlim(icMAG['srDateMAG'][0], icMAG['srDateMAG'][len(icMAG['srDateMAG'])-1])
    # ----		
    #V_avg = icSWA['V'].rolling(window=100).mean()
    ax3 = plt.subplot(7, 1,3)
    #plt.plot(icSWA['srDateSWA'], icSWA['V'], linewidth='0.5')
    #plt.plot(icSWA['srDateSWA'], V_avg, linewidth='0.5')
    plt.ylabel(r'$V_{th} [km/s]$')
    plt.axvline(pd.to_datetime(ICME_ti), color='black', linestyle='-', lw=2)
    plt.axvline(pd.to_datetime(ME_ti), color='black', linestyle='--', lw=2)
    plt.axvline(pd.to_datetime(ME_tf), color='black', linestyle='--', lw=2)
    plt.axvspan(ME_ti, ME_tf, alpha=0.2)
    ax3.axes.xaxis.set_ticklabels([])
    plt.xlim(icMAG['srDateMAG'][0], icMAG['srDateMAG'][len(icMAG['srDateMAG'])-1])
    # ----
    #Np_avg = icSWA['Np'].rolling(window=100).mean()
    ax4 = plt.subplot(7, 1,4)
    #plt.plot(icSWA['srDateSWA'], icSWA['Np'], linewidth='0.5')
    #plt.plot(icSWA['srDateSWA'], Np_avg, linewidth='0.5')
    plt.ylabel(r'$N_{p}[\#/cm^{3}]$')
    plt.axvline(pd.to_datetime(ICME_ti), color='black', linestyle='-', lw=2)
    plt.axvline(pd.to_datetime(ME_ti), color='black', linestyle='--', lw=2)
    plt.axvline(pd.to_datetime(ME_tf), color='black', linestyle='--', lw=2)
    plt.axvspan(ME_ti, ME_tf, alpha=0.2)
    ax4.axes.xaxis.set_ticklabels([])
    plt.xlim(icMAG['srDateMAG'][0], icMAG['srDateMAG'][len(icMAG['srDateMAG'])-1])
    # ----
    #Tp_avg = icSWA['Tp'].rolling(window=100).mean() 
    ax5 = plt.subplot(7, 1,5)
    #plt.plot(icSWA['srDateSWA'], icSWA['Tp'], linewidth='0.5')
    #plt.plot(icSWA['srDateSWA'], Tp_avg, linewidth='0.5')
    plt.ylabel(r'$T [K]$')
    plt.axvline(pd.to_datetime(ICME_ti), color='black', linestyle='-', lw=2)
    plt.axvline(pd.to_datetime(ME_ti), color='black', linestyle='--', lw=2)
    plt.axvline(pd.to_datetime(ME_tf), color='black', linestyle='--', lw=2)
    plt.axvspan(ME_ti, ME_tf, alpha=0.2)
    ax5.axes.xaxis.set_ticklabels([])
    plt.xlim(icMAG['srDateMAG'][0], icMAG['srDateMAG'][len(icMAG['srDateMAG'])-1])
    # ----
    #betta_avg = icMAG['betta'].rolling(window=100).mean()
    ax6 = plt.subplot(7, 1,6)
    #icMAG['srDateMAG'] = pd.to_datetime(icMAG['srDateMAG'], format='%Y-%m-%d %H:%M:%S')
    #plt.plot(icMAG['srDateMAG'], icMAG['betta'], linewidth='0.5')
    #plt.plot(icMAG['srDateMAG'], betta_avg, linewidth='0.5')
    #plt.ylim(0,3)
    plt.ylabel(r'$\beta$')
    plt.xlabel("mm-dd hh")
    plt.axvline(pd.to_datetime(ICME_ti), color='black', linestyle='-', lw=2)
    plt.axvline(pd.to_datetime(ME_ti), color='black', linestyle='--', lw=2)
    plt.axvline(pd.to_datetime(ME_tf), color='black', linestyle='--', lw=2)
    plt.axvspan(ME_ti, ME_tf, alpha=0.2)
    plt.xlim(icMAG['srDateMAG'][0], icMAG['srDateMAG'][len(icMAG['srDateMAG'])-1])
    plt.tick_params(rotation=45)
    
    #plt.savefig('/Users/cperezal/Documents/SoLo/plots/PSP/'+dt['name_file'][i][24:32]+'.png')


    # ============================================================================
    #   analysis
    # ============================================================================
    ME = icMAG[icMAG.srDateMAG.ge(pd.to_datetime(ME_ti)) & icMAG.srDateMAG.le(pd.to_datetime(ME_tf))]

    punto_inicial = parse_time(ME_ti).plot_date
    punto_final = parse_time(ME_tf).plot_date
    
    icMAG = icMAG.dropna() #para eliminar NaN y se pueda aplicar la funcion curve_fit
    
    interval = icMAG[icMAG.srDateMAGVal.ge(punto_inicial ) & icMAG.srDateMAGVal.le(punto_final)]
    
    num_points = len(interval)
    xdata = linspace(1.1, 10.1, num_points) 
    ydata = interval['B']
    
    # logarithmic function
    def func(x, p1,p2):
      return p1*np.log(x)+p2

    popt, pcov = curve_fit(func, xdata, ydata)      #aqui hay un error
    # curve params
    p1 = popt[0]        # -----> alpha
    p2 = popt[1]        # -----> slope
    
    # ============================================================================
    # Results
    #
    print('ICME: ', dt['ICME_ti'][i], 'alpha: ', p1)
    print('ME_ti: ', dt['ME_ti'][i])
    print('ME_tf: ', dt['ME_tf'][i])





