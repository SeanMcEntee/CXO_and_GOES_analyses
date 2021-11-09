# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:29:47 2021

@author: dmw1n18
"""
from astropy.io import fits as pyfits
import astropy.units as u
from astropy.time import Time                   #convert between different time coordinates
from astropy.time import TimeDelta              #add/subtract time intervals 

import matplotlib
import matplotlib.pyplot as plt

import pandas as pd
import os
import numpy as np
from math import ceil

obsIDs = ['1862', '20002']
# obsIDs = ['1862', '2519', '15669', '18608', '20000', '18678', '22146']
'''obsIDs = []
dirs = os.listdir('/Users/mcentees/Desktop/Chandra')
for x in dirs:
    if os.path.isdir(os.path.join('/Users/mcentees/Desktop/Chandra', x)):
        obsIDs.append(x)

obsIDs.remove('ArLac')
obsIDs.remove('G21')
obsIDs.remove('Horizons_files')
obsIDs.remove('JGR_light_curves')
obsIDs.remove('18303')
obsIDs.remove('23368')
obsIDs.remove('23370')
obsIDs.remove('23371')
obsIDs.remove('23372')
obsIDs.remove('Rayleigh')

obsIDs_int = [int(x) for x in obsIDs]
obsIDs = [str(x) for x in np.sort(obsIDs_int)]'''
# Accounting for different filepaths of ObsIDs that originlly had SAMP values and others that did not.
df = pd.read_csv('/Users/mcentees/Desktop/Chandra/ObsIDs_with_samp.txt', header=None, delimiter='\t')
samp_ids = np.array(df.iloc[:,0])

df2 = pd.read_csv('/Users/mcentees/Desktop/Chandra/no_samp.txt', header=None, delimiter='\t')
no_samp_ids = np.array(df2.iloc[:,0])

cts_D_cxo = []
ct_rate_cxo = []

int_flux_goes = []
med_flux_goes = []

# Reading in excel file
catalogue = pd.read_excel('/Users/mcentees/Desktop/Chandra/catalogue_PI_w_angles.xlsx')
J_E_dist = np.array(catalogue['J-E Distance (AU)'])
esj_ang = np.array([catalogue['E-S-J Angle (deg)'][0], catalogue['E-S-J Angle (deg)'][19]])
esj_lt_60 = []
J_E_dist_new = []

count = 0
for obsID in obsIDs:
    if int(obsID) in samp_ids:
        folder_path = '/Users/mcentees/Desktop/Chandra/' + str(obsID) + '/primary'
    else:
        folder_path = '/Users/mcentees/Desktop/Chandra/' + str(obsID) + '/repro'
    cor_evt_location = []
    for file in os.listdir(str(folder_path)):
        if file.startswith("hrcf") and file.endswith("pytest_evt2.fits"):
            cor_evt_location.append(os.path.join(str(folder_path), file))
    hdulist = pyfits.open(cor_evt_location[0], dtype=float)
    img_head = hdulist[1].header
    img_events = hdulist['EVENTS'].data
    date_evt = img_head['DATE-OBS'] # Start date of observation
    date_end = img_head['DATE-END'] # End date of observation
    tstart_evt = img_head['TSTART']
    tstop_evt = img_head['TSTOP']
    hdulist.close()
    exp_evt = tstop_evt - tstart_evt 
    exp_ks = np.round(exp_evt/1.e3, 2)

    # Panel plot with Jup disk lc in first panel and goes lc in second panel with times overlapping


    # Reading in data within PI filter 
    PI_file = pd.read_csv(str(folder_path) + f'/{obsID}_photonlist_PI_filter.txt')
    lat = np.array(PI_file['lat (deg)'].tolist())
    if esj_ang[count] < 60.0:
        cts_D_cxo.append(len(lat))
        ct_rate_cxo.append(len(lat)/exp_ks)
        J_E_dist_new.append(J_E_dist[count])
    
    lon = np.array(PI_file['SIII_lon (deg)'].tolist())
    time = PI_file['# t(s)'].tolist()
    time_min = [(time[i] - tstart_evt)/60 for i in range(len(time))]
    tstop_min = (tstop_evt - tstart_evt)/60
    
    bin_arr = np.arange(ceil(tstop_min + 5), step=5)
    counts, bins = np.histogram(time_min, bins=bin_arr)
    centre = (bins[:-1] + bins[1:])/2
    fig2 = plt.figure(figsize=(8,5))
    ax2 = fig2.add_subplot(211)
    loc, label =  plt.xticks()
    ax2.plot(centre, counts, 'k', zorder=0, label='5-min bin')
    '''from astropy.convolution import convolve, Box1DKernel
    boxcar = convolve(counts, Box1DKernel(11))
    ax2.plot(centre, boxcar, color='cyan', zorder = 10, label = '11-min Boxcar')'''
    major_xticks_time = np.arange(0, ceil(tstop_min + 6), step=100)
    minor_xticks_time = np.arange(0, ceil(tstop_min + 6), step=10)


    major_yticks_cnts = np.arange(0, 12, step=1)
    minor_yticks_cnts = np.arange(0 , 11, step=0.1)
    # major_yticks_cnts = np.arange(0, max(counts) + 1, step=1)
    # minor_yticks_cnts = np.arange(0 , max(counts), step=0.1)
    ax2.set_xticks(major_xticks_time)
    ax2.set_xticks(minor_xticks_time, minor=True)
    # ax2.set_yticks(major_yticks_cnts)
    # ax2.set_yticks(minor_yticks_cnts, minor=True)
    ax2.tick_params(which='both', direction='in', top=True, right=True)
    ax2.set_xticks([])
    ax2.legend(loc='upper right')
    plt.ylabel('Counts')
    plt.xlim(-5, ceil(tstop_min) + 8)
    #plt.ylim(-0.2, 11.2)
    plt.ylim(-1, 25)
    plt.title(f'ObsID: {obsID} - Exp time: {exp_ks} ks - Earth-Sun-Jupiter angle: {np.round(esj_ang[count], 2)} deg')
    # plt.savefig(str(folder_path) + f'/{obsID}_lightcurve_w_PI_filter.png', dpi=500)
    # goes data that aligns with cxo data
    ax3 = fig2.add_subplot(212)


    # Reading in GOES data
    GOES_file = pd.read_csv(str(folder_path) + f'/{obsID}_GOES_XRS.txt')
    t_goes = np.array(GOES_file['# time (min)'].tolist())
    flux_goes = np.array(GOES_file['flux (W/m^2)'].tolist())
    if esj_ang[count] < 60.0:
        med_flux_goes.append(np.nanmedian(flux_goes))
    

    # Rebinning data from 3-sec to 60-sec bins in order to compare lcs
    df_goes = pd.DataFrame({'time': t_goes, 'flux': flux_goes})
    df_goes['bin_i'] = np.digitize(np.array(df_goes['time']), bin_arr)
    binned_data = df_goes.groupby('bin_i').agg({'flux': np.sum})

    ax3.plot(centre, binned_data, 'r')
    ax3.set_xticks(major_xticks_time)
    ax3.set_xticks(minor_xticks_time, minor=True)

    ax3.set_yscale('log')
    ax3.set_ylim(2e-6, 2e-3)
    # ax3.set_yticks(major_yticks_cnts)
    # ax3.set_yticks(minor_xticks_time, minor=True)
    ax3.tick_params(which='both', direction='in', top=True, right=True)
    plt.xlabel('Time (min)'); plt.ylabel(r'Flux ($\mathrm{W m^{-2}}$)')
    plt.xlim(-5, ceil(tstop_min) + 8)
    plt.subplots_adjust(hspace=0)
    # plt.savefig(folder_path + f'/{obsID}_panel_plot_lcs.png', dpi=500)

    # Correlation plot
    fig3, ax4  = plt.subplots(figsize=(8,5))
    ax4.scatter(binned_data, counts, color='k', s=5)
    from scipy.stats import pearsonr
    corr, _ = pearsonr(np.array(binned_data.iloc[:,0]), counts)
    plt.title(f'ObsID: {obsID} - Earth-Sun-Jupiter angle: {np.round(esj_ang[count], 2)} deg - Pearsons Correlation: {np.round(corr, 2)}')
    plt.ylabel('Counts'); plt.xlabel(r'Flux ($\mathrm{W m^{-2}}$)')
    count += 1
    # plt.savefig(folder_path + f'/{obsID}_corr_plot.png', dpi=500)

'''plt.figure()
plt.scatter(int_flux_goes, cts_D_cxo, color='k', s=5)
from scipy.stats import pearsonr
corr, _ = pearsonr(int_flux_goes, cts_D_cxo)
plt.title(f'Total Counts vs Total Integrated Flux \nPearsons Correlation: {np.round(corr, 2)}')
plt.ylabel('Counts'); plt.xlabel(r'Flux ($\mathrm{W m^{-2}}$)')'''

'''plt.figure()
plt.scatter(med_flux_goes, ct_rate_cxo/np.array(J_E_dist_new), color='k', s=5)
from scipy.stats import pearsonr
corr, _ = pearsonr(med_flux_goes, ct_rate_cxo/np.array(J_E_dist_new))
plt.title(f'Normalised Count Rate vs Median Flux \nPearsons Correlation: {np.round(corr, 2)}')
plt.ylabel(r'Normalised count rate $\left(\frac{c}{ks\;AU}\right)$'); plt.xlabel(r'Flux ($\mathrm{W m^{-2}}$)')
plt.tight_layout()
plt.ylim(3.9, 12.1)
# plt.xlim(-0.1e-6, 2.1e-6)
plt.xscale('log')'''

# plt.figure()
# plt.scatter(med_flux_goes, ct_rate_cxo, color='k', s=5)
# from scipy.stats import pearsonr
# corr, _ = pearsonr(med_flux_goes, ct_rate_cxo)
# plt.title(f'Count Rate vs Median Flux \nPearsons Correlation: {np.round(corr, 2)}')
# plt.ylabel('Count rate (c/ks)'); plt.xlabel(r'Flux ($\mathrm{W m^{-2}}$)')
# plt.xscale('log')
plt.show()


