#relevant packages 
import numpy as np
import pandas as pd
from astropy.io import fits as pyfits
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Ellipse
from matplotlib.patches import Arc 
from matplotlib.patches import Rectangle 
from math import ceil
import os
from scipy import interpolate
from astropy.time import Time

# Reading in sunspot data
sun_data = pd.read_csv("/Users/mcentees/Desktop/Chandra/SN_m_tot_V2.0.csv", header=None, sep=";")

t = np.array(sun_data.iloc[:, 2])
ssp = np.array(sun_data.iloc[:, 3])

# Need arrays for sunspot plot disk count rate panel
date_obs_dec = []
disk_cr = []

# obsIDs = ['1862', '2519', '15669', '18608', '20000', '18678', '22146']
# obsIDs = ['15671']
# obsIDs = ['1862', '2519', '20002']
obsIDs = ['1862']

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
obsIDs.remove('1862')
obsIDs.remove('2519')
obsIDs.remove('15669')
obsIDs.remove('18608')
obsIDs.remove('20000')
obsIDs.remove('18678')
obsIDs.remove('22146')
obsIDs.remove('23368')
obsIDs.remove('23370')
obsIDs.remove('23371')
obsIDs.remove('23372')'''


obsIDs_int = [int(x) for x in obsIDs]
obsIDs = [str(x) for x in np.sort(obsIDs_int)]
# Accounting for different filepaths of ObsIDs that originlly had SAMP values and others that did not.
df = pd.read_csv('/Users/mcentees/Desktop/Chandra/ObsIDs_with_samp.txt', header=None, delimiter='\t')
samp_ids = np.array(df.iloc[:,0])

df2 = pd.read_csv('/Users/mcentees/Desktop/Chandra/no_samp.txt', header=None, delimiter='\t')
no_samp_ids = np.array(df2.iloc[:,0])

for obsID in obsIDs:
    # Read in fits file from correct path
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
    obs_id = img_head['OBS_ID'] # ObsID
    date = img_head['DATE-OBS'] # Start date of observation
    tstart = img_head['TSTART']
    tstop = img_head['TSTOP']
    hdulist.close()
    duration = tstop - tstart
    exp_time_ks = duration/1.e3

    # Reading in original go_chandra_data
    go_chandra = pd.read_csv(folder_path + f'/{obsID}_photonlist_full_obs_w_samp.txt')
    x_gc = np.array(go_chandra.iloc[:,1])
    y_gc = np.array(go_chandra.iloc[:,2])

    # Reading in ellipse data
    ellipse_data = pd.read_csv(folder_path + f'/{obsID}_photonlist_full_obs_ellipse.txt')
    x_ell = np.array(ellipse_data.iloc[:,1])
    y_ell = np.array(ellipse_data.iloc[:,2])

    # reading in amplifier signals
    av1_jup = np.array(ellipse_data['av1'])
    av2_jup = np.array(ellipse_data['av2'])
    av3_jup = np.array(ellipse_data['av3'])

    au1_jup = np.array(ellipse_data['au1'])
    au2_jup = np.array(ellipse_data['au2'])
    au3_jup = np.array(ellipse_data['au3'])
    amp_sf_jup = np.array(ellipse_data['amp_sf'])
    sumamp_jup = av1_jup + av2_jup + av3_jup + au1_jup + au2_jup + au3_jup 
    samp_jup = (sumamp_jup * (2. ** (amp_sf_jup - 1.0)))/148.0


    # Reading in orignial photon results
    full_image= pd.read_csv(folder_path + f'/{obsID}_all_photons.txt', header=None, delimiter=' ')
    bg_reg = full_image[(full_image[0] > -100) & (full_image[0] < 100) & (full_image[1] > -100) & (full_image[1] < 100)]
    x = np.array(bg_reg.iloc[:,0]) 
    y = np.array(bg_reg.iloc[:,1])

    # reading in amplifier signals for whole region
    av1_bg = np.array(bg_reg.iloc[:,8])
    av2_bg = np.array(bg_reg.iloc[:,9])
    av3_bg = np.array(bg_reg.iloc[:,10])

    au1_bg = np.array(bg_reg.iloc[:,11])
    au2_bg = np.array(bg_reg.iloc[:,12])
    au3_bg = np.array(bg_reg.iloc[:,13])
    amp_sf_bg = np.array(bg_reg.iloc[:,7])
    sumamp_bg = av1_bg + av2_bg + av3_bg + au1_bg + au2_bg + au3_bg
    samp_bg = (sumamp_bg * (2. ** (amp_sf_bg - 1.0)))/148.0

    # Reading in JPL data
    JPL_data = pd.read_csv(folder_path + f'/{obsID}_JPL_ellipse_vals.txt')
    ang_diam = JPL_data.iloc[0,0]
    tilt_ang_deg = JPL_data.iloc[0,1]
    tilt_ang_rad = np.deg2rad(tilt_ang_deg) # Converting degrees to radians

    ecc = 0.3543163789650 
    R_eq = (ang_diam/2.)/np.cos(tilt_ang_rad) # Equatorial radius of Jupiter
    R_p = R_eq * np.sqrt(1 - ecc**2) # Polar radius of Jupiter
    area_jup = np.pi * R_eq * R_p # area of jupiter in units of arcsec**2

    R_eq_bg = 1.5 * (ang_diam/2.)/np.cos(tilt_ang_rad) # Equatorial radius of Jupiter
    R_p_bg = R_eq_bg * np.sqrt(1 - ecc**2) # Polar radius of Jupiter
    area_bg = (200. * 200) - (np.pi * R_eq_bg * R_p_bg)

    ell = Ellipse((0,0), R_eq * 2., R_p * 2., tilt_ang_deg, facecolor='None', edgecolor='m', lw=1, zorder=50) # ellipse outline
    ell_bg = Ellipse((0,0), R_eq_bg * 2., R_p_bg * 2., tilt_ang_deg, facecolor='None', edgecolor='r', lw=1, zorder=50) # ellipse outline
    circ = Ellipse((0,0), 60., 60., facecolor='None', edgecolor='k', lw=1, zorder=50) # original go_chandra outline

    box = Rectangle((-100, -100), 200, 200, facecolor='None', edgecolor='r', lw=1, zorder=50)

    x_bg = []
    y_bg = []
    samp_bg_reg = []
    for i in range(len(x)):
        if (x[i] * np.cos(tilt_ang_rad) + y[i] * np.sin(tilt_ang_rad)) ** 2./(R_eq_bg ** 2.) + (x[i] * np.sin(tilt_ang_rad) - y[i] * np.cos(tilt_ang_rad)) ** 2./(R_p_bg ** 2.) > 1.0:
            x_bg.append(x[i])
            y_bg.append(y[i])
            samp_bg_reg.append(samp_bg[i])

    #PI calculation 
    samp_bg_reg = np.array(samp_bg_reg)

    date_dec = Time(date).decimalyear
    date_obs_dec.append(date_dec)

    g= 1.0418475 + 0.020125799 * (date_dec - 2000.) + 0.010877227 * (date_dec - 2000.) ** 2. + - 0.0014310146 * (date_dec - 2000.) ** 3. + 5.8426766e-05 * (date_dec - 2000.) ** 4. # multiplicative scaling factor
    PI_jup = g * samp_jup 
    PI_bg = g * samp_bg_reg 


    bins_PI = np.arange(0, 1100)
    counts_PI_jup, bin_PI_jup = np.histogram(PI_jup, bins=bins_PI)
    centre_PI = (bin_PI_jup[:-1] + bin_PI_jup[1:])/2

    counts_PI_bg, bin_PI_bg = np.histogram(PI_bg, bins=bins_PI)

    norm_counts_jup = counts_PI_jup/area_jup
    norm_counts_bg = counts_PI_bg/area_bg

    '''# Linear gain correction factor 
    gC = 1.8822853 + 7.0183144 * (date_dec - 2000) + -0.70972911 * (date_dec - 2000) ** 2 + 0.022186861 * (date_dec - 2000) ** 3
    gS = 1.0109385 + 0.022563058 * (date_dec - 2000) + 0.0061835047 * (date_dec - 2000) ** 2 +-0.00088251078 * (date_dec - 2000) ** 3 + 3.7804516e-05 * (date_dec - 2000)**4
    
    PI_lin_jup = gC + gS * samp_jup
    PI_lin_bg = gC + gS * samp_bg

    counts_PI_lin_jup = np.histogram(PI_lin_jup, bins=bins_PI)[0]
    counts_PI_lin_bg = np.histogram(PI_lin_bg, bins=bins_PI)[0]'''
    
    # Schematic diagram showing location Jup and bg regions
    fig, ax = plt.subplots(figsize=(7,7))
    ax.add_artist(ell)
    ax.add_artist(ell_bg)
    ax.add_artist(circ)
    ax.add_artist(box)
    ax.plot(x, y, 'w.', markersize=1, zorder=0) #, alpha=0.3)
    ax.plot(x_gc, y_gc,'k.', markersize=1, zorder=5)
    ax.plot(x_bg, y_bg, 'r.', markersize=1, zorder=20)
    ax.plot(x_ell, y_ell, 'm.', markersize=1, zorder=10)
    ax.set_xlim(-105, 105)
    ax.set_ylim(-105, 105)
    ax.set_xlabel('x (arcsec)'); ax.set_ylabel('y (arcsec)')
    ax.set_title(f'ObsID: {obs_id} \n{date}')
    plt.tight_layout()
    # plt.savefig(folder_path + f'/{obsID}_bg_plus_ell.png', dpi=500)

    # cumulative distributions for background and source
    bgcumul = np.cumsum(counts_PI_bg)
    bgcdf = bgcumul / max(bgcumul)

    jupcumul = np.cumsum(counts_PI_jup)
    jupcdf = jupcumul / max(jupcumul)

    xfrac = 0.05
    pimin = np.interp(xfrac/2, jupcdf, centre_PI)
    pimax = np.interp(1-xfrac/2, jupcdf, centre_PI)
    lo = int(pimin)
    hi = int(pimax) + 1
    
    print('\n' + '-' * 40)
    print(f'ObsID: {obsID} \nYear: {date[0:4]}')
    print(f'g = {g}')
    print(f'PI range: {lo} to {hi}')

    blo = np.interp(lo, centre_PI, bgcdf)
    bhi = 1.0 - np.interp(hi, centre_PI, bgcdf)
    bfrac = bhi + blo

    # bg_reg = full_image[(full_image[0] > -100) & (full_image[0] < 100) & (full_image[1] > -100) & (full_image[1] < 100)]
    df_filtered = ellipse_data.drop(['PHA','samp', 'sumamps', 'pi'], axis=1)
    df_filtered.insert(15, "PI", PI_jup)
    df_filtered['lat (deg)'] = df_filtered['lat (deg)'] - 90.
    df_filtered_new = df_filtered[(df_filtered["lat (deg)"] < 45) & (df_filtered["lat (deg)"] > -55) & (df_filtered["PI"] > lo) & (df_filtered["PI"] < hi)]
    df_filtered_Jup_full = df_filtered[(df_filtered["PI"] > lo) & (df_filtered["PI"] < hi)]
    disk_cr.append(len(df_filtered_new)/exp_time_ks)
    # df_filtered_new.to_csv(folder_path + f'/{obsID}_photonlist_PI_filter.txt', index=False)
    # df_filtered_Jup_full.to_csv(folder_path + f'/{obsID}_photonlist_PI_filter_Jup_full.txt', index=False)

    print(f'Fraction of background excluded is {np.round(bfrac, 3)}')

    # PI filter on Jupiter cumulative distribution function
    '''fig3, ax3 = plt.subplots()
    ax3.set_xlabel('PI', fontsize=12)
    ax3.set_ylabel(r'$\Sigma$ (counts $\leq$ PI) / $\Sigma$ (counts)', fontsize=12)
    ax3.set_title(f'ObsID: {obs_id} \n{date}', fontsize=14)
    ax3.plot(centre_PI, jupcdf, 'r', zorder=10, label='Jupiter')
    ax3.plot(centre_PI, bgcdf, 'k', zorder=5, label='Background')
    # ax3.vlines(lo, ymin=-0.1, ymax=1.1, color='darkgrey', zorder=0, label=f'PI range: {lo} to {hi}')
    # ax3.vlines(hi, ymin=-0.1, ymax=1.1, color='darkgrey', zorder=0)
    ax3.set_ylim(-0.05, 1.05)
    ax3.set_xlim(-20, 1020)
    ax3.fill_between([lo, hi], -0.05, 1.05, color='darkgrey', zorder=0, label=f'PI range: {lo} to {hi}')
    ax3.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(folder_path + f'/{obsID}_PI_filter.png', dpi=500)'''

    # PI spectra
    '''fig2, ax2 = plt.subplots(figsize=(8,6))
    ax2.plot(centre_PI, norm_counts_jup, 'r', label='Jupiter')
    ax2.plot(centre_PI, norm_counts_bg, 'k', label='Background')
    # ax2.fill_between([lo, hi], -0.001, max(norm_counts_jup) + .001, color='darkgrey', zorder=0, label=f'PI range: {lo} to {hi}')
    ax2.set_xlabel('PI', fontsize=14); ax2.set_ylabel(r'$\mathrm{counts/arcsec}^2$', fontsize=14)
    plt.xticks(fontsize=14); plt.yticks(fontsize=14)
    ax2.set_title(f'ObsID: {obs_id} \n{date}', fontsize=16)
    # ax2.set_xlim(0,hi + 30)
    ax2.set_xlim(0, 1000)
    # ax2.set_ylim(-.001, max(norm_counts_jup) + .001)
    plt.tight_layout()
    # plt.savefig(folder_path + f'/{obsID}_PI_spectrum.png', dpi=500)
    # plt.savefig(folder_path + f'/{obsID}_PI_spectrum_w_filter.png', dpi=500)'''

    # Normalised SAMP plot for powerpoint slides
    '''bins_SAMP = np.arange(0, 361)
    fig6, ax6 = plt.subplots(figsize=(8,6))

    counts_SAMP_jup, bin_arr = np.histogram(samp_jup, bins=bins_SAMP)
    counts_SAMP_bg = np.histogram(samp_bg_reg, bins=bins_SAMP)[0]
    centre_SAMP = (bin_arr[:-1] + bin_arr[1:])/2
    ax6.plot(centre_SAMP, counts_SAMP_jup/area_jup, 'r', label='Jupiter')
    ax6.plot(centre_SAMP, counts_SAMP_bg/area_bg, 'k', label='Background')
    ax6.set_xlabel('SAMP'); ax6.set_ylabel(r'$\mathrm{counts/arcsec}^2$')
    ax6.set_title(f'ObsID: {obs_id} \n{date}')
    plt.tight_layout()
    plt.savefig(folder_path + f'/{obsID}_norm_SAMP.png', dpi=500)'''

# Sunspot plot comparison with filtered disk count rates
'''fig4 = plt.figure(figsize=(8,6))
ax4 = fig4.add_subplot(211)
loc, label =  plt.xticks()
ax4.plot(t, ssp, 'k', zorder=10)
major_xticks = np.arange(2000, 2021, 5)
minor_xticks = np.arange(1999, 2022, 1)
plt.ylabel('Sunspot Number')
plt.xlim(1999, max(t)); plt.ylim(-10, 360)
plt.title('Chandra X-ray Observations')
ax4.set_xticks(major_xticks)
ax4.set_xticks(minor_xticks, minor=True)

major_yticks = np.arange(0, 351, 50)
minor_yticks = np.arange(0, 351, 10)
ax4.set_yticks(major_yticks)
ax4.set_yticks(minor_yticks, minor=True)

ax4.tick_params(which='both', direction='in', top=True, right=True)

count = 0
for x in date_obs_dec:
    if (count == 0) or (count == 19):
        ax4.vlines(x, ymin=-10, ymax=300, color='b', zorder=0, label='HRC-I')
    else:
        ax4.vlines(x, ymin=-10, ymax=300, color='b', zorder=0, label='HRC-I')
    count += 1

ax5 = fig4.add_subplot(212)
hrc_cr = ax5.scatter(date_obs_dec, disk_cr, color='b', marker='^', zorder=0)
# ax5.scatter([date_obs_dec[0], date_obs_dec[19]], [disk_cr[0], disk_cr[19]], color='r', marker= '^', zorder=100)
plt.xlim(1999, max(t))
ax5.set_xticks(major_xticks)
ax5.set_xticks(minor_xticks, minor=True)
major_yticks = np.arange(20, 51, 5)
minor_yticks = np.arange(20, 51, 1)
ax5.set_yticks(major_yticks)
ax5.set_yticks(minor_yticks, minor=True)
ax5.tick_params(which='both', direction='in', top=True, right=False)
plt.ylim(19, 51)
plt.xlabel('Year')
plt.ylabel('Counts/ks')
plt.subplots_adjust(hspace=0)'''

plt.show()
