#!/usr/bin/env python3
# -*- eoding: utf-8 -*-

#relevant packages 
import numpy as np
import pandas as pd
from astropy.io import fits as pyfits
from astropy.time import Time
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from math import ceil
import os
from astropy.time import Time, TimeDelta
import astropy.units as u

# obsIDs = ['23368', '23370', '23371', '23372']
obsIDs = ['1862']
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

obsIDs_int = [int(x) for x in obsIDs]
obsIDs = [str(x) for x in np.sort(obsIDs_int)]'''
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
    date = img_head['DATE-OBS'] # Start date of observation
    date_end = img_head['DATE-END'] # End date of observation
    tstart = img_head['TSTART']
    tstop = img_head['TSTOP']
    hdulist.close()
    tstop_min = (tstop - tstart)/60

    cxo_tstart = Time(date, format='isot')
    cxo_tend = Time(date_end, format='isot')
    tstart_iso = cxo_tstart.to_value('iso', subfmt='date_hm')
    tend_iso = cxo_tend.to_value('iso', subfmt='date_hm')

    # Reading in data within PI filter for all of Jupiter 
    PI_file = pd.read_csv(str(folder_path) + f'/{obsID}_photonlist_PI_filter_Jup_full.txt')
    lat = np.array(PI_file['lat (deg)'].tolist()) 
    lon = np.array(PI_file['SIII_lon (deg)'].tolist())
    time = PI_file['# t(s)'].tolist()
    time_min = [(time[i] - tstart)/60 for i in range(len(time))]
    exp_t = tstop - tstart
    exp_t_ks = np.round(exp_t/1.e3, 2)


    # Making data periodic
    lon2 = lon + 360.0 
    lon_neg = lon  - 360.0

    lat_final = np.append(lat, (lat, lat))
    lon_final = np.append(lon_neg, (lon, lon2))

    # Plotting point maps for each obsID
    point_map = True 
    if point_map:
        fig = plt.figure(figsize=(8,5))
        ax = fig.add_subplot(111)
        major_xticks = np.arange(350, -1, step=-50)
        minor_xticks = np.arange(360, -1, step=-10)

        major_yticks = np.arange(-55, 46, step=10)
        minor_yticks = np.arange(-55, 46, step=5)
        plt.title(f'ObsID: {obsID} \n{tstart_iso} - {tend_iso} ({exp_t_ks} ks)')
        plt.scatter(lon_final, lat_final, s=2, color='k', zorder=20) #, alpha=0.3)
        # plt.xlim(362, -2)
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.set_yticks(major_yticks)
        ax.set_yticks(minor_yticks, minor=True)
        ax.tick_params(which='both', direction='in', top=True, right=True)
        # plt.ylim(-57, 47)
        plt.xlabel('S3 Longitude (deg)'); plt.ylabel('S3 Latitude (deg)')
        # plt.tight_layout()
        # plt.savefig(str(folder_path) + f'/{obsID}_point_plot_w_PI_filter.png', dpi=500)
        # Adding Voronoi regions to plot
        voronoi = True
        if voronoi:
            points = np.column_stack((lon_final, lat_final))
            from scipy.spatial import Voronoi, voronoi_plot_2d
            vor = Voronoi(points)
            # fig, ax3 = plt.figure(figsize=(8, 5))
            fig_new = voronoi_plot_2d(vor, ax, show_vertices=False, show_points=False, point_size=2, color='k')
            plt.xlim(362, -2)
            # plt.xlim(722, -362)
            # plt.ylim(-55, 45)
            plt.ylim(-60, 50)
            # plt.ylim(-157, 147)

            def PolyArea(x,y):
                # if (len(np.where(y_poly > 45)[0]) > 0) or (len(np.where(y_poly < -55)[0])):
                #     return 0.0
                # else:
                return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
            area_list = []
            norm = mpl.colors.Normalize(vmin=0, vmax=150, clip=True)
            # norm = mpl.colors.LogNorm(vmin=1, vmax=1e5, clip=True)
            mapper = cm.ScalarMappable(norm=norm, cmap=cm.viridis)# cmap=cm.viridis)

            # for r in range(len(vor.point_region)):
            #     region = vor.regions[vor.point_region[r]]
            
            # fixing regions so only calculate area within bounds
            for i, point in enumerate(vor.points):
                if (point[0] >= 0) and (point[0] < 360) and (point[1] > -55) and (point[1] < 45):
                    region = vor.regions[vor.point_region[i]]
                    x_poly = np.array([vor.vertices[i][0] for i in region])
                    y_poly = np.array([vor.vertices[i][1] for i in region])

                    polygon = [vor.vertices[i] for i in region]

                    # Find area of regions
                    area_poly = PolyArea(x_poly, y_poly)
                    area_list.append(area_poly)

                    plt.fill(*zip(*polygon), color=mapper.to_rgba(area_poly), alpha=1.0, zorder=0)
            cbar = fig.colorbar(mapper)
            cbar.set_label(r'Area ($\mathrm{deg}^2$)', rotation=270, labelpad=30)
            print(max(area_list))
            print(len(area_list))
            plt.tight_layout()
            # plt.savefig(str(folder_path) + f'/{obsID}_VR_tess.png', dpi=500)

            area_hist_plot = False 
            if area_hist_plot:
                area_bins = np.arange(401)
                counts, bins = np.histogram(area_list, bins=area_bins)
                norm_counts = counts/len(area_list) * 100.
                # print(max(counts))
                centre = (bins[:-1] + bins[1:])/2
                fig_area, ax_area = plt.subplots(figsize=(8,5))
                ax_area.plot(centre, norm_counts, 'k', label = f'N = {len(area_list)} polygons')
                ax_area.set_xlabel('Area'); ax_area.set_ylabel('# Polygons as % of Total')
                major_xticks = np.arange(0, 401, 50)
                minor_xticks = np.arange(0, 401, 10)

                major_yticks = np.arange(0, 9, 1)
                minor_yticks = np.arange(0, 8.1, 0.1)

                ax_area.set_xticks(major_xticks)
                ax_area.set_xticks(minor_xticks, minor=True)

                ax_area.tick_params(which='both', direction='in', top=True, right=True)

                ax_area.set_yticks(major_yticks)
                ax_area.set_yticks(minor_yticks, minor=True)
                ax_area.set_xlim(-2, 402);  ax_area.set_ylim(-.2, 8)
                ax_area.legend(loc='upper right', handlelength=0)
                plt.title(f'ObsID: {obsID} \n{tstart_iso} - {tend_iso} ({exp_t_ks} ks)')
                plt.tight_layout()
                # plt.savefig(str(folder_path) + f'/{obsID}_area_hist.png', dpi=500)

    # Plotting lightcurves for each obsID
    lightcurve = True 
    if lightcurve:
        bin_arr = np.arange(ceil(tstop_min + 1))
        counts, bins = np.histogram(time_min, bins=bin_arr)
        centre = (bins[:-1] + bins[1:])/2
        fig2, ax2 = plt.subplots(figsize=(8,5))
        ax2.plot(centre, counts, 'k', zorder=0, label='60-sec bin')
        from astropy.convolution import convolve, Box1DKernel
        boxcar = convolve(counts, Box1DKernel(11))
        ax2.plot(centre, boxcar, color='cyan', zorder = 10, label = '11-min Boxcar')
        # major_xticks_time = np.arange(0, ceil(tstop_min + 2), step=100)
        # minor_xticks_time = np.arange(0, ceil(tstop_min + 2), step=10)


        # major_yticks_cnts = np.arange(0, 12, step=1)
        # minor_yticks_cnts = np.arange(0 , 11, step=0.1)
        # major_yticks_cnts = np.arange(0, max(counts) + 1, step=1)
        # minor_yticks_cnts = np.arange(0 , max(counts), step=0.1)
        # ax2.set_xticks(major_xticks_time)
        # ax2.set_xticks(minor_xticks_time, minor=True)
        # ax2.set_yticks(major_yticks_cnts)
        # ax2.set_yticks(minor_yticks_cnts, minor=True)
        ax2.tick_params(which='both', direction='in', top=True, right=True)
        ax2.legend(loc='upper right')
        plt.xlabel('Time (min)'); plt.ylabel('Counts/min')
        plt.xlim(-5, ceil(tstop_min) + 5)
        # plt.ylim(-0.2, 11.2)
        # plt.ylim(-0.2, max(counts) + 0.2)
        ax2.set_title(f'ObsID: {obsID} \n{date}')
        plt.tight_layout()
        # plt.savefig(str(folder_path) + f'/{obsID}_lightcurve_w_PI_filter.png', dpi=500)
plt.show()
