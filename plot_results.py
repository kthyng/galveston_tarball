'''
Plot specific dates/times with data
'''

import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from glob import glob
import netCDF4 as netCDF
from datetime import datetime, timedelta
from pyproj import Proj
from scipy import ndimage
import matplotlib as mpl
import cmocean.cm as cmo


land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face')
dwh = [-88.386944, 28.736667]

pc = ccrs.PlateCarree()
merc = ccrs.Mercator(central_longitude=-85.0)
utm = ccrs.UTM(zone=15)

cmap = cmo.thermal


# oilnames = ['Light_TB', 'Light', 'Moderate', 'Mod_TB', 'Neg_TB', 'Very_Light', 'Heavy_TB', 'Heavy']
# removed tarball names since those are all older
# oilnames = ['Very_Light', 'Light', 'Moderate', 'Heavy']
# oilcolors = ['0.8', '0.6', '0.4', '0.2']
oilprops = {'Very_Light': {'color': '0.8', 'zorder': 5, 'ms': 18},
            'Light': {'color': '0.6', 'zorder': 6, 'ms': 16},
            'Moderate': {'color': '0.4', 'zorder': 7, 'ms': 14},
            'Heavy': {'color': '0.2', 'zorder': 8, 'ms': 12}}

lastdate = datetime(2010, 7, 14, 0, 0)  # last date possible in sims
firstdate = datetime(2010, 6, 29, 0, 0)  # last date possible in sims

# loop through dates and plot one per day to start
dates = [datetime(2010, 6, 27, 0, 0)]
# dates = [datetime(2010, 5, 2, 0, 0) + timedelta(seconds=86400*i) for i in range(72)]
Files = glob('data/erma/layer_*/2010????_??')
for i in range(len(Files)):
    try:
        if '.zip' in Files[i]:
            Files.pop(i)
    except:  # will run out of indices to use
        pass
for date in dates:
    print(date)
    # Set up figure and map of area
    fig = plt.figure(figsize=(12,6))
    ax = plt.axes(projection=merc)
    gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabels_top = False  # turn off labels where you don't want them
    gl.ylabels_right = False
    ax.set_extent([-95.5, -88, 27.6, 30.15], pc)
    ax.add_feature(land_10m, facecolor='#D8ECA7')
    ax.plot(dwh[0], dwh[1], 'k^', ms=10, transform=pc, label='Deepwater Horizon site', zorder=11)


    # read in coastal data
    fnamecoast = 'data/erma/layer_19872/18611/SCAT_Max_OilZones_Louisiana_2014_09_30Sorted.shp'
    # coast = cfeature.ShapelyFeature(shpreader.Reader(fnamecoast).geometries(), utm, facecolor='none')

    # ## plot coastal oiling
    # reader = shpreader.Reader(fnamecoast)
    # records = reader.records()
    # oiled = []
    # for record in records:
    #     rdate = record.attributes['max_date']
    #     oilcode = record.attributes['OilCatCode']
    #     west, south, east, north = record.bounds  # bounding box in UTM coords
    #     lonplot, latplot = pc.transform_points(utm, np.array([[west]]), np.array([[north]]))[0][0][0:2]
    #     if oilcode in oilprops.keys() and rdate[:8] == date.strftime('%Y%m%d'):
    #         # print(oilcode, rdate)
    #         color = oilprops[oilcode]['color']
    #         zorder = oilprops[oilcode]['zorder']
    #         ms = oilprops[oilcode]['ms']
    #         coast = cfeature.ShapelyFeature(record.geometry, utm, facecolor='none')
    #         ax.plot(lonplot, latplot, marker='o', ms=ms, mfc='none', mec=color, mew=3.5, zorder=zorder, transform=pc)
    #         # # this looks poorly plotted since it's a small area. Instead plot a marker.
    #         # ax.add_feature(coast, facecolor=color, edgecolor=color, linewidth=15, zorder=15)


    ## read in model output from date
    # can't use tracpy because of a cartopy problem. recreating code needed here.
    inputs = {'proj': 'lcc', 'ellps': 'clrk66', 'datum': 'NAD27',
              'lat_1': 22.5, 'lat_2': 31.0, 'lat_0': 30, 'lon_0': -94,
              'x_0': 0, 'y_0': 0}
    proj = Proj(**inputs)
    grid = netCDF.Dataset('../grid.nc')
    dFiles = glob('tracks/*.nc')
    # save over loops to plot afterward
    lonpsave = []; latpsave = []
    timecolors = []; dtimes = []
    for File in dFiles:
        print(File)
        d = netCDF.Dataset(File)
        xg = d['xg'][:]; yg = d['yg'][:]
        tp = d['tp'][0,:]
        d.close()

        drifterdates = netCDF.num2date(tp[:], 'seconds since 1970-01-01')
        idate = np.where(drifterdates == date)[0][0]

        # change from grid space to ll space
        lonp = ndimage.map_coordinates(grid['lon_rho'][:].T, np.array([xg[:,idate].flatten()+.5,
                                                               yg[:,idate].flatten()+.5]),
                                     order=1, mode='nearest',
                                     cval=0.).reshape(xg[:,idate].shape)
        latp = ndimage.map_coordinates(grid['lat_rho'][:].T, np.array([xg[:,idate].flatten()+.5,
                                                               yg[:,idate].flatten()+.5]),
                                     order=1, mode='nearest',
                                     cval=0.).reshape(yg[:,idate].shape)
        if lonp.size > 1:
            print('using drifters')
            ind = np.isnan(xg[:,idate])
            lonp[ind] = np.nan
            latp[ind] = np.nan

            lonpsave.append(lonp)
            latpsave.append(latp)



        # lonp, latp, _ = tracpy.tools.interpolate2d(xg[:,idate], yg[:,idate], grid, 'm_ij2ll')
        # ax.plot(lonp, latp, 'o', color='darkcyan', alpha=0.3, transform=pc, mec=None, mew=0, ms=6, zorder=9);

            # time for a drifter is between idate and tracks earliest time (2010/4/20)
            dtimes.append((drifterdates[0] - date).days)
            # fractional colors for each circle in order
            # fracs = dtimes/(vmax-vmin)
            # cmap.N = 1
            # timecolors = cmap(norm(dtime))

        # # Add on vulnerability of coastline
        # # Need to plot the coast boxes as patches and color them according to vulnerability level
        # # http://matplotlib.org/1.2.1/examples/pylab_examples/hist_colormapped.html
        # # we need to normalize the data to 0..1 for the full
        # # range of the colormap
        # fracs = ndbox[i,j,:].astype(float)/ndbox[:,j,:].max() # max across years
        # norm = colors.Normalize(fracs.min(), fracs.max())
        #
        # # Save patches together
        # patches = []
        # for path in pathsxy:
        #     patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, edgecolor=None, zorder=5))
        #
        # # assign shades of colormap to the patches according to values, and plot
        # for thisfrac, thispatch in zip(fracs, patches):
        #     color = cmo.matter(norm(thisfrac))
        #     thispatch.set_facecolor(color)
        #     ax.add_patch(thispatch)

    vmin = (firstdate-date).days
    vmax = (lastdate-date).days
    # norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    mappable = ax.scatter(np.asarray(lonpsave), np.asarray(latpsave), s=50,
                          c=np.asarray(dtimes)[:,np.newaxis], alpha=0.3,
                          cmap=cmap, transform=pc, zorder=9, vmin=vmin, vmax=vmax);#,
                          # norm=norm);#, mec=None, mew=0, ms=6);

    # colorbar
    cax = fig.add_axes([0.7, 0.05, 0.32, 0.018]) #colorbar axes
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label('Days to reach Galveston from date', fontsize=13, color='0.2')
    # cb.ax.tick_params(labelsize=12, length=2, color='0.2', labelcolor='0.2')
    # cb.set_ticks(ticks)

    ##


    ## plot oil slicks (boxes)
    # locations of files matching date. Could include thick and thin oil data and zip files.
    inds = np.where([date.strftime('%Y%m%d') in File for File in Files])[0]
    # import pdb; pdb.set_trace()
    filesdate = [Files[i] for i in inds]

    for filedate in filesdate:  # loop over thick and thin oil files
        fnameoiling = filedate + '/Daily_Grid_' + filedate.split('/')[-1]
        print(fnameoiling)
        oiling = cfeature.ShapelyFeature(shpreader.Reader(fnameoiling).geometries(), pc, facecolor='none')

        if 'TN' in fnameoiling:
            color = oilprops['Light']['color']
            lw = .4
        elif 'TK' in fnameoiling:
            color = oilprops['Heavy']['color']
            lw = .7
        ax.add_feature(oiling, edgecolor=color, linewidth=lw, zorder=10)  # add oiling data
    # ax.add_feature(oiling, edgecolor='r', linewidth=0.5)  # add oiling data
    # ax.add_feature(coast, edgecolor='r', linewidth=0.5)  # add coast oiling data
    ##

    # import pdb; pdb.set_trace()
    ax.text(-92.6, 29.9, date.strftime('%B %d, %Y %H:%M'), transform=pc, fontsize=16)
    ax.legend(loc=2, numpoints=1, frameon=False)

    # import pdb; pdb.set_trace()
    fig.savefig('figures/comp/' + date.strftime('%Y%m%d') + '.png', bbox_inches='tight', dpi=300)
    plt.close(fig)
    # plt.show()
