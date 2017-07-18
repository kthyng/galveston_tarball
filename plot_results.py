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
import tracpy
from datetime import datetime


land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face')
dwh = [-88.386944, 28.736667]

pc = ccrs.PlateCarree()
merc = ccrs.Mercator(central_longitude=-85.0)
utm = ccrs.UTM(zone=15)


# read in data
fnameoiling = 'data/erma/layer_19763/20100520_TN/Daily_Grid_20100520_TN.shp'
oiling = cfeature.ShapelyFeature(shpreader.Reader(fnameoiling).geometries(), pc, facecolor='none')
fnamecoast = 'data/erma/layer_19872/18611/SCAT_Max_OilZones_Louisiana_2014_09_30Sorted.shp'
coast = cfeature.ShapelyFeature(shpreader.Reader(fnamecoast).geometries(), utm, facecolor='none')


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
ax.plot(dwh[0], dwh[1], 'k^', ms=10, transform=pc, label='Deepwater Horizon site')



#

# oilnames = ['Light_TB', 'Light', 'Moderate', 'Mod_TB', 'Neg_TB', 'Very_Light', 'Heavy_TB', 'Heavy']
# removed tarball names since those are all older
# oilnames = ['Very_Light', 'Light', 'Moderate', 'Heavy']
# oilcolors = ['0.8', '0.6', '0.4', '0.2']
oilprops = {'Very_Light': {'color': '0.8', 'zorder': 5},
            'Light': {'color': '0.6', 'zorder': 6},
            'Moderate': {'color': '0.4', 'zorder': 7},
            'Heavy': {'color': '0.2', 'zorder': 8}}
reader = shpreader.Reader(fnamecoast)
records = reader.records()
oiled = []
for record in records:
    date = record.attributes['max_date']
    oilcode = record.attributes['OilCatCode']
    # west, south, east, north = record.bounds  # bounding box in UTM coords
    # pc.transform_points(utm, np.array([[west, east]]), np.array([[north, south]]))[0][:,:2]
    if oilcode in oilprops.keys() and date[:8] == '20100521':
        print(oilcode, date)
        # if oilcode == 'Very_Light':
        #     oilcolor = oilcolors[0]
        color = oilprops[oilcode]['color']
        zorder = oilprops[oilcode]['zorder']
        coast = cfeature.ShapelyFeature(record.geometry, utm, facecolor='none')
        ax.add_feature(coast, facecolor=color, edgecolor=color, linewidth=5)
        # ax.add_geometries([record.geometry], utm, facecolor=color, edgecolor=color, zorder=zorder, linewidth=10)

        # oiled.append(record)
color = oilprops['Light']['color']
ax.add_feature(oiling, edgecolor=color, linewidth=0.2)  # add oiling data
# ax.add_feature(oiling, edgecolor='r', linewidth=0.5)  # add oiling data
# ax.add_feature(coast, edgecolor='r', linewidth=0.5)  # add coast oiling data


# read in model output from May 20, 2010
proj = tracpy.tools.make_proj('nwgom-pyproj')
grid = tracpy.inout.readgrid(loc, proj)
Files = glob('tracks/*.nc')
for File in Files:
    d = netCDF.Dataset(File)
    xg = d['xg'][:]; yg = d['yg'][:]
    tp = d['tp'][0,:]
    d.close()

    dates = netCDF.num2date(tp[:], 'seconds since 1970-01-01')
    idate = np.where(dates == datetime(2010, 5, 20, 0, 0))[0][0]
    # change from grid space to ll space
    lonp, latp, _ = tracpy.tools.interpolate2d(xg[:,idate], yg[:,idate], grid, 'm_ij2ll')
    ax.plot(lonp[:].T, latp[:].T, 'darkcyan.', alpha=0.3, transform=pc);
ax.text(-97.7, 29.3, date.strftime('%B %d, %Y %H:%M'), transform=pc, fontsize=16)
ax.legend(loc=2, numpoints=1, frameon=False)

plt.show()
