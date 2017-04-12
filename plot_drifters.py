'''
Plot drifters backward from Galveston Bay
'''

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cartopy
crs = cartopy.crs
from glob import glob
# import xarray as xr
import netCDF4 as netCDF
import os
import tracpy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker


mpl.rcParams.update({'font.size': 12})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg'
Files = glob('tracks/*.nc')
dwh = [-88.386944, 28.736667]
import cartopy.feature as cfeature
land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face')
proj = tracpy.tools.make_proj('nwgom-pyproj')
grid = tracpy.inout.readgrid(loc, proj)

for File in Files:

    d = netCDF.Dataset(File)
    xg = d['xg'][:]; yg = d['yg'][:]
    tp = d['tp'][0,:]
    # d = xr.open_dataset(File)
    # lonp = d['lonp'].data; latp = d['latp'].data
    d.close()

    # change from grid space to ll space
    lonp, latp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2ll')

    dates = netCDF.num2date(tp[:], 'seconds since 1970-01-01')
    startdate = File.split('/')[-1][:10]
    if not os.path.exists('figures/' + startdate):
        os.mkdir('figures/' + startdate)

    for i, date in enumerate(dates[::5]):
        figname = 'figures/' + startdate + '/'  + str(i).zfill(3) + '.png'
        # bfigname = 'backward-' + str(len(dates[::5])-i).zfill(3) + '.png'  # to show backwards
        # figname = 'figures/' + startdate + '/' + date.isoformat()
        if os.path.exists(figname):
            continue

        fig = plt.figure(figsize=(12,6))
        ax = plt.axes(projection=crs.Mercator(central_longitude=-85.0))
        gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabels_top = False  # turn off labels where you don't want them
        gl.ylabels_right = False
        ax.set_extent([-98, -88, 26, 30.15], crs.PlateCarree())
        ax.add_feature(land_10m, facecolor='0.8')
        ax.plot(dwh[0], dwh[1], 'k^', ms=10, transform=crs.PlateCarree(), label='Deepwater Horizon site')
        ax.plot(lonp[:,i].T, latp[:,i].T, 'k.', lw=0.3, alpha=0.7, transform=crs.PlateCarree());
        ax.text(-97.7, 29.3, date.strftime('%B %d, %Y %H:%M'), transform=crs.PlateCarree(), fontsize=16)
        ax.legend(loc=2, numpoints=1, frameon=False)
        fig.savefig(figname, bbox_inches='tight')
        # os.symlink(figname, bfigname)
        plt.close(fig)
