'''
Plot some drifter tracks and location histograms, as determined by input indices.
'''

import numpy as np
import matplotlib.pyplot as plt
import glob
import netCDF4 as netCDF
import tracpy
import tracpy.plotting
import matplotlib as mpl
import pdb
import op
from matplotlib import ticker, colors, cm
from matplotlib.mlab import find
import os


mpl.rcParams.update({'font.size': 14})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
currents_filename = list(np.sort(glob.glob('/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_????.nc')))
if 'grid' not in locals():
    grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)

howplot = 'log'

cmap = 'YlGn'

shelf_depth = 100

fname = 'calcs/S.npz'

# Simulations to plot
Files = glob.glob('tracks/2010-*gc.nc')

S = np.empty((len(Files), grid['xpsi'].shape[0], grid['xpsi'].shape[1])) # initialize

for i, File in enumerate(Files):

    datestr = File.split('/')[-1].split('gc')[0]

    fig = plt.figure(figsize=(6.8375, 6.6125))
    fig.subplots_adjust(left=0.04, bottom=0.15, right=1.0, top=0.96, wspace=0.07, hspace=0.04)
    ax = fig.add_subplot(111)
    tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
    ax.set_title('Simulation starting ' + datestr)

    if not os.path.exists(fname):

        d = netCDF.Dataset(File)
        U = d.variables['U'][:]; V = d.variables['V'][:]
        d.close()
        S[i,:,:] = np.sqrt(op.resize(U, 1)**2 + op.resize(V, 0)**2)
        # S[i,:,:] = S[i,:,:] + Stemp

        # locator = ticker.MaxNLocator(11)
        # locator.create_dummy_axis()
        # locator.set_bounds(0, 1) 
        # levels = locator()
        # extend = 'max'
        # H = H/Hmax
        Smax = 1.

    else:
        d = np.load(fname); S = d['S']; d.close()
        Smax = S.max()

    if howplot=='log':
        lev_exp = np.linspace(np.log10(0.00001),np.log10(1.0),100)
        levs = np.power(10, lev_exp)
    elif howplot=='linear':
        levs = np.linspace(0,1,100)
    mappable = ax.contourf(grid['xpsi'], grid['ypsi'], S[i,:,:]/Smax, cmap=cmap, levels=levs, norm=colors.LogNorm())#, extend=extend)

    # lev_exp = np.linspace(1, np.ceil(np.log10(S.max())-1),50)
    # # lev_exp = np.arange(1e-13, np.ceil(np.log10(S.max())+1))
    # levs = np.power(10, lev_exp)

    # mappable = ax.contourf(grid['xpsi'], grid['ypsi'], S[i,:,:], cmap=cmap, levels=levs, norm=colors.LogNorm())#, extend=extend)
    # # mappable = ax.contourf(grid['xpsi'], grid['ypsi'], S[i,:,:]/Smax, cmap=cmap, levels=np.linspace(0,0.25), extend='max')
    # ax.contour(grid['xr'], grid['yr'], grid['h'], [shelf_depth], colors='0.1', linewidth=3)
    # # pdb.set_trace()

    # plot starting drifter locations for this region 
    loc_shelftransport = '/home/kthyng/projects/shelf_transport/'
    dconn = np.load(loc_shelftransport + 'calcs/galvestonpts.npz')
    lon = dconn['lon']; lat = dconn['lat']
    xp, yp = grid['basemap'](lon,lat)
    dconn.close()
    ax.plot(xp, yp, 'k', alpha=0.5, lw=2)

    # Plot DWH site
    xdwh, ydwh = grid['basemap'](-88.386944, 28.736667)
    ax.plot(xdwh, ydwh, 'ok', ms=15)

    # Horizontal colorbar below plot
    if i == (len(Files)-1):
        cax = fig.add_axes([0.4, 0.25, 0.5, 0.02]) #colorbar axes
        cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
        cb.set_label('Backward transport')
        levscb = np.hstack((levs[::20],levs[-1]))
        cb.set_ticks(levscb)
        cb.set_ticklabels(["%1.4f" % lev for lev in levscb])
        cb.ax.tick_params(labelsize=10)

    fig.savefig('figures/transport/' + datestr + '-log.png', bbox_inches='tight')
    plt.close()

np.savez(fname, S=S, x=grid['xpsi'], y=grid['ypsi'])


# Plot all combined
Ssum = S.sum(axis=0)
Ssummax = Ssum.max()

fig = plt.figure(figsize=(6.8375, 6.6125))
fig.subplots_adjust(left=0.04, bottom=0.15, right=1.0, top=0.96, wspace=0.07, hspace=0.04)
ax = fig.add_subplot(111)
tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
ax.set_title('All simulations')

if howplot=='log':
    lev_exp = np.linspace(np.log10(0.00001),np.log10(1.0),100)
    levs = np.power(10, lev_exp)
elif howplot=='linear':
    levs = np.linspace(0,1,100)
mappable = ax.contourf(grid['xpsi'], grid['ypsi'], Ssum/Ssummax, cmap=cmap, levels=levs, norm=colors.LogNorm())#, extend=extend)

# plot starting drifter locations for this region 
ax.plot(xp, yp, 'k', alpha=0.5, lw=2)

# Plot DWH site
xdwh, ydwh = grid['basemap'](-88.386944, 28.736667)
ax.plot(xdwh, ydwh, '^k', ms=15)
ax.plot(486538, 221471, '^k', ms=15)
ax.text(0.475, 0.2175, 'Deepwater Horizon site', transform=ax.transAxes)

# Horizontal colorbar below plot
if i == (len(Files)-1):
    cax = fig.add_axes([0.4, 0.25, 0.5, 0.02]) #colorbar axes
    cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label('Backward transport')
    levscb = np.hstack((levs[::20],levs[-1]))
    cb.set_ticks(levscb)
    cb.set_ticklabels(["%1.4f" % lev for lev in levscb])
    cb.ax.tick_params(labelsize=10)

fig.savefig('figures/transport/all-log.png', bbox_inches='tight')
