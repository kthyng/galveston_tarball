'''
Compare winds year to year for DWH time period in TX.
Just plot June and July based on results shown in paper.
'''

# pip install windrose
import pandas as pd
import tabs
import seaborn as sns
import numpy as np
from windrose import WindroseAxes
import cmocean.cm as cmo
import os


# Read in data
# 42035 missing 2013
# 8766072 missing 2004-2008
# N missing 2010
# V missing 2004, 2006, 2012
fname = 'V.csv'
if os.path.exists(fname):
    Vall = pd.read_csv(fname, index_col=0, parse_dates=True)
else:
    Vall = tabs.read('V', '2004', '2015-1-1')
    Vall.to_csv('V.csv')

fname = 'N.csv'
if os.path.exists(fname):
    Nall = pd.read_csv(fname, index_col=0, parse_dates=True)
else:
    Nall = tabs.read('N', '2004', '2015-1-1')
    Nall.to_csv('N.csv')


# only take subset of months each year
# V = Vall[(Vall.index.month == 5) + (Vall.index.month == 6) + (Vall.index.month == 7) + (Vall.index.month == 8)]
# N = Nall[(Nall.index.month == 5) + (Nall.index.month == 6) + (Nall.index.month == 7) + (Nall.index.month == 8)]

# V = Vall[(Vall.index.month == 7) + (Vall.index.month == 8) + (Vall.index.month == 9)]
# N = Nall[(Nall.index.month == 7) + (Nall.index.month == 8) + (Nall.index.month == 9)]

# V = Vall[(Vall.index.month == 6) + (Vall.index.month == 7)]
# N = Nall[(Nall.index.month == 6) + (Nall.index.month == 7)]

V = Vall[(Vall.index.month == 7) + (Vall.index.month == 8)]
N = Nall[(Nall.index.month == 7) + (Nall.index.month == 8)]
# V [2005, 2007, 2010, 2011, 2013, 2014]
# N [2004, 2005, 2006, 2007, 2008, 2009, 2011, 2012, 2013, 2014]
# Use N but for 2010 use V


# df = N.copy(); buoy = 'N'

# Make a subplot for each available year of data
# Use N and V together to get years that match connectivity plot
years = []; rmax = 0; fmax = 0; dfs = []; buoys = []
for year in np.arange(2004, 2015):
    if year in [2004, 2005, 2006, 2007, 2008, 2009, 2011, 2012, 2013, 2014]:
        df = N.copy(); buoy = 'N'
    elif year == 2010:
        df = V.copy(); buoy = 'V'
    dfs.append(df); buoys.append(buoy)
    dfyear = df[df.index.year == year]
    theta = dfyear[buoy + ': Dir from [deg T] (wind)'].copy()
    r = dfyear[buoy + ': Speed [m/s] (wind)'].copy()
    # at least half of data should be present and not nan
    if theta.count() > theta.size/2:
        # save years where this is true
        years.append(year)
        ind = ~(np.isnan(theta) | np.isnan(r))
        ax = WindroseAxes.from_ax()
        # normed=True does frequency as a percentage
        ax.bar(theta[ind].values, r[ind].values, opening=0.8, normed=True, edgecolor='white', cmap=cmo.speed)
        rmax = max(rmax, r.max())  # max speed
        table = ax._info['table']
        wd_freq = np.sum(table, axis=0)
        fmax = max(fmax, wd_freq.max())  # max frequency
        plt.close(plt.gcf())

fmax = np.round(fmax, -1)
fmax = 40
step = 10
rmax = np.round(rmax, 0)
rmax = 12

# len(years) == 12
fig, axes = plt.subplots(4,3, figsize=(6,9),
                         subplot_kw=dict(projection="windrose"))
fig.subplots_adjust(left=0, right=1, hspace=0.04, wspace=0.0, bottom=0.04, top=0.95)
for i, (ax, year, df, buoy) in enumerate(zip(axes.flatten(), years, dfs, buoys)):
    dfyear = df[df.index.year == year]
    theta = dfyear[buoy + ': Dir from [deg T] (wind)'].copy()
    r = dfyear[buoy + ': Speed [m/s] (wind)'].copy()
    # at least half of data should be present and not nan
    if theta.count() > theta.size/2:

        ind = ~(np.isnan(theta) | np.isnan(r))
        # print(year, theta.count(), r.count())
        # WindroseAxes.from_ax(ax=ax, fig=fig)
        ax.bar(theta[ind].values, r[ind].values, normed=True, opening=0.7,
               edgecolor='0.5', cmap=cmo.matter, bins=np.arange(0, rmax, 3),
               nsector=12)
        ax.set_yticks(np.arange(10, fmax, step=step))
        ax.set_yticklabels(np.arange(10, fmax, step=step, dtype=int))
        ax.text(0.15, 0.2, year, transform=ax.transAxes, fontsize=14,
                bbox=dict(facecolor='w', alpha=0.75, edgecolor='w'))
        ax.set_ylim(0,fmax)
        # turn off compass labels except first subplot
        if i > 0:
            ax.set_xticklabels('')
ax.set_legend(loc=(1.25,0.25))

axes[-1,-1].axis('off')  # turn off last subplot
axes[0,0].set_zorder(10)  # set zorder for first subplot to be on top

# # clean up final subplot
# axes[-1,-1].set_yticks(np.arange(10, fmax, step=step))
# axes[-1,-1].set_yticklabels(np.arange(10, fmax, step=step, dtype=int))
# axes[-1,-1].set_ylim(0,fmax)
#

# fig.savefig('figures/winds.png', bbox_inches='tight', dpi=100)
fig.savefig('figures/winds.pdf', bbox_inches='tight', dpi=100)
# fig.savefig('figures/winds_high.png', bbox_inches='tight', dpi=300)

# TODOS:
# be clear about direction to/from
# Use may wind also. August? Comparing wind w oil locations first and bring in river from shelf transport paper. Then additionally connectivity plot. Add to methods and results. Compare n and v winds when both available to make sure consistent
# Tabs package should label non tabs buoy speed and fit as currents or wind too. Ports are currents

# #
# ax = plt.subplot(111, polar=True)
# # ax.set_ylim(0,0.01)
#
# ax.quiver(0, 0, np.cos(theta) - np.sin(theta), np.sin(theta) + np.cos(theta))
#
# ax.quiver(0,0, theta, gb['V: Speed [cm/s]'])
#
# ax.scatter(x=gb.mean()['V: Speed [cm/s]'].values, y=gb.mean()['V: Dir [deg T]'].values)
# ax.set_theta_zero_location('N')
# ax.set_theta_direction(-1)
