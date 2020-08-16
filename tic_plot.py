# -*- coding: utf-8 -*-
"""
Helpers to plot the lightcurve of a TESS subject, given a
LigthCurveFileCollection
"""

import matplotlib.pyplot as plt
import matplotlib as matplotlib
from matplotlib.ticker import (FormatStrFormatter, AutoMinorLocator)
import numpy as np
import pandas as pd

def lcf_of_sector(lcf_coll, sectorNum):
    for lcf in lcf_coll:
        if lcf.header()['SECTOR'] == sectorNum:
            return lcf
    return None

def lcfs_of_sectors(*args):
    lcf_coll = args[0]
    sectorNums = args[1:]
    res = []
    for lcf in lcf_coll:
        if lcf.header()['SECTOR'] in sectorNums:
            res.append(lcf)
    return res

# Plot the flux changes (not flux themselves) to get a sense of the rate of changes, not too helpful yet.
def plot_lcf_flux_delta(lcf, ax, xmin=None, xmax=None, moving_avg_window='30min'):

    # possible input arguments
    lc = lcf.PDCSAP_FLUX.normalize(unit='percent')

    # Basic scatter of the observation
#    ax = lc.scatter(ax=ax)

    # convert to dataframe to add moving average
    df = lc.to_pandas(columns=['time', 'flux'])
    df['time_ts'] = df['time'].apply(lambda x: pd.Timestamp(x, unit='D'))
    # the timestamp above is good for relative time.
    # if we want the timestamp to reflect the actual time, we need to convert the BTJD in time to timetamp, e.g.
    # df['time_ts'] = df['time'].apply(lambda x: pd.Timestamp(astropy.time.Time(x + 2457000, format='jd', scale='tdb').datetime.timestamp(), unit='s'))
    df['flux_mavg'] = df.rolling(moving_avg_window, on='time_ts')['flux'].mean()
#    ax.plot(lc.time, df['flux_mavg'], c='black', label=f"Moving average ({moving_avg_window})")

    df['flux_delta'] = df.rolling(moving_avg_window, on='time_ts')['flux_mavg'].apply(lambda vary: vary[-1] - vary[0], raw=True)
    ax.plot(lc.time, df['flux_delta'], c='blue', label=f"Flux delta ({moving_avg_window})")

    ax.set_xlim(xmin, xmax)

    return ax

def lcf_fig():
    return plt.figure(figsize=(15, 6))

def flux_near(lc, time):
    if time is None or lc is None:
        return None
    else:
        idx = (np.abs(lc.time - time)).argmin()
        return lc.flux[idx]

def flux_mavg_near(df, time):
    if time is None or df is None:
        return None
    else:
        idx = (np.abs(df['time'].values - time)).argmin()
        # must use df.iloc[idx]['flux_mavg'], rather than df['flux_mavg'][idx]
        # because dataframe from lightkurve is indexed by time (rather than regular 0-based index)
        # df.iloc[] ensures we can still access the value by 0-based index
        return df.iloc[idx]['flux_mavg']

def as_4decimal(float_num):
    if (float_num is None):
        return None
    elif (isinstance(float_num, tuple) or isinstance(float_num, list)):
        return [float('{0:.4f}'.format(n)) for n in float_num]
    else:
        return float('{0:.4f}'.format(float_num))

def add_flux_moving_average(lc, moving_avg_window):
    df = lc.to_pandas(columns=['time', 'flux'])
    df['time_ts'] = df['time'].apply(lambda x: pd.Timestamp(x, unit='D'))
    # the timestamp above is good for relative time.
    # if we want the timestamp to reflect the actual time, we need to convert the BTJD in time to timetamp, e.g.
    # df['time_ts'] = df['time'].apply(lambda x: pd.Timestamp(astropy.time.Time(x + 2457000, format='jd', scale='tdb').datetime.timestamp(), unit='s'))
    df['flux_mavg'] = df.rolling(moving_avg_window, on='time_ts')['flux'].mean()
    return df

def add_relative_time(lc, lcf):
    t_start = lcf.header()['TSTART']
    lc.time_rel = lc.time - t_start
    return lc.time_rel

def mask_gap(x, y, min_x_diff):
    """
    Help to plot graphs with gaps in the data, so that straight line won't be draw to fill the gap.
    Return a masked y that can be passed to pyplot.plot() that can show the gap.
    """
    x_diff = np.diff(x, prepend=-min_x_diff)
    return np.ma.masked_where(x_diff > min_x_diff, y)

def plot_n_annotate_lcf(lcf, ax, xmin=None, xmax=None, t0=None, t_start=None, t_end=None, moving_avg_window='30min', t0mark_ymax = 0.3, lc_tweak_fn=None, ax_tweak_fn=None):
    if lcf == None:
        print("Warning: lcf is None. Plot skipped")
        return

    # possible input arguments

    lc = lcf.PDCSAP_FLUX.normalize(unit='percent')

    if lc_tweak_fn is not None:
        lc = lc_tweak_fn(lc)

    # Basic scatter of the observation
    ax = lc.scatter(ax=ax)

    # convert to dataframe to add moving average
    if moving_avg_window is not None:
        df = add_flux_moving_average(lc, moving_avg_window)
        ax.plot(lc.time, df['flux_mavg'], c='black', label=f"Moving average ({moving_avg_window})")
    else:
        df = add_flux_moving_average(lc, '2min') # still needed for some subsequent calc, but don't plot it


    # annotate the graph
    lcfh = lcf.header()
    if xmin is None and t_start is not None:
        xmin = t_start - 0.5
    if xmax is None and t_end is not None:
        xmax = t_end + 0.5
    ax.set_xlim(xmin, xmax)
    if t_start is not None:
        ax.axvline(t_start)
    if t_end is not None:
        ax.axvline(t_end)
    if t0 is not None:
        ax.axvline(t0, ymin=0, ymax=t0mark_ymax, color='black', linewidth=3, linestyle='--', label=f"t0 ~= {t0}")

    transit_duration_msg = ''
    if t_start is not None and t_end is not None:
        transit_duration_msg = f'\ntransit duration ~= {as_4decimal(24 * (t_end - t_start))}h'
    flux_t0 = flux_mavg_near(df, t0)
    flux_dip = None
    if flux_t0 is not None:
        flux_begin = max(flux_mavg_near(df, t_start), flux_mavg_near(df, t_end))
        flux_dip = flux_begin - flux_t0
    ax.set_title(f"{lc.label}, sector {lcfh['SECTOR']} \nflux@t0 ~= {as_4decimal(flux_t0)}%, dip ~= {as_4decimal(flux_dip)}%{transit_duration_msg}", {'fontsize': 24})
    ax.legend()
    ax.xaxis.label.set_size(18)
    ax.yaxis.label.set_size(18)

    # to avoid occasional formating in scentific notations
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    if ax_tweak_fn is not None:
        ax_tweak_fn(ax)

    return ax

# Do the actual plots
def plot_all(lcf_coll, moving_avg_window=None, lc_tweak_fn=None, ax_fn=None, use_relative_time=False, ax_tweak_fn=None):
    # choice 1: use the built-in plot method
#    ax_all = plt.figure(figsize=(30, 15)).gca()
#     lcf_coll.PDCSAP_FLUX.plot(ax=ax_all) # Or lcf_coll.SAP_FLUX.plot()

    # choice 2: stitch lightcurves of the collection together, and then use more flexible methods, e.g., scatter
    #   Note: pass lambda x: x to stitch() so that the code won't normalize the flux value sector by sector
#     lc_all = lcf_coll.PDCSAP_FLUX.stitch(lambda x: x)
#     lc_all.scatter(ax=ax_all, normalize=True)

    # choice 3: plot the lightcurve sector by sector: each sector has its own color
#     for i in range(0, len(lcf_coll)):
#         lcf_coll[i].PDCSAP_FLUX.scatter(ax=ax_all)

#     ax_all.set_title(f"TIC {lcf_coll[0].PDCSAP_FLUX.label}, sectors {list(map(lambda lcf: lcf.header()['SECTOR'], lcf_coll))}")
#     return ax_all

    # choice 4: plot the lightcurve sector by sector: each sector in its own graph
    for i in range(0, len(lcf_coll)):
        if ax_fn is None:
            ax = lcf_fig().gca()
        else:
            ax = ax_fn()
        lcf = lcf_coll[i]
        lc = lcf.PDCSAP_FLUX
        lc = lc.normalize(unit='percent')
        if lc_tweak_fn is not None:
            lc = lc_tweak_fn(lc)

        # temporarily change time to a relative one if specified
        if use_relative_time:
            add_relative_time(lc, lcf)
            lc.time_orig = lc.time
            lc.time = lc.time_rel

        lc.scatter(ax=ax)

        # convert to dataframe to add moving average
        if moving_avg_window is not None:
            df = add_flux_moving_average(lc, moving_avg_window)
            # mask_gap: if there is a gap larger than 2 hours,
            # show the gap rather than trying to fill the gap with a straight line.
            ax.plot(lc.time, mask_gap(lc.time, df['flux_mavg'], 2/24), c='black', label=f"Moving average ({moving_avg_window})")

        title_extras = ''
        if lc_tweak_fn is not None:
            title_extras = '\nLC tweaked, e.g., outliers removed'

        ax.set_title(f"{lcf_coll[0].PDCSAP_FLUX.label}, sectors {lcf_coll[i].header()['SECTOR']}{title_extras}", {'fontsize': 36})
#        ax.set_title(f"{lcf_coll[0].PDCSAP_FLUX.label}, sectors N/A - Kepler")
#         ax.legend()
        if use_relative_time:
            ax.xaxis.set_label_text('Time - relative')
            # restore original time after plot is done
            lc.time = lc.time_orig

        # to avoid occasional formating in scentific notations
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        # minor tick, 1 day interval in practice
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis='x', which='minor', length=4)
        # ax.xaxis.grid(True, which='minor') # too noisy to be there by default

        ax.xaxis.label.set_size(fontsize=18)
        ax.yaxis.label.set_size(fontsize=18)
        if ax_tweak_fn is not None:
            ax_tweak_fn(ax)
    return None


def scatter_centroids(lcf, fig=None, highlight_time_range=None, time_range=None):
    '''
    Scatter centroids, and highlight the specific time range
    '''

    if fig is None:
        fig = plt.figure(figsize=(12,12))

    lc = lcf.PDCSAP_FLUX.normalize(unit='percent')
    sector = lcf.header()['SECTOR']

    df = lc.to_pandas(columns=['time', 'flux', 'centroid_row', 'centroid_col'])
    if time_range is not None:
        df = df[(df.time >= time_range[0]) & (df.time <= time_range[1])]
        if (len(df) < 1):
            raise Exception(f'Zoomed lightcurve has no observation. time_range={time_range}')

    fig.gca().yaxis.set_major_formatter(FormatStrFormatter('%.3f')) # avoid scientific notations
    fig.gca().scatter(df.centroid_col.values, df.centroid_row.values, c='blue', label=f'TIC {lc.targetid}')

    if highlight_time_range is not None:
        df_highlight = df[(df.time >= highlight_time_range[0]) & (df.time <= highlight_time_range[1])]
        if (len(df_highlight) < 1):
            print('WARNING: scatter_centroids() no observations in highlight_time_range')
        fig.gca().scatter(df_highlight.centroid_col.values, df_highlight.centroid_row.values, c='red', label='highlights')

    title = f'TIC {lc.targetid} Centroids, sector {sector}'
    if time_range is not None:
        title += f'\n{as_4decimal(time_range)}'
    if highlight_time_range is not None:
        title += f'\nHighlights:{as_4decimal(highlight_time_range)}'
    fig.gca().set_title(title)
    fig.legend()
    return fig


import matplotlib.animation as animation

def _update_anim(n, ax, lc, label, num_centroids_to_show, use_relative_time, c):
    ax.cla()
    # fix the x/y scale to ensure it doesn't change over the animation
    ax.set_xlim(np.nanmin(lc.centroid_col), np.nanmax(lc.centroid_col))
    ax.set_ylim(np.nanmin(lc.centroid_row), np.nanmax(lc.centroid_row))

    # avoid scientific notation for y-axis
    # x-axis might need scentific notation so that the labels won't get too cramped with long decimals
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    if use_relative_time:
        time_label = f"{as_4decimal(lc.time_rel[n])} ({as_4decimal(lc.time[n])})"
    else:
        time_label = f"{as_4decimal(lc.time[n])}"


    if num_centroids_to_show is None:
        col = lc.centroid_col[:n]
        row = lc.centroid_row[:n]
    else:
        n_start = max(0, n - num_centroids_to_show)
        col = lc.centroid_col[n_start:n]
        row = lc.centroid_row[n_start:n]

    ax.set_title(f'TIC {lc.targetid} Centroids, {label}\nday: {time_label}')
    ax.scatter(col, row, c=c)

def animate_centroids(lcf, fig=None, frames=None, num_obs_per_frame=240, interval=250, use_relative_time=False, time_range=None, accumulative=True, c=None, display=True):
    '''
    Animate centroids to visualize changes over time.

    '''
    lc = lcf.PDCSAP_FLUX
    label = f"sector {lcf.header()['SECTOR']}"

    # Zoom to a particular time range if specified
    if time_range is not None:
        # use pandas to zoom to a particular time_range
        df = lcf.PDCSAP_FLUX.normalize(unit='percent').to_pandas(columns=['time', 'flux', 'centroid_row', 'centroid_col'])
        df = df[(df.time >= time_range[0]) & (df.time <= time_range[1])]
        if (len(df) < 1):
            raise Exception(f'Zoomed lightcurve has no observation. time_range={time_range}')

        lc_z = lambda: None # zoomed-in lightcurve-like object for the purpose of animation
        setattr(lc_z, 'time', df.time.values)
        setattr(lc_z, 'flux', df.flux.values)
        setattr(lc_z, 'centroid_row', df.centroid_row.values)
        setattr(lc_z, 'centroid_col', df.centroid_col.values)
        setattr(lc_z, 'targetid', lc.targetid)
        lc = lc_z

    if fig is None:
        fig = plt.figure(figsize=(12,12))
    if frames is None:
        num_obs = len(lc.centroid_row)
        num_frames = int(num_obs / num_obs_per_frame) # default 240 is about every 8 hours, given 2-minute intervals
        ary_n = np.linspace(1, num_obs, num=num_frames, endpoint=False)
        ary_n[0] = np.ceil(ary_n[1] /2)
        ary_n = list(map(lambda n: int(n), ary_n))
    else:
        ary_n = frames
        num_obs_per_frame = frames[-1] - frames[-2] # assume the spacing of input is linear

    num_centroids_to_show = num_obs_per_frame
    if accumulative:
        num_centroids_to_show = None
#     print(f'Steps: {ary_n}')
    if use_relative_time:
        add_relative_time(lc, lcf)
    anim = animation.FuncAnimation(fig, _update_anim, frames=ary_n
                                   , fargs=(fig.gca(), lc, label, num_centroids_to_show, use_relative_time, c)
                                   , interval=interval, blit=False)
    if display:
        # for inline display in jupyter
        try:
            from IPython.display import HTML
            return HTML(anim.to_jshtml(default_mode='once'))
        except ImportError:
            print('WARNING: animate_centroids() - inline display not possible Not in IPython envrionment.')
            return anim

