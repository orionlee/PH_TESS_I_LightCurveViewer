# -*- coding: utf-8 -*-
"""
Helpers to plot the lightcurve of a TESS subject, given a
LightCurveCollection
"""

# so that TransitTimeSpec can be referenced in type annotation in the class itself
# see: https://stackoverflow.com/a/49872353
from __future__ import annotations


import inspect
import warnings
from pathlib import Path
import re
from types import SimpleNamespace

from memoization import cached
import xmltodict

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, AutoMinorLocator
import matplotlib.animation as animation
import numpy as np
import pandas as pd

from astropy.io import fits
from astropy import units as u
from astropy.table import Table
from astroquery.exceptions import NoResultsWarning
from astroquery.mast import Observations

import IPython
from IPython.display import display, HTML, Audio
from ipywidgets import interactive, interactive_output, fixed
import ipywidgets as widgets

from lightkurve import LightCurveCollection, LightkurveWarning
from lightkurve.utils import TessQualityFlags
from lightkurve_ext import of_sectors
import lightkurve_ext as lke


def parse_dvs_filename(filename):
    # e.g.: tess2020267090513-s0030-s0030-0000000142087638-01-00394_dvs.pdf
    match = re.match(r"^tess\d+-(s\d+-s\d+)-(\d+)-(\d+)-.+_dvs[.]pdf", filename)
    if not match:
        return {}
    sector_range, tic_id_padded, tce_num_padded = (
        match.group(1),
        match.group(2),
        match.group(3),
    )
    tic_id = re.sub(r"^0+", "", tic_id_padded)
    tce_num = re.sub(r"^0+", "", tce_num_padded)
    # sufficient to identify one for a given TIC, less visually busy
    tce_id_short = f"{sector_range}:TCE{tce_num}"

    # tce_id is the format used on ExoMAT, e.g,  TIC142087638S0030S0030TCE1
    tce_id = f"""TIC{tic_id}{re.sub("-", "", sector_range.upper())}TCE{tce_num}"""

    return dict(
        tce_id=tce_id,
        tce_id_short=tce_id_short,
        sector_range=sector_range,
        tic_id=tic_id,
        tce_num=tce_num,
    )


def parse_dvr_filename(filename):
    match = re.match(r"^tess\d+-(s\d+-s\d+)-(\d+)-.+_dvr[.](pdf|xml)", filename)
    if not match:
        return {}
    sector_range, tic_id_padded, file_type = (
        match.group(1),
        match.group(2),
        match.group(3),
    )
    tic_id = re.sub(r"^0+", "", tic_id_padded)

    return dict(sector_range=sector_range, tic_id=tic_id, file_type=file_type)


@cached
def get_dv_products_of_tic(tic_id, productSubGroupDescription, download_dir=None):
    # Based on:
    # - https://outerspace.stsci.edu/display/TESS/7.0+-+Tips+and+Tricks+to+Getting+TESS+Data+At+MAST
    # https://github.com/spacetelescope/notebooks/blob/master/notebooks/MAST/TESS/beginner_astroquery_dv/beginner_astroquery_dv.ipynb

    # Note: for TESS, tic_id (the number without TIC) is what an exact match works
    # Kepler / K2 ids will need some additional processing for exact match to work.
    exact_target_name = tic_id
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=NoResultsWarning, message=".*No products to download.*")
        obs_wanted = Observations.query_criteria(
            target_name=exact_target_name,
            dataproduct_type="timeseries",
            obs_collection="TESS",
        )
        data_products = Observations.get_product_list(obs_wanted)
        return Observations.filter_products(data_products, productSubGroupDescription=productSubGroupDescription)


@cached
def parse_dvr_xml(file_path):
    def as_list(data):
        """Wrap an item as a list, if it's not one.
        Useful for handling dict from XML where elements might be one or multiple elements"""
        if type(data) is list:
            return data
        else:
            return [data]

    def param_value(model_params_dict, param_name):
        param_dict = model_params_dict.get(param_name)
        if param_dict is None:
            return None
        val_str = param_dict.get("@value")
        if val_str is None:
            return None
        return float(val_str)

    # the body
    with open(file_path, "r") as f:
        dvr_xml_str = f.read()
    parsed = xmltodict.parse(dvr_xml_str)

    planets_dict = {}

    e_pr_list = as_list(parsed["dv:dvTargetResults"]["dv:planetResults"])
    for e_pr in e_pr_list:
        e_afit = e_pr["dv:allTransitsFit"]
        planet_num = e_afit["@planetNumber"]

        params_dict = {}  # a temporary structure to access params internally
        for mp in e_afit["dv:modelParameters"]["dv:modelParameter"]:
            params_dict[mp["@name"]] = mp

        # TODO: add other DV fitting parameters, odd/even test, centroid, etc.
        # use the underlying xml attribute names, even thought it breaks the convention
        a_planet_dict = dict(
            planetNumber=planet_num,
            transitEpochBtjd=param_value(params_dict, "transitEpochBtjd"),
            planetRadiusEarthRadii=param_value(params_dict, "planetRadiusEarthRadii"),
            transitDurationHours=param_value(params_dict, "transitDurationHours"),
            orbitalPeriodDays=param_value(params_dict, "orbitalPeriodDays"),
            transitDepthPpm=param_value(params_dict, "transitDepthPpm"),
            minImpactParameter=param_value(params_dict, "minImpactParameter"),
        )

        planets_dict[planet_num] = a_planet_dict

    return planets_dict


def get_tce_infos_of_tic(tic_id, download_dir=None):
    def filter_by_dataURI_suffix(products, suffix):
        # Helper to filter products into summary, full report, full report xml using suffix.
        # It replaces the logic to filter by "description" column, as description is sometimes unreliable
        # E.g., for the TCE for TIC 43843023 sector 5, the dvr xml has incorrect description
        # so that the entry is treated as a dvr pdf
        return products[np.char.endswith(products["dataURI"], suffix)]

    products_wanted = get_dv_products_of_tic(tic_id, ["DVS", "DVR"], download_dir=download_dir)

    res = []
    # basic info
    for p in filter_by_dataURI_suffix(products_wanted, "_dvs.pdf"):
        tce_info = parse_dvs_filename(p["productFilename"])
        entry = dict(
            obsID=p["obsID"],
            tic_id=tce_info.get("tic_id"),
            sector_range=tce_info.get("sector_range"),
            tce_num=tce_info.get("tce_num"),
            tce_id=tce_info.get("tce_id"),
            tce_id_short=tce_info.get("tce_id_short"),
            dvs_dataURI=p["dataURI"],
        )
        res.append(entry)

    # DVR pdf link
    for p in filter_by_dataURI_suffix(products_wanted, "_dvr.pdf"):
        # find TCEs for the same observation (sometimes there are multiple TCEs for the same observation)
        for entry in [e for e in res if e["obsID"] == p["obsID"]]:
            entry["dvr_dataURI"] = p["dataURI"]

    products_dvr_xml = filter_by_dataURI_suffix(products_wanted, "_dvr.xml")
    manifest = Observations.download_products(products_dvr_xml, download_dir=download_dir)
    if manifest is None:
        return res
    for m in manifest:
        dvr_xml_local_path = m["Local Path"]

        dvr_info = parse_dvr_filename(Path(dvr_xml_local_path).name)
        for entry in [e for e in res if e["tic_id"] == dvr_info["tic_id"] and e["sector_range"] == dvr_info["sector_range"]]:
            entry["dvr_xml_local_path"] = dvr_xml_local_path

        planets_dict = parse_dvr_xml(dvr_xml_local_path)
        for a_planet_dict in planets_dict.values():
            for entry in [
                e
                for e in res
                if e["tic_id"] == dvr_info["tic_id"]
                and e["sector_range"] == dvr_info["sector_range"]
                and e["tce_num"] == a_planet_dict["planetNumber"]
            ]:
                entry["planet"] = a_planet_dict

    return res


def get_tic_meta_in_html(lc, a_subject_id=None, download_dir=None):
    # This function does not do the actual display,
    # so that the caller can call it in background
    # and display it whereever it's needed
    def link(link_text, url):
        return f"""<a href="{url}" target="_blank">{link_text}</a>"""

    def prop(prop_name, prop_value):
        return f"""    <tr><td>{prop_name}</td><td>{prop_value}</td></tr>\n"""

    def row(*args):
        return "<tr>" + "".join(f"<td>{v}</td>" for v in args) + "</tr>"

    # main logic
    m = lc.meta
    tic_id = str(m.get("TICID"))

    def safe_m_get(key, default_val):
        # in some meta, the key exists but the value is None
        # this helper handles it
        res = m.get(key, default_val)
        return res if res is not None else default_val

    html = f"""
<h3>TIC {tic_id}</h3>
"""
    html += "&emsp;" + link("ExoFOP", f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic_id}")
    html += "\n&emsp;|&emsp;"
    html += link(
        "PHT Talk",
        f"https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/talk/search?query={tic_id}",
    )
    if a_subject_id is not None:
        # note, a TIC can have multiple subjects, here is just one of them.
        html += "\n , a subject: "
        html += link(
            a_subject_id,
            f"https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/talk/subjects/{a_subject_id}",
        )
        # show the sector number (here we assume a_subject_id does correspond the the sector)
        # the sector is useful to be included so that users can easily locate the TCE matching the sector.
        html += f' (sector {safe_m_get("SECTOR", "")})'
    html += "<br>\n"
    html += "<table>\n"
    html += prop("R<sub>S</sub> (in R<sub>☉</sub>)", f'{safe_m_get("RADIUS", 0):.3f}')
    html += prop("Magnitude (TESS)", f'{safe_m_get("TESSMAG", 0):.2f}')
    html += prop("T_eff (in K)", safe_m_get("TEFF", 0))
    html += "</table>\n"

    # TODO: For TCE, query MAST download / parse results (the _dvr.xml), tho show
    # - basic planet parameters and orbital info
    # - red flags in vetting report
    # see: https://archive.stsci.edu/missions-and-data/tess/data-products
    tce_info_list = get_tce_infos_of_tic(tic_id, download_dir=download_dir)
    if len(tce_info_list) < 1:
        return html

    header = [
        ("TCE", ""),
        ("Reports", ""),
        ("R<sub>p</sub>", "R<sub>j</sub>"),
        ("Epoch", "BTJD"),
        ("Duration", "hr"),
        ("Period", "day"),
        ("Depth", "%"),
        ("Impact P.", "<i>b</i>"),
        ("Codes", ""),
    ]
    html += """<br>TCEs: <table>
<thead>"""
    html += "<tr>"
    html += " ".join([f"<th>{h[0]}</th>" for h in header])
    html += "</tr>\n"
    html += "<tr>"
    html += " ".join([f"<th>{h[1]}</th>" for h in header])
    html += "</tr>\n"
    html += """
</thead>
<tbody>
"""
    R_EARTH_TO_R_JUPITER = 6378.1 / 71492
    for info in tce_info_list:
        exomast_url = f'https://exo.mast.stsci.edu/exomast_planet.html?planet={info.get("tce_id")}'
        dvs_url = f'https://exo.mast.stsci.edu/api/v0.1/Download/file?uri={info.get("dvs_dataURI")}'
        dvr_url = f'https://exo.mast.stsci.edu/api/v0.1/Download/file?uri={info.get("dvr_dataURI")}'
        p_i = info.get("planet", {})
        html += row(
            link(info.get("tce_id_short"), exomast_url),
            f"""{link("dvs", dvs_url)},&emsp;{link("full", dvr_url)}""",
            f'{p_i.get("planetRadiusEarthRadii", 0) * R_EARTH_TO_R_JUPITER:.3f}',
            f'{p_i.get("transitEpochBtjd", 0):.4f}',
            f'{p_i.get("transitDurationHours", 0):.4f}',
            f'{p_i.get("orbitalPeriodDays", 0):.6f}',
            f'{p_i.get("transitDepthPpm", 0) / 10000:.4f}',
            f'{p_i.get("minImpactParameter", 0):.2f}',
            # code fragments to so that users can easily use a TCE as an entry in transit_specs
            f"""\
<input type="text" style="margin-left: 3ch; font-size: 90%; color: #666; width: 10ch;"
       value='epoch={p_i.get("transitEpochBtjd", 0):.4f}, duration_hr={p_i.get("transitDurationHours", 0):.4f}, \
period={p_i.get("orbitalPeriodDays", 0):.6f}, label="{info.get("tce_id_short")}",'>""",
        )
        html += "\n"

    html += "</tbody></table>\n"

    # TODO: check if there is a TOI?!

    return html


def beep():
    """Emits a beep sound. It works only in IPython / Jupyter environment only"""
    # a beep to remind the users that the data has been downloaded
    # css tweak to hide beep
    display(
        HTML(
            """<script>
function tweakCSS() {
  if (document.getElementById("hide-beep-css")) {
      return;
  }
  document.head.insertAdjacentHTML('beforeend', `<style id="hide-beep-css" type="text/css">
  #beep { /* hide the audio control for the beep, generated from tplt.beep() */
    width: 1px;
    height: 1px;
  }
</style>`);
}
tweakCSS();
</script>
"""
        )
    )
    # the actual beep
    beep_url = "https://upload.wikimedia.org/wikipedia/commons/f/fb/NEC_PC-9801VX_ITF_beep_sound.ogg"
    if int(re.sub(r"[.].+", "", IPython.__version__)) < 7:
        # compatibility with older older IPython (e.g., google colab)
        audio = Audio(url=beep_url, autoplay=True, embed=True)
    else:
        audio = Audio(url=beep_url, autoplay=True, embed=True, element_id="beep")
    display(audio)


def _normalize_to_percent_quiet(lc):
    # Some product are in normalized flux, e.g., as 1, we still want to normalize them to percentage
    # for consistency
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=LightkurveWarning, message=".*in relative units.*")
        return lc.normalize(unit="percent")


# Plot the flux changes (not flux themselves) to get a sense of the rate of changes, not too helpful yet.
def plot_lcf_flux_delta(lcf, ax, xmin=None, xmax=None, moving_avg_window="30min"):

    # possible input arguments
    lc = _normalize_to_percent_quiet(lcf)

    # Basic scatter of the observation
    #    ax = lc.scatter(ax=ax)

    # convert to dataframe to add moving average
    df = lc.to_pandas()
    df["time_ts"] = [pd.Timestamp(x, unit="D") for x in df.index]
    # the timestamp above is good for relative time.
    # if we want the timestamp to reflect the actual time, we need to convert the BTJD in time to timetamp, e.g.
    # pd.Timestamp(astropy.time.Time(x + 2457000, format='jd', scale='tdb').datetime.timestamp(), unit='s')
    df["flux_mavg"] = df.rolling(moving_avg_window, on="time_ts")["flux"].mean()
    #    ax.plot(lc.time.value, df['flux_mavg'], c='black', label=f"Moving average ({moving_avg_window})")

    df["flux_delta"] = df.rolling(moving_avg_window, on="time_ts")["flux_mavg"].apply(
        lambda vary: vary[-1] - vary[0], raw=True
    )
    ax.plot(
        lc.time.value,
        df["flux_delta"],
        c="blue",
        label=f"Flux delta ({moving_avg_window})",
    )

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
        idx = (np.abs(df.index.values - time)).argmin()
        # must use df.iloc[idx]['flux_mavg'], rather than df['flux_mavg'][idx]
        # because dataframe from lightkurve is indexed by time (rather than regular 0-based index)
        # df.iloc[] ensures we can still access the value by 0-based index
        return df.iloc[idx]["flux_mavg"]


def _to_unitless(n):
    if hasattr(n, "value"):
        return n.value
    else:
        return n


def as_4decimal(float_num):
    if float_num is None:
        return None
    elif isinstance(float_num, tuple) or isinstance(float_num, list):
        return [float("{0:.4f}".format(_to_unitless(n))) for n in float_num]
    else:
        return float("{0:.4f}".format(_to_unitless(float_num)))


def add_flux_moving_average(lc, moving_avg_window):
    df = lc.to_pandas()
    begin_t = df.index[0]
    df["time_ts"] = [pd.Timestamp(t - begin_t, unit="D") for t in df.index]
    # the timestamp above is good for relative time.
    # 1. we subtract the time with the timestamp because for some products, e.g., CDIPS, the time value itself
    #    is so large that creating pd.Timestamp with it causes Overflow error
    # 2. if we want the timestamp to reflect the actual time, we need to convert the BTJD in time to timetamp, e.g.
    #      pd.Timestamp(astropy.time.Time(x + 2457000, format='jd', scale='tdb').datetime.timestamp(), unit='s')
    df["flux_mavg"] = df.rolling(moving_avg_window, on="time_ts")["flux"].mean()
    return df


def add_relative_time(lc, lcf):
    t_start = lcf.meta.get("TSTART")
    if t_start is None:
        return False
    lc["time_rel"] = lc.time - t_start
    return True


def mask_gap(x, y, min_x_diff):
    """
    Help to plot graphs with gaps in the data, so that straight line won't be draw to fill the gap.
    Return a masked y that can be passed to pyplot.plot() that can show the gap.
    """
    # handle case that x is a astropy Time object, rather than simple float array
    x = _to_unitless(x)

    x_diff = np.diff(x, prepend=-min_x_diff)
    return np.ma.masked_where(x_diff > min_x_diff, y)


def normalize_percent(lc):
    """
    A syntactic surgar for lambda for normalize as percentage.
    Useful when calling ``lc.fold()``, ``tpf.interact()``, etc.
    """
    return lc.normalize(unit="percent")


def _to_lc_with_flux(lc, flux_col):
    """Return a Lightcurve object with the named column as the flux column"""

    # analogous lkv1's way: lc = getattr(lcf, flux_col)
    res = lc.copy()
    res["flux"] = lc[flux_col.lower()]  # e.g., PDCSAP_FLUX (how we do in lkv1) will be lowerecased
    return res


_cache_plot_n_annotate_lcf = dict(lcf=None, flux_col=None, lc=None)


def plot_n_annotate_lcf(
    lcf,
    ax,
    flux_col="flux",
    xmin=None,
    xmax=None,
    t0=None,
    t_start=None,
    t_end=None,
    moving_avg_window="30min",
    t0mark_ymax=0.3,
    set_title=True,
    show_r_obj_estimate=True,
    title_fontsize=18,
    lc_tweak_fn=None,
    ax_tweak_fn=None,
):
    if lcf is None:
        print("Warning: lcf is None. Plot skipped")
        return

    # cache lc to speed up plots repeatedly over the same lcf
    global _cache_plot_n_annotate_lcf
    if lcf is _cache_plot_n_annotate_lcf["lcf"] and flux_col == _cache_plot_n_annotate_lcf["flux_col"]:
        lc = _cache_plot_n_annotate_lcf["lc"]
    else:
        lc = _normalize_to_percent_quiet(_to_lc_with_flux(lcf, flux_col))
        _cache_plot_n_annotate_lcf["lcf"] = lcf
        _cache_plot_n_annotate_lcf["flux_col"] = flux_col
        _cache_plot_n_annotate_lcf["lc"] = lc

    if lc_tweak_fn is not None:
        lc = lc_tweak_fn(lc)

    if xmin is None and t_start is not None:
        xmin = t_start - 0.5
    if xmax is None and t_end is not None:
        xmax = t_end + 0.5

    # implement xmin / xmax by limiting the LC itself, rather than using ax.set_xlim after the plot
    # - the Y-scale will then automatically scaled to the specified time range, rather than over entire lightcurve
    # - make plotting faster (fewer data points)
    if xmin is not None:
        lc = lc[lc.time.value >= xmin]
    if xmax is not None:
        lc = lc[lc.time.value <= xmax]

    lcfh = lcf.meta

    # Basic scatter of the observation
    if "long" == lke.estimate_cadence_type(lc):
        # long cadence has more spare data, use a larger "x" to represent them
        # "x" is also useful to distinguish it from moving average,
        # which will likely overlap with the points given the sparse data
        ax = lc.scatter(ax=ax, s=36, marker="x")
    else:
        ax = lc.scatter(ax=ax)

    if len(lc) < 1:
        print(
            (
                "Warning: specified (xmin, xmax) is out of the range of the lightcurve "
                f"{lc.label} sector {lcfh['SECTOR']}. Nothing to plot"
            )
        )
        return ax

    # convert to dataframe to add moving average
    if moving_avg_window is not None:
        df = add_flux_moving_average(lc, moving_avg_window)
        # mask_gap: if there is a gap larger than 2 hours,
        # show the gap rather than trying to fill the gap with a straight line.
        ax.plot(
            lc.time.value,
            mask_gap(lc.time, df["flux_mavg"], 2 / 24),
            c="black",
            label=f"Moving average ({moving_avg_window})",
        )
    else:
        df = add_flux_moving_average(lc, "10min")  # still needed for some subsequent calc, but don't plot it

    # annotate the graph
    if t_start is not None:
        ax.axvline(t_start)
    if t_end is not None:
        ax.axvline(t_end)
    if t0 is not None:
        t_lc_start = lcf.meta.get("TSTART", None)
        t0_rel_text = ""
        if t_lc_start is not None:
            t0_rel = t0 - t_lc_start
            t0_rel_text = f" ({as_4decimal(t0_rel)})"
        ax.axvline(
            t0,
            ymin=0,
            ymax=t0mark_ymax,
            color="black",
            linewidth=3,
            linestyle="--",
            label=f"t0 ~= {t0}{t0_rel_text}",
        )

    if set_title:
        title_text = f"{lc.label}, sector {lcfh['SECTOR']}"
        if lc.author is not None and lc.author != "SPOC":
            title_text += f", by {lc.author}"
        if t0 is not None:
            transit_duration_msg = ""
            if t_start is not None and t_end is not None:
                transit_duration_msg = f"\ntransit duration ~= {as_4decimal(24 * (t_end - t_start))}h"
            flux_t0 = flux_mavg_near(df, t0)
            if flux_t0 is not None:
                flux_begin = max(flux_mavg_near(df, t_start), flux_mavg_near(df, t_end))
                flux_dip = flux_begin - flux_t0
                r_obj_msg = ""
                r_obj = lke.estimate_object_radius_in_r_jupiter(lc, flux_dip / 100)  # convert flux_dip in percent to fractions
                if show_r_obj_estimate and r_obj is not None:
                    r_obj_msg = f", R_p ~= {r_obj:0.2f} R_j"
                title_text += (
                    f" \nflux@$t_0$ ~= {as_4decimal(flux_t0)}%, "
                    f"dip ~= {as_4decimal(flux_dip)}%{r_obj_msg}{transit_duration_msg}"
                )
        ax.set_title(title_text, {"fontsize": title_fontsize})
    ax.legend()
    ax.xaxis.label.set_size(18)
    ax.yaxis.label.set_size(18)

    # to avoid occasional formating in scientific notations
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="x", which="minor", length=4)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="y", which="minor", length=4)

    if ax_tweak_fn is not None:
        ax_tweak_fn(ax)

    return ax


def plot_transit(lcf, ax, t0, duration, surround_time, **kwargs):
    return plot_n_annotate_lcf(
        lcf,
        ax=ax,
        t0=t0 if duration > 0 else None,
        t_start=t0 - duration / 2 if duration > 0 else None,
        t_end=t0 + duration / 2 if duration > 0 else None,
        xmin=t0 - (duration + surround_time) / 2,
        xmax=t0 + (duration + surround_time) / 2,
        **kwargs,
    )


def plot_transits(lcf_coll, transit_specs, ax_fn=lambda: lcf_fig().gca(), **kwargs):
    """Helper to plot transits zoomed-in."""
    axs = []
    for spec in transit_specs:
        for lcf in of_sectors(lcf_coll, spec["sector"]):  # in case we have multiple lcf per sector
            #  process the supplied spec and apply defaults
            t0 = spec["epoch"]
            duration = spec["duration_hr"] / 24
            period = spec["period"]
            steps_to_show = spec["steps_to_show"]
            surround_time = spec.get("surround_time", 1.5)  # a hardcoded last resort default

            # TODO: warn if period is 0, but steps to show is not [0]

            for i in steps_to_show:
                cur_t0 = t0 + period * i
                ax = plot_n_annotate_lcf(
                    lcf,
                    ax=ax_fn(),
                    t0=cur_t0,
                    t_start=cur_t0 - duration / 2,
                    t_end=cur_t0 + duration / 2,
                    xmin=cur_t0 - (duration + surround_time) / 2,
                    xmax=cur_t0 + (duration + surround_time) / 2,
                    **kwargs,
                )
                axs.append(ax)
    return axs


def print_data_range(lcf_coll):
    """Print the data range for the given LightCurveCollection

    For each LightCurveFile:
    * sector start/stop time
    * first / last observation time
    * camera used
    """
    html = '<pre style="line-height: 1.1;">\n'
    html += "<summary>Sectors: " + str(list(map(lambda lc: lc.meta.get("SECTOR"), lcf_coll))) + f" ({len(lcf_coll)})" + "\n"
    html += "Observation period range / data range:" + "\n"
    html += "<details>"
    for lc in lcf_coll:
        html += f"  Sector {lc.meta.get('SECTOR')}: {lc.meta.get('TSTART')} - {lc.meta.get('TSTOP')}" + "\n"
        html += f"   (cam {lc.meta.get('CAMERA')})   {lc.time.min()} - {lc.time.max()}" + "\n"
    html += "</details></summary></pre>"
    display(HTML(html))


# Do the actual plots
def plot_all(
    lcf_coll,
    flux_col="flux",
    moving_avg_window=None,
    lc_tweak_fn=None,
    ax_fn=None,
    use_relative_time=False,
    mark_quality_issues=True,
    mark_momentum_dumps=True,
    set_title=True,
    ax_tweak_fn=None,
):
    """Plot the given LightCurveFile collection, one graph for each LightCurve

    Returns
    -------
    axs : the list of plots in `matplotlib.Axes`
    """
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

    #     ax_all.set_title((f"TIC {lcf_coll[0].PDCSAP_FLUX.label}, "
    #                       f"sectors {list(map(lambda lcf: lcf.meta.get('SECTOR'), lcf_coll))}"))
    #     return ax_all

    # choice 4: plot the lightcurve sector by sector: each sector in its own graph
    axs = []
    for i in range(0, len(lcf_coll)):
        if ax_fn is None:
            ax = lcf_fig().gca()
        else:
            ax = ax_fn()
        lcf = lcf_coll[i]
        lc = _to_lc_with_flux(lcf, flux_col)

        lc = _normalize_to_percent_quiet(lc)
        if lc_tweak_fn is not None:
            lc = lc_tweak_fn(lc)

        # temporarily change time to a relative one if specified
        if use_relative_time:
            rel_time_added = add_relative_time(lc, lcf)
            if rel_time_added:
                lc["time_orig"] = lc.time
                lc.time = lc.time_rel
            else:
                # the file has no observation start time, so we cannot add it
                use_relative_time = False

        # tweak label to include sector if any
        sector = lcf_coll[i].meta.get("SECTOR", None)
        label_long = lc.label
        if sector is not None:
            lc.label += f", s.{sector}"
            label_long += f", sector {sector}"
        if lc.author is not None and lc.author != "SPOC":
            label_long += f", by {lc.author}"

        if "long" == lke.estimate_cadence_type(lc):
            # long cadence has more spare data, use a larger "x" to represent them
            # "x" is also useful to distinguish it from moving average,
            # which will likely overlap with the points given the sparse data
            ax = lc.scatter(ax=ax, s=16, marker="x")
        else:
            ax = lc.scatter(ax=ax)

        # convert to dataframe to add moving average
        if moving_avg_window is not None:
            df = add_flux_moving_average(lc, moving_avg_window)
            # mask_gap: if there is a gap larger than 2 hours,
            # show the gap rather than trying to fill the gap with a straight line.
            ax.plot(
                lc.time.value,
                mask_gap(lc.time, df["flux_mavg"], 2 / 24),
                c="black",
                label=f"Moving average ({moving_avg_window})",
            )

        title_extras = ""
        if lc_tweak_fn is not None:
            title_extras = "\nLC tweaked, e.g., outliers removed"

        if set_title:
            ax.set_title(f"{label_long} {title_extras}", {"fontsize": 36})
        if use_relative_time:
            ax.xaxis.set_label_text("Time - relative")
            # restore original time after plot is done
            lc.time = lc.time_orig
        else:
            t_start = lc.meta.get("TSTART")
            if t_start is not None:
                ax.xaxis.set_label_text(ax.xaxis.label.get_text() + f", TSTART={t_start:0.2f}")

        # to avoid occasional formating in scientific notations
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

        # minor tick, 1 day interval in practice
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis="x", which="minor", length=4)
        # ax.xaxis.grid(True, which='minor') # too noisy to be there by default

        ax.xaxis.label.set_size(fontsize=18)
        ax.yaxis.label.set_size(fontsize=18)
        if ax_tweak_fn is not None:
            ax_tweak_fn(ax)

        # mark quality issue is applied after ax_tweak_fn, in case users use ax_tweak_fn and change the graph's ylim
        if mark_quality_issues:
            # the time where flux might have potential issues, using the suggested starting quality flag mask
            time = lc.time if not use_relative_time else lc.time_rel
            time_w_quality_issues = time[lke.create_quality_issues_mask(lc)]
            if len(time_w_quality_issues) > 0:
                # add marks as vertical lines at bottom 10% of the plot
                # Note: ax.vlines's ymin/ymax refers to the data. To specify them relative to y-axis

                # I have to 1) use transform, and
                #           2) tell the plot not to auto-scale Y-axis
                #               (if auto-scaled is done, it will treat the line's coodrinate as data)
                # somehow it doesn't work all the time. it could crop the y axis such that
                # only the vlines are visible
                #                 ax.set_autoscaley_on(False)
                #                 ax.vlines(time_w_quality_issues, ymin=0, ymax=0.1, transform=ax.get_xaxis_transform()
                #                           , color='red', linewidth=1, linestyle='--', label="potential quality issue")

                # back to visually less appealing one (that vline doesn't start from the bottom
                ybottom, ytop = ax.get_ylim()
                ax.vlines(
                    time_w_quality_issues.value,
                    ymin=ybottom,
                    ymax=ybottom + 0.1 * (ytop - ybottom),
                    color="red",
                    linewidth=1,
                    linestyle="--",
                    label="potential quality issue",
                )

        if mark_momentum_dumps:
            # Note: momentum_dump signals are by default masked out in LightCurve objects.
            # To access times marked as such, I need to access the raw LightCurveFile directly.
            with fits.open(lcf.filename) as hdu:
                if "TIME" in hdu[1].columns.names:
                    time = hdu[1].data["TIME"]
                    if use_relative_time:
                        t_start = lcf.meta.get("TSTART")
                        time = time - t_start
                    mom_dumps_mask = np.bitwise_and(hdu[1].data["QUALITY"], TessQualityFlags.Desat) >= 1
                    time_mom_dumps = time[mom_dumps_mask]
                    if len(time_mom_dumps) > 0:
                        ybottom, ytop = ax.get_ylim()
                        ax.vlines(
                            time_mom_dumps,
                            ymin=ybottom,
                            ymax=ybottom + 0.15 * (ytop - ybottom),
                            color="red",
                            linewidth=1,
                            linestyle="-.",
                            label="Momentum dumps",
                        )
                else:
                    # case the file has no TIME column, typically non SPOC-produced ones, e.g., CDIPS,
                    # the logic of finding momentum dump would not apply to such files anyway.
                    pass

        ax.legend()
        axs.append(ax)
    return axs


_lcf_4_plot_interactive = None


def _update_plot_lcf_interactive(figsize, flux_col, xrange, moving_avg_window, ymin, ymax, widget_out2):
    # use global to accept lct
    global _lcf_4_plot_interactive
    lcf = _lcf_4_plot_interactive

    ax = plt.figure(figsize=figsize).gca()
    plot_n_annotate_lcf(
        lcf,
        ax,
        flux_col=flux_col,
        xmin=xrange[0],
        xmax=xrange[1],
        moving_avg_window=moving_avg_window,
    )
    codes_text = f"ax.set_xlim({xrange[0]}, {xrange[1]})"
    ymin_to_use = ymin if ymin >= 0 else None
    ymax_to_use = ymax if ymax >= 0 else None
    if (ymin_to_use is not None) or (ymax_to_use is not None):
        ax.set_ylim(ymin_to_use, ymax_to_use)
        codes_text += f"\n\nax.set_ylim({ymin_to_use}, {ymax_to_use})"

    widget_out2.clear_output()
    with widget_out2:
        print(codes_text)

    return None


def plot_lcf_interactive(lcf, figsize=(15, 8), flux_col="flux"):
    desc_style = {"description_width": "25ch"}
    slider_style = {"description_width": "25ch"}
    slider_layout = {"width": "100ch"}
    t_start = lcf.meta.get("TSTART")
    t_stop = lcf.meta.get("TSTOP")
    # Add a second output for textual
    widget_out2 = widgets.Output()

    # pass lcf with a global rather than the slow fixed(lcf) with lkv2
    #
    # import warnings
    # with warnings.catch_warnings():
    #     # lkv2 workaround: to suppress astropy table warning, stating that the semantics of == will be changed in the future.
    #     warnings.filterwarnings("ignore", category=FutureWarning)
    #     fixed_lcf = fixed(lcf)
    global _lcf_4_plot_interactive
    _lcf_4_plot_interactive = lcf

    w = interactive(
        _update_plot_lcf_interactive,
        figsize=fixed(figsize),
        # lcf = fixed_lcf,
        flux_col=fixed(flux_col),
        xrange=widgets.FloatRangeSlider(
            min=t_start,
            max=t_stop,
            step=0.1,
            value=(t_start, t_stop),
            description="Time",
            continuous_update=False,
            readout_format=".1f",
            layout=slider_layout,
            style=slider_style,
        ),
        moving_avg_window=widgets.Dropdown(
            options=[
                ("None", None),
                ("10 min", "20min"),
                ("20 min", "20min"),
                ("30 min", "30min"),
                ("1 hour", "1h"),
                ("2 hours", "2h"),
                ("4 hours", "4h"),
            ],
            value="30min",
            description="Moving average window",
            style=desc_style,
        ),
        ymin=widgets.FloatText(value=-1, description="Flux min, -1 for default", style=desc_style),
        ymax=widgets.FloatText(value=-1, description="Flux max, -1 for default", style=desc_style),
        widget_out2=fixed(widget_out2),
    )
    w.layout.border = "1px solid lightgray"
    w.layout.padding = "1em 0px"

    widget_out2.layout.padding = "1em"
    w.children = w.children + (widget_out2,)

    display(w)
    return w


_lcf_4_plot_transit_interactive = None


def _update_plot_transit_interactive(
    figsize,
    flux_col,
    t0,
    duration_hr,
    period,
    step,
    surround_time,
    moving_avg_window,
    t0mark_ymax,
    ymin,
    ymax,
    widget_out2,
):
    # a clumsy way to pass lcf, without using fixed(lcf), which is very slow in lkv2
    global _cache_plot_n_annotate_lcf
    lcf = _lcf_4_plot_transit_interactive

    ax = plt.figure(figsize=figsize).gca()
    codes_text = "# Snippets to generate the plot"
    moving_avg_window_for_codes = "None" if moving_avg_window is None else f"'{moving_avg_window}'"
    if t0 < 0:
        plot_n_annotate_lcf(lcf, ax, flux_col=flux_col, moving_avg_window=moving_avg_window)
        codes_text += f"\nplot_n_annotate_lcf(lcf, ax, moving_avg_window={moving_avg_window_for_codes})"
    else:
        t0_to_use = t0 + step * period
        plot_transit(
            lcf,
            ax,
            t0_to_use,
            duration_hr / 24,
            surround_time,
            flux_col=flux_col,
            moving_avg_window=moving_avg_window,
            t0mark_ymax=t0mark_ymax,
        )
        codes_text += f"""
#   transit parameters - t0: BTJD {t0}, duration: {duration_hr} hours, period: {period} days
plot_transit(lcf, ax, {t0_to_use}, {duration_hr} / 24, {surround_time}, \
moving_avg_window={moving_avg_window_for_codes}, t0mark_ymax={t0mark_ymax})

# transit_specs for calling plot_transits()
transit_specs = TransitTimeSpecList(
    dict(epoch={t0}, duration_hr={duration_hr}, period={period}, label="dip",
         sector={lcf.meta.get('SECTOR')}, steps_to_show=[{step}],
        ),
    defaults=dict(surround_time={surround_time})
    )
"""

    ymin_to_use = ymin if ymin >= 0 else None
    ymax_to_use = ymax if ymax >= 0 else None
    if (ymin_to_use is not None) or (ymax_to_use is not None):
        ax.set_ylim(ymin_to_use, ymax_to_use)
        codes_text += f"""
# Zoom in on flux
ax.set_ylim({ymin_to_use}, {ymax_to_use})
"""
    widget_out2.clear_output()
    with widget_out2:
        print(codes_text)
    return None


def plot_transit_interactive(lcf, figsize=(15, 8), flux_col="flux"):
    desc_style = {"description_width": "25ch"}

    # Add a second output for textual
    widget_out2 = widgets.Output()

    t0 = widgets.FloatText(
        value=-1,
        step=0.01,
        description=r"$t_{epoch}$, -1 for unspecified",
        style=desc_style,
    )
    duration_hr = widgets.FloatText(value=1, step=0.01, description="duration (hours)", style=desc_style)
    period = widgets.FloatText(value=999, step=0.01, description="period (days)", style=desc_style)
    step = widgets.IntText(value=0, description=r"cycle (0 for transit at $t_{epoch}$)", style=desc_style)
    surround_time = widgets.FloatText(value=7, step=0.5, description="padding (days)", style=desc_style)
    moving_avg_window = widgets.Dropdown(
        options=[
            ("None", None),
            ("10 min", "10min"),
            ("20 min", "20min"),
            ("30 min", "30min"),
            ("1 hour", "1h"),
            ("2 hours", "2h"),
            ("4 hours", "4h"),
        ],
        value="20min",
        description="moving average window",
        style=desc_style,
    )
    ymin = widgets.FloatText(value=-1, step=0.1, description="flux min, -1 for default", style=desc_style)
    ymax = widgets.FloatText(value=-1, step=0.1, description="flux max, -1 for default", style=desc_style)
    t0mark_ymax = widgets.BoundedFloatText(
        value=0.05,
        step=0.05,
        min=0.0,
        max=1.0,
        description=r"$t_{epoch}$ mark height",
        style=desc_style,
    )
    VB = widgets.VBox
    HB = widgets.HBox
    ui = VB(
        [
            HB([t0, duration_hr, period]),
            HB([step, surround_time, moving_avg_window]),
            HB([ymin, ymax, t0mark_ymax]),
        ]
    )

    # pass lcf via a global, as fixed(lcf) is very slow with lkv2
    #
    # import warnings
    # with warnings.catch_warnings():
    #     # lkv2 workaround: to suppress astropy table warning, stating that the semantics of == will be changed in the future.
    #     warnings.filterwarnings("ignore", category=FutureWarning)
    #     fixed_lcf = fixed(lcf)
    global _lcf_4_plot_transit_interactive
    _lcf_4_plot_transit_interactive = lcf

    w = interactive_output(
        _update_plot_transit_interactive,
        dict(
            figsize=fixed(figsize),
            # lcf=fixed_lcf,
            flux_col=fixed(flux_col),
            t0=t0,
            duration_hr=duration_hr,
            period=period,
            step=step,
            surround_time=surround_time,
            moving_avg_window=moving_avg_window,
            t0mark_ymax=t0mark_ymax,
            ymin=ymin,
            ymax=ymax,
            widget_out2=fixed(widget_out2),
        ),
    )

    w.layout.border = "1px solid lightgray"
    w.layout.padding = "1em 0px"

    widget_out2.layout.padding = "1em"

    display(ui, w, widget_out2)
    return w


def plot_flux_sap_flux_comparison(lc, sap_col="sap_flux", ax=None, offset=None, **kwargs):
    """Plot flux (typically PDCSAP_FLUX) and sap_flux together,
    to spot any anomaly in processed lightcurve."""
    lc_sap = lc.copy()
    lc_sap["flux"] = lc[sap_col]
    if sap_col + "_err" in lc.colnames:
        lc_sap["flux_err"] = lc[sap_col + "_err"]
    else:  # some products, e.g., QLP, does not give err
        # Hit a bug - ValueError: TessLightCurve object is invalid - expected 'time' as the first columns but found 'time'
        # lc_sap.remove_column('flux_err')
        # zero out the column as a workaround
        lc_sap["flux_err"] = np.zeros_like(lc_sap["flux_err"])

    if offset is None:
        # auto offset: move lc_sap curve so that
        # - its median is at about 10 percentile of main flux
        # - move farther down by a factor of the 40% amplitude of the main flux, so that it is most below the main flux
        #    without much gap, but not overlapping too much either.
        offset = (
            np.nanpercentile(lc.flux.value, 10)
            - np.nanmedian(lc_sap.flux.value)
            - (np.nanmedian(lc.flux.value) - np.nanpercentile(lc.flux.value, 10)) * 3
        )

    lc_sap.label += f" {sap_col} + {offset:.0f}"

    ax = LightCurveCollection([lc, lc_sap]).plot(ax=ax, offset=offset, **kwargs)
    ax.set_title(f"{lc.label}, sector {lc.sector} - flux vs {sap_col}")
    return ax


class TransitTimeSpec(dict):
    def __init__(
        self,
        epoch: float = None,
        period: float = None,
        duration_hr: float = None,
        sector: int = None,
        steps_to_show: list = None,
        surround_time: float = None,
        label: str = None,
        defaults: TransitTimeSpec = None,
    ):
        # core parameters
        self["epoch"] = epoch
        self["period"] = period
        self["duration_hr"] = duration_hr

        # used for plotting
        self["sector"] = sector
        self["steps_to_show"] = steps_to_show
        self["surround_time"] = surround_time
        self["label"] = label

        if defaults is None:
            defaults = {}
        self._defaults = defaults  # put it as a custom attribute

    def __getitem__(self, key):
        res = super().get(key)
        if res is None:
            res = self._defaults.get(key)
        return res

    def get(self, key, default=None):
        res = self.__getitem__(key)
        if res is None:
            res = default
        return res


class TransitTimeSpecList(list):
    def __init__(self, *tt_spec_dict_list, defaults={}):
        self._defaults = TransitTimeSpec(**defaults)
        for tt_spec_dict in tt_spec_dict_list:
            self.append(TransitTimeSpec(**tt_spec_dict, defaults=self._defaults))

    def _spec_property_values(self, property_name):
        return np.array([tt[property_name] for tt in self])

    #
    # The following properties return the specific transit parameters
    # in an array. Together they can be used to create a mask
    # for the transits using ``LightCurve.create_transit_mask()``
    #

    @property
    def epoch(self):
        return self._spec_property_values("epoch")

    @property
    def period(self):
        return self._spec_property_values("period")

    @property
    def duration_hr(self):
        return self._spec_property_values("duration_hr")

    @property
    def duration(self):
        return self.duration_hr / 24

    @property
    def label(self):
        return self._spec_property_values("label")

    def to_table(self, columns=("label", "epoch", "duration_hr", "period")):
        """Convert the specs to an ``astropy.Table``"""
        data = [getattr(self, col) for col in columns]
        return Table(data, names=columns)


def mark_transit_times(
    lc, tt_specs, axvline_kwargs_specs=None, skip_no_transit_plot=False, lc_plot_func_name="scatter", ax=None
):
    """Plot the given LC, and mark the transit times based on `tt_specs`."""
    tt_list = [lke.get_transit_times_in_lc(lc, a_spec["epoch"], a_spec["period"]) for a_spec in tt_specs]

    # skip if no transit found
    # (tt_list is a list of list, so it needs to be flattend for counting)
    if skip_no_transit_plot and len(np.array(tt_list, dtype=object).flatten()) < 1:
        print(f"{lc._repr_simple_()} is skipped - no matching transits.")
        return None, None

    # base plot
    #
    if ax is None:
        #         ax = plt.figure(figsize=(30, 10)).gca()
        ax = plt.figure(figsize=(15, 5)).gca()
    ax = getattr(lc, lc_plot_func_name)(ax=ax, color="black", label=f"{lc.label} s.{getattr(lc, 'sector', 'N/A')}")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="x", which="minor", length=4)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="y", which="minor", length=4)

    # pre-process axvline_kwargs
    #
    if axvline_kwargs_specs is None:
        axvline_kwargs_specs = [dict(label="dip", linestyle="--", color="red")]

    # use the label in tt_specs if not specified in axvline_kwargs
    for (a_spec, an_axvline_kwargs, idx_0_based) in zip(tt_specs, axvline_kwargs_specs, range(len(axvline_kwargs_specs))):
        if an_axvline_kwargs.get("label") is None:
            an_axvline_kwargs["label"] = a_spec.get("label", f"dip {idx_0_based + 1}")

    # Mark transit times on the base plot
    #
    # a hack: mark the first line for each tt set, then set legend
    # so that each tt set will have 1 legend
    # if we simply set legend at the end, each dip will have its own legend!
    for (transit_times, axvline_kwargs) in zip(tt_list, axvline_kwargs_specs):
        if len(transit_times) > 0 and axvline_kwargs is not None:
            ax.axvline(transit_times[0], 0, 0.1, **axvline_kwargs)
    ax.legend()

    for (transit_times, axvline_kwargs) in zip(tt_list, axvline_kwargs_specs):
        if axvline_kwargs is not None:
            for tt in transit_times:
                ax.axvline(tt, 0, 0.1, **axvline_kwargs)

    return ax, tt_list


def scatter_centroids(
    lcf,
    fig=None,
    highlight_time_range=None,
    time_range=None,
    c="blue",
    c_highlight="red",
):
    """
    Scatter centroids, and highlight the specific time range
    """

    if fig is None:
        fig = plt.figure(figsize=(12, 12))

    lc = _normalize_to_percent_quiet(lcf)
    sector = lcf.meta.get("SECTOR")

    if time_range is not None:
        lc = lc.truncate(time_range[0], time_range[1])

    fig.gca().yaxis.set_major_formatter(FormatStrFormatter("%.3f"))  # avoid scientific notations
    fig.gca().scatter(lc.centroid_col.value, lc.centroid_row.value, c=c, label=f"TIC {lc.targetid}")

    if highlight_time_range is not None:
        lc_highlight = lc.truncate(highlight_time_range[0], highlight_time_range[1])
        if len(lc_highlight) < 1:
            print("WARNING: scatter_centroids() no observations in highlight_time_range")
        fig.gca().scatter(
            lc_highlight.centroid_col.value,
            lc_highlight.centroid_row.value,
            c=c_highlight,
            label="highlights",
        )

    title = f"TIC {lc.targetid} Centroids, sector {sector}"
    if time_range is not None:
        title += f"\n{as_4decimal(time_range)}"
    if highlight_time_range is not None:
        title += f"\nHighlights:{as_4decimal(highlight_time_range)}"
    fig.gca().set_title(title)
    fig.legend()
    return fig


def _update_anim(n, ax, lc, label, num_centroids_to_show, use_relative_time, c):
    ax.cla()
    # fix the x/y scale to ensure it doesn't change over the animation
    c_col, c_row = _to_unitless(lc.centroid_col), _to_unitless(lc.centroid_row)
    ax.set_xlim(np.nanmin(c_col), np.nanmax(c_col))
    ax.set_ylim(np.nanmin(c_row), np.nanmax(c_row))

    # avoid scientific notation for y-axis
    # x-axis might need scientific notation so that the labels won't get too cramped with long decimals
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))

    if num_centroids_to_show is None:
        col = lc.centroid_col[:n]
        row = lc.centroid_row[:n]
        time_label = f"{as_4decimal(lc.time[n])}"
        if use_relative_time:
            time_label = time_label + f" ({as_4decimal(lc.time_rel[n])})"
    else:
        n_start = max(0, n - num_centroids_to_show)
        col = lc.centroid_col[n_start:n]
        row = lc.centroid_row[n_start:n]
        time_label = f"{as_4decimal(lc.time[n_start])} - {as_4decimal(lc.time[n])}"
        if use_relative_time:
            time_label = time_label + f" ({as_4decimal(lc.time_rel[n_start])} - {as_4decimal(lc.time_rel[n])})"

    ax.set_title(f"TIC {lc.targetid} Centroids, {label}\nday: {time_label}")
    ax.scatter(col, row, c=c)


def animate_centroids(
    lcf,
    fig=None,
    frames=None,
    num_obs_per_frame=240,
    interval=250,
    use_relative_time=False,
    time_range=None,
    accumulative=True,
    c=None,
    display=True,
):
    """
    Animate centroids to visualize changes over time.

    """
    lc = lcf
    label = f"sector {lcf.meta.get('SECTOR')}"

    # Zoom to a particular time range if specified
    if time_range is not None:
        # use pandas to zoom to a particular time_range
        df = _normalize_to_percent_quiet(lc).to_pandas(columns=["time", "flux", "centroid_row", "centroid_col"])
        df = df[(df.time >= time_range[0]) & (df.time <= time_range[1])]
        if len(df) < 1:
            raise Exception(f"Zoomed lightcurve has no observation. time_range={time_range}")

        lc_z = SimpleNamespace()  # zoomed-in lightcurve-like object for the purpose of animation
        setattr(lc_z, "time", df.time.values)
        setattr(lc_z, "flux", df.flux.values)
        setattr(lc_z, "centroid_row", df.centroid_row.values)
        setattr(lc_z, "centroid_col", df.centroid_col.values)
        setattr(lc_z, "targetid", lc.targetid)
        lc = lc_z

    if fig is None:
        fig = plt.figure(figsize=(12, 12))
    if frames is None:
        num_obs = len(lc.centroid_row)
        num_frames = int(num_obs / num_obs_per_frame)  # default 240 is about every 8 hours, given 2-minute intervals
        ary_n = np.linspace(1, num_obs, num=num_frames, endpoint=False)
        ary_n[0] = np.ceil(ary_n[1] / 2)
        ary_n = list(map(lambda n: int(n), ary_n))
    else:
        ary_n = frames
        num_obs_per_frame = frames[-1] - frames[-2]  # assume the spacing of input is linear

    num_centroids_to_show = num_obs_per_frame
    if accumulative:
        num_centroids_to_show = None
    #     print(f'Steps: {ary_n}')
    if use_relative_time:
        rel_time_added = add_relative_time(lc, lcf)
        if not rel_time_added:
            use_relative_time = False
    anim = animation.FuncAnimation(
        fig,
        _update_anim,
        frames=ary_n,
        fargs=(fig.gca(), lc, label, num_centroids_to_show, use_relative_time, c),
        interval=interval,
        blit=False,
    )
    if display:
        # for inline display in jupyter
        try:
            from IPython.display import HTML
            from IPython.display import display as iDisplay

            return iDisplay(HTML(anim.to_jshtml(default_mode="once")))
        except ImportError:
            print("WARNING: animate_centroids() - inline display not possible Not in IPython environment.")
            return anim


def markTimes(ax, times, **kwargs):
    """Helper to mark specifics time as vertical lines on a plot"""
    axvline_kwargs = kwargs.copy()
    # apply defaults
    axvline_kwargs.setdefault("c", "gray")
    axvline_kwargs.setdefault("linewidth", 1)
    axvline_kwargs.setdefault("linestyle", "--")
    for t in times:
        ax.axvline(t, **axvline_kwargs)


def fold_and_plot_odd_even(lc, period, epoch_time, figsize=(12, 6), title_extra=""):
    lc_folded = lc.fold(period=period, epoch_time=epoch_time, epoch_phase=0)

    ax = plt.figure(figsize=figsize).gca()
    lc_f_odd = lc_folded[lc_folded.odd_mask]
    lc_f_odd.scatter(ax=ax, c="r", label="odd", marker=".", s=4)
    lc_f_even = lc_folded[lc_folded.even_mask]
    lc_f_even.scatter(ax=ax, c="b", label="even", marker="x", s=4)

    pct01_odd = np.nanpercentile(lc_f_odd.flux, 0.1)
    pct01_even = np.nanpercentile(lc_f_even.flux, 0.1)

    ax.axhline(
        pct01_odd * 100,
        c="r",
        linestyle="--",
        label=f"odd 0.1 pctile {pct01_odd:0.4f}",
    )
    ax.axhline(
        pct01_even * 100,
        c="b",
        linestyle="dotted",
        label=f"even 0.1 pctile {pct01_even:0.4f}",
    )

    ax.legend()
    plt.title(f"{lc.label} folded {title_extra}")

    print("odd  0.1 percentile: ", pct01_odd)
    print("even 0.1 percentile: ", pct01_even)
    return ax, lc_folded


def fold_2x_periods_and_plot(lc, period, epoch_time, figsize=(12, 6), title_extra=""):
    lc_folded = lc.fold(period=period * 2, epoch_time=epoch_time, epoch_phase=period / 2)

    ax = plt.figure(figsize=figsize).gca()
    lc_folded.scatter(ax=ax)

    ax.legend()
    plt.title(f"{lc.label} folded at 2X periods {title_extra}")

    return ax, lc_folded


#
# TargetPixelFile helpers
#


def show_tpf_orientation(tpf):
    """ "Helper to visualize the TPF's orientation in the sky. Requires IPython.
    Long arm is north, short arm with arrow is east.
    """
    coord_bottom_left = tpf.wcs.pixel_to_world(0, 0)
    coord_upper_right = tpf.wcs.pixel_to_world(tpf.shape[2] - 1, tpf.shape[1] - 1)
    coord_upper_left = tpf.wcs.pixel_to_world(0, tpf.shape[2] - 1)
    deg_from_north = coord_bottom_left.position_angle(coord_upper_left).to(u.deg).value

    display(
        HTML(
            f"""<div style="position: relative; margin-left: 16px;height: 64px;">
    <div title="Long arm: North; Short arm with arrow: East"
         style="float: left; max-width: 64px;font-size: 32px;margin: 16px;\
transform: rotate({-deg_from_north}deg);transform-origin: left; cursor:pointer;">↳</div>
        <div style="font-family: monospace;">Upper right offset from bottom left - <br>
        RA: {(coord_upper_right.ra - coord_bottom_left.ra).to(u.arcmin):0.6},
        Dec: {(coord_upper_right.dec - coord_bottom_left.dec).to(u.arcmin):0.6}
        </div>
    </div>"""
        )
    )


def interact_sky(tpf, notebook_url="localhost:8888", aperture_mask="empty", magnitude_limit=18):
    """tpf.interact_sky wrapper to handle different lightkurve versions."""
    if "aperture_mask" in inspect.getfullargspec(tpf.interact_sky).args:
        # case using a pre-release lightkurve that supports aperture_mask
        return tpf.interact_sky(notebook_url=notebook_url, aperture_mask=aperture_mask, magnitude_limit=magnitude_limit)
    else:
        # using release lightkurve that not yet supports aperture_mask
        return tpf.interact_sky(notebook_url=notebook_url, magnitude_limit=magnitude_limit)


def show_nearby_tic_summary_form():
    """Display a form that create a 1-line summary of a nearby TIC from the first three rows of the selected TIC table"""
    display(
        HTML(
            r"""
First 3 rows of the TIC info table:<br>
<textarea id="inStarInfo" style="width: 40ch; height: 6em;" placeholder="TIC \t12345678\nTESS Mag \t10.123\nSeparation ...">
</textarea><br>
<button id="ctlStarInfo">Create Nearby TIC summary</button>
<input id="outStarInfo" style="width: 40ch;" value="" readonly>

<script>

function convertToMultiLinePlaceholder(elem) { // create multiline placeholder
  elem.placeholder = elem.placeholder.replace(/\\n/g, '\n');
  elem.placeholder = elem.placeholder.replace(/\\t/g, '\t');
}
convertToMultiLinePlaceholder(document.querySelector('#inStarInfo'));


function createNearbyTicSummary(text) {
const toCells = (line) => {
    return line.split("\t");
}

const lines = text.split(/[\r\n]+/);


const ticId = toCells(lines[0])[1].replace(/^\s*(\d+).*$/, '$1');
const tessMag = toCells(lines[1])[1];
const separation = toCells(lines[2])[1];

return `TIC ${ticId} (TESS magnitude ${tessMag}, ${separation} arcsec away)`;
}

document.querySelector('#ctlStarInfo').onclick = (evt) => {
const summary = createNearbyTicSummary(document.querySelector('#inStarInfo').value);
document.querySelector('#outStarInfo').value = summary;
};
</script>
"""
        )
    )
