#
# Helpers to download TESS-specific non-lightcurve data: TOIs, TCEs, etc.
#

from collections.abc import Sequence
import pathlib
import re
from types import SimpleNamespace
import warnings

from memoization import cached
from retry import retry

import numpy as np
import pandas as pd
from pandas.io.formats.style import Styler

import astropy
from astropy import coordinates as coord
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from astroquery.utils import TableList

import asyncio_compat

import download_utils
import lightkurve as lk
import lightkurve_ext as lke

# Ues to resolve data files relative the the module (used by MomentumDumpsAccessor)
_MODULE_PATH_ = pathlib.Path(__file__).parent.resolve()

#
# Misc constants
#

R_earth = 6371000  # radius of the Earth [m]
R_jup = 69911000  # radius of Jupiter [m]
BTJD_REF = 2457000

#
# Generic CSV Download helper
#


def _get_csv(url, filename, download_dir, cache_policy_func, **kwargs):
    local_filename = download_utils.download_file(
        url, filename=filename, download_dir=download_dir, cache_policy_func=cache_policy_func
    )
    return pd.read_csv(local_filename, **kwargs)


def _single_row(df):
    if len(df) > 0:
        return df.iloc[0]
    else:
        return None


#
# TOIs / CTOIs
#


class TOIAccessor:
    Headers = SimpleNamespace(
        TIC="TIC ID",
        TOI="TOI",
        MASTER_PRIORITY="Master",
        EPOCH_BJD="Epoch (BJD)",
        EPOCH_BTJD="Epoch (BTJD)",  # derived
        PERIOD="Period (days)",
        DURATION_HR="Duration (hours)",
        DEPTH_PPM="Depth (ppm)",
        DEPTH_PCT="Depth (percent)",  # derived
        PLANET_RADIUS_E="Planet Radius (R_Earth)",
        PLANET_RADIUS_J="Planet Radius (R_Jupiter)",  # derived
        TESS_DISPOSITION="TESS Disposition",
        TFOPWG_DISPOSITION="TFOPWG Disposition",
        COMMENTS="Comments",
    )

    # TODO: in-memory cache (with @cached) needs to be redone to properly support cache_policy_func
    @classmethod
    def get_all_tois(cls, download_dir=None, cache_policy_func=None):
        url = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
        filename = "tess_tois.csv"
        res = _get_csv(url, filename, download_dir, cache_policy_func=cache_policy_func, dtype={cls.Headers.TOI: str})
        # add derived columns
        res[cls.Headers.EPOCH_BTJD] = res[cls.Headers.EPOCH_BJD] - BTJD_REF
        res[cls.Headers.PLANET_RADIUS_J] = res[cls.Headers.PLANET_RADIUS_E] * R_earth / R_jup
        res[cls.Headers.DEPTH_PCT] = res[cls.Headers.DEPTH_PPM] / 10000
        return res

    def __init__(self, download_dir=None, cache_policy_func=None):
        self._all = self.get_all_tois(download_dir=download_dir, cache_policy_func=cache_policy_func)

    def all(self):
        return self._all

    def of_toi(self, toi):
        tois = self._all[self._all[self.Headers.TOI] == str(toi)]
        return _single_row(tois)

    def of_tic(self, tic):
        return self._all[self._all[self.Headers.TIC] == int(tic)]

    def of_tics(self, tics):
        return self._all[np.isin(self._all[self.Headers.TIC], tics)]


class CTOIAccessor:
    Headers = SimpleNamespace(
        TIC="TIC ID",
        CTOI="CTOI",
        TOI="Promoted to TOI",
        EPOCH_BJD="Transit Epoch (BJD)",
        EPOCH_BTJD="Transit Epoch (BTJD)",  # derived
        PERIOD="Period (days)",
        DURATION_HR="Duration (hrs)",
        DEPTH_PPM="Depth ppm",
        DEPTH_PCT="Depth percent",  # derived
        PLANET_RADIUS_E="Planet Radius (R_Earth)",
        PLANET_RADIUS_J="Planet Radius (R_Jupiter)",  # derived
        COMMENTS="Notes",
    )

    @classmethod
    def get_all_ctois(cls, download_dir=None, cache_policy_func=None):
        url = "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv"
        filename = "tess_ctois.csv"
        res = _get_csv(
            url,
            filename,
            download_dir,
            cache_policy_func=cache_policy_func,
            dtype={cls.Headers.CTOI: str, cls.Headers.TOI: str},
        )
        # add derived columns
        res[cls.Headers.EPOCH_BTJD] = res[cls.Headers.EPOCH_BJD] - BTJD_REF
        res[cls.Headers.PLANET_RADIUS_J] = res[cls.Headers.PLANET_RADIUS_E] * R_earth / R_jup
        res[cls.Headers.DEPTH_PCT] = res[cls.Headers.DEPTH_PPM] / 10000
        return res

    def __init__(self, download_dir=None, cache_policy_func=None):
        self._all = self.get_all_ctois(download_dir=download_dir, cache_policy_func=cache_policy_func)

    def all(self):
        return self._all

    def of_ctoi(self, ctoi):
        ctois = self._all[self._all[self.Headers.CTOI] == str(ctoi)]
        return _single_row(ctois)

    def of_tic(self, tic):
        return self._all[self._all[self.Headers.TIC] == int(tic)]


def add_transit_as_codes_column_to_df(df, headers, label_value_func):
    h = headers
    # string interpolation does not work. So use old-school concatenation

    # for single transit TOI/CTOIs, period returned is often nan or 0
    # to make the codes (used in transit_specs) usable later on
    # we substitute it with a large period
    def handle_nan_or_zero(per):
        if np.isnan(per) or per == 0.0:
            return 9999.9
        else:
            return per

    # somehow `period = pd.Series([handle_nan_or_zero(p) for p in df[h.PERIOD]])`
    # does not work. I temporarily created a new column as a workaround
    df["_period_nan_fixed"] = [handle_nan_or_zero(p) for p in df[h.PERIOD]]

    df["Codes"] = (
        "epoch="
        + df[h.EPOCH_BTJD].map("{:.4f}".format)
        + ", duration_hr="
        + df[h.DURATION_HR].map("{:.4f}".format)
        + ", period="
        + df["_period_nan_fixed"].map("{:.6f}".format)
        + ', label="'
        + label_value_func(df)
        + '", transit_depth_percent='
        + df[h.DEPTH_PCT].map("{:.4f}".format)
        + ","
    )
    df.drop(columns=["_period_nan_fixed"], inplace=True)  # drop the temp column
    return df


def _get_tois_in_html(tic, download_dir=None):
    h = TOIAccessor.Headers
    # Consider cache TOIAccessor in some module global (keyed by download_dir) to avoid
    # repeated loading/parsing the underlying TOI csv
    tois = TOIAccessor(download_dir=download_dir).of_tic(tic)
    if len(tois) < 1:
        return "<p>No TOIs.</p>"

    add_transit_as_codes_column_to_df(tois, h, label_value_func=lambda df: "TOI " + df[h.TOI])
    report_view = tois[
        [
            h.TOI,
            h.MASTER_PRIORITY,
            h.TFOPWG_DISPOSITION,
            h.PLANET_RADIUS_J,
            h.EPOCH_BTJD,
            h.DURATION_HR,
            h.PERIOD,
            h.DEPTH_PCT,
            h.COMMENTS,
            "Codes",
        ]
    ]
    # tweak output styling
    styler = Styler(report_view, cell_ids=False)  # avoid unnecessary long cell ids
    styler.hide(axis="index")
    styler.format(
        formatter={
            (h.PLANET_RADIUS_J): "{:.3f}",
            (h.EPOCH_BTJD, h.DURATION_HR): "{:.4f}",
            (h.PERIOD): "{:.6f}",
            (h.DEPTH_PCT): "{:.4f}",
        }
    )
    styler.set_table_styles(
        [
            # make the TOI table align (roughly) with the TCE table
            {"selector": "td.col0", "props": [("padding-left", "10px")]},
        ]
    )

    html = styler._repr_html_()

    # make the headers to make them more compact
    html = html.replace(h.MASTER_PRIORITY, "Master<br>priority", 1)
    html = html.replace(h.TFOPWG_DISPOSITION, "TFOPWG<br>Dispo.", 1)
    html = html.replace(h.PLANET_RADIUS_J, "R<sub>p</sub><br>R<sub>j</sub>", 1)
    html = html.replace(h.EPOCH_BTJD, "Epoch<br>BTJD", 1)
    html = html.replace(h.DURATION_HR, "Duration<br>hr", 1)
    html = html.replace(h.PERIOD, "Period<br>day", 1)
    html = html.replace(h.DEPTH_PCT, "Depth<br>%", 1)

    # render nan as -- (as nan is really no value in our case)
    # - styler.format()'s na_rep option seems to fix some but not all, so we do it ourselves
    # - replace the pattern of <td class="..." >nan</td>
    html = html.replace(">nan</td>", ">--</td>")

    # turn Codes column into html input element (easier to be selected)
    html = re.sub(
        r"<td([^>]+)>(epoch=.+,)</td>",
        r"""<td\1><input  type="text" style="margin-left: 3ch; font-size: 90%; color: #666; width: 10ch;" onclick="this.select();" readonly="" value='\2'></td>""",
        html,
    )

    return html


def _get_ctois_in_html(tic, download_dir=None):
    # TODO: lots of codes similar to _get_tois_in_html(). factor them out

    h = CTOIAccessor.Headers
    # Consider cache TOIAccessor in some module global (keyed by download_dir) to avoid
    # repeated loading/parsing the underlying TOI csv
    ctois = CTOIAccessor(download_dir=download_dir).of_tic(tic)
    if len(ctois) < 1:
        return "<p>No CTOIs.</p>"

    add_transit_as_codes_column_to_df(ctois, h, label_value_func=lambda df: "CTOI " + df[h.CTOI])
    report_view = ctois[
        [
            h.CTOI,
            h.TOI,
            h.PLANET_RADIUS_J,
            h.EPOCH_BTJD,
            h.DURATION_HR,
            h.PERIOD,
            h.DEPTH_PCT,
            h.COMMENTS,
            "Codes",
        ]
    ]
    # tweak output styling
    styler = Styler(report_view, cell_ids=False)  # avoid unnecessary long cell ids
    styler.hide(axis="index")
    styler.format(
        formatter={
            (h.PLANET_RADIUS_J): "{:.3f}",
            (h.EPOCH_BTJD, h.DURATION_HR): "{:.4f}",
            (h.PERIOD): "{:.6f}",
            (h.DEPTH_PCT): "{:.4f}",
        }
    )
    styler.set_table_styles(
        [
            # make the CTOI table align (roughly) with the TCE table
            {"selector": "td.col0", "props": [("padding-left", "20px")]},
            # min-width to ensure TOI column, often no value, are wide enough to hold typical TOI value
            # so as to make alignment more consistent
            {"selector": "td.col1", "props": [("min-width", "80px")]},
        ]
    )

    html = styler._repr_html_()

    # make the headers to make them more compact
    html = html.replace(h.TOI, "TOI?", 1)
    html = html.replace(h.PLANET_RADIUS_J, "R<sub>p</sub><br>R<sub>j</sub>", 1)
    html = html.replace(h.EPOCH_BTJD, "Epoch<br>BTJD", 1)
    html = html.replace(h.DURATION_HR, "Duration<br>hr", 1)
    html = html.replace(h.PERIOD, "Period<br>day", 1)
    html = html.replace(h.DEPTH_PCT, "Depth<br>%", 1)

    # render nan as -- (as nan is really no value in our case)
    # - styler.format()'s na_rep option seems to fix some but not all, so we do it ourselves
    # - replace the pattern of <td class="..." >nan</td>
    html = html.replace(">nan</td>", ">--</td>")

    # turn Codes column into html input element (easier to be selected)
    html = re.sub(
        r"<td([^>]+)>(epoch=.+,)</td>",
        r"""<td\1><input  type="text" style="margin-left: 3ch; font-size: 90%; color: #666; width: 10ch;" onclick="this.select();" readonly="" value='\2'></td>""",
        html,
    )

    return html


def get_tic_meta_in_html(
    lc_or_tic, a_subject_id=None, download_dir=None, tce_filter_func=None, include_transit_model_stellar_density=False
):
    # tess_dv_fast.py is at https://github.com/orionlee/tess_dv_fast
    # copy it over (or include it in sys.path)
    import tess_dv_fast

    # This function does not do the actual display,
    # so that the caller can call it in background
    # and display it wherever it's needed
    def link(link_text, url):
        return f"""<a href="{url}" target="_blank">{link_text}</a>"""

    def prop(prop_name, prop_value):
        return f"""    <tr><td>{prop_name}</td><td>{prop_value}</td></tr>\n"""

    # main logic
    if isinstance(lc_or_tic, lk.LightCurve):
        tic_id = str(lc_or_tic.meta.get("TICID"))
    elif isinstance(lc_or_tic, (str, int)):
        tic_id = lc_or_tic
    else:
        raise TypeError("lc_or_tic must be either a LightCurve object or a tic id (int/str)")

    m = _to_stellar_meta(lc_or_tic)

    def safe_m_get(key, default_val):
        # in some meta, the key exists but the value is None
        # this helper handles it
        res = getattr(m, key, default_val)
        return res if res is not None else default_val

    html = f"""
<div id="tic_metadata_ctr">
<div id="tic_metadata_body">
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
        html += f' (sector {safe_m_get("sector", "")})'
    html += "<br>\n"

    # stellar parameters
    html += "<table>\n"
    s_radius = safe_m_get("radius", -1)
    s_mass = safe_m_get("mass", -1)
    html += prop("R<sub>S</sub> (in R<sub>☉</sub>)", f"{s_radius:.3f}")
    html += prop("M<sub>S</sub> (in M<sub>☉</sub>)", f"{s_mass:.3f}")
    if s_radius > 0 and s_mass > 0:
        s_rho = lke.estimate_rho(s_mass, s_radius, return_unit=u.g / u.cm**3).value
    else:
        s_rho = -1
    if include_transit_model_stellar_density:
        html += prop("rho<sub>S</sub> (in g/cm<sup>3</sup>)", f"{s_rho:.3f}")
    html += prop("Magnitude (TESS)", f'{safe_m_get("tess_mag", -1):.2f}')
    html += prop("T_eff (in K)", safe_m_get("teff", -1))
    html += "</table>\n"

    html += "<p>TCEs:</p>"

    # TODO: not yet implemented in tess_dv_fast
    #   tce_filter_func=tce_filter_func,
    # Note: tess_dv_fast cannot support
    #   include_transit_model_stellar_density=include_transit_model_stellar_density,
    df_tces = tess_dv_fast.get_tce_infos_of_tic(tic_id)
    html += tess_dv_fast.display_tce_infos(df_tces, return_as="html", no_tce_html="<p>No TCE.</p>")

    # TOIs/CTOIs
    html += "<p>TOIs / CTOIs:</p>"
    html += _get_tois_in_html(tic_id, download_dir=download_dir)
    html += _get_ctois_in_html(tic_id, download_dir=download_dir)

    html += """
</div> <!-- id="tic_metadata_body" -->
</div> <!--  id="tic_metadata_ctr" -->
"""
    return html


#
# TESS Momentum dump accessor
#


def get_momentum_dump_times(lcf):
    """Get the momentum dump times from the given lightcurve.
    It is usually one from a sector.
    The output can be added to data/tess_mom_dumps.txt for further plotting usage.
    """
    # Note: momentum_dump signals are by default masked out in LightCurve objects.
    # To access times marked as such, I need to access the raw LightCurveFile directly.
    filename = lcf.meta.get("FILENAME", None)
    if filename is None:
        warnings.warn("get_momentum_dump_times(): No-Op, because there is the LightCurve object has no backing FITS file.")
        return np.array([])
    with fits.open(filename) as hdu:
        if "TIME" not in hdu[1].columns.names:
            # case the file has no TIME column, typically non SPOC-produced ones, e.g., CDIPS,
            # the logic of finding momentum dump would not apply to such files anyway.
            return np.array([])

        # normal flow
        time = hdu[1].data["TIME"]
        mom_dumps_mask = np.bitwise_and(hdu[1].data["QUALITY"], lk.utils.TessQualityFlags.Desat) >= 1
        time_mom_dumps = time[mom_dumps_mask]
        return time_mom_dumps


class MomentumDumpsAccessor:
    _mom_dumps_tab = None

    @classmethod
    def _load_mom_dumps_from_file(cls, cache_policy=download_utils.CachePolicy.ALWAYS_USE):
        local_filename = download_utils.download_file(
            "https://tess.mit.edu/public/files/Table_of_momentum_dumps.csv",
            filename="tess_mom_dumps.csv",
            download_dir=f"{_MODULE_PATH_}/data/",
            cache_policy_func=cache_policy,
        )

        # I use np.genfromtxt rather than pandas, as filtering from numpy array filtering is easier for the use case
        raw_table = np.genfromtxt(local_filename, delimiter=",", comments="#")
        # Time in BTJD in 2nd column, I don't care about the rest.
        cls._mom_dumps_tab = raw_table[:, 1]

    @classmethod
    def refresh(cls):
        cls._load_mom_dumps_from_file(cache_policy=download_utils.CachePolicy.ALWAYS_REJECT)

    @classmethod
    def get_all(cls, refresh=False):
        if refresh or cls._mom_dumps_tab is None:
            cls.refresh()
        return cls._mom_dumps_tab

    @classmethod
    def get_in_range(cls, lc_or_tpf=None, start=None, end=None, refresh=False):
        times = cls.get_all(refresh=refresh)

        if lc_or_tpf is not None:
            start, end = (
                lc_or_tpf.time.value.min(),
                lc_or_tpf.time.value.max() + 1e-6,
            )

        if start is not None:
            times = times[times >= start]
        if end is not None:
            times = times[times < end]
        return times

    @classmethod
    def exclude_around(cls, lc_or_tpf=None, window_before=15 / 60 / 24, window_after=15 / 60 / 24):
        """Exclude cadences of the given LC / TPF around momentum dumps.
        Useful to exclude data points that are often skewed.
        """

        def compress_as_exclude_ranges(mom_dumps, window_before, window_after):
            """Transform the list of momentum dumps to a list of range of time to exclude.
            The function also compress the momentum dump list, by consolidating
            multiple nearby timestamps to a single range.
            This is done as an performance optimization, to reduce the number of actual lc/tpf truncation needed.
            (In practice, it cuts the time or processing a typical TESS 2-minute cadence tpf from a few seconds to ~500ms)
            """
            if len(mom_dumps) < 1:
                return []
            res = []
            cur_range = [mom_dumps[0] - window_before, mom_dumps[0] + window_after]
            for t in mom_dumps[1:]:
                if t <= cur_range[1]:
                    cur_range[1] = t + window_after
                else:
                    res.append(cur_range)
                    cur_range = [t - window_before, t + window_after]
            res.append(cur_range)
            return res

        mom_dumps = cls.get_in_range(lc_or_tpf)
        exclude_ranges = compress_as_exclude_ranges(mom_dumps, window_before, window_after)

        res = lc_or_tpf
        for an_exclude in exclude_ranges:
            t = res.time.value
            res = res[(t < an_exclude[0]) | (t >= an_exclude[1])]

        return res


class WTVResultAccessor:
    @classmethod
    def get_all(cls, wtv_csv_path, add_sectors_summary=True, start_sector=1, end_sector_inclusive=69):
        res = Table.read(wtv_csv_path)
        if not add_sectors_summary:
            return res

        def to_sectors_summary(row):
            summary = ""
            for sector in np.arange(start_sector, end_sector_inclusive + 1):
                if row[f"S{sector}"] > 0:
                    summary = f"{summary} {sector},"
            return summary

        summary_ary = []
        for row in res:
            summary_ary.append(to_sectors_summary(row))
        res["Sectors"] = summary_ary

        return res


#
# TGLC lightcurve search based on downloaded target list CSV files locally
#   https://archive.stsci.edu/hlsp/tglc
#


def search_tess_point(tic):
    from tess_stars2px import tess_stars2px_function_entry

    rs = catalog_info_of_tics_in_vizier(tic)
    ra, dec = rs[0]["RAJ2000"], rs[0]["DEJ2000"]

    # note: the search below does not use tic,
    # it's just there for convenience to keep track of targets
    outID, outEclipLong, outEclipLat, outSec, outCam, outCcd, outColPix, outRowPix, scinfo = tess_stars2px_function_entry(
        tic, ra, dec
    )

    return Table(
        data=dict(
            tic=outID,
            sector=outSec,
            camera=outCam,
            ccd=outCcd,
            column=outColPix,
            row=outRowPix,
        )
    )


def _search_tglc_lightcurve_csv(tic, csv_dir=".", grep_cmd="grep -H"):
    def parse_line(line):
        # eg: s0005.csv:5573000755560118784,393294857,91.82816617709965,-40.32806441931107
        f = re.split("[:,]", line)
        sector = int(re.search(r"\d+", f[0])[0])
        # return (sector, TIC, gaia_source)
        return (sector, int(f[2]), int(f[1]))

    def parse_grep_out(lines):
        return [parse_line(l) for l in lines.splitlines()]

    import subprocess

    if "rg" not in grep_cmd:
        cmdline = rf"{grep_cmd}  ,{tic}, s*.csv"
    else:
        # special case for ripgrep, glob does not work on Windows
        # https://github.com/BurntSushi/ripgrep/issues/234
        cmdline = rf"{grep_cmd}  ,{tic}, ."
    res = subprocess.run(cmdline, cwd=csv_dir, capture_output=True, text=True)
    if res.returncode == 0:
        return parse_grep_out(res.stdout)
    elif res.returncode == 1:  # rg signals no match found
        return []

    # case unexpected error
    raise Exception(f"Error in search TGLC lcs: {res.stderr} . Returncode: {res.returncode}", res)


async def search_tglc_lightcurve(tic, csv_dir=".", grep_cmd="grep -H"):
    tess_point_task = asyncio_compat.create_background_task(search_tess_point, tic)
    csv_out = _search_tglc_lightcurve_csv(tic, csv_dir=csv_dir, grep_cmd=grep_cmd)
    tess_point_out = await tess_point_task

    out = []
    for r in csv_out:
        sector, tic, gaia_source = r
        tp = tess_point_out[tess_point_out["sector"] == sector][0]
        camera, ccd = tp["camera"], tp["ccd"]

        # eg:
        #   'https://archive.stsci.edu/hlsps/tglc/s0045/cam1-ccd1/0033/2814/1948/3267/
        #    hlsp_tglc_tess_ffi_gaiaid-3328141948326716800-s0045-cam1-ccd1_tess_v1_llc.fits'
        #
        #   'https://archive.stsci.edu/hlsps/tglc/s0039/cam3-ccd4/0046/2303/6865/3737/
        #    hlsp_tglc_tess_ffi_gaiaid-4623036865373793408-s0039-cam3-ccd4_tess_v1_llc.fits'

        gp = f"{gaia_source:021}"

        download_url = (
            f"https://archive.stsci.edu/hlsps/tglc/s{sector:04}/cam{camera}-ccd{ccd}/"
            f"{gp[0:4]}/{gp[4:8]}/{gp[8:12]}/{gp[12:16]}/"
            f"hlsp_tglc_tess_ffi_gaiaid-{gaia_source}-s{sector:04}-cam{camera}-ccd{ccd}_tess_v1_llc.fits"
        )
        out.append(
            dict(
                tic=tic,
                sector=sector,
                url=download_url,
            )
        )

    # convert to astropy Table
    if len(out) > 0:
        tab = Table(rows=out)
    else:
        tab = Table(data=dict(tic=[], sector=[], url=[]))
    tab.sort(keys=["tic", "sector"])
    return tab


#
# TESS Flux - Magnitude Conversion
#


def tess_flux_to_mag(flux):
    """Convert flux from TESS observation to magnitude."""
    # Based on https://heasarc.gsfc.nasa.gov/docs/tess/observing-technical.html#saturation
    # TESS CCDs produce 15000 e/s for magnitude 10 light source

    if isinstance(flux, u.Quantity):
        flux_raw = flux.to(u.electron / u.second).value
    else:
        flux_raw = flux

    # np.log10 does not work on Quantity, unless it's dimensionless
    res_raw = 10 - 2.5 * np.log10((flux_raw / 15000))

    if isinstance(flux, u.Quantity):
        return res_raw * u.mag
    else:
        return res_raw


def mag_to_tess_flux(mag):
    """Convert magnitude to flux in TESS observation."""
    if isinstance(mag, u.Quantity):
        mag_raw = mag.to(u.mag).value
    else:
        mag_raw = mag

    flux_raw = (10 ** ((10 - mag_raw) / 2.5)) * 15000
    if isinstance(mag, u.Quantity):
        return flux_raw * ((u.electron / u.second))
    else:
        return flux_raw


def calc_flux_range(lcf_coll, flux_column="flux", accepted_authors=["SPOC", "TESS-SPOC"], stitched_lc_corrector=lambda lc: lc):
    """Derive flux range (in % and magnitude) from normalized lightcurve with mean TESS mag from TIC as the base"""

    # the default SPOC, TESS-SPOC is to avoid inconsistency between SPOC and QLP
    lcf_coll_filtered = lke.select(lcf_coll, lambda lc: lc.author in accepted_authors)
    lc = lke.stitch(
        lcf_coll_filtered,
        corrector_func=lambda lc: (
            lc.select_flux(flux_column).remove_nans()
            # normalize on per-sector basis, it seems TESS calibration across sectors is not necessarily consistent
            .normalize(unit="percent")
        ),
    )

    # optionally let caller tweak teh stitched LC, e.g., excluding some cadences, say, if flares are to be ignored.
    lc = stitched_lc_corrector(lc)

    flux_range_pct = np.asarray([lc.flux.max(), lc.flux.min()])
    base_mag = lc.meta.get("TESSMAG")
    flux_range_mag = lke.normalized_flux_val_to_mag(flux_range_pct, base_mag=base_mag)

    return SimpleNamespace(
        flux_range_pct=flux_range_pct,
        flux_range_mag=flux_range_mag,
        lc_stitched=lc,
        lcf_coll_filtered=lcf_coll_filtered,
    )


#
# Misc. TESS specific utilities
#


def display_crowdsap(lc):
    """Display `CROWDSAP` of a TESS SPOC lightcurve. Highlight it if it could be too low.

    Warning that the field might be crowded using CROWDSAP header,
    when CROWDSAP < 0.8 .

    From section 4.1.2 of the paper TOI catalog from TESS primary mission
    https://arxiv.org/pdf/2103.12538.pdf
    """
    from IPython.display import display, HTML

    if lc is not None and lc.meta.get("CROWDSAP") is not None:
        display(
            HTML(
                f"""Fraction of flux in aperture attributed to the target, <span style="font-family: monospace;">CROWDSAP</span>:
        <span style="background-color: {'red' if lc.meta.get("CROWDSAP") < 0.8 else 'transparent'}; padding: 2px;">{lc.meta.get("CROWDSAP")}</span>
        &emsp;<span style="font-family: monospace;">FLFRCSAP</span>: {lc.meta.get('FLFRCSAP')}
        &emsp;<a href="https://heasarc.gsfc.nasa.gov/docs/tess/UnderstandingCrowding.html" target="_crowdsap_tutorial">(Help)</a>
        """
            )
        )


def btjd_to_hjd_utc(time_val, position):
    t_btjd = Time(time_val, format="btjd", scale="tdb")
    t_bjd = t_btjd.copy("jd")

    if isinstance(position, coord.SkyCoord):
        sky_coord = position
    elif isinstance(position, str):
        ra, dec = position.split(",")
        sky_coord = coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")
    elif isinstance(position, dict):
        ra, dec = position["ra"], position["dec"]
        sky_coord = coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")
    else:
        raise TypeError(f"position, of type {type(position)} is not in supported types.")

    return lke.to_hjd_utc(t_bjd, sky_coord).value


#
# TIC Metadata in catalogs (TIC / Gaia)
#


MAST_QUERY_NUM_RETRIES = 4


# OPEN: consider to cache result in disk with diskcache
@cached
@retry(IOError, tries=MAST_QUERY_NUM_RETRIES, delay=0.5, backoff=2, jitter=(0, 0.5))
def catalog_info_of_tics(tic):
    """Return the info of a TIC in the TIC catalog"""
    from astroquery.mast import Catalogs

    return Catalogs.query_criteria(catalog="Tic", ID=tic)


def catalog_info_of_tics_in_vizier(tic, columns=["*", "Tmag"]):
    """Return the info of a TIC in the TIC catalog (from Vizier)"""
    # (Vizier tends to be speedier. The result is cached by astroquery on disk)
    from astroquery.vizier import Vizier

    Vizier.ROW_LIMIT = -1
    return Vizier(catalog="IV/39", columns=columns).query_constraints(TIC=tic)[0]


def _to_stellar_meta(target):
    if hasattr(target, "meta"):  # case LC, TPF, etc.
        meta = target.meta
        ra, dec, equinox = meta.get("RA"), meta.get("DEC"), meta.get("EQUINOX")
        pmra, pmdec = meta.get("PMRA"), meta.get("PMDEC")
        tess_mag = meta.get("TESSMAG")
        tic = meta.get("TICID")
        if tic is not None:
            label = f"{tic}"
        else:
            label = f"[{ra:4f} {dec:4f}]"
        return SimpleNamespace(
            # for Gaia DR3 query use case
            ra=ra,
            dec=dec,
            equinox=equinox,
            pmra=pmra,
            pmdec=pmdec,
            e_pmra=np.nan,
            e_pmdec=np.nan,
            tess_mag=tess_mag,
            label=label,
            # additional attributes for get_tic_meta_in_html() use case
            sector=meta.get("SECTOR"),
            radius=meta.get("RADIUS"),
            teff=meta.get("TEFF"),
        )

    # case target is a tic id
    if isinstance(target, (int, str)):
        result = catalog_info_of_tics(target)
        if len(result) < 1:
            return None
        row = result[0]
        ra, dec, equinox = row["ra"], row["dec"], 2000
        pmra, pmdec = row["pmRA"], row["pmDEC"]
        e_pmra, e_pmdec = row["e_pmRA"], row["e_pmDEC"]
        tess_mag = row["Tmag"]
        tic = row["ID"]
        if tic is not None:
            label = f"{tic}"
        else:
            label = f"[{ra:4f} {dec:4f}]"
        gaiadr2_id = row["GAIA"]  # useful to crossmatch with Gaia data
        return SimpleNamespace(
            ra=ra,
            dec=dec,
            equinox=equinox,
            pmra=pmra,
            pmdec=pmdec,
            e_pmra=e_pmra,
            e_pmdec=e_pmdec,
            tess_mag=tess_mag,
            label=label,
            gaiadr2_id=gaiadr2_id,
            # additional attributes for get_tic_meta_in_html() use case
            radius=row["rad"],
            mass=row["mass"],
            teff=row["Teff"],
        )
    raise TypeError(f"target, of type {type(target)} is not supported")


def decode_gaiadr3_nss_flag(nss_flag):
    """Decode NSS (NON_SINGLE_STAR) flag in Gaia DR3 Main.
    Reference:
    https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html#p344
    """
    flags = []
    for mask, nss_type in [
        (0b1, "AB"),  # astrometric binary
        (0b10, "SB"),  # spectroscopic binary
        (0b100, "EB"),  # eclipsing binary
    ]:
        if nss_flag & mask > 0:
            flags.append(nss_type)
    return flags


def decode_gaiadr3_nss_solution_flag(sol_flag):
    """Decode the flag for NSS solution in Gaia DR3 NSS 2 body orbit tables.
    Not to be confused with the flag in Gaia DR3 Main table.

    Reference:
    https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_non--single_stars_tables/ssec_dm_nss_two_body_orbit.html#p155
    """

    bits_meaning = {
        # for AB
        0: "AB_No_solution_searched",
        1: "AB_No_stochastic_solution_searched",
        2: "AB_Failure_to_compute_a_stochastic_solution",
        6: "AB_RV_available",
        7: "AB_RV_used_for_perspective_acceleration_correction",
        # for SB
        8: "SB_BAD_UNCHECKED_NUMBER_OF_TRANSITS",
        9: "SB_NO_MORE_VARIABLE_AFTER_FILTERING",
        10: "SB_BAD_CHECKED_NUMBER_OF_TRANSITS",
        11: "SB2_REDIRECTED_TO_SB1_CHAIN_NOT_ENOUGH_COUPLE_MEASURES",
        12: "SB2_REDIRECTED_TO_SB1_CHAIN_PERIODS_NOT_COHERENT",
        13: "SB_NO_SIGNIFICANT_PERIODS_CAN_BE_FOUND",
        14: "SB_REFINED_SOLUTION_DOES_NOT_CONVERGE",
        15: "SB_REFINED_SOLUTION_SINGULAR_VARIANCE_COVARIANCE_MATRIX",
        16: "SB_CIRCULAR_SOLUTION_SINGULAR_VARIANCE_COVARIANCE_MATRIX",
        17: "SB_TREND_SOLUTION_SINGULAR_VARIANCE_COVARIANCE_MATRIX",
        18: "SB_REFINED_SOLUTION_NEGATIVE_DIAGONAL_OF_VARIANCE_COVARIANCE_MATRIX",
        19: "SB_CIRCULAR_SOLUTION_NEGATIVE_DIAGONAL_OF_VARIANCE_COVARIANCE_MATRIX",
        20: "SB_TREND_SOLUTION_NEGATIVE_DIAGONAL_OF_VARIANCE_COVARIANCE_MATRIX",
        21: "SB_CIRCULAR_SOLUTION_DOES_NOT_CONVERGE",
        22: "SB_LUCY_TEST_APPLIED",
        23: "SB_TREND_SOLUTION_NOT_APPLIED",
        24: "SB_SOLUTION_OUTSIDE_E_LOGP_ENVELOP",
        25: "SB_PERIOD_FOUND_IN_CU7_PERIODICITY",
        26: "SB_FORTUITOUS_SB2",
        # for EB
        32: "EB_No_variance-covariance_matrix",
        # for Combined solutions
        48: "CO_NOCOMBINATION_FOUND",
        49: "CO_BAD_GOF_COMBINATION",
        50: "CO_WRONG_COMPONENT_COMBINATION",
        51: "CO_SB2_TREATED_AS_SB1",
        52: "CO_STOCHA_TO_ORBITAL",
        53: "CO_STOCHA_TO_MULTIPLE",
        54: "CO_ORBITALALTERNATIVE_TO_ORBITAL",
        55: "CO_TRIPLE_COMBINATION",
        56: "CO_TREND_COMBINATION",
        57: "CO_DU434_INPUT_USED",
    }

    flags = []
    for bit, meaning in bits_meaning.items():
        if sol_flag & 2**bit > 0:
            flags.append(meaning)
    return flags


def _is_all_finite(num_or_num_list):
    def is_one_finite(num):
        return num is not None and np.isfinite(num)

    if not isinstance(num_or_num_list, (list, tuple, np.ndarray)):
        num_or_num_list = [num_or_num_list]
    return np.all([is_one_finite(n) for n in num_or_num_list])


def linkify_gaiadr3_result_html(result: Table, max_num_rows=999):
    """Return Gaia DR3 Vizier Search Result in HtML, with links back to Vizier"""

    # the linkify depends on Source column.
    # so simply return html if Source is not there
    if "Source" not in result.colnames:
        with astropy.conf.set_temp("max_lines", max_num_rows):
            return result._repr_html_()

    # the result table content is to be tweaked for linkify implementation
    # so we use a copy.
    result = result.copy()

    # Include the source for variable and NSS when applicable
    # They will be used for reformatting as link.
    if "VarFlag" in result.colnames:
        c = np.char.add(result["VarFlag"], result["Source"].astype(str))
        c[result["VarFlag"] != "VARIABLE"] = "NOT_AVAILABLE"
        result["VarFlag"] = c
    if "NSS" in result.colnames:
        c = np.char.add(result["NSS"].astype(str), result["Source"].astype(str))
        c = np.char.add(np.full_like(c, "NSS"), c)
        c[result["NSS"] == 0] = "0"
        result["NSS"] = c

    with astropy.conf.set_temp("max_lines", max_num_rows):
        html = result._repr_html_()

    # linkify Gaia DR3 Source, Variable and NSS
    for id in result["Source"]:
        # the long URL includes both Gaia DR3 Main and Astrophysical, with frequently used columns included.
        gaiadr3_url = f"https://vizier.cds.unistra.fr/viz-bin/VizieR-4?-ref=VIZ6578bb1b54eda&-to=-4b&-from=-4&-this=-4&%2F%2Fsource=I%2F355%2Fgaiadr3&%2F%2Ftables=I%2F355%2Fgaiadr3&%2F%2Ftables=I%2F355%2Fparamp&-out.max=50&%2F%2FCDSportal=http%3A%2F%2Fcdsportal.u-strasbg.fr%2FStoreVizierData.html&-out.form=HTML+Table&%2F%2Foutaddvalue=default&-order=I&-oc.form=sexa&-out.src=I%2F355%2Fgaiadr3%2CI%2F355%2Fparamp&-nav=cat%3AI%2F355%26tab%3A%7BI%2F355%2Fgaiadr3%7D%26tab%3A%7BI%2F355%2Fparamp%7D%26key%3Asource%3DI%2F355%2Fgaiadr3%26HTTPPRM%3A&-c=&-c.eq=J2000&-c.r=++2&-c.u=arcmin&-c.geom=r&-source=&-x.rs=10&-source=I%2F355%2Fgaiadr3+I%2F355%2Fparamp&-out.orig=standard&-out=RA_ICRS&-out=DE_ICRS&-out=Source&Source={id}&-out=Plx&-out=PM&-out=pmRA&-out=pmDE&-out=sepsi&-out=IPDfmp&-out=RUWE&-out=Dup&-out=Gmag&-out=BPmag&-out=RPmag&-out=BP-RP&-out=RV&-out=e_RV&-out=VarFlag&-out=NSS&-out=XPcont&-out=XPsamp&-out=RVS&-out=EpochPh&-out=EpochRV&-out=MCMCGSP&-out=MCMCMSC&-out=Teff&-out=logg&-out=%5BFe%2FH%5D&-out=Dist&-out=A0&-out=HIP&-out=PS1&-out=SDSS13&-out=SKYM2&-out=TYC2&-out=URAT1&-out=AllWISE&-out=APASS9&-out=GSC23&-out=RAVE5&-out=2MASS&-out=RAVE6&-out=RAJ2000&-out=DEJ2000&-out=Pstar&-out=PWD&-out=Pbin&-out=ABP&-out=ARP&-out=GMAG&-out=Rad&-out=SpType-ELS&-out=Rad-Flame&-out=Lum-Flame&-out=Mass-Flame&-out=Age-Flame&-out=Flags-Flame&-out=Evol&-out=z-Flame&-meta.ucd=0&-meta=0&-usenav=1&-bmark=GET"
        html = html.replace(
            f">{id}<",
            f"><a target='vizier_gaia_dr3' href='{gaiadr3_url}'>{id}</a><",
        )

        gaiadr3_var_url = f"https://vizier.cds.unistra.fr/viz-bin/VizieR-4?-ref=VIZ65ac1f481b91d6&-to=-4b&-from=-3&-this=-4&%2F%2Fsource=%2BI%2F358%2Fvarisum%2BI%2F358%2Fvclassre%2BI%2F358%2Fveb%2BI%2F358%2Fvcc%2BI%2F358%2Fvst&%2F%2Fc=06%3A59%3A36.3+%2B23%3A28%3A51.14&%2F%2Ftables=I%2F358%2Fvarisum&%2F%2Ftables=I%2F358%2Fvclassre&%2F%2Ftables=I%2F358%2Fvcc&%2F%2Ftables=I%2F358%2Fveb&%2F%2Ftables=I%2F358%2Fvst&-out.max=50&%2F%2FCDSportal=http%3A%2F%2Fcdsportal.u-strasbg.fr%2FStoreVizierData.html&-out.form=HTML+Table&-out.add=_r&%2F%2Foutaddvalue=default&-sort=_r&-order=I&-oc.form=sexa&-out.src=I%2F358%2Fvarisum%2CI%2F358%2Fvclassre%2CI%2F358%2Fveb%2CI%2F358%2Fvcc%2CI%2F358%2Fvst&-nav=cat%3AI%2F358%26tab%3A%7BI%2F358%2Fvarisum%7D%26tab%3A%7BI%2F358%2Fvclassre%7D%26tab%3A%7BI%2F358%2Fvcc%7D%26tab%3A%7BI%2F358%2Fveb%7D%26tab%3A%7BI%2F358%2Fvst%7D%26key%3Asource%3D%2BI%2F358%2Fvarisum%2BI%2F358%2Fvclassre%2BI%2F358%2Fveb%2BI%2F358%2Fvcc%2BI%2F358%2Fvst%26key%3Ac%3D06%3A59%3A36.3+%2B23%3A28%3A51.14%26pos%3A06%3A59%3A36.3+%2B23%3A28%3A51.14%28+60+arcsec%29%26HTTPPRM%3A&-c=&-c.eq=J2000&-c.r=+60&-c.u=arcsec&-c.geom=r&-source=&-x.rs=10&-source=I%2F358%2Fvarisum+I%2F358%2Fvclassre+I%2F358%2Fveb+I%2F358%2Fvcc+I%2F358%2Fvst&-out.orig=standard&-out=Source&Source={id}&-out=RA_ICRS&-out=DE_ICRS&-out=TimeG&-out=DurG&-out=Gmagmean&-out=TimeBP&-out=DurBP&-out=BPmagmean&-out=TimeRP&-out=DurRP&-out=RPmagmean&-out=VCR&-out=VRRLyr&-out=VCep&-out=VPN&-out=VST&-out=VLPV&-out=VEB&-out=VRM&-out=VMSO&-out=VAGN&-out=Vmicro&-out=VCC&-out=SolID&-out=Classifier&-out=Class&-out=ClassSc&-out=Rank&-out=TimeRef&-out=Freq&-out=magModRef&-out=PhaseGauss1&-out=sigPhaseGauss1&-out=DepthGauss1&-out=PhaseGauss2&-out=sigPhaseGauss2&-out=DepthGauss2&-out=AmpCHP&-out=PhaseCHP&-out=ModelType&-out=Nparam&-out=rchi2&-out=PhaseE1&-out=DurE1&-out=DepthE1&-out=PhaseE2&-out=DurE2&-out=DepthE2&-out=Per&-out=T0G&-out=T0BP&-out=T0RP&-out=HG0&-out=HG1&-out=HG2&-out=HG3&-out=HG4&-out=HG5&-out=HBP0&-out=HBP1&-out=HBP2&-out=HBP3&-out=HBP4&-out=HBP5&-out=HRP0&-out=HRP1&-out=HRP2&-out=HRP3&-out=HRP4&-out=HRP5&-out=Gmodmean&-out=BPmodmean&-out=RPmodmean&-out=Mratiomin&-out=alpha&-out=Ampl&-out=NfoVTrans&-out=FoVAbbemean&-out=NTimeScale&-out=TimeScale&-out=Variogram&-meta.ucd=2&-meta=1&-meta.foot=1&-usenav=1&-bmark=GET"
        html = html.replace(
            f">VARIABLE{id}<",
            (
                f"><a target='vizier_gaia_dr3_var' href='{gaiadr3_var_url}' "
                "style='background-color: rgba(255, 255, 0, 0.5); font-weight: bold;'>VARIABLE</a><"
            ),
        )

        gaiadr3_nss_url = f"https://vizier.cds.unistra.fr/viz-bin/VizieR-4?-ref=VIZ65a1a2351812e4&-to=-4b&-from=-3&-this=-4&%2F%2Fsource=I%2F357&%2F%2Ftables=I%2F357%2Ftboasb1c&%2F%2Ftables=I%2F357%2Ftboeb&%2F%2Ftables=I%2F357%2Ftboes&%2F%2Ftables=I%2F357%2Ftbooc&%2F%2Ftables=I%2F357%2Ftbooac&%2F%2Ftables=I%2F357%2Ftbooavc&%2F%2Ftables=I%2F357%2Ftbootsc&%2F%2Ftables=I%2F357%2Ftbootsvc&%2F%2Ftables=I%2F357%2Ftbosb1&%2F%2Ftables=I%2F357%2Ftbosb1c&%2F%2Ftables=I%2F357%2Ftbosb2&%2F%2Ftables=I%2F357%2Ftbosb2c&%2F%2Ftables=I%2F357%2Facc7&%2F%2Ftables=I%2F357%2Facc9&%2F%2Ftables=I%2F357%2Flinspec1&%2F%2Ftables=I%2F357%2Flinspec2&%2F%2Ftables=I%2F357%2Fvimfl&-out.max=50&%2F%2FCDSportal=http%3A%2F%2Fcdsportal.u-strasbg.fr%2FStoreVizierData.html&-out.form=HTML+Table&-out.add=_r&%2F%2Foutaddvalue=default&-sort=_r&-order=I&-oc.form=sexa&-out.src=I%2F357%2Ftboasb1c%2CI%2F357%2Ftboeb%2CI%2F357%2Ftboes%2CI%2F357%2Ftbooc%2CI%2F357%2Ftbooac%2CI%2F357%2Ftbooavc%2CI%2F357%2Ftbootsc%2CI%2F357%2Ftbootsvc%2CI%2F357%2Ftbosb1%2CI%2F357%2Ftbosb1c%2CI%2F357%2Ftbosb2%2CI%2F357%2Ftbosb2c%2CI%2F357%2Facc7%2CI%2F357%2Facc9%2CI%2F357%2Flinspec1%2CI%2F357%2Flinspec2%2CI%2F357%2Fvimfl&-nav=cat%3AI%2F357%26tab%3A%7BI%2F357%2Ftboasb1c%7D%26tab%3A%7BI%2F357%2Ftboeb%7D%26tab%3A%7BI%2F357%2Ftboes%7D%26tab%3A%7BI%2F357%2Ftbooc%7D%26tab%3A%7BI%2F357%2Ftbooac%7D%26tab%3A%7BI%2F357%2Ftbooavc%7D%26tab%3A%7BI%2F357%2Ftbootsc%7D%26tab%3A%7BI%2F357%2Ftbootsvc%7D%26tab%3A%7BI%2F357%2Ftbosb1%7D%26tab%3A%7BI%2F357%2Ftbosb1c%7D%26tab%3A%7BI%2F357%2Ftbosb2%7D%26tab%3A%7BI%2F357%2Ftbosb2c%7D%26tab%3A%7BI%2F357%2Facc7%7D%26tab%3A%7BI%2F357%2Facc9%7D%26tab%3A%7BI%2F357%2Flinspec1%7D%26tab%3A%7BI%2F357%2Flinspec2%7D%26tab%3A%7BI%2F357%2Fvimfl%7D%26key%3Asource%3DI%2F357%26HTTPPRM%3A&-c=&-c.eq=J2000&-c.r=++2&-c.u=arcmin&-c.geom=r&-source=&-x.rs=10&-source=I%2F357%2Ftboasb1c+I%2F357%2Ftboeb+I%2F357%2Ftboes+I%2F357%2Ftbooc+I%2F357%2Ftbooac+I%2F357%2Ftbooavc+I%2F357%2Ftbootsc+I%2F357%2Ftbootsvc+I%2F357%2Ftbosb1+I%2F357%2Ftbosb1c+I%2F357%2Ftbosb2+I%2F357%2Ftbosb2c+I%2F357%2Facc7+I%2F357%2Facc9+I%2F357%2Flinspec1+I%2F357%2Flinspec2+I%2F357%2Fvimfl&-out.orig=standard&-out=Source&Source={id}&-out=NSSmodel&-out=RA_ICRS&-out=DE_ICRS&-out=Plx&-out=pmRA&-out=pmDE&-out=ATI&-out=BTI&-out=FTI&-out=GTI&-out=CTI&-out=HTI&-out=Per&-out=Tperi&-out=ecc&-out=Vcm&-out=Flags&-out=_RA.icrs&-out=_DE.icrs&-out=ffactp&-out=ffacts&-out=inc&-out=Tratio&-out=Teclp&-out=Tecls&-out=Durp&-out=Durs&-out=K1&-out=MassRatio&-out=K2&-out=dpmRA&-out=dpmDE&-out=ddpmRA&-out=ddpmDE&-out=Velmean&-out=dVel%2Fdt&-out=dVel%2Fdt2&-out=RAVIM&-out=DEVIM&-meta.ucd=2&-meta=1&-meta.foot=1&-usenav=1&-bmark=GET"
        html = re.sub(
            rf">NSS(\d){id}<",
            (
                rf"><a target='vizier_gaia_dr3_nss' href='{gaiadr3_nss_url}' "
                rf"style='background-color: rgba(255, 255, 0, 0.5); font-weight: bold;'>&ensp;\1&ensp;</a><"
            ),
            html,
        )
    return html


def search_gaiadr3_of_tics(
    targets,
    radius_arcsec=15,
    magnitude_range=2.5,
    pm_error_factor=None,  # e.g., 3
    pm_range_fraction=0.25,
    pm_range_minimum=1.0,
    warn_if_all_filtered=True,
    calc_separation_from_first_row=False,
    compact_columns=True,
    also_return_html=True,
    also_return_astrophysical=False,  # defaulted to False for backward compatibility
    verbose_html=True,
    include_nss_summary_in_html=True,
):
    """Locate the lightcurve target's correspond entry in GaiaDR3.
    The match is by an heuristics based on coordinate and magnitude.

    Parameters
    ----------
    target : int, LightCurve, TargetPixelFile, or a list of them
        targets to be searched. Either the TIC, or LightCurve/TargetPixelFile of a TIC.

    pm_error_factor, pm_range_fraction, pm_range_minimum
        range of proper motion to include in the search result,
        with the pmRA / pmDEC of the target as the reference.
        `pm_error_factor` is used if e_pmRA, e_pmDEC is present (case the target is from TIC catalog).
        The pmRA range will be `e_pmRA` * `pm_error_factor` (ditto for pmDEC)
        `pm_range_fraction` is used if e_pmRA, e_pmDEC is not present.
        The pmRA range will be `pmRA` * `pm_range_fraction` (ditto for pmDEC)
        In all cases if `pm_range_minimum` is defined, the range will have a minimum of
        `+/- pm_range_minimum`.
        It is useful to handle the case the derived range is very small and overly restrictive.
    """

    # OPEN:
    # Consider alternative by crossmatching Gaia DR2 of the TIC (available on MAST) with Gaia EDR3
    # https://gea.esac.esa.int/archive/documentation/GEDR3/Catalogue_consolidation/chap_cu9dr2xm/sec_cu9dr2xm_adql_queries/sec_cu9dr2xm_closest_edr3_neighbour_to_each_dr2_source_10m.html
    # some suggestion: limit cross match result by comparing GMag (difference < 0.1 was suggested in some doc)
    # Other Gaia crossmatch tips:
    # https://www.cosmos.esa.int/web/gaia-users/archive/combine-with-other-data

    if isinstance(targets, (Sequence, np.ndarray, lk.collections.Collection)):
        add_target_as_col = True
    else:
        add_target_as_col = False
        targets = [targets]

    # _paramp: Gaia DR3 Astrophysical, "I/355/paramp"
    result_list, result_paramp_list = [], []

    targets = np.asarray([_to_stellar_meta(t) for t in targets])
    targets = targets[targets != None]

    for t in targets:
        gaiadr2_id = getattr(t, "gaiadr2_id", np.nan)
        if magnitude_range is not None:
            lower_limit, upper_limit = t.tess_mag - magnitude_range, t.tess_mag + magnitude_range
        else:
            lower_limit, upper_limit = None, None

        pmra_lower, pmra_upper, pmdec_lower, pmdec_upper = None, None, None, None
        if _is_all_finite([pm_error_factor, t.pmra, t.e_pmra]):
            pmra_range = t.e_pmra * pm_error_factor
            if pm_range_minimum is not None:
                pmra_range = max(pmra_range, pm_range_minimum)
            pmra_lower, pmra_upper = t.pmra - pmra_range, t.pmra + pmra_range
            pmdec_range = t.e_pmdec * pm_error_factor
            if pm_range_minimum is not None:
                pmdec_range = max(pmdec_range, pm_range_minimum)
            pmdec_lower, pmdec_upper = t.pmdec - pmdec_range, t.pmdec + pmdec_range
        elif _is_all_finite([pm_range_fraction, t.pmra]):
            pmra_range = np.abs(t.pmra) * pm_range_fraction
            if pm_range_minimum is not None:
                pmra_range = max(pmra_range, pm_range_minimum)
            pmra_lower, pmra_upper = t.pmra - pmra_range, t.pmra + pmra_range
            pmdec_range = np.abs(t.pmdec) * pm_range_fraction
            if pm_range_minimum is not None:
                pmdec_range = max(pmdec_range, pm_range_minimum)
            pmdec_lower, pmdec_upper = t.pmdec - pmdec_range, t.pmdec + pmdec_range

        # print("DBG pm range for filter -  ra:", t.pmra, pmra_lower, pmra_upper, "dec: ", t.pmdec, pmdec_lower, pmdec_upper)
        a_result, a_result_paramp = lke.search_nearby(
            t.ra,
            t.dec,
            equinox=f"J{t.equinox}",
            radius_arcsec=radius_arcsec,
            magnitude_limit_column="RPmag",
            magnitude_lower_limit=lower_limit,
            magnitude_upper_limit=upper_limit,
            pmra=t.pmra,
            pmdec=t.pmdec,
            pmra_lower=pmra_lower,
            pmra_upper=pmra_upper,
            pmdec_lower=pmdec_lower,
            pmdec_upper=pmdec_upper,
            calc_separation_from_first_row=calc_separation_from_first_row,
            include_gaiadr3_astrophysical=True,
            warn_if_all_filtered=warn_if_all_filtered,
        )

        if a_result is not None:
            a_result["target"] = [t.label for i in range(0, len(a_result))]
            a_result["target_gaia_dr2_source"] = [gaiadr2_id for i in range(0, len(a_result))]

            result_list.append(a_result)
            result_paramp_list.append(a_result_paramp)

    with warnings.catch_warnings():
        # Avoid spurious "MergeConflictWarning: Cannot merge meta key 'null' types <class 'float'>
        #  and <class 'float'>, choosing null=nan [astropy.utils.metadata]"
        result = astropy.table.vstack(result_list) if len(result_list) > 0 else Table()
        result_paramp = astropy.table.vstack(result_paramp_list) if len(result_paramp_list) > 0 else Table()

    if len(result) < 1:
        if also_return_html:
            return None, None, ""
        else:
            return (
                None,
                None,
            )

    # flag entries (that could indicate binary systems, etc.)
    flag_column_values = []
    for row in result:
        flag = ""
        # RUWE > 1.4 cutoff source:
        # https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_ruwe.html
        if row["RUWE"] > 1.4:
            flag += "!"
        # astrometric excess noise sig > 2 cutoff source
        # https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html
        # https://web.archive.org/web/20211121142803/https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html
        if row["sepsi"] > 2:
            flag += "!"
        # e_RV > 1.5 heuristics suggested by mhuten that seems to be reliable
        if row["e_RV"] > 1.5:
            flag += "!"
        if str(row["Source"]) == str(row["target_gaia_dr2_source"]):  # use str to avoid str / int type complication
            flag += " ✓"  # signify Gaia DR3 ID match with TIC's Gaia DR2 ID
        flag_column_values.append(flag)
    result.add_column(flag_column_values, name="flag")
    result_all_columns = result
    result_paramp_all_columns = result_paramp

    if compact_columns:  # select the most useful columns
        # prefer RA/DEC in Epoch 2000 as they can be compared more easily with those those in TESS LC metadata
        result = result[
            "target",
            "flag",
            "_r",  # angular separation
            "_p",  # positional angle
            "Source",
            "RPmag",  # prioritize RPmag, as its pass band is close to TESS
            "Gmag",
            "BPmag",
            "BP-RP",
            "Vmag",  # calculated locally from Gmag and BP-RP
            "Teff",  # "Tefftemp" in Gaia DR2 / Gaia EDR3
            "RUWE",
            "sepsi",
            "epsi",
            "NSS",  # Gaia DR3: 1 if there is entry in Non-Single Star tables
            "Plx",
            "pmRA",
            "pmDE",
            "VarFlag",  # Gaia DR3: variability
            "RV",  # Gaia DR3
            "e_RV",  # Gaia DR3, e_RV > 1.5 km/s also possibly signifies non single star
            "IPDfmp",  # Gaia DR3, Percent of successful-IPD windows with more than one peak. High => possibly visually double
            "Dup",  # Gaia DR3: if there are multiple source/Gaia DR3 entries for the same target
            "RAJ2000",
            "DEJ2000",
            "EpochPh",  # Gaia DR3: 1 if epoch photometry is available
            "EpochRV",  # Gaia DR3: 1 if epoch RV is available
        ]
        if not add_target_as_col:
            result.remove_column("target")
        result_paramp = result_paramp[
            # separation "_r" is included, to make it easier to cross-reference
            # a row with the main result above
            "_r",
            "Source",
            "Pstar",
            "Pbin",
            "Teff",
            "logg",
            "__Fe_H_",
            "Dist",
            "GMAG",
            "Rad",
            "SpType-ELS",
            "ClassELS",
            "Rad-Flame",
            "Mass-Flame",
            "Lum-Flame",
            "Age-Flame",
            "Evol",
            "Flags-Flame",
        ]

    if also_return_html:
        html = ""
        if verbose_html:
            for t in targets:
                stellar_summary = (
                    f"TIC {t.label} - TESS mag: {t.tess_mag} ; coordinate: {t.ra}, {t.dec} ; " f"PM: {t.pmra}, {t.pmdec} ."
                )
                gaiadr2_id = getattr(t, "gaiadr2_id", None)
                if gaiadr2_id is not None:
                    stellar_summary += f" Gaia DR2 {gaiadr2_id}"
                html += f"<pre>{stellar_summary}</pre>"

        html += linkify_gaiadr3_result_html(result)
        # remove the Table length=n message
        result_len = len(result)
        html = html.replace(f"<i>Table length={result_len}</i>", "")

        html += '<div style="padding-top: 6px; font-size: 120%;">Gaia DR3 Astrophysical parameters:<div>'
        with astropy.conf.set_temp("max_lines", 999):
            html += result_paramp._repr_html_()
        # remove the Table length=n message
        result_paramp_len = len(result_paramp)
        html = html.replace(f"<i>Table length={result_paramp_len}</i>", "")

        if include_nss_summary_in_html:
            with warnings.catch_warnings():
                # ignore the warning
                # "Format strings passed to MaskedConstant are ignored, but in future may error or produce different behavior"
                warnings.filterwarnings(
                    "ignore", category=FutureWarning, message=".*Format strings passed to MaskedConstant are ignored,.*"
                )
                for i in range(0, len(result_all_columns)):
                    html_inner = (
                        f"RUWE: {result_all_columns[i]['RUWE']}, "
                        f"astrometric excess noise significance: {result_all_columns[i]['sepsi']:.3f}, "
                        f"e_RV: {result_all_columns[i]['e_RV']:.2f} km/s"
                        f"; ipd_frac_multi_peak: {result_all_columns[i]['IPDfmp']}%"
                    )
                    nss_flag = result_all_columns[i]["NSS"]
                    if nss_flag > 0:
                        html_inner += f"; NSS: {nss_flag} ({decode_gaiadr3_nss_flag(nss_flag)})"
                    html += f"<pre>{html_inner}</pre>"

        if verbose_html:
            html += """
<br>flag - &emsp; !: non single star proxy indicator (RUWE, sepsi, e_RV) ; &emsp; ✓: Gaia source matched with the one in TIC
<br>Reference:
    <a target ="_doc_gaiadr3_vizier" href="https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3">Column description on Vizier</a> &nbsp; | &nbsp;
    <a target="_doc_gaiadr3_datamodel" href="https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html">data model doc on ESA</a>
    &ensp;( <a target="_doc_gaiadr3_datamodel_astrophysical" href="https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_astrophysical_parameter_tables/ssec_dm_astrophysical_parameters.html">astrophysical params</a>,
    <a target="_doc_evolstage_basti" href="https://iopscience.iop.org/article/10.3847/1538-4357/aab158#apjaab158s4" style="font-size: 90%;">evolstage from BASTI</a>
    )
"""
    else:
        html = None

    return_list = [result_all_columns, result]
    if also_return_astrophysical:
        return_list += [result_paramp_all_columns, result_paramp]
    if html is not None:
        return_list.append(html)

    return return_list


from astroquery.vizier import Vizier


def _to_tic_str_list(targets):
    def to_tic_str(tic_id_or_obj):
        if isinstance(tic_id_or_obj, str):
            return tic_id_or_obj
        if isinstance(tic_id_or_obj, int):
            return str(tic_id_or_obj)
        if isinstance(tic_id_or_obj, (lk.LightCurve, lk.TessTargetPixelFile)):
            return tic_id_or_obj.meta.get("TICID")
        raise TypeError(f"Unsupported type: {type(tic_id_or_obj)}")

    if not isinstance(targets, (Sequence, np.ndarray, lk.collections.Collection)):
        targets = [targets]

    # Vizier requires a list of TIC in string
    targets = [to_tic_str(t) for t in targets]
    targets = [t for t in targets if t is not None]
    return targets


def search_tesseb_of_tics(
    targets,
    compact_columns=True,
    compact_preferred_model="pf",  # "pf" or "2g"
    also_return_html=True,
):
    """Locate the lightcurve target's correspond entry in TESS EB.
    (The static version in Vizier for sectors 1 - 26)

    Parameters
    ----------
    target : int, LightCurve, TargetPixelFile, or a list of them
        targets to be searched. Either the TIC, or LightCurve/TargetPixelFile of a TIC.

    """

    targets = _to_tic_str_list(targets)

    columns = ["*", "Sectors", "UpDate"]  # UpDate: the "date_modified" column
    result_list = Vizier(catalog="J/ApJS/258/16/tess-ebs", columns=columns).query_constraints(TIC=targets)

    if len(result_list) < 1:
        return None, None, "None found"
    result = result_list[0]  # there is only 1 table in the catalog

    result.rename_column("_tab1_10", "BJD0")  # somehow the column name for BJD0 is incorrect

    # add convenience columns

    result["Epochp"] = result["BJD0"]

    for m in ["pf", "2g"]:
        result[f"Durationp-{m}"] = result["Per"] * result[f"Wp-{m}"] * 24  # assuming "Per" is in unit day
        result[f"Durationp-{m}"].unit = u.hour

        result[f"Epochs-{m}"] = result["BJD0"] + result["Per"] * (result[f"Phis-{m}"] - result[f"Phip-{m}"])
        result[f"Epochs-{m}"].unit = result["BJD0"].unit

        result[f"Durations-{m}"] = result["Per"] * result[f"Ws-{m}"] * 24  # assuming "Per" is in unit day
        result[f"Durations-{m}"].unit = u.hour

    result_all_columns = result

    if compact_columns:
        _cpm = compact_preferred_model  # shortern it for convenience
        result = result[
            "TIC",
            "m_TIC",
            "Morph",
            "Per",
            "Epochp",
            f"Durationp-{_cpm}",
            f"Dp-{_cpm}",
            f"Epochs-{_cpm}",
            f"Durations-{_cpm}",
            f"Ds-{_cpm}",
            "Sectors",
        ]

    if also_return_html:
        rs = result.copy()  # create a copy so that I can add / change columns optimized for display
        # add link to live TESS EB
        rs["TESSebs"] = ["!TESSEB-" + tic for tic in rs["TIC"]]  # to be converted to link in html

        _cpm = compact_preferred_model  # shortern it for convenience
        rs["Code"] = [
            f"""epoch={epoch_p}, duration_hr={dur_p}, period={per}, label="primary",   epoch={epoch_s}, duration_hr={dur_s}, period={per}, label="secondary","""
            for per, epoch_p, dur_p, epoch_s, dur_s in zip(
                rs["Per"], rs["Epochp"], rs[f"Durationp-{_cpm}"], rs[f"Epochs-{_cpm}"], rs[f"Durations-{_cpm}"]
            )
        ]

        html = rs._repr_html_()

        # linkify TESSebs
        for id in rs["TIC"]:
            html = html.replace(
                f">!TESSEB-{id}<",
                f"><a target='tesseb' href='http://tessebs.villanova.edu/{id}'>{id}</a><",
            )

        # make Code column as HTML <input> for ease of copy
        html = re.sub(
            r"<td([^>]*)>(epoch=.+,)\s*</td>",
            r"""<td\1><input  type="text" style="margin-left: 3ch; font-size: 90%; color: #666; width: 10ch;" onclick="this.select();" readonly="" value='\2'></td>""",
            html,
        )

        # remove the Table length=n message
        result_len = len(result)
        html = html.replace(f"<i>Table length={result_len}</i>", "")

        return result_all_columns, result, html
    else:
        return result_all_columns, result


def search_kelt_fp_of_tics(targets, kelt_fp_csv_path="data/kelt_fps.csv"):
    if not isinstance(targets, (Sequence, np.ndarray)):
        targets = [targets]

    targets = [int(t) for t in targets]

    df = pd.read_csv(kelt_fp_csv_path)
    res = df[df["TIC"].isin(targets)]
    return res


def download_kelt_fps_as_csv(kelt_fp_csv_path="data/kelt_fps.csv"):
    v = Vizier()
    v.ROW_LIMIT = -1
    # KELT transit false positive catalog for TESS (Collins+, 2018)
    #   https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/AJ/156/234
    tb_kelt_fp = v.get_catalogs("J/AJ/156/234")[0]
    df_kelt_fp = tb_kelt_fp.to_pandas()
    df_kelt_fp.to_csv(kelt_fp_csv_path, index=False)
    return df_kelt_fp


def search_multi_star_systems_of_tics(targets, also_display=True):
    """Search various published multi-star systems with TESS data used."""
    if also_display:
        from IPython.display import display, HTML

    def display_tab_list(tab_list: TableList):
        if not also_display:
            return

        for tab_name in tab_list.keys():
            display(
                HTML(
                    f"""
<a href="https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source={tab_name}" target="_blank">{tab_name}</a>
"""
                )
            )
            display(tab_list[tab_name])

    def do_display(tab_list_dict: dict):
        if not also_display:
            return

        for k in tab_list_dict:
            display(HTML(f"<b>{k}</b>:"))
            display_tab_list(tab_list_dict[k])
        # additional sources not available at Vizier (requires manual search)
        display(
            HTML(
                """
<b>Kosotv+ 2023 (101 eclipsing quadruple stars in TESS FFI)</b>:<br>
Search <a href="https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mnras/527/2/10.1093_mnras_stad2947/2/stad2947_supplemental_file.pdf?Expires=1745801097&Signature=3Qhfwr6UR0T80alKbwOBBWMqqwSRx4Vtdhu~3o~JX1Q9bBgX~8zfKqaX4HLuM2iyoJfP4YawlfeT8t6meses5mUHJI7MFqbRxcKobsCZvJYiQQz50hhqJ2K4GP-ujT4qACaWx4KxxISoNwRjEdAe82lD1issmkSGLCR7MOUF4u9tNBdN0RqDOtTeWBIBsn9CZ~wWMi-MU-JTpr-i0Qkebd31b0oGWlh-coiB-6L06AxTaE4SJXBNlxmf8ipH7CUaGaGNurff-IPq8UEKOdRc-8zvb4aDZXVbjbG5wDb69M2HJXiGoFgQUf91uEpib7ZewBxB6VP3ID4lj7~9cpyU3Q__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA" target="_blank">table 1 pdf</a>&emsp;|&emsp;
<a href="https://academic.oup.com/mnras/article/527/2/3995/7284401#supplementary-data" target="_blank">paper</a>

<br>

<b>Zasche+ 2023 (7 2+2 doubly eclipsing quadruples)</b>:<br>
Search <a href="https://www.aanda.org/articles/aa/full_html/2023/07/aa46848-23/T1.html" target="_blank">table 1 in html</a>&emsp;|&emsp;
<a href="https://www.aanda.org/articles/aa/full_html/2023/07/aa46848-23/aa46848-23.html" target="_blank">paper</a>

<br>

<b>Zasche+ 2024 (8 new 2+2 doubly eclipsing quadruples)</b>:<br>
Search <a href="https://www.aanda.org/articles/aa/full_html/2024/07/aa50400-24/T1.html" target="_blank">table 1 in html</a>&emsp;|&emsp;
<a href="https://www.aanda.org/articles/aa/full_html/2024/07/aa50400-24/aa50400-24.html" target="_blank">paper</a>

<br>
<b>Eisner+ 2020 (Planet Hunters TESS II)</b>:<br>
Search <a href="https://academic.oup.com/view-large/225349534" target="_blank">table 2 in html</a>&emsp;|&emsp;
<a href="https://academic.oup.com/mnras/article/501/4/4669/6027708" target="_blank">paper</a>
"""
            )
        )

    def prepend_to_fixed_len(id, prefix, str_len):
        """Format a number to a fixed length string with spaces after prefix.
        E.g., 'TIC   50445777'
        Use case: J/A+A/664/A96 uses a fixed str14
        """
        id = str(id)
        num_spaces = str_len - len(id) - len(prefix)
        if num_spaces < 0:
            num_spaces = 0
        return prefix + " " * num_spaces + id

    targets = _to_tic_str_list(targets)

    res = {}

    result_list = Vizier(catalog=["J/A+A/664/A96/table1", "J/A+A/664/A96/tablea1"], columns=["*"]).query_constraints(
        TIC=[prepend_to_fixed_len(t, "TIC", 14) for t in targets]  # TICs formatted as TIC nnnnn
    )
    res["Zasche+, 2022 (Multiply eclipsing candidates from TESS)"] = result_list

    result_list = Vizier(catalog="J/ApJS/259/66/table1", columns=["*"]).query_constraints(TIC=targets)
    res["Kostov+, 2022 (97 eclipsing quadruple stars in TESS FFI)"] = result_list

    result_list = Vizier(catalog="J/A+A/683/A158/table1", columns=["*"]).query_constraints(TIC=targets)
    res["Zasche+, 2024 (6  new eccentric eclipsing systems with a third body.)"] = result_list

    # coordinate search:
    tic_rs = catalog_info_of_tics(targets)
    # set the coordinates to the columns expected by Vizier search
    tic_rs["_RAJ2000"] = tic_rs["ra"]
    tic_rs["_RAJ2000"].unit = u.deg
    tic_rs["_DEJ2000"] = tic_rs["dec"]
    tic_rs["_DEJ2000"].unit = u.deg

    result_list = Vizier(catalog="J/A+A/682/A164/table1", columns=["*", "+_r"]).query_region(tic_rs, radius=30 * u.arcsec)
    res["Vaessen+, 2024 (Double eclipsing binaries from ZTF)"] = result_list

    do_display(res)

    return res
