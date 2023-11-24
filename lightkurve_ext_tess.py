#
# Helpers to download TESS-specific non-lightcurve data: TOIs, TCEs, etc.
#

from collections.abc import Sequence
import re
from types import SimpleNamespace
import warnings

import numpy as np
import pandas as pd
from pandas.io.formats.style import Styler

import astropy
from astropy import coordinates as coord
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import astropy.units as u

import download_utils
import lightkurve as lk
import lightkurve_ext as lke
import tess_dv

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
        + '",'
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


def get_tic_meta_in_html(lc_or_tic, a_subject_id=None, download_dir=None, tce_filter_func=None):
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
    html += prop("R<sub>S</sub> (in R<sub>☉</sub>)", f'{safe_m_get("radius", -1):.3f}')
    html += prop("M<sub>S</sub> (in M<sub>☉</sub>)", f'{safe_m_get("mass", -1):.3f}')
    html += prop("Magnitude (TESS)", f'{safe_m_get("tess_mag", -1):.2f}')
    html += prop("T_eff (in K)", safe_m_get("teff", -1))
    html += "</table>\n"

    html += "<p>TCEs:</p>"
    html += tess_dv._get_tces_in_html(tic_id, download_dir=download_dir, tce_filter_func=tce_filter_func)

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
    def _load_mom_dumps_from_file(cls):
        # data/tess_mom_dumps.txt is a tab-delimited list of all TESS momentum dumps
        # (the same one used by LATTE).
        # It can be generated by reading a sample of 2-minute cadence lightcurve FITS files (1 for each sector),
        # looking for times where quality flag bit 6 (value 32) is 1.

        # cls._mom_dumps_tab = pd.read_csv("data/tess_mom_dumps.txt", sep="\t")
        # I use np.genfromtxt rather than pandas, as filtering from numpy array filtering is easier for the use case
        cls._mom_dumps_tab = np.genfromtxt("data/tess_mom_dumps.txt", delimiter="\t", names=True)

    @classmethod
    def refresh(cls):
        cls._load_mom_dumps_from_file()

    @classmethod
    def get_all(cls, refresh=False):
        if refresh or cls._mom_dumps_tab is None:
            cls.refresh()
        return cls._mom_dumps_tab

    @classmethod
    def get_in_range(cls, lc_or_tpf=None, start=None, end=None, refresh=False):
        times = cls.get_all(refresh=refresh)["time"]

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

    ra, dec = position.split(",")
    sky_coord = coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")

    return lke.to_hjd_utc(t_bjd, sky_coord).value


#
# TIC Metadata in catalogs (TIC / Gaia)
#


def catalog_info_of_tics(tic):
    """Return the info of a TIC in the TIC catalog"""
    from astroquery.mast import Catalogs

    return Catalogs.query_criteria(catalog="Tic", ID=tic)


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
    """Decode NSS (NON_SINGLE_STAR) flag in Gaia DR3.
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


def search_gaiadr3_of_tics(
    targets,
    radius_arcsec=15,
    magnitude_range=2.5,
    pm_range_fraction=0.25,
    warn_if_all_filtered=True,
    compact_columns=True,
    also_return_html=True,
    verbose_html=True,
    include_nss_summary_in_html=True,
):
    """Locate the lightcurve target's correspond entry in GaiaDR3.
    The match is by an heuristics based on coordinate and magnitude.

    Parameters
    ----------
    target : int, LightCurve, TargetPixelFile, or a list of them
        targets to be searched. Either the TIC, or LightCurve/TargetPixelFile of a TIC.

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

    result_list = []

    targets = np.asarray([_to_stellar_meta(t) for t in targets])
    targets = targets[targets != None]

    for t in targets:
        gaiadr2_id = getattr(t, "gaiadr2_id", np.nan)
        if magnitude_range is not None:
            lower_limit, upper_limit = t.tess_mag - magnitude_range, t.tess_mag + magnitude_range
        else:
            lower_limit, upper_limit = None, None

        a_result = lke.search_nearby(
            t.ra,
            t.dec,
            equinox=f"J{t.equinox}",
            radius_arcsec=radius_arcsec,
            magnitude_limit_column="RPmag",
            magnitude_lower_limit=lower_limit,
            magnitude_upper_limit=upper_limit,
            pmra=t.pmra,
            pmdec=t.pmdec,
            pm_range_fraction=pm_range_fraction,
            warn_if_all_filtered=warn_if_all_filtered,
        )

        if a_result is not None:
            a_result["target"] = [t.label for i in range(0, len(a_result))]
            a_result["target_gaia_dr2_source"] = [gaiadr2_id for i in range(0, len(a_result))]

            result_list.append(a_result)

    with warnings.catch_warnings():
        # Avoid spurious "MergeConflictWarning: Cannot merge meta key 'null' types <class 'float'>
        #  and <class 'float'>, choosing null=nan [astropy.utils.metadata]"
        result = astropy.table.vstack(result_list) if len(result_list) > 0 else []

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

    if compact_columns:  # select the most useful columns
        # prefer RA/DEC in Epoch 2000 as they can be compared more easily with those those in TESS LC metadata
        result = result[
            "target",
            "flag",
            "separation",
            "RAJ2000",
            "DEJ2000",
            "RPmag",  # prioritize RPmag, as its pass band is close to TESS
            "Gmag",
            "BPmag",
            "BP-RP",
            "Teff",  # "Tefftemp" in Gaia DR2 / Gaia EDR3
            "RUWE",
            "sepsi",
            "epsi",
            "NSS",  # Gaia DR3: 1 if there is entry in Non-Single Star tables
            "Plx",
            "pmRA",
            "pmDE",
            "GRVSmag",  # Gaia DR3 Gmag from Radial Velocity spectrometer
            "VarFlag",  # Gaia DR3: variability
            "EpochPh",  # Gaia DR3: 1 if epoch photometry is available
            "RV",  # Gaia DR3
            "e_RV",  # Gaia DR3, e_RV > 1.5 km/s also possibly signifies non single star
            "EpochRV",  # Gaia DR3
            "Dup",  # Gaia DR3: if there are multiple source/Gaia DR3 entries for the same target
            "Source",
        ]
        if not add_target_as_col:
            result.remove_column("target")

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

            with astropy.conf.set_temp("max_lines", 250):
                html = html + result._repr_html_()

        # linkify Gaia DR3 ID
        for id in result["Source"]:
            html = html.replace(
                f">{id}<",
                f"><a target='vizier_gaia_dr3' href='https://vizier.u-strasbg.fr/viz-bin/VizieR-S?Gaia%20DR3%20{id}'>{id}</a><",
            )

        # remove the Table length=n message
        result_len = len(result)
        html = html.replace(f"<i>Table length={result_len}</i>", "")

        if include_nss_summary_in_html:
            with warnings.catch_warnings():
                # ignore "Format strings passed to MaskedConstant are ignored, but in future may error or produce different behavior"
                warnings.filterwarnings(
                    "ignore", category=FutureWarning, message=".*Format strings passed to MaskedConstant are ignored,.*"
                )
                for i in range(0, len(result_all_columns)):
                    html_inner = f"RUWE: {result_all_columns[i]['RUWE']}, astrometric excess noise significance: {result_all_columns[i]['sepsi']:.3f}, e_RV: {result_all_columns[i]['e_RV']:.2f} km/s"
                    nss_flag = result_all_columns[i]["NSS"]
                    if nss_flag > 0:
                        html_inner += f"; NSS: {nss_flag} ({decode_gaiadr3_nss_flag(nss_flag)})"
                    html += f"<pre>{html_inner}</pre>"

        if verbose_html:
            html += """
<br>flag - &emsp; !: non single star proxy indicator (RUWE, sepsi, e_RV) ; &emsp; ✓: Gaia source matched with the one in TIC
<br>Reference:
    <a target ="_doc_gaiadr3_vizier" href="https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3">Column descption on Vizier</a> &nbsp; | &nbsp;
    <a target="_doc_gaiadr3_datamodel" href="https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html">data model doc on ESA</a>
"""

        return result_all_columns, result, html
    else:
        return result_all_columns, result


from astroquery.vizier import Vizier


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
