#
# Helpers to download TESS-specific non-lightcurve data: TOIs, TCEs, etc.
#

import os
from pathlib import Path
import re
import shutil
import time
from types import SimpleNamespace
import warnings

from memoization import cached
import requests
import numpy as np
import pandas as pd
from pandas.io.formats.style import Styler

# for accessing / parsing TCEs from MAST
from astroquery.exceptions import NoResultsWarning
from astroquery.mast import Observations
import xmltodict

#
# Misc constatants
#

R_earth = 6371000  # radius of the Earth [m]
R_jup = 69911000  # radius of Jupiter [m]
BTJD_REF = 2457000

#
# Generic file download, and CSV helper
#


def _policy_always_use(url, filename):
    return True


def _policy_always_reject(url, filename):
    return False


def _create_policy_ttl_in_seconds(ttl_in_seconds):
    def _policy_ttl(url, filename):
        try:
            time_since_last_modified = time.time() - os.path.getmtime(filename)
            if time_since_last_modified <= ttl_in_seconds:
                return True
            else:
                return False
        except Exception as e:
            warnings.warn(
                f"Unexpected error in determining if local file should be used. Local file is thus not used. Error: {e}",
            )
            return False

    return _policy_ttl


def _create_policy_ttl_in_days(ttl_in_days):
    return _create_policy_ttl_in_seconds(ttl_in_days * 86400)


LocalFileUsePolicy = SimpleNamespace(
    ALWAYS_USE=_policy_always_use,
    ALWAYS_REJECT=_policy_always_reject,
    TTL_IN_SECONDS=_create_policy_ttl_in_seconds,
    TTL_IN_DAYS=_create_policy_ttl_in_days,
)


def _create_local_filename(url, filename, download_dir):
    if filename is not None:
        local_filename = filename
    else:
        local_filename = url.split("/")[-1]
        local_filename = re.sub(r"\?.*$", "", local_filename)

    return os.path.join(download_dir, local_filename)


def _download_file(url, filename=None, download_dir=""):
    local_filename = _create_local_filename(url, filename, download_dir)

    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        # write to a temporary file. If successful, make it the real local file
        # it is to preven interrupted download leaving a partial file
        local_filename_temp = f"{local_filename}.download"
        with open(local_filename_temp, "wb") as out_file:
            shutil.copyfileobj(response.raw, out_file)
        os.replace(local_filename_temp, local_filename)
    return local_filename


def _download_file_if_needed(url, filename=None, download_dir="", use_localfile_func=None):
    local_filename = _create_local_filename(url, filename, download_dir)
    if os.path.isfile(local_filename):
        if use_localfile_func is None or use_localfile_func(url, local_filename):
            return local_filename

    return _download_file(url, filename, download_dir)


def _get_csv(url, filename, download_dir, use_localfile_func, **kwargs):
    local_filename = _download_file_if_needed(
        url, filename=filename, download_dir=download_dir, use_localfile_func=use_localfile_func
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
        EPOCH_BTJD="Epoch (BTJD)",
        EPOCH_BJD="Epoch (BJD)",
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

    # TODO: in-memory cache (with @cached) needs to be redone to properly support use_localfile_func
    @classmethod
    def get_all_tois(cls, download_dir="", use_localfile_func=None):
        url = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
        filename = "tess_tois.csv"
        res = _get_csv(url, filename, download_dir, use_localfile_func=use_localfile_func, dtype={cls.Headers.TOI: str})
        # add dervied columns
        res[cls.Headers.EPOCH_BTJD] = res[cls.Headers.EPOCH_BJD] - BTJD_REF
        res[cls.Headers.PLANET_RADIUS_J] = res[cls.Headers.PLANET_RADIUS_E] * R_earth / R_jup
        res[cls.Headers.DEPTH_PCT] = res[cls.Headers.DEPTH_PPM] / 10000
        return res

    def __init__(self, download_dir="", use_localfile_func=None):
        self._all = self.get_all_tois(download_dir=download_dir, use_localfile_func=use_localfile_func)

    def all(self):
        return self._all

    def of_toi(self, toi):
        tois = self._all[self._all[self.Headers.TOI] == str(toi)]
        return _single_row(tois)

    def of_tic(self, tic):
        return self._all[self._all[self.Headers.TIC] == int(tic)]


class CTOIAccessor:

    Headers = SimpleNamespace(
        TIC="TIC ID",
        CTOI="CTOI",
        TOI="Promoted to TOI",
        EPOCH="Epoch (BJD)",
        PERIOD="Period (days)",
        DURATION_HR="Duration (hrs)",
        DEPTH_PPM="Depth ppm",
        DEPTH_PCT="Depth percent",  # derived
        PLANET_RADIUS_E="Radius (R_Earth)",
        PLANET_RADIUS_J="Radius (R_Jupiter)",  # derived
        COMMENTS="Notes",
    )

    @classmethod
    def get_all_ctois(cls, download_dir="", use_localfile_func=None):
        url = "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv"
        filename = "tess_ctois.csv"
        res = _get_csv(
            url,
            filename,
            download_dir,
            use_localfile_func=use_localfile_func,
            dtype={cls.Headers.CTOI: str, cls.Headers.TOI: str},
        )
        # add dervied columns
        res[cls.Headers.PLANET_RADIUS_J] = res[cls.Headers.PLANET_RADIUS_E] * R_earth / R_jup
        res[cls.Headers.DEPTH_PCT] = res[cls.Headers.DEPTH_PPM] / 10000
        return res

    def __init__(self, download_dir="", use_localfile_func=None):
        self._all = self.get_all_ctois(download_dir=download_dir, use_localfile_func=use_localfile_func)

    def all(self):
        return self._all

    def of_ctoi(self, ctoi):
        ctois = self._all[self._all[self.Headers.CTOI] == str(ctoi)]
        return _single_row(ctois)

    def of_tic(self, tic):
        return self._all[self._all[self.Headers.TIC] == int(tic)]


#
# Download and parse TCEs
#


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


#
# Top-level report logic, which provides
# - stellar metadata (in LightCurve file)
# - TCEs
# - TOIs
#


def _tce_info_to_html(tce_info_list):
    if len(tce_info_list) < 1:
        return ""

    def link(link_text, url):
        return f"""<a href="{url}" target="_blank">{link_text}</a>"""

    def row(*args):
        return "<tr>" + "".join(f"<td>{v}</td>" for v in args) + "</tr>"

    html = ""
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
    onclick="this.select();" readonly
    value='epoch={p_i.get("transitEpochBtjd", 0):.4f}, duration_hr={p_i.get("transitDurationHours", 0):.4f}, \
period={p_i.get("orbitalPeriodDays", 0):.6f}, label="{info.get("tce_id_short")}",'>""",
        )
        html += "\n"

    html += "</tbody></table>\n"
    return html


def add_codes_column_to_toi_df(df, headers):
    h = headers
    # string interpolation does not work. So use old-school concatenation
    df["Codes"] = (
        "epoch="
        + df[h.EPOCH_BTJD].map("{:.4f}".format)
        + ", duration_hr="
        + df[h.DURATION_HR].map("{:.4f}".format)
        + ", period="
        + df[h.PERIOD].map("{:.6f}".format)
        + ', label="TOI '
        + df[h.TOI].astype(str)
        + '",'
    )
    return df


def _get_tois_in_html(tic, download_dir=None):
    h = TOIAccessor.Headers
    # Consider cache TOIAccessor in some module global (keyed by download_dir) to avoid
    # repeated loading/parsing the underlying TOI csv
    tois = TOIAccessor(download_dir=download_dir).of_tic(tic)
    if len(tois) < 1:
        return "<p>No TOIs.</p>"

    add_codes_column_to_toi_df(tois, h)
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
    styler.hide_index()
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

    # turn Codes column into html input element (easier to be selected)
    html = re.sub(
        r"<td([^>]+)>(epoch=.+,)</td>",
        r"""<td\1><input  type="text" style="margin-left: 3ch; font-size: 90%; color: #666; width: 10ch;" onclick="this.select();" readonly="" value='\2'></td>""",
        html,
    )

    return html


def get_tic_meta_in_html(lc, a_subject_id=None, download_dir=None):
    # This function does not do the actual display,
    # so that the caller can call it in background
    # and display it whereever it's needed
    def link(link_text, url):
        return f"""<a href="{url}" target="_blank">{link_text}</a>"""

    def prop(prop_name, prop_value):
        return f"""    <tr><td>{prop_name}</td><td>{prop_value}</td></tr>\n"""

    # main logic
    m = lc.meta
    tic_id = str(m.get("TICID"))

    def safe_m_get(key, default_val):
        # in some meta, the key exists but the value is None
        # this helper handles it
        res = m.get(key, default_val)
        return res if res is not None else default_val

    html = f"""
<div id="tic_metadata_ctr">
<div id="tic_metadata_ctl">
    <span id="float_expand_toggle" title="Toggle whether the metadata is shown or not"></span>
    <span id="float_fixed_toggle" title="Toggle whether the metadata is shown in a floating box or a regular cell"></span>
</div>
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
        html += f' (sector {safe_m_get("SECTOR", "")})'
    html += "<br>\n"
    html += "<table>\n"
    html += prop("R<sub>S</sub> (in R<sub>☉</sub>)", f'{safe_m_get("RADIUS", 0):.3f}')
    html += prop("Magnitude (TESS)", f'{safe_m_get("TESSMAG", 0):.2f}')
    html += prop("T_eff (in K)", safe_m_get("TEFF", 0))
    html += "</table>\n"

    # For TCE, query MAST download / parse results (the _dvr.xml), then show:
    # - basic planet parameters and orbital info
    # - TODO: red flags in vetting report
    # see: https://archive.stsci.edu/missions-and-data/tess/data-products
    tce_info_list = get_tce_infos_of_tic(tic_id, download_dir=download_dir)
    html += _tce_info_to_html(tce_info_list)

    # TOIs/CTOIs
    html += "<p>TOIs / CTOIs:</p>"
    html += _get_tois_in_html(tic_id, download_dir=download_dir)

    html += """
</div> <!-- id="tic_metadata_body" -->
</div> <!--  id="tic_metadata_ctr" -->
<style id="tic_metadata_out_style">
    #tic_metadata_ctr.float {
        position: fixed;
        bottom: 12px;
        right: 36px;
        z-index: 999;
        background-color: rgba(255, 255, 0, 0.3);
        padding: 6px;
        max-height: 75vh;  /* ensure for TIC with large number of TCEs, the floating box won't fill up too much space */
        overflow-y: scroll;
    }

    #tic_metadata_ctl {
        margin-left: 12em; /* make it roughly to the left of TIC heading */
    }
    #tic_metadata_ctr.float #tic_metadata_ctl {
        margin-left: 0;
        float: right;
    }

    #float_fixed_toggle {
        cursor: pointer;
        padding: 6px;
        font-size: 16px;
        font-weight: normal;
    }
    #float_fixed_toggle:before {
        content: "[To float >]";
    }
    #tic_metadata_ctr.float #float_fixed_toggle:before {
        content: "[X]";
    }

    #float_expand_toggle {
        cursor: pointer;
        padding: 6px;
        font-size: 16px;
        font-weight: normal;
        margin-left: 10px;
    }

    #tic_metadata_ctr.float #float_expand_toggle:before {
        content: "<<";
    }

    #tic_metadata_ctr.float.expand #float_expand_toggle:before {
        content: ">>";
    }

    #tic_metadata_ctr.float #tic_metadata_body {
        display: none;
    }

    #tic_metadata_ctr.float.expand #tic_metadata_body {
        display: block;
    }

</style>
<script>
    document.getElementById("float_fixed_toggle").onclick = function(evt) {
        const ctr = document.getElementById("tic_metadata_ctr");
        ctr.classList.toggle("float");
        if (ctr.classList.contains("float")) {
            ctr.classList.add("expand");
        }
};
    document.getElementById("float_expand_toggle").onclick = function(evt) {
        document.getElementById("tic_metadata_ctr").classList.toggle("expand");
};
</script>
"""
    return html