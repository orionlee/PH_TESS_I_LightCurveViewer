from collections import OrderedDict
import logging
from pathlib import Path
import re
import warnings

from memoization import cached
import numpy as np
import pandas as pd
import xmltodict

from astroquery.exceptions import NoResultsWarning
from astroquery.mast import Observations


#
# Search, download and parse TCEs
#


def parse_dvs_filename(filename):
    # e.g.: tess2020267090513-s0030-s0030-0000000142087638-01-00394_dvs.pdf
    # name pattern reference: https://archive.stsci.edu/missions-and-data/tess/data-products.html#name_schema
    match = re.match(r"^tess\d+-s(\d+)-s(\d+)-(\d+)-(\d+)-(\d+)_dvs[.]pdf", filename)
    if not match:
        return {}
    sector_start, sector_stop, tic_id_padded, tce_num_padded, pipeline_run_padded = (
        match.group(1),
        match.group(2),
        match.group(3),
        match.group(4),
        match.group(5),
    )
    sector_range = f"s{sector_start}-s{sector_stop}"
    tic_id = int(tic_id_padded)  # it's a 0-padded string
    tce_num = int(tce_num_padded)  # it's a 0-padded string
    # sufficient to identify one for a given TIC, less visually busy
    tce_id_short = f"{sector_range}:TCE{tce_num}"

    # tce_id is the format used on ExoMAT, e.g,  TIC142087638S0030S0030TCE1
    # TODO: unclear on the format of the sector_range portion.
    #       In those cases, there is dash "-" in the range, e.g, instead of
    #       S0039S0039, it ie S0039-S0039
    tce_id = f"""TIC{tic_id}{re.sub("-", "", sector_range.upper())}TCE{tce_num}"""

    # pipeline_run: occasionally, a TIC in a sector has multiple dvs, they are differentiated by pipeline_run
    pipeline_run = int(pipeline_run_padded)  # it's a 0-padded string

    # convert sector start/stop to int
    sector_start, sector_stop = int(sector_start), int(sector_stop)

    return dict(
        tce_id=tce_id,
        tce_id_short=tce_id_short,
        sector_range=sector_range,
        sector_range_start=sector_start,
        sector_range_stop=sector_stop,
        sector_range_span=1 + (sector_stop - sector_start),
        tic_id=tic_id,
        tce_num=tce_num,
        pipeline_run=pipeline_run,
    )


def parse_dvr_filename(filename):
    match = re.match(r"^tess\d+-(s\d+-s\d+)-(\d+)-(\d+)_dvr[.](pdf|xml)", filename)
    if not match:
        return {}
    sector_range, tic_id_padded, pipeline_run_padded, file_type = (
        match.group(1),
        match.group(2),
        match.group(3),
        match.group(4),
    )
    tic_id = int(tic_id_padded)  # it's a 0-padded string
    pipeline_run = int(pipeline_run_padded)  # it's a 0-padded string

    return dict(sector_range=sector_range, tic_id=tic_id, file_type=file_type, pipeline_run=pipeline_run)


def get_dv_products_of_tic(tic_id, productSubGroupDescription):
    # Based on:
    # - https://outerspace.stsci.edu/display/TESS/7.0+-+Tips+and+Tricks+to+Getting+TESS+Data+At+MAST
    # https://github.com/spacetelescope/notebooks/blob/master/notebooks/MAST/TESS/beginner_astroquery_dv/beginner_astroquery_dv.ipynb

    # Note: for TESS, tic_id (the number without TIC) is what an exact match works
    # Kepler / K2 ids will need some additional processing for exact match to work.
    exact_target_name = tic_id
    with warnings.catch_warnings():
        # to filter WARNING: NoResultsWarning: Query returned no results. [astroquery.mast.discovery_portal]
        warnings.filterwarnings("ignore", category=NoResultsWarning, message=".*Query returned no results.*")
        obs_wanted = Observations.query_criteria(
            target_name=exact_target_name,
            dataproduct_type="timeseries",
            obs_collection="TESS",
        )
        data_products = Observations.get_product_list(obs_wanted)
        res = Observations.filter_products(data_products, productSubGroupDescription=productSubGroupDescription)
        # somehow astype() does not work, it raised "ValueError: invalid literal for int() with base 10: 'N/A'"
        # from numpy, even though there is no missing value in obsID column
        # res["obsID"] = res["obsID"].astype(int)
        res["obsID"] = [int(v) for v in res["obsID"]]
        return res


def parse_dvr_xml(file_path):
    def as_list(data):
        """Wrap an item as a list, if it's not one.
        Useful for handling dict from XML where elements might be one or multiple elements"""
        if type(data) is list:
            return data
        else:
            return [data]

    def param_value(model_params_dict, param_name, attr="value"):
        param_dict = model_params_dict.get(param_name)
        if param_dict is None:
            return None
        val_str = param_dict.get(f"@{attr}")
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
        #  planet / transit parameters
        #
        e_afit = e_pr["dv:allTransitsFit"]
        planet_num = int(e_afit["@planetNumber"])

        params_dict = {}  # a temporary structure to access planet params internally
        for mp in e_afit["dv:modelParameters"]["dv:modelParameter"]:
            params_dict[mp["@name"]] = mp

        # use the underlying xml attribute names, even thought it breaks  Python convention
        a_planet_dict = dict(
            planetNumber=planet_num,
            transitEpochBtjd=param_value(params_dict, "transitEpochBtjd"),
            planetRadiusEarthRadii=param_value(params_dict, "planetRadiusEarthRadii"),
            transitDurationHours=param_value(params_dict, "transitDurationHours"),
            orbitalPeriodDays=param_value(params_dict, "orbitalPeriodDays"),
            transitDepthPpm=param_value(params_dict, "transitDepthPpm"),
            minImpactParameter=param_value(params_dict, "minImpactParameter"),
        )

        # centroid offsets, under <dv:centroidResults><dv:differenceImageMotionResults> element
        #
        # the TicOffset-rm in pdf
        e_centroid_tic = e_pr["dv:centroidResults"]["dv:differenceImageMotionResults"]["dv:msTicCentroidOffsets"]
        meanSkyOffsetTic = param_value(e_centroid_tic, "dv:meanSkyOffset")
        meanSkyOffsetErrTic = param_value(e_centroid_tic, "dv:meanSkyOffset", attr="uncertainty")
        meanSkyOffsetSigTic = meanSkyOffsetTic / meanSkyOffsetErrTic  # in sigma

        # the OotOffset-rm in pdf
        e_centroid_oot = e_pr["dv:centroidResults"]["dv:differenceImageMotionResults"]["dv:msControlCentroidOffsets"]
        meanSkyOffsetOot = param_value(e_centroid_oot, "dv:meanSkyOffset")
        meanSkyOffsetErrOot = param_value(e_centroid_oot, "dv:meanSkyOffset", attr="uncertainty")
        meanSkyOffsetSigOot = meanSkyOffsetOot / meanSkyOffsetErrOot  # in sigma

        a_planet_dict.update(
            # for centroid offset, the primary interest is SkyOffset and its sigma.
            # Ra/Dec offset is added for potential usage in the future, e.g., if the centroid RA/Dec is needed.
            dict(
                meanSkyOffsetTic=meanSkyOffsetTic,
                meanSkyOffsetErrTic=meanSkyOffsetErrTic,
                meanSkyOffsetSigTic=meanSkyOffsetSigTic,
                meanRaOffsetTic=param_value(e_centroid_tic, "dv:meanRaOffset"),
                meanRaOffsetErrTic=param_value(e_centroid_tic, "dv:meanRaOffset", attr="uncertainty"),
                meanDecOffsetTic=param_value(e_centroid_tic, "dv:meanDecOffset"),
                meanDecOffsetErrTic=param_value(e_centroid_tic, "dv:meanDecOffset", attr="uncertainty"),
                meanSkyOffsetOot=meanSkyOffsetOot,
                meanSkyOffsetErrOot=meanSkyOffsetErrOot,
                meanSkyOffsetSigOot=meanSkyOffsetSigOot,
                meanRaOffsetOot=param_value(e_centroid_oot, "dv:meanRaOffset"),
                meanRaOffsetErrOot=param_value(e_centroid_oot, "dv:meanRaOffset", attr="uncertainty"),
                meanDecOffsetOot=param_value(e_centroid_oot, "dv:meanDecOffset"),
                meanDecOffsetErrOot=param_value(e_centroid_oot, "dv:meanDecOffset", attr="uncertainty"),
            )
        )

        # TODO: add other DV fitting parameters: odd/even, dv:weakSecondary , suspectedEclipsingBinary, etc.

        planets_dict[planet_num] = a_planet_dict

    return planets_dict


def get_tce_minimal_infos_of_tic(tic_id, also_return_dvr_xml_table=True):
    """Get the list of TCEs of the given TIC with minimal information: links to the underlying pdfs / XMLs, etc.
    No download is done.
    """
    # recent astroquery emits INFO messages by default, e.g., in downloading cached dvr xml files
    # suppress them
    logging.getLogger("astroquery").setLevel(logging.WARNING)

    def filter_by_dataURI_suffix(products, suffix):
        # Helper to filter products into summary, full report, full report xml using suffix.
        # It replaces the logic to filter by "description" column, as description is sometimes unreliable
        # E.g., for the TCE for TIC 43843023 sector 5, the dvr xml has incorrect description
        # so that the entry is treated as a dvr pdf
        return products[np.char.endswith(products["dataURI"], suffix)]

    def get_existing_entry_of_obsID(entries, obsID, tce_num=None):
        # help to identify duplicates for a given TIC + sector
        # (obsID uniquely identifies a TIC + sector)
        # Specify `tce_num` to only remove duplicates for a specific planet
        # (for cases that multiple planet candidates for a given TIC + sector)
        for e in entries:
            if e["obsID"] == obsID:
                if tce_num is None or e["tce_num"] == tce_num:
                    return e
                # else the same obsID but different tce_num, i.e., different planet candidate
                # continue to search
        return None

    def remove_duplicates(tab):
        # The input table should be is assumed to be for 1 type of products, e.g., dvr_xml
        def get_pipeline_run(filename):
            return parse_dvr_filename(filename)["pipeline_run"]

        tab = tab.copy()
        tab["priority_key"] = [-get_pipeline_run(f) for f in tab["productFilename"]]
        tab.sort(["obsID", "priority_key"], reverse=False)

        # create an empty table for results, with the same set of columns
        res_t = tab[np.zeros(len(tab), dtype=bool)].copy()

        # for each obsID, select the one with the largest pipeline_run (implemented as the smallest priority_key)
        uniq_obsIDs = list(OrderedDict.fromkeys(tab["obsID"]))
        for oid in uniq_obsIDs:
            res_t.add_row(tab[tab["obsID"] == oid][0])
        return res_t

    products_wanted = get_dv_products_of_tic(tic_id, ["DVS", "DVR", "DVM"])

    res = []
    # basic info
    for p in filter_by_dataURI_suffix(products_wanted, "_dvs.pdf"):
        tce_info = parse_dvs_filename(p["productFilename"])
        existing_entry = get_existing_entry_of_obsID(res, p["obsID"], tce_num=tce_info.get("tce_num"))
        if existing_entry is None:
            entry = dict(
                obsID=p["obsID"],
                tic_id=tce_info.get("tic_id"),
                sector_range=tce_info.get("sector_range"),
                sector_range_start=tce_info.get("sector_range_start"),
                sector_range_stop=tce_info.get("sector_range_stop"),
                sector_range_span=tce_info.get("sector_range_span"),
                tce_num=tce_info.get("tce_num"),
                tce_id=tce_info.get("tce_id"),
                tce_id_short=tce_info.get("tce_id_short"),
                pipeline_run=tce_info.get("pipeline_run"),
                dvs_dataURI=p["dataURI"],
            )
            res.append(entry)
        else:
            # case Multiple DVS for a TCE. We use a heuristics to retain the one with the largest pipeline_run number.
            if existing_entry.get("pipeline_run") < tce_info.get("pipeline_run"):
                warnings.warn(
                    f"""get_tce_minimal_infos_of_tic(): Multiple DVS for {existing_entry["tce_id_short"]}. """
                    f"""Discard pipeline_run {existing_entry.get("pipeline_run")}."""
                )
                existing_entry["pipeline_run"] = tce_info["pipeline_run"]
                existing_entry["dvs_dataURI"] = p["dataURI"]
            else:
                warnings.warn(
                    f"""get_tce_minimal_infos_of_tic(): Multiple DVS for {existing_entry["tce_id_short"]}. """
                    f"""Discard pipeline_run {tce_info.get("pipeline_run")}."""
                )

    # DVM pdf link
    for p in filter_by_dataURI_suffix(products_wanted, "_dvm.pdf"):
        # find TCEs for the same observation (sometimes there are multiple TCEs for the same observation)
        # Consideration for multiple DVM for a given TCE: we should explicitly pick one, similar to what is done for DVS.
        # we skip it in practice here, because given we choose the one with the largest pipeline run,
        # it generally is the last one in the list, i.e., it will be the one written to entry["dvm_dataURI"]
        for entry in [e for e in res if e["obsID"] == p["obsID"]]:
            entry["dvm_dataURI"] = p["dataURI"]

    # DVR pdf link
    for p in filter_by_dataURI_suffix(products_wanted, "_dvr.pdf"):
        # find TCEs for the same observation (sometimes there are multiple TCEs for the same observation)
        # Consideration for multiple DVR for a given TCE: we should explicitly pick one, similar to what is done for DVS.
        # we skip it in practice here, because given we choose the one with the largest pipeline run,
        # it generally is the last one in the list, i.e., it will be the one written to entry["dvr_dataURI"]
        for entry in [e for e in res if e["obsID"] == p["obsID"]]:
            entry["dvr_dataURI"] = p["dataURI"]

    products_dvr_xml = filter_by_dataURI_suffix(products_wanted, "_dvr.xml")
    # for single TCE+sector with multiple pipeline runs, only retain the one with the largest pipeline_run
    # (analogous to what's done for dvs above)
    products_dvr_xml = remove_duplicates(products_dvr_xml)

    if also_return_dvr_xml_table:
        return res, products_dvr_xml
    else:
        return res


def add_info_from_tce_xml(tce_infos, products_dvr_xml, download_dir=None):
    """Download and prase relevant TCE XML files, and add extracted parameters to the given TCE info list."""

    res = tce_infos

    # recent astroquery emits INFO messages by default, e.g., in downloading cached dvr xml files
    # suppress them
    logging.getLogger("astroquery").setLevel(logging.WARNING)

    # remove product rows that are not needed (not specified in tce_infos)
    # to avoid unnecessary download/ parsing
    tce_info_obsIDs = [r["obsID"] for r in tce_infos]
    products_dvr_xml = products_dvr_xml[np.isin(products_dvr_xml["obsID"], tce_info_obsIDs)]

    with warnings.catch_warnings():
        # filter WARNING: NoResultsWarning: No products to download. [astroquery.mast.observations]
        warnings.filterwarnings("ignore", category=NoResultsWarning, message=".*No products to download.*")
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


@cached
def get_tce_infos_of_tic(tic_id, tce_filter_func=None, download_dir=None):
    """Get the list of TCES of a given TIC, including downloading and extract some parameters from the correspond DVR XMLs."""
    res, products_dvr_xml = get_tce_minimal_infos_of_tic(tic_id, also_return_dvr_xml_table=True)
    if tce_filter_func is not None:
        res = tce_filter_func(res)
    return add_info_from_tce_xml(res, products_dvr_xml, download_dir=download_dir)


def filter_top_2_tces_for_eb(tce_infos):
    """Given a list of TCEs, try to get top 2 for eclipsing binary use case."""
    # convert tce_infos to a DataFrame for ease of filtering
    if len(tce_infos) < 1:
        return []

    df = pd.DataFrame(tce_infos)
    num_tics = df["tic_id"].nunique()
    if num_tics > 1:
        raise ValueError(f"The tce_infos list should be for a TIC. There are {num_tics}.")

    df.sort_values(["sector_range_span", "obsID", "tce_num"], ascending=[False, True, True], inplace=True)
    if len(df) > 1 and df["obsID"].iloc[1] == df["obsID"].iloc[0]:
        # we always return the top one, we return the second one only if it is from the same TIC/sector (i.e., obsID)
        # it might capture shallow eclipses, if the top TCE capture the deep eclipses; or vice versa
        df = df[:2]
    else:
        df = df[:1]

    return df.to_dict(orient="records")  # convert the filtered result back to a list of dict


#
# Top-level report logic, which provides
# - stellar metadata (in LightCurve file)
# - TCEs
# - TOIs
#


def _tce_info_to_html(tce_info_list):
    if len(tce_info_list) < 1:
        return "No TCEs."

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
        ("TicOffset", "σ"),
        ("OotOffset", "σ"),
        ("Codes", ""),
    ]
    html += """<table>
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
        dvs_url = f'https://mast.stsci.edu/api/v0.1/Download/file?uri={info.get("dvs_dataURI")}'
        dvm_url = f'https://mast.stsci.edu/api/v0.1/Download/file?uri={info.get("dvm_dataURI")}'
        dvr_url = f'https://mast.stsci.edu/api/v0.1/Download/file?uri={info.get("dvr_dataURI")}'
        p_i = info.get("planet", {})
        html += row(
            link(info.get("tce_id_short"), exomast_url),
            f"""{link("dvs", dvs_url)},&emsp;{link("mini", dvm_url)},&emsp;{link("full", dvr_url)}""",
            f'{p_i.get("planetRadiusEarthRadii", 0) * R_EARTH_TO_R_JUPITER:.3f}',
            f'{p_i.get("transitEpochBtjd", 0):.4f}',
            f'{p_i.get("transitDurationHours", 0):.4f}',
            f'{p_i.get("orbitalPeriodDays", 0):.6f}',
            f'{p_i.get("transitDepthPpm", 0) / 10000:.4f}',
            f'{p_i.get("minImpactParameter", 0):.2f}',
            f'{p_i.get("meanSkyOffsetSigTic", -1):.2f}',
            f'{p_i.get("meanSkyOffsetSigOot", -1):.2f}',
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


def _get_tces_in_html(tic, download_dir=None, tce_filter_func=None):
    # For TCEs, query MAST download / parse results (the _dvr.xml), then show:
    # - basic planet parameters and orbital info
    # - TODO: red flags in vetting report
    # see: https://archive.stsci.edu/missions-and-data/tess/data-products
    tce_info_list = get_tce_infos_of_tic(tic, download_dir=download_dir, tce_filter_func=tce_filter_func)
    return _tce_info_to_html(tce_info_list)
