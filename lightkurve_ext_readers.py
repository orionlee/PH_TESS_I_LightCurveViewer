#
# Misc. lightkurve readers for different formats
#

import re
import urllib

from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy import units as u

import numpy as np

import lightkurve as lk


def read_asas_sn_csv(url=None, asas_sn_uuid=None):
    """Read ASAS-SN lightcurve CSV data."""
    # csv url format: https://asas-sn.osu.edu/variables/4bccdf4d-f98c-5697-ad69-4199e5cf3905.csv
    if url is None:
        if asas_sn_uuid is not None:
            url = f"https://asas-sn.osu.edu/variables/{asas_sn_uuid}.csv"
        else:
            raise ValueError("either url or asas_sn_uuid must be supplied")
    else:
        match = re.match(r".*asas-sn.osu.edu/([^/]+([^/]+)?)/([^/]+)[.]csv", url)
        if match is not None:
            asas_sn_uuid = match[2]

    # print("DBG1: url=", url)
    tab = Table.read(
        url,
        format="ascii.csv",
        comment="#",  # SkyPatrol v2 tweak: it has comments begin with "#" (metadata like info)
    )
    # print("DBG2: colnames=", tab.colnames)

    # SkyPatrol v2 tweak: lower case the column name, the logic does not hurt v1 data
    for c in tab.colnames:
        tab.rename_column(c, c.lower())
    # SkyPatrol v2 tweak: rename columns to be v1 like
    for colname, colname_new in [("jd", "hjd"), ("flux error", "flux_err"), ("mag error", "mag_err")]:
        if colname in tab.colnames:
            tab.rename_column(colname, colname_new)

    # make columns follow Lightkurve convention
    if "flux(mjy)" in tab.colnames:
        # column flux(mJy) shows up in csv from path "/variables/" or "/sky-patrol/coordinate/"
        # for csv from path "/photometry/" , the column is simply flux
        tab.rename_column("flux(mjy)", "flux")

    for c in tab.colnames:
        tab.rename_column(c, re.sub(" ", "_", c.lower()))

    tab.sort(keys="hjd")

    tab["mag"].unit = u.mag
    tab["mag_err"].unit = u.mag
    # OPEN: unit of flux is actually mJy. I don't know what it means though
    tab["flux"].unit = u.dimensionless_unscaled
    tab["flux_err"].unit = u.dimensionless_unscaled

    # SkyPatrol csv contains rows with mag / flux as 99.99, presumably the underlying data is bad
    tab = tab[tab["mag"] < 99.99]  # for SkyPatrol v1
    tab = tab[tab["mag_err"] < 99.99]  # for SkyPatrol v2

    lc = lk.LightCurve(time=Time(tab["hjd"], format="jd", scale="utc"), data=tab)

    lc.meta["FILEURL"] = url

    if asas_sn_uuid is None:  # deduce ASAS-SN uuid from url
        uuid_match = re.search(r"[a-zA-Z0-9]{8}-[a-zA-Z0-9]{4}-[a-zA-Z0-9]{4}-[a-zA-Z0-9]{4}-[a-zA-Z0-9]{8}", url)
        if uuid_match is not None:
            asas_sn_uuid = uuid_match[0].lower()

    if asas_sn_uuid is not None:
        lc.meta["ASAS-SN_UUID"] = asas_sn_uuid
        label = f"ASAS-SN {asas_sn_uuid}"
        if "camera" in lc.colnames:
            cameras = np.unique(lc["camera"])
            if len(cameras) > 1:
                label += f" , cameras: {', '.join(cameras)}"
        lc.meta["LABEL"] = label

    return lc


def read_superwasp_dr1_data(superwasp_id):
    """Read SuperWASP DR1 data, in the format (FITS + CSV) available at
    https://wasp.cerit-sc.cz/
    """
    from urllib.parse import quote_plus

    fits_url = f"https://wasp.cerit-sc.cz/download?object={quote_plus(superwasp_id)}"
    # camera data is not in FITs but can be found in CSV
    csv_url = f"https://wasp.cerit-sc.cz/csv?object={quote_plus(superwasp_id)}"

    with fits.open(fits_url) as hdul:
        jd_ref = hdul[0].header.get("JD_REF")
        # unsure if the HJD here is in UTC
        time = Time(hdul[1].data["TMID"] / 86400 + jd_ref, format="jd", scale="utc")

        lc = lk.LightCurve(time=time)
        e_s = u.electron / u.s
        # TAMFLUX2 is corrected flux. See: https://exoplanetarchive.ipac.caltech.edu/docs/SuperWASPProcessing.html
        flux_column = "TAMFLUX2"

        # the unit of flux in micro-vega is basically electron per sec
        # https://www.zooniverse.org/projects/ajnorton/superwasp-variable-stars/about/faq
        lc["flux"] = hdul[1].data[flux_column] * e_s
        lc["flux_err"] = hdul[1].data[f"{flux_column}_ERR"] * e_s
        lc.meta["FLUX_ORIGIN"] = flux_column

        lc["flux2"] = hdul[1].data["FLUX2"] * e_s
        lc["flux2_err"] = hdul[1].data["FLUX2_ERR"] * e_s
        lc["tamflux2"] = hdul[1].data["TAMFLUX2"] * e_s
        lc["tamflux2_err"] = hdul[1].data["TAMFLUX2_ERR"] * e_s

        lc["ccdx"] = hdul[1].data["CCDX"]
        lc["ccdy"] = hdul[1].data["CCDY"]

        # need the following columns to make interact_bls() works
        lc["quality"] = hdul[1].data["FLAG"]
        lc["cadenceno"] = hdul[1].data["TMID"]  # a fake cadenceno

        # misc. that might be useful
        lc["imageid"] = hdul[1].data["IMAGEID"]

        lc.meta.update(hdul[0].header)
        lc.meta["LABEL"] = lc.meta["OBJNAME"]  # useful in plots

        lc.meta["FILEURL"] = fits_url

        # get camera data from csv data
        if csv_url is not None:
            csv = Table.read(csv_url, format="ascii")
            lc["camera"] = csv["camera"]
            lc.meta["FILEURL2"] = csv_url

        return lc


def read_astroimagej_tbl(
    url=None,
    time_column="BJD_TDB",
    time_format="jd",
    time_scale="tdb",
    flux_column="rel_flux_T1",
    flux_err_column="rel_flux_T1",
):
    """Read AstroImageJ measurement export .tbl file as LightCurve"""
    import pandas as pd

    # use pandas rather than Table.read(), because
    # Table.read() does not handle the header line correctly
    # AstroImageJ tbl's first column is a 1-based index that " " as column name
    # Table.read() does not recognize it as a column in the header
    df = pd.read_csv(url, sep="\t")

    df.rename(columns={" ": "recordno"}, inplace=True)

    if time_column not in df.columns:
        raise ValueError(f"time column {time_column} not found")

    time = Time(df[time_column], format=time_format, scale=time_scale)
    df.drop(columns=[time_column], inplace=True)

    flux = df[flux_column]
    flux_err = df[flux_err_column] if flux_err_column is not None else None

    lc = lk.LightCurve(
        time=time,
        flux=flux,
        flux_err=flux_err,
        data=Table.from_pandas(
            df,
        ),
    )
    lc.meta.update(
        {
            "FILEURL": url,
            "FLUX_ORIGIN": flux_column,
        }
    )

    return lc


def read_ztf_csv(
    url=None,
    time_column="hjd",
    time_format="jd",
    time_scale="utc",
    flux_column="mag",
    flux_err_column="magerr",
    mask_func=lambda lc: lc["catflags"] != 0,
):
    """Return ZTF Archive lightcurve files in IPAC Table csv.

    Parameters
    ----------
    mask_func : function, optional
        a function that returns a boolean mask given a `Lightcurve` object
        of the data. Cadences with `True` will be masked out.
        Pass `None` to disable masking.
        The default is to exclude cadences where `catflags` is not 0, the
        guideline for VSX submission.
        https://www.aavso.org/vsx/index.php?view=about.notice
    """
    # Note: First tried to read ZTF's ipac table .tbl, but got
    #   TypeError: converter type does not match column type
    # due to: https://github.com/astropy/astropy/issues/15989

    def get_required_column(tab, colname):
        if colname not in tab.colnames:
            raise ValueError(f"Required column {colname} is not found")
        return tab[colname]

    tab = Table.read(
        url,
        format="ascii.csv",
        converters={
            "oid": np.int64,
            "expid": np.int64,
            "filefracday": np.int64,
        },
    )

    time = get_required_column(tab, time_column)
    time = Time(time, format=time_format, scale=time_scale)
    flux = get_required_column(tab, flux_column)
    flux_err = get_required_column(tab, flux_err_column)

    lc = lk.LightCurve(
        time=time,
        flux=flux,
        flux_err=flux_err,
        data=tab,
    )

    # assign units
    for col in ["flux", "flux_err", flux_column, flux_err_column, "limitmag"]:
        if col in lc.colnames:
            lc[col] = lc[col] * u.mag

    lc.meta.update(
        {
            "FILEURL": url,
            "FLUX_ORIGIN": flux_column,
            "TIME_ORIGIN": time_column,
        }
    )

    if mask_func is not None:
        mask = mask_func(lc)
        lc = lc[~mask]

    return lc


def read_hipparcos_data(url):
    """Read Hipparcos and Tycho Epoch Photometry hosted on Vizier.
    https://cdsarc.cds.unistra.fr/viz-bin/VizieR?-meta&-meta.ucd&-source=I/239
    """

    def get_hip_id(url):
        """Extract HIP ID from Vizier's data URL"""
        # e.g.,
        #  https://cdsarc.cds.unistra.fr/viz-bin/nph-Plot/Vgraph/txt?I%2f239%2f.%2f17076&0&P=0&-Y&mag&-y&-&-&-
        #  https://cdsarc.cds.unistra.fr/viz-bin/nph-Plot/Vgraph/txt?I/239/17076
        matches = re.match(".*txt[?]([^&]+)", url)
        if matches is None:
            return None
        param = urllib.parse.unquote(matches[1])
        if param.startswith("I/239/"):  # I/239/ or I/239/.
            try:
                hip_id = int(re.sub("^I/239/([.]/)?", "", param))
                return hip_id
            except Exception:
                # Unexpected pattern
                return None
        else:
            # unrecognized pattern
            return None

    tab = Table.read(url, format="ascii")
    time = tab["col1"]  # JD-2440000
    # Time in TT scale
    # reference: EAS doc "The Hipparcos and Tycho Catalogues", sections 1.2.3 and 1.2.6
    # https://www.cosmos.esa.int/documents/532822/552851/vol1_all.pdf/99adf6e3-6893-4824-8fc2-8d3c9cbba2b5
    time = Time(time + 2440000, format="jd", scale="tt")
    flux = tab["col2"]  # mag
    flux = flux * u.mag
    flux_err = tab["col3"]  # (error)
    flux_err = flux_err * u.mag

    lc = lk.LightCurve(time=time, flux=flux, flux_err=flux_err)
    lc.meta["FILEURL"] = url

    hip_id = get_hip_id(url)
    if hip_id is not None:
        # mimic convention used in TESS SPOC Lightcurve objects
        lc.meta["TARGETID"] = hip_id
        lc.meta["OBJECT"] = f"HIP {hip_id}"
        lc.meta["LABEL"] = f"HIP {hip_id}"
    return lc


def read_asas3(url, flux_column="mag_3", grade_mask=["A", "B", "C"]):
    # https://www.astrouw.edu.pl/cgi-asas/asas_cgi_get_data?083112-3906.8,asas3

    # OPEN: ASAS3 data is actually a series of data, each series is a dataset,
    # possibly with slightly different mean magnitude, target positions, etc.
    # for now we just treat them as one single dataset
    headers = [
        "HJD",
        "MAG_3",
        "MAG_0",
        "MAG_1",
        "MAG_2",
        "MAG_4",
        "MER_3",
        "MER_0",
        "MER_1",
        "MER_2",
        "MER_4",
        "GRADE",
        "FRAME",
    ]
    headers = [h.lower() for h in headers]  # Lightkurve convention is to use lower case for column names
    tab = Table.read(url, format="ascii", comment="#", names=headers)

    for n in range(0, 5):
        # rename mag error columns to follow Lightkurve convention
        tab.rename_column(f"mer_{n}", f"mag_{n}_err")

        # c, ce = f"mag_{n}", f"mag_{n}_err"  # shorthand for the current mag / mag err column
        # convert the magic 99.999 (not measured) to nan
        tab[f"mag_{n}"][tab[f"mag_{n}"] == 99.999] = np.nan
        tab[f"mag_{n}_err"][tab[f"mag_{n}_err"] == 99.999] = np.nan

        # add units
        tab[f"mag_{n}"] = tab[f"mag_{n}"] * u.mag
        tab[f"mag_{n}_err"] = tab[f"mag_{n}_err"] * u.mag

    time = Time(tab["hjd"] + 2450000, format="jd", scale="utc")
    tab.remove_column("hjd")
    lc = lk.LightCurve(time=time, flux=tab[flux_column.lower()], flux_err=tab[f"{flux_column.lower()}_err"], data=tab)

    lc = lc[np.isin(lc["grade"], grade_mask)]

    lc.meta["FLUX_ORIGIN"] = flux_column.lower()
    lc.meta["FILEURL"] = url

    # deduce ASAS id, if the URL is the standard form from ASAS3 web site
    asas_id_match = re.search(r"https://www.astrouw.edu.pl/cgi-asas/asas_cgi_get_data[?]([0-9.+-]+)", url)
    if asas_id_match is not None:
        asas_id = asas_id_match[1]
    else:
        asas_id = None
    if asas_id is not None:
        lc.meta["LABEL"] = f"ASAS {asas_id}"
    return lc
