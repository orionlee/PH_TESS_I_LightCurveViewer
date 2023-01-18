#
# Misc. lightkurve readers for different formats
#

import re

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
    tab = Table.read(url, format="ascii")
    # print("DBG2: colnames=", tab.colnames)

    # make columns follow Lightkurve convention
    if "flux(mJy)" in tab.colnames:
        # column flux(mJy) shows up in csv from path "/variables/" or "/sky-patrol/coordinate/"
        # for csv from path "/photometry/" , the column is simply flux
        tab.rename_column("flux(mJy)", "flux")

    for c in tab.colnames:
        tab.rename_column(c, re.sub(" ", "_", c.lower()))

    tab.sort(keys="hjd")

    tab["mag"].unit = u.mag
    tab["mag_err"].unit = u.mag
    # OPEN: unit of flux is actually mJy. I don't know what it means though
    tab["flux"].unit = u.dimensionless_unscaled
    tab["flux_err"].unit = u.dimensionless_unscaled

    # SkyPatrol csv contains rows with mag / flux as 99.99, presumably the underlying data is bad
    tab = tab[tab["mag"] < 99.99]

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
