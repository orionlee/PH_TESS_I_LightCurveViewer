#
# Helpers to download TESS-specific non-lightcurve data
#

import os
import re
import shutil
import time
from types import SimpleNamespace
import warnings

import pandas as pd
import requests


R_earth = 6371000  # radius of the Earth [m]
R_jup = 69911000  # radius of Jupiter [m]


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
            # TODO: add a parameter, force_download_fn=None, to let caller specify criteria to re download
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


class TOIAccessor:

    Headers = SimpleNamespace(
        TIC="TIC ID",
        TOI="TOI",
        EPOCH="Epoch (BJD)",
        PERIOD="Period (days)",
        DURATION_HR="Duration (hours)",
        DEPTH_PPM="Depth (ppm)",
        DEPTH_PCT="Depth (percent)",  # derived
        PLANET_RADIUS_E="Planet Radius (R_Earth)",
        PLANET_RADIUS_J="Planet Radius (R_Jupiter)",  # derived
        COMMENTS="Comments",
    )

    # TODO: in-memory cache (with @cached) needs to be redone to properly support use_localfile_func
    @classmethod
    def get_all_tois(cls, download_dir="", use_localfile_func=None):
        url = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
        filename = "tess_tois.csv"
        res = _get_csv(url, filename, download_dir, use_localfile_func=use_localfile_func, dtype={cls.Headers.TOI: str})
        # add dervied columns
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
