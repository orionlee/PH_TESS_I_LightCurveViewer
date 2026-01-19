#
# Generic file download utilities that support file-based cache
#

import os
import re
import shutil
import time
import warnings
from types import SimpleNamespace

import requests


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


CachePolicy = SimpleNamespace(
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


def _do_download_file(url, filename=None, download_dir=None):
    if download_dir is None:
        download_dir = ""
    os.makedirs(download_dir, exist_ok=True)

    local_filename = _create_local_filename(url, filename, download_dir)

    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        # write to a temporary file. If successful, make it the real local file
        # it is to prevent interrupted download leaving a partial file
        local_filename_temp = f"{local_filename}.download"
        with open(local_filename_temp, "wb") as out_file:
            shutil.copyfileobj(response.raw, out_file)
        os.replace(local_filename_temp, local_filename)
    return local_filename


def download_file(
    url,
    filename=None,
    download_dir=None,
    cache_policy_func=None,
    return_is_cache_used=False,
):
    if download_dir is None:
        download_dir = ""

    is_cache_used = False
    local_filename = _create_local_filename(url, filename, download_dir)
    if os.path.isfile(local_filename):
        if cache_policy_func is None or cache_policy_func(url, local_filename):
            is_cache_used = True

    if not is_cache_used:
        _do_download_file(url, filename, download_dir)

    if return_is_cache_used:
        return local_filename, is_cache_used
    else:
        return local_filename
