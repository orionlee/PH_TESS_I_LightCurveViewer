import os
from pathlib import Path
import re
import sqlite3
import shutil

import numpy as np
import pandas as pd

try:
    from astropy.table import Column

    _HAS_ASTROPY = False
except Exception:
    _HAS_ASTROPY = False

# for tic parameter in get_tce_infos_of_tic(), case a list of TICs
#
# Note: use list, tuple explicitly, instead of collection.abc.Sequence,
# because types such as str also implements Sequence
if _HAS_ASTROPY:
    _ARRAY_LIKE_TYPES = (list, tuple, set, np.ndarray, pd.Series, Column)
else:
    _ARRAY_LIKE_TYPES = (list, tuple, set, np.ndarray, pd.Series)


DATA_BASE_DIR = f"{Path(__file__).parent}/data/tess_dv_fast"

TCESTATS_FILENAME = f"tess_tcestats.csv"
TCESTATS_DBNAME = f"tess_tcestats.db"

# csv source: https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html
"""javascript
urls = Array.from(document.querySelectorAll('table#TABLE_4 a')).map(a => a.href)
urls.reverse()
urls.map(v => `    "${v}",`).join("\n")
"""
sources_tcestats_single_sector = [
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018206190142-s0001-s0001_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018235142541-s0002-s0002_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018263124740-s0003-s0003_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018292093539-s0004-s0004_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018319112538-s0005-s0005_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018349182737-s0006-s0006_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019008025936-s0007-s0007_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019033200935-s0008-s0008_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019059170935-s0009-s0009_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019085221934-s0010-s0010_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019113062933-s0011-s0011_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019141104532-s0012-s0012_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019170095531-s0013-s0013_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019199201929-s0014-s0014_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019227203528-s0015-s0015_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019255032927-s0016-s0016_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019281041526-s0017-s0017_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019307033525-s0018-s0018_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019332134924-s0019-s0019_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019358235523-s0020-s0020_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020021221522-s0021-s0021_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020050191121-s0022-s0022_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020079142120-s0023-s0023_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020107065519-s0024-s0024_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020135030118-s0025-s0025_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020161181517-s0026-s0026_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020187183116-s0027-s0027_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020213081515-s0028-s0028_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020239173514-s0029-s0029_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020267090513-s0030-s0030_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020296001112-s0031-s0031_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020325171311-s0032-s0032_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2020353052510-s0033-s0033_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021014055109-s0034-s0034_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021040113508-s0035-s0035_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021066093107-s0036-s0036_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021092173506-s0037-s0037_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021119082105-s0038-s0038_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021147062104-s0039-s0039_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021176033103-s0040-s0040_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021205113501-s0041-s0041_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021233042500-s0042-s0042_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021259155059-s0043-s0043_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021285162058-s0044-s0044_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021311000057-s0045-s0045_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021337012456-s0046-s0046_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021365070455-s0047-s0047_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022028101454-s0048-s0048_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022057231053-s0049-s0049_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022085182052-s0050-s0050_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022113103451-s0051-s0051_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022139030450-s0052-s0052_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022164114449-s0053-s0053_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022190092448-s0054-s0054_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022217141447-s0055-s0055_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022245180045-s0056-s0056_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022273202044-s0057-s0057_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022302194443-s0058-s0058_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022330181042-s0059-s0059_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2022357093041-s0060-s0060_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023018070040-s0061-s0061_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023043222439-s0062-s0062_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023069204038-s0063-s0063_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023096143437-s0064-s0064_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023124053436-s0065-s0065_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023153040035-s0066-s0066_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023182031434-s0067-s0067_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023210023032-s0068-s0068_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023237202031-s0069-s0069_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023263202030-s0070-s0070_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023289124029-s0071-s0071_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023315161428-s0072-s0072_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2023341070027-s0073-s0073_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2024003083426-s0074-s0074_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2024030064025-s0075-s0075_dvr-tcestats.csv",
]

# Similar javascript codes to scrape the url, except the initial CSS selector is: "table#TABLE_5 a"
sources_tcestats_multi_sector = [
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018206190142-s0001-s0002_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018206190142-s0001-s0003_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018206190142-s0001-s0006_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018206190142-s0001-s0009_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018206190142-s0001-s0013_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019199201929-s0014-s0016_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019199201929-s0014-s0019_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019199201929-s0014-s0023_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019199201929-s0014-s0026_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018206190142-s0001-s0036_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018206190142-s0001-s0039_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019199201929-s0014-s0041_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021233042500-s0042-s0043_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018206190142-s0001-s0046_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2021233042500-s0042-s0046_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019199201929-s0014-s0050_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019199201929-s0014-s0055_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2019199201929-s0014-s0060_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018206190142-s0001-s0065_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018206190142-s0001-s0069_dvr-tcestats.csv",
    "https://archive.stsci.edu/missions/tess/catalogs/tce/tess2018235142541-s0002-s0072_dvr-tcestats.csv",
]

# dv products bulk download source: https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html
"""javascript
urls = Array.from(document.querySelectorAll('table a')).map(a => a.href)
urls = urls.filter(u => u.match(/_sector_.+_dv.sh$/))
urls.reverse()
urls.map(v => `    "${v}",`).join("\n")
"""
sources_dv_sh_single_sector = [
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_1_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_2_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_3_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_4_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_5_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_6_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_7_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_8_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_9_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_10_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_11_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_12_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_13_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_14_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_15_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_16_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_17_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_18_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_19_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_20_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_21_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_22_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_23_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_24_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_25_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_26_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_27_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_28_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_29_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_30_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_31_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_32_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_33_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_34_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_35_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_36_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_37_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_38_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_39_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_40_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_41_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_42_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_43_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_44_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_45_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_46_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_47_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_48_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_49_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_50_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_51_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_52_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_53_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_54_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_55_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_56_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_57_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_58_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_59_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_60_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_61_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_62_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_63_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_64_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_65_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_66_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_67_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_68_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_69_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_70_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_71_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_72_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_73_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_74_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_75_dv.sh",
]

# urls = urls.filter(u => u.match(/_multisector_.+_dv.sh$/))
sources_dv_sh_multi_sector = [
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0001-s0002_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0001-s0003_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0001-s0006_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0001-s0009_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0001-s0013_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0014-s0016_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0014-s0019_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0014-s0023_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0014-s0026_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0001-s0036_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0001-s0039_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0014-s0041_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0042-s0043_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0001-s0046_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0042-s0046_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0014-s0050_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0014-s0055_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0014-s0060_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0001-s0065_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0001-s0069_dv.sh",
    "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_multisector_s0002-s0072_dv.sh",
]


def _filename(url):
    match = re.search("[^/]+$", url)
    if match is not None:
        return match[0]
    else:
        raise ValueError(f"Failed to extract filename from url: {url}")


def _add_exomast_id(df):
    # id to construct exomast URL, e.g., TIC232646881S0073S0073TCE1
    # it can serves as a unique ID for the TCE across sectors too.
    df.insert(
        loc=0,
        column="exomast_id",
        value=(
            "TIC"
            + df["ticid"].astype(str)
            + df["sectors"].str.upper().str.replace("-", "")
            + "TCE"
            + df["tce_plnt_num"].astype(str)
        ),
    )


def _get_dv_products_of_sectors(sectors):
    # sectors: the value of sectors column in csv, e.g., s0002-s0072
    sector_start, sector_end = sectors.split("-")
    if sector_start == sector_end:
        # case single sector
        match = re.search(r"[1-9]\d+$", sectors)  # 2+ digits sector
        if match is not None:
            sector_str_in_filename = match[0]
        else:
            match = re.search("[1-9]$", sectors)  # single sector
            if match is not None:
                sector_str_in_filename = match[0]
            else:
                raise ValueError(f"Cannot prase sector from string {sectors}")
        script_name = f"{DATA_BASE_DIR}/tesscurl_sector_{sector_str_in_filename}_dv.sh"
    else:
        # case multi-sector
        script_name = f"{DATA_BASE_DIR}/tesscurl_multisector_{sectors}_dv.sh"

    filename = pd.read_csv(
        script_name,
        comment="#",
        sep=" ",
        names=["curl", "C", "-", "L", "o", "filename", "url"],
        usecols=["filename"],
    )["filename"]

    dvs_filename = filename[filename.str.contains("_dvs.pdf")]
    dvs = dvs_filename.str.extract(rf"tess\d+-{sectors}-0+?(?P<ticid>[1-9]\d+)-(?P<tce_plnt_num>\d\d)-")
    dvs.ticid = dvs.ticid.astype("int64")
    dvs.tce_plnt_num = dvs.tce_plnt_num.astype(int)
    dvs["dvs"] = dvs_filename
    # in case multiple runs for the same TIC-sector, use the last one only
    dvs = dvs.groupby(by=["ticid", "tce_plnt_num"], as_index=False).last()

    def get_products_of_type(colname, filename_suffix):
        product_filename = filename[filename.str.contains(filename_suffix)]
        df = product_filename.str.extract(rf"tess\d+-{sectors}-0+?(?P<ticid>[1-9]\d+)-")
        df.ticid = df.ticid.astype("int64")
        df[colname] = product_filename
        # in case multiple runs for the same TIC-sector, use the last one only
        df = df.groupby(by="ticid", as_index=False).last()
        return df

    dvm = get_products_of_type("dvm", "_dvm.pdf")
    dvr = get_products_of_type("dvr", "_dvr.pdf")
    dvr_xml = get_products_of_type("dvr_xml", "_dvr.xml")
    dvt = get_products_of_type("dvt", "_dvt.fits")

    res = dvs
    for df in [dvm, dvr, dvr_xml, dvt]:
        res = pd.merge(res, df, on="ticid", how="left", validate="many_to_one")

    res["sectors"] = sectors
    _add_exomast_id(res)
    res.drop(["ticid", "tce_plnt_num", "sectors"], axis="columns", inplace=True)
    return res


def _append_to_tcestats_csv(filepath, sectors_val, dest):
    print(f"DEBUG appending to master tcestats csv from: {filepath}")

    # base tcestats csv of a sector
    df = pd.read_csv(filepath, comment="#")

    # replace sectors column with a uniform format, e.g., s0002-s0002
    # that describes the actual range for both single sector and multi sector csvs
    df["sectors"] = sectors_val

    # uniquely identify a TCE across all sectors
    _add_exomast_id(df)

    # include filenames of dvs, dvm, etc.
    df_filenames = _get_dv_products_of_sectors(sectors_val)
    df = pd.merge(df, df_filenames, on="exomast_id", validate="one_to_one")

    # only write header when the file is first created
    write_header = not os.path.isfile(dest)
    df.to_csv(dest, index=False, header=write_header, mode="a")


def download_all_data(minimal_db=False):
    """Download all relevant data locally."""
    import download_utils

    # TODO: scrape the listing pages to get the list of URLs, instead of hardcoded lists above

    # dv products download scripts (for urls to the products)
    # - they need to be first downloaded: as creating master csv below relies on the scripts
    for url in sources_dv_sh_single_sector + sources_dv_sh_multi_sector:
        filename = _filename(url)
        filepath, is_cache_used = download_utils.download_file(
            url, filename=filename, download_dir=DATA_BASE_DIR, return_is_cache_used=True
        )
        if not is_cache_used:
            print(f"DEBUG Downloaded to {filepath} from: {url}")

    # for tce stats csv files, download and merge them to a single csv
    # - first write to a temporary master csv. Once done, overwrite the existing master (if any)
    dest_csv_tmp = f"{DATA_BASE_DIR}/{TCESTATS_FILENAME}.tmp"
    dest_csv = f"{DATA_BASE_DIR}/{TCESTATS_FILENAME}"
    Path(dest_csv_tmp).unlink(missing_ok=True)
    for url in sources_tcestats_single_sector + sources_tcestats_multi_sector:
        filename = _filename(url)

        filepath, is_cache_used = download_utils.download_file(
            url, filename=filename, download_dir=DATA_BASE_DIR, return_is_cache_used=True
        )
        if not is_cache_used:
            print(f"DEBUG Downloaded to {filepath} from: {url}")

        sectors_val = re.search(r"s\d{4}-s\d{4}", url)[0]  # e.g., s0002-s0002
        _append_to_tcestats_csv(filepath, sectors_val, dest_csv_tmp)
    shutil.move(dest_csv_tmp, dest_csv)

    # convert the master csv into a sqlite db for speedier query by ticid
    print(f"DEBUG Convert master tcestats csv to sqlite db, minimal_db={minimal_db}...")
    _export_tcestats_as_db(minimal_db)


# minimal list of columns in db in order to support display_tce_infos
# the resulting db is about 30% of the full db
_MIN_DB_COLS = [
    "exomast_id",
    "ticid",
    "tce_prad",  # for deriving: "Rp in Rj"
    "tce_time0bt",  #  "Epoch"
    "tce_period",  # "Period"
    "tce_duration",  # "Duration"
    "tce_impact",  # "Impact b"
    "tce_depth",  # for deriving "Depth"
    # for deriving "TicOffset"
    "tce_ditco_msky",
    "tce_ditco_msky_err",
    # for deriving "OotOffset",
    "tce_dicco_msky",
    "tce_dicco_msky_err",
    "dvs",
    "dvm",
    "dvr",
    "tce_sectors",  # for deriving "tce_num_sectors"
    # the following is not needed strictly speaking, but
    # they are components of a TCE's identifier, so it might be handy.
    "tce_plnt_num",
    "sectors",
]


def _export_tcestats_as_db(minimal_db=False):
    db_path_tmp = f"{DATA_BASE_DIR}/{TCESTATS_DBNAME}.tmp"
    db_path = f"{DATA_BASE_DIR}/{TCESTATS_DBNAME}"

    usecols = None if not minimal_db else _MIN_DB_COLS
    df = read_tcestats_csv(
        usecols=usecols,
    )

    Path(db_path_tmp).unlink(missing_ok=True)
    con = sqlite3.connect(db_path_tmp)
    try:  # use try / finally instead of with ... because sqlite3 context manager does not close the connection
        df.to_sql("tess_tcestats", con, if_exists="replace", index=False)

        # nice-to-have, but not critical
        sql_index = "create index tess_tcestats_ticid on tess_tcestats(ticid);"
        cursor = con.cursor()
        cursor.execute(sql_index)
        cursor.close()
        con.commit()
    finally:
        con.close()

    shutil.move(db_path_tmp, db_path)


def read_tcestats_csv(**kwargs):
    # for ~230k rows of TCE stats data, it took 4-10secs, taking up 200+Mb memory.
    csv_path = f"{DATA_BASE_DIR}/{TCESTATS_FILENAME}"
    return pd.read_csv(csv_path, comment="#", dtype={"tce_sectors": str}, **kwargs)


def _query_tcestats_from_db(sql, **kwargs):
    db_uri = f"file:{DATA_BASE_DIR}/{TCESTATS_DBNAME}?mode=ro"  # read-only
    with sqlite3.connect(db_uri, uri=True) as con:
        return pd.read_sql(sql, con, **kwargs)


def _get_tcestats_of_tic_from_db(tic):
    # OPEN: support optional columns parameter?
    # - double quote the column names in the constructed SQL
    #   would be sufficient wo avoid SQL injection
    # - need to make _add_helpful_columns_to_tcestats()
    #   handle missing columns though
    if isinstance(tic, (int, float, str)):
        return _query_tcestats_from_db(
            "select * from tess_tcestats where ticid = ?",
            params=[tic],
        )
    elif isinstance(tic, _ARRAY_LIKE_TYPES):
        #
        # convert to sqlite driver acceptable type
        # 1. list, and
        # 2. of type `int`, e.g., using `np.int32` would result in no match
        #    (probably requires registering some conversion in the driver)
        tic = [int(v) for v in tic]
        in_params_place_holder = ",".join(["?" for i in range(len(tic))])
        return _query_tcestats_from_db(
            f"select * from tess_tcestats where ticid in ({in_params_place_holder})",
            params=tic,
        )
    else:
        raise TypeError(f"tic must be a scalar or array-like. Actual type: {type(tic).__name__}")


def get_tce_infos_of_tic(tic, tce_filter_func=None):
    df = _get_tcestats_of_tic_from_db(tic)
    _add_helpful_columns_to_tcestats(df)
    # sort the result to the standard form
    # so that it is predictable for tce_filter_func
    df = df.sort_values(by=["ticid", "tce_num_sectors", "exomast_id"], ascending=[True, False, True])
    if tce_filter_func is not None and len(df) > 0:
        df = tce_filter_func(df)

    return df


R_EARTH_TO_R_JUPITER = 6378.1 / 71492


def _add_helpful_columns_to_tcestats(df):
    # convert the bit pattern in tce_sectors column to number of sectors a TCE covers
    df["tce_num_sectors"] = df["tce_sectors"].str.count("1")
    df["tce_prad_jup"] = df["tce_prad"] * R_EARTH_TO_R_JUPITER
    df["tce_depth_pct"] = df["tce_depth"] / 10000
    df["tce_ditco_msky_sig"] = df["tce_ditco_msky"] / df["tce_ditco_msky_err"]  # TicOffset sig
    df["tce_dicco_msky_sig"] = df["tce_dicco_msky"] / df["tce_dicco_msky_err"]  # OotOffset sig
    # Note: model's stellar density, `starDensitySolarDensity` in dvr xml, is not available in csv


def to_product_url(filename):
    """Convert the product filenames in columns such as dvs, dvr, etc., to URL to MAST server"""
    return f"https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/{filename}"


def display_tce_infos(df, return_as=None, no_tce_html=None):
    df["Codes"] = (
        "epoch=" + df["tce_time0bt"].astype(str) + ", "
        "duration_hr=" + df["tce_duration"].astype(str) + ", "
        "period=" + df["tce_period"].astype(str) + ", "
        "label=" + '"' + df["exomast_id"].str.replace(r"TIC\d+", "", regex=True).str.lower() + '", '
        "transit_depth_percent=" + df["tce_depth_pct"].map("{:.4f}".format) + ","
    )

    df = df.rename(
        columns={
            "tce_prad_jup": "Rp",
            "tce_time0bt": "Epoch",
            "tce_period": "Period",
            "tce_duration": "Duration",
            "tce_impact": "Impact b",
            "tce_depth_pct": "Depth",
        }
    )

    df["TicOffset"] = df["tce_ditco_msky"].astype(str) + "|" + df["tce_ditco_msky_sig"].astype(str)
    df["OotOffset"] = df["tce_dicco_msky"].astype(str) + "|" + df["tce_dicco_msky_sig"].astype(str)

    display_columns = [
        "exomast_id",
        # "ticid",
        # "tce_plnt_num",
        # "sectors",
        # "tce_num_sectors",
        "dvs",
        "dvm",
        "dvr",
        "Rp",
        "Epoch",
        "Duration",
        "Period",
        "Depth",
        "Impact b",
        "TicOffset",
        "OotOffset",
        "Codes",
    ]

    if len(df["ticid"].unique()) > 1:
        # case multiple TICs in the result
        # prepend ticid to the columns to be displayed to differentiate between them
        display_columns = ["ticid"] + display_columns

    def format_exomast_id(id):
        short_name = re.sub(r"TIC\d+", "", id).lower()
        return f'<a target="_exomast" href="https://exo.mast.stsci.edu/exomast_planet.html?planet={id}">{short_name}</a>'

    def format_offset_n_sigma(val_sigma_str):
        # format TicOffset / OotOffset, to value (sigma)
        try:
            val, sigma = val_sigma_str.split("|")
            val = float(val)
            sigma = float(sigma)
            if sigma >= 3:
                sigma_style = ' style="color: red; font-weight: bold;"'
            else:
                sigma_style = ""
            return f"{val:.0f} <span{sigma_style}>({sigma:.1f})<span>"
        except Exception:
            # in case something unexpected, fallback to raw str to avoid exception in display
            return val_sigma_str

    def format_codes(codes):
        return f"""\
<input type="text" style="margin-left: 3ch; font-size: 90%; color: #666; width: 10ch;"
    onclick="this.select();" readonly value='{codes}'>"""

    format_specs = {
        "exomast_id": format_exomast_id,
        "dvs": lambda f: f'<a target="_blank" href="{to_product_url(f)}">dvs</a>',
        "dvm": lambda f: f'<a target="_blank" href="{to_product_url(f)}">dvm</a>',
        "dvr": lambda f: f'<a target="_blank" href="{to_product_url(f)}">dvr</a>',
        "Rp": "{:.3f}",
        "Epoch": "{:.2f}",  # the csv has 2 digits precision
        "Duration": "{:.4f}",
        "Period": "{:.6f}",
        "Depth": "{:.4f}",
        "Impact b": "{:.2f}",
        "TicOffset": format_offset_n_sigma,
        "OotOffset": format_offset_n_sigma,
        "Codes": format_codes,
    }

    with pd.option_context("display.max_colwidth", None, "display.max_rows", 999, "display.max_columns", 99):
        styler = df[display_columns].style.format(format_specs).hide(axis="index")
        # hack to add units to the header
        html = styler.to_html()
        html = html.replace(">Rp</th>", ">R<sub>p</sub><br>R<sub>j</sub></th>", 1)
        html = html.replace(">Epoch</th>", ">Epoch<br>BTJD</th>", 1)
        html = html.replace(">Duration</th>", ">Duration<br>hr</th>", 1)
        html = html.replace(">Period</th>", ">Period<br>day</th>", 1)
        html = html.replace(">Depth</th>", ">Depth<br>%</th>", 1)
        html = html.replace(">TicOffset</th>", '>TicOffset<br>" (σ)</th>', 1)
        html = html.replace(">OotOffset</th>", '>OotOffset<br>" (σ)</th>', 1)
        if len(df) < 1 and no_tce_html is not None:
            html = no_tce_html
        if return_as is None:
            from IPython.display import display, HTML

            return display(HTML(html))
        elif return_as == "html":
            return html


# update the DB / master csv from command line
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Update master TESS TCE data")
    parser.add_argument("--update", dest="update", action="store_true", help="Update master csv and sqlite db.")
    parser.add_argument(
        "--minimal_db",
        dest="minimal_db",
        action="store_true",
        default=False,
        help="make the sqlite db minimal for typical use cases / webapp usage",
    )

    args = parser.parse_args()
    if not args.update:
        print("--update must be specified")
        parser.print_help()
        parser.exit()

    print(f"Downloading data to create master csv, minimal_db={args.minimal_db}")
    download_all_data(minimal_db=args.minimal_db)
