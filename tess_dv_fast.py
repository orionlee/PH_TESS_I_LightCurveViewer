import os
import re
import sqlite3
import numpy as np
import pandas as pd

import download_utils

DATA_BASE_DIR = "data/tess_dv_fast"

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


def _append_to_tcestats_csv(filepath, sectors_val):
    print(f"DEBUG appending to master tcestats csv from: {filepath}")
    dest = f"{DATA_BASE_DIR}/{TCESTATS_FILENAME}"

    df = pd.read_csv(filepath, comment="#")
    # replace sectors column with a uniform format, e.g., s0002-s0002
    # that describes the actual range for both single sector and multi sector csvs
    df["sectors"] = sectors_val

    # only write header when the file is first created
    write_header = not os.path.isfile(dest)
    df.to_csv(dest, index=False, header=write_header, mode="a")


def download_all_data():
    """Download all relevant data locally."""
    # TODO: scrape the listing pages to get the list of URLs, instead of hardcoded lists above
    # for tce stats csv files
    for url in sources_tcestats_single_sector + sources_tcestats_multi_sector:
        filename = _filename(url)
        filepath = f"{DATA_BASE_DIR}/{filename}"
        if os.path.isfile(filepath):
            continue

        print(f"DEBUG Download to {filepath} from: {url}")
        download_utils.download_file(url, filename=filename, download_dir=DATA_BASE_DIR)
        sectors_val = re.search(r"s\d{4}-s\d{4}", url)[0]  # e.g., s0002-s0002
        _append_to_tcestats_csv(filepath, sectors_val)

    # dv products download scripts (for urls to the products)
    for url in sources_dv_sh_single_sector + sources_dv_sh_multi_sector:
        filename = _filename(url)
        filepath = f"{DATA_BASE_DIR}/{filename}"
        if os.path.isfile(filepath):
            continue

        print(f"DEBUG Download to {filepath} from: {url}")
        download_utils.download_file(url, filename=filename, download_dir=DATA_BASE_DIR)

    # convert the csv into a sqlite db for speedier query by ticid
    print(f"DEBUG Convert master tcestats csv to sqlite db...")
    _export_tcestats_as_db()


def _export_tcestats_as_db():
    db_path = f"{DATA_BASE_DIR}/{TCESTATS_DBNAME}"
    df = read_tcestats_csv()
    with sqlite3.connect(db_path) as con:
        df.to_sql("tess_tcestats", con, if_exists="replace")

        # nice-to-have, but not critical
        sql_index = "create index tess_tcestats_ticid on tess_tcestats(ticid);"
        cursor = con.cursor()
        cursor.execute(sql_index)
        cursor.close()


def read_tcestats_csv():
    # for ~230k rows of TCE stats data, it took 4-10secs, taking up 200+Mb memory.
    csv_path = f"{DATA_BASE_DIR}/{TCESTATS_FILENAME}"
    return pd.read_csv(csv_path, comment="#", dtype={"tce_sectors": str})


def _query_tcestats_from_db(sql):
    db_path = f"{DATA_BASE_DIR}/{TCESTATS_DBNAME}"
    with sqlite3.connect(db_path) as con:
        return pd.read_sql(sql, con)


def _get_tcestats_of_tic_from_db(tic):
    return _query_tcestats_from_db(f"select * from tess_tcestats where ticid = {tic}")


# TODO: support tce_filter_func parameter, e.g.,
def get_tce_infos_of_tic(tic):
    df = _get_tcestats_of_tic_from_db(tic)
    df = _add_helpful_columns_to_tcestats(df)

    return df


R_EARTH_TO_R_JUPITER = 6378.1 / 71492


def _add_helpful_columns_to_tcestats(df):
    # 1. add convenience columns

    # id to construct exomast URL, e.g., TIC232646881S0073S0073TCE1
    # it can serves as a unique ID for the TCE too.
    df["exomast_id"] = (
        "TIC"
        + df["ticid"].astype(str)
        + df["sectors"].str.upper().str.replace("-", "")
        + "TCE"
        + df["tce_plnt_num"].astype(str)
    )
    # convert the bit pattern in tce_sectors column to number of sectors a TCE covers
    df["tce_num_sectors"] = df["tce_sectors"].str.count("1")
    df["tce_prad_jup"] = df["tce_prad"] * R_EARTH_TO_R_JUPITER
    df["tce_depth_pct"] = df["tce_depth"] / 10000
    df["tce_ditco_msky_sig"] = df["tce_ditco_msky"] / df["tce_ditco_msky_err"]  # TicOffset sig
    df["tce_dicco_msky_sig"] = df["tce_dicco_msky"] / df["tce_dicco_msky_err"]  # OotOffset sig
    # Note: model's stellar density, `starDensitySolarDensity` in dvr xml, is not available in csv

    # 2. add links to data products
    report_cols = dict(
        # dtype=object: for arbitrary string length
        dvs=np.full_like(df.index, "", dtype=object),
        dvm=np.full_like(df.index, "", dtype=object),
        dvr=np.full_like(df.index, "", dtype=object),
        dvr_xml=np.full_like(df.index, "", dtype=object),
        dvt_fits=np.full_like(df.index, "", dtype=object),
    )

    for i in range(len(df)):
        tic = df.iloc[i]["ticid"]
        sectors = df.iloc[i]["sectors"]
        planet_num = df.iloc[i]["tce_plnt_num"]
        report_dict = _get_dv_products(tic, sectors, planet_num)
        for type in report_cols.keys():
            report_cols[type][i] = report_dict[type]

    for type in report_cols.keys():
        df[type] = report_cols[type]

    return df


def _get_dv_products(tic, sectors, planet_num):
    # sectors: the value of sectors column in csv, e.g., s0002-s0072
    sector_start, sector_end = sectors.split("-")
    if sector_start == sector_end:
        # case single sector
        match = re.search("[1-9]\d+$", sectors)  # 2+ digits sector
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

    df = pd.read_csv(
        script_name,
        comment="#",
        sep=" ",
        names=["curl", "C", "-", "L", "o", "filename", "url"],
        usecols=["filename"],
    )

    df = df[df["filename"].str.contains(f"0{tic}-", regex=False)]
    a_res = {}

    if len(df) < 1:
        return a_res

    def to_url(filename_vals):
        if len(filename_vals) > 0:
            # use the last one, in case there are values, corresponding to multiple runs for a single TCE
            # (probably should not happen within a bulk download script, but just in case)
            filename = filename_vals.iat[-1]
            return f"https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/{filename}"
        else:
            return ""

    n = planet_num
    vals = df[df["filename"].str.contains(f"0{tic}-{n:02d}-\d+_dvs.pdf", regex=True)]["filename"]
    a_res["dvs"] = to_url(vals)
    vals = df[df["filename"].str.endswith("_dvm.pdf")]["filename"]
    a_res["dvm"] = to_url(vals)
    vals = df[df["filename"].str.endswith("_dvr.pdf")]["filename"]
    a_res["dvr"] = to_url(vals)
    vals = df[df["filename"].str.endswith("_dvr.xml")]["filename"]
    a_res["dvr_xml"] = to_url(vals)
    vals = df[df["filename"].str.endswith("_dvt.fits")]["filename"]
    a_res["dvt_fits"] = to_url(vals)
    return a_res


def display_tce_infos(df, return_as=None):
    from IPython.display import display, HTML

    df = df.sort_values(by=["tce_num_sectors", "exomast_id"], ascending=[False, True])

    # TODO: set the precision for period, duration_hr, etc.
    df["Codes"] = (
        "epoch=" + df["tce_time0bt"].astype(str) + ", "
        "duration_hr=" + df["tce_duration"].astype(str) + ", "
        "period=" + df["tce_period"].astype(str) + ", "
        "label=" + '"' + df["exomast_id"].str.replace(r"TIC\d+", "", regex=True).str.lower() + '",'
    )

    df = df.rename(
        columns={
            "tce_prad_jup": "Rp",
            "tce_time0bt": "Epoch",
            "tce_period": "Period",
            "tce_duration": "Duration",
            "tce_impact": "Impact b",
            "tce_depth_pct": "Depth",
            "tce_ditco_msky_sig": "TicOffset",
            "tce_dicco_msky_sig": "OotOffset",
        }
    )

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

    def format_exomast_id(id):
        short_name = re.sub(r"TIC\d+", "", id).lower()
        return f'<a target="_exomast" href="https://exo.mast.stsci.edu/exomast_planet.html?planet={id}">{short_name}</a>'

    def format_codes(codes):
        return f"""\
<input type="text" style="margin-left: 3ch; font-size: 90%; color: #666; width: 10ch;"
    onclick="this.select();" readonly value='{codes}'>"""

    format_specs = {
        "exomast_id": format_exomast_id,
        "dvs": lambda url: f'<a target="_blank" href="{url}">dvs</a>',
        "dvm": lambda url: f'<a target="_blank" href="{url}">dvm</a>',
        "dvr": lambda url: f'<a target="_blank" href="{url}">dvr</a>',
        "Rp": "{:.3f}",
        "Epoch": "{:.2f}",  # the csv has 2 digits precision
        "Duration": "{:.4f}",
        "Period": "{:.6f}",
        "Depth": "{:.4f}",
        "Impact b": "{:.2f}",
        "TicOffset": "{:.2f}",
        "OotOffset": "{:.2f}",
        "Codes": format_codes,
    }

    with pd.option_context("display.max_colwidth", None, "display.max_rows", 999, "display.max_columns", 99):
        styler = df[display_columns].style.format(format_specs).hide(axis="index")
        # hack to add units to the header
        html = styler.to_html()
        html = html.replace(">Rp</th>", "title='Rjup'>Rp</th>")
        html = html.replace(">Epoch</th>", "title='BTJD'>Epoch</th>")
        html = html.replace(">Duration</th>", "title='hr'>Duration</th>")
        html = html.replace(">Period</th>", "title='d'>Period</th>")
        html = html.replace(">Depth</th>", "title='percent'>Depth</th>")
        html = html.replace(">TicOffset</th>", "title='sigma'>TicOffset</th>")
        html = html.replace(">OotOffset</th>", "title='sigma'>OotOffset</th>")
        if return_as is None:
            return display(HTML(html))
        elif return_as == "html":
            return html
