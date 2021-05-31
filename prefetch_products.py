# A standalone script to batch prefetch products used in vetting

import asyncio
from datetime import datetime
import lightkurve as lk
import lightkurve_ext as lke
import tic_plot as tplt


# For vetting EB Candidates in pilot
# the list has subject id as well. It is ignored.
tic_sector_list_str = """\
1159914, 13, 36036102
10467747, 26, 48951863
22877676, 23, 44572924
36728572, 4, 30175638
38461578, 13, 36046252
43843023, 5, 30959993
49701947, 6, 31339817
76073981, 26, 48934888
80254085, 13, 36028870
85711393, 14, 36962309
86629981, 13, 36028759
140801898, 11, 34472719
146439804, 5, 30960461
149600376, 13, 36026492
150429807, 9, 41437360
159208074, 23, 44580309
160991022, 13, 36039086
165721304, 22, 43293255
167303776, 6, 31325090
167524561, 13, 36024151
167554898, 9, 41544056
167656297, 9, 41425746
177308817, 13, 36022114
178171080, 4, 30253517
198537349, 22, 43283518
199688472, 26, 48951484
207385593, 23, 44564164
229689154, 26, 48942603
233060434, 25, 48227121
235598526, 26, 48950926
237941053, 7, 31992482
248990782, 5, 30974885
254278612, 13, 36039811
255921197, 26, 48948447
260352096, 13, 36021377
260659986, 13, 36030526
261205462, 13, 36035541
261749149, 13, 36043805
278863150, 6, 31323969
291340620, 6, 31334910
293345927, 10, 34107629
299778628, 13, 36047821
320205468, 13, 36049500
340065168, 13, 36026668
341615749, 26, 48951817
342061072, 26, 48943266
348319604, 26, 48944619
348898049, 11, 34471723
349057976, 13, 36049092
349309442, 13, 36030675
349523518, 5, 30896790
363149298, 9, 41419270
367838785, 26, 48944782
368260545, 26, 48943897
372419506, 13, 36049659
388508344, 25, 48227275
397058826, 23, 44568773
422284103, 23, 44563169
423027012, 26, 48942511
464380201, 13, 36030879
468276605, 25, 48224154
"""


async def prefetch_products(tic_sector_list, max_num_sectors_to_download, download_dir=None):
    def info(msg):
        ts_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"{ts_str} \t{msg}", flush=True)

    def parse_tic_sector_list(text):
        lines = text.split("\n")
        lines_tokens = [a_line.split(",") for a_line in lines if a_line.strip() != "" and not a_line.startswith("#")]
        return [(int(row_tokens[0].strip()), int(row_tokens[1].strip())) for row_tokens in lines_tokens]

    if isinstance(tic_sector_list, str):
        tic_sector_list = parse_tic_sector_list(tic_sector_list)

    info("BEGIN prefetch_products()...")

    for tic, sector in tic_sector_list:
        info(f"Prefetching TIC {tic}, sector {sector}...")

        # emulate "Enter TIC" cell logic
        # TODO: for now it is semi hard-coded to be the one used by EB_Vetting.ipynb
        def limit_sr_to_download(sr):
            # use primary mission data only for this vetting
            sr = lk.SearchResult(sr.table[sr.table["sequence_number"] <= 26])
            if max_num_sectors_to_download is None:
                return sr
            return lke.of_sector_n_around(sr, sector, num_additions=max_num_sectors_to_download - 1)

        # OPEN: redirect stdout and IPython.display.display to null
        # Problems:
        # - IPython display: don't know how to do it.
        # - stdout: using "contextlib.redirect_stdout(open(os.devnull, "w")):" would work

        # Used to try to get LCs in background , but I got intermittent
        #   ValueError: I/O operation on closed file.
        # apparently from some of the messages about lcf_coll,
        # so it's run in frontground again
        lcf_coll, sr, sr_unfiltered = lke.download_lightcurves_of_tic_with_priority(
            tic,
            download_filter_func=limit_sr_to_download,
            download_dir=download_dir,
        )

        tpf_task = lke.create_download_tpf_task(
            f"TIC{tic}", sector=sector, exptime="short", author="SPOC", mission="TESS", download_dir=download_dir
        )

        tce_res = tplt.get_tce_infos_of_tic(tic, download_dir=download_dir)

        tpf_coll, sr_tpf = await tpf_task

        result_msg = f"Downloaded - num. LCs: {len(lcf_coll)}, num. TPFs: {len(tpf_coll)}, num. TCEs: {len(tce_res)}"

        info(f"result: {result_msg}\n")

    info("END  prefetch_products()")


#
# The main logic
#

# BEGIN copied from "Enter TIC" cell in EB_Vetting.ipynb
max_num_sectors_to_download = 3

if hasattr(lk.search, "sr_cache"):
    lk.search.sr_cache.cache_dir = "data"
    lk.search.sr_cache.expire_seconds = 86400 * 7
# END copied from "Enter TIC" cell

asyncio.run(
    prefetch_products(tic_sector_list_str, max_num_sectors_to_download=max_num_sectors_to_download, download_dir="data")
)
