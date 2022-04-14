#
# Lightkurve Extension to interface with Pyaneti Modeling Library
# - https://github.com/oscaribv/pyaneti
#

from collections.abc import Iterable
import os
from os import path
from pathlib import Path
import re
import warnings

from astropy.table import Table
from astropy.time import Time
import numpy as np
import lightkurve as lk


class PyanetiEnv:
    """Define the directories used for 1 modeling session"""

    def __init__(self, home_dir, alias, sector):
        self.home_dir = Path(home_dir)
        self.alias = alias
        self.sector = sector

    @property
    def base_in_dir(self):
        """The base directory for all modeling input of the Pyaneti installation"""
        return Path(self.home_dir, "inpy")

    @property
    def base_out_dir(self):
        """The base directory for all modeling output of the Pyaneti installation"""
        return Path(self.home_dir, "outpy")

    @property
    def target_in_dir(self):
        return Path(self.base_in_dir, self.alias)

    @property
    def target_out_dir(self):
        return Path(self.base_out_dir, f"{self.alias}_out")

    @property
    def lc_dat_filename(self):
        sector_str = to_sector_str(self.sector)
        return f"{self.alias}_lc_s{sector_str}.dat"

    @property
    def lc_dat_filepath(self):
        return Path(self.target_in_dir, self.lc_dat_filename)


#
# Prepare and export `LightCurve` to Pyaneti input data
#


def to_sector_str(sector):
    def format(a_sector):
        return f"{a_sector:02d}"

    if isinstance(sector, Iterable):
        return "_".join([format(i) for i in sector])
    else:
        format(sector)


def _create_dir_if_needed(path):
    basedir = os.path.dirname(path)
    if not os.path.isdir(basedir):
        os.makedirs(basedir)


def truncate_lc_to_around_transits(lc, transit_specs):
    def calc_duration_to_use(spec):
        """Calc a duration for the purpose of masking,
        by adding transit duration with an additional `surround_time`,
        with a default of (`max(transit duration * 2, 1 day)`
        """
        duration = spec["duration_hr"] / 24
        surround_time = spec.get("surround_time", None)
        if surround_time is None:
            surround_time = max(duration * 3, 1)  # duration + 2 * duration
        return surround_time

    def calc_period_to_use(spec):
        period = spec.get("period", None)
        if period is not None and period > 0:
            return period
        # else period is not really specified,
        # for the use case that a single dip is observed with noted
        # use an arbitrary large period as a filler
        return 99999999

    lc = lc.remove_nans().normalize()

    period = [calc_period_to_use(t) for t in transit_specs]
    duration = [calc_duration_to_use(t) for t in transit_specs]
    transit_time = [t["epoch"] for t in transit_specs]
    # a mask to include the transits and their surrounding (out of transit observations)
    mask = lc.create_transit_mask(period=period, duration=duration, transit_time=transit_time)

    return lc[mask]


def to_pyaneti_dat(lc, transit_specs, pyaneti_env, return_processed_lc=False):
    "Output lc data to a file readable by Pyaneti, with lc pre-processed to be suitable for Pyaneti modeling"
    if transit_specs is not None:
        lc = truncate_lc_to_around_transits(lc, transit_specs)

    out_path = pyaneti_env.lc_dat_filepath

    # finally write to output
    #  lc["time", "flux", "flux_err"]  # somehow "time" is causing problem
    lc1 = type(lc)(time=lc.time.copy(), flux=lc.flux, flux_err=lc.flux_err)
    _create_dir_if_needed(out_path)
    lc1.write(out_path, format="ascii.commented_header", overwrite=True)
    if return_processed_lc:
        return out_path, lc
    else:
        return out_path


def catalog_info_TIC(tic_id):
    """Takes TIC_ID, returns stellar information from online catalog using Vizier"""
    if type(tic_id) is not int:
        raise TypeError('tic_id must be of type "int"')
    try:
        from astroquery.mast import Catalogs
    except:
        raise ImportError("Package astroquery required but failed to import")

    result_tab = Catalogs.query_criteria(catalog="Tic", ID=tic_id)
    return {c: result_tab[0][c] for c in result_tab[0].colnames}


def get_limb_darkening_params(tic_meta):
    # derived from:
    # https://github.com/hippke/tls/blob/v1.0.31/transitleastsquares/catalog.py
    logg, Teff, = (
        tic_meta["logg"],
        tic_meta["Teff"],
    )
    if logg is None:
        logg = 4
        warnings.warn("No logg in metadata. Proceeding with logg=4")

    if Teff is None:
        Teff = 6000
        warnings.warn("No Teff in metadata Proceeding with Teff=6000")

    ld = np.genfromtxt(
        path.join("catalogs", "ld_claret_tess.csv"),
        skip_header=1,
        delimiter=",",
        dtype="f8, int32, f8, f8",
        names=["logg", "Teff", "a", "b"],
    )

    """
        - Take Teff from star catalog and find nearest entry in LD catalog
        - Same for logg, but only for the Teff values returned before
        - Return  best-match LD
    """
    nearest_Teff = ld["Teff"][(np.abs(ld["Teff"] - Teff)).argmin()]
    idx_all_Teffs = np.where(ld["Teff"] == nearest_Teff)
    relevant_lds = np.copy(ld[idx_all_Teffs])
    idx_nearest = np.abs(relevant_lds["logg"] - logg).argmin()
    q1 = relevant_lds["a"][idx_nearest]
    q2 = relevant_lds["b"][idx_nearest]

    # Provide a rough guess on error for q1/q2
    # it's a rough heuristics that tried to be conservation (i.e, erred to be larger than actual)
    e_q1 = np.ceil(q1 * 0.33 * 100) / 100
    e_q2 = np.ceil(q2 * 0.33 * 100) / 100

    return dict(q1=q1, q2=q2, e_q1=e_q1, e_q2=e_q2)


def _round_n_decimals(num, num_decimals):
    factor = np.power(10, num_decimals)
    return np.round(num * factor) / factor


def estimate_planet_radius_in_r_star(r_star, depth):
    """Return a back of envelope estimate of a planet object's radius,
    based on the simple model of a planet with circular orbit,
    transiting across the center of the host star (impact parameter `b` = 0)
    """
    if r_star is None or r_star < 0 or depth is None or depth <= 0:
        return None  # cannot estimate

    R_JUPITER_IN_R_SUN = 71492 / 695700

    r_planet = np.sqrt(r_star * r_star * depth)

    r_planet_in_r_star = r_planet / r_star

    # Provide some rough min / max estimate
    r_planet_in_r_star_min = 0
    # a rough guess for max: 2 times of the above estimate, capped to the size of 2.5 R_jupiter
    max_r_planet_in_r_star = 2.5 * R_JUPITER_IN_R_SUN / r_star
    r_planet_in_r_star_max = min(r_planet_in_r_star * 2, max_r_planet_in_r_star)

    # somehow the number in full precision causes pyaneti to behave strangely
    # (getting invalid numeric number in calculation, results in `nan` in `T_full`, etc.
    # capped at 4 decimal, taking a cue from pyaneti output
    return dict(
        r_planet_in_r_star=_round_n_decimals(r_planet_in_r_star, 4),
        r_planet_in_r_star_min=_round_n_decimals(r_planet_in_r_star_min, 4),
        r_planet_in_r_star_max=_round_n_decimals(r_planet_in_r_star_max, 4),
    )


def estimate_orbital_distance_in_r_star(tic_meta):
    # TODO: possibly use Kepler third law for better constraints
    return dict(a_min=2.0, a_max=30.0)


def create_input_fit(
    template_name,
    tic,
    alias,
    pti_env,
    lc,
    transit_specs,
    meta,
    q1_q2,
    r_planet_dict,
    a_planet_dict,
    write_to_file=True,
    return_content=False,
):
    """Output parts of Pyaneti `input_fit.py` based on the specification included"""
    template = Path("pyaneti_templates", f"input_{template_name}.py").read_text()

    pyaneti_target_in_dir = pti_env.target_in_dir
    lc_pyaneti_dat_filename = pti_env.lc_dat_filename

    mapping = meta.copy()
    mapping.update(q1_q2)
    mapping.update(r_planet_dict)
    mapping.update(a_planet_dict)
    mapping["tic"] = tic
    mapping["alias"] = alias
    mapping["fname_tr"] = lc_pyaneti_dat_filename

    rho = mapping.get("rho")
    e_rho = mapping.get("e_rho")
    if rho is not None and e_rho is not None:
        mapping["rho_min"] = rho - e_rho
        mapping["rho_max"] = rho + e_rho

    # map transit_specs to epoch_min/max, period_min/max
    # TODO: handle multiple transit_specs
    epoch_error = transit_specs[0].get("epoch_error", None)
    if epoch_error is None:
        epoch_error = transit_specs[0]["duration_hr"] * 0.05 / 24  # default to a +/- 5% of the duration
    mapping["epoch_min"] = transit_specs[0]["epoch"] - epoch_error
    mapping["epoch_max"] = transit_specs[0]["epoch"] + epoch_error

    period_error = transit_specs[0].get("period_error", None)
    if period_error is None:
        period_error = transit_specs[0]["period"] * 0.5  # deault to +/- 50% of the period
    mapping["period_min"] = transit_specs[0]["period"] - period_error
    mapping["period_max"] = transit_specs[0]["period"] + period_error

    lc_time_label = lc.time.format.upper()
    if lc.time.format == "btjd":
        lc_time_label = "BJD - 2457000 (BTJD days)"
    elif lc.time.format == "bkjd":
        lc_time_label = "BJD - 2454833 (BKJD days)"
    mapping["lc_time_label"] = lc_time_label

    result = template
    for key, value in mapping.items():
        result = result.replace("{" + key + "}", str(value))

    if re.search(r"{[^}]+}", result):
        warnings.warn("create_input_fit(): the created `input_fit.py` still has values not yet defined.")

    input_fit_filepath = Path(pyaneti_target_in_dir, "input_fit.py")
    if write_to_file:
        input_fit_filepath.write_text(result)

    if return_content:
        return input_fit_filepath, result
    else:
        return input_fit_filepath


def display_pyaneti_input_py_location(input_fit_filepath):
    from IPython.display import display, HTML

    display(
        HTML(
            f"""
    <a href="{input_fit_filepath}" target="_input_fit_py">{input_fit_filepath}</a>
    """
        )
    )


def display_pyaneti_instructions(pti_env):
    from IPython.display import display, Markdown

    display(
        Markdown(
            f"""
Run `Pyaneti` to do the modeling:
```
cd {pti_env.home_dir}
python pyaneti.py  {pti_env.alias}
```
    """
        )
    )


#
# Read Pyaneti model output, lightcurve files, etc.
#


def display_model(
    pti_env,
    show_params=True,
    show_posterior=True,
    show_correlations=False,
    show_transits=True,
    show_lightcurve=True,
    show_chains=False,
):
    from IPython.display import display, Image, HTML

    target_out_dir = pti_env.target_out_dir
    alias = pti_env.alias
    display(HTML(f"<h3>Model for {alias}</h3>"))

    if show_params:
        file_params = Path(target_out_dir, f"{alias}_params.dat")
        url_params = f"file:///{file_params}".replace("\\", "/")
        file_init = Path(target_out_dir, f"{alias}_init.dat")
        url_init = f"file:///{file_init}".replace("\\", "/")
        display(
            HTML(
                f"""<ul>
    <li><a target="_params" href="{url_params}">Model params</a> - {file_params}</li>
    <li><a target="_init" href="{url_init}">Init params</a> - {file_init}</li>
</ul>
(Copy the link to open in a new tab if clicking them does not work.)
"""
            )
        )

    if show_posterior:
        display(Image(Path(target_out_dir, f"{alias}_posterior.png")))
    if show_correlations:
        display(Image(Path(target_out_dir, f"{alias}_correlations.png")))
    if show_transits:
        display(Image(Path(target_out_dir, f"{alias}b_tr.png")))  # TODO: handle multiple planets
    if show_lightcurve:
        display(Image(Path(target_out_dir, f"{alias}_lightcurve.png")))
    if show_chains:
        display(Image(Path(target_out_dir, f"{alias}_chains.png")))


def read_pyaneti_lc_dat(filename, time_format="btjd", time_converter_func=None):
    """Read Pyaneti lightcurve files as `LightCurve` objects, e.g.,
    `inpy/.../<starname>.dat`, `outpy/.../<starname>-trdata_lightcurve.txt`."""
    # format="ascii.commented_header" does not work
    tab = Table.read(filename, format="ascii")
    if len(tab.colnames) >= 3:
        (n_time, n_flux, n_flux_err, *rest) = tab.colnames
        flux = tab[n_flux]
        flux_err = tab[n_flux_err]
    else:
        n_time, n_flux = tab.colnames  # for trmodel files
        flux = tab[n_flux]
        flux_err = np.zeros_like(flux)

    if time_converter_func is not None:
        time = time_converter_func(tab[n_time])
    else:
        time = Time(tab[n_time], format=time_format)

    lc_cls = lk.LightCurve
    if time.format == "btjd":
        lc_cls = lk.TessLightCurve
    return lc_cls(time=time, flux=flux, flux_err=flux_err)
