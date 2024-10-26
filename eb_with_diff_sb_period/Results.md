<style>
body { counter-reset: h1counter h2counter h3counter h4counter h5counter h6counter; }

h1 { counter-reset: h2counter; }
h2 { counter-reset: h3counter; }
h3 { counter-reset: h4counter; }
h4 { counter-reset: h5counter; }
h5 { counter-reset: h6counter; }
h6 {}

h2:before {
    counter-increment: h2counter;
    content: counter(h2counter) ".\0000a0\0000a0";
}

h3:before {
    counter-increment: h3counter;
    content: counter(h2counter) "." counter(h3counter) ".\0000a0\0000a0";
}

h4:before {
    counter-increment: h4counter;
    content: counter(h2counter) "." counter(h3counter) "." counter(h4counter) ".\0000a0\0000a0";
}

h5:before {
    counter-increment: h5counter;
    content: counter(h2counter) "." counter(h3counter) "." counter(h4counter) "." counter(h5counter) ".\0000a0\0000a0";
}

h6:before {
    counter-increment: h6counter;
    content: counter(h2counter) "." counter(h3counter) "." counter(h4counter) "." counter(h5counter) "." counter(h6counter) ".\0000a0\0000a0";
}
</style>

<!--
regex for images
^!(.+)$
<img src="results_assets/$1" width=429 height=323 alt="$1">
-->

# Identifying Multi Star System Candidates using VSX, Gaia DR3 NSS and TESS

## Results

- Identified 51 targets of which the orbital period from eclipses (VSX / TESS) is different from the period from spectroscopy (Gaia DR3 NSS)
- ETV analysis was attempted for the targets. The breakdown:

| Disposition                               | Num. of TICs |
| ----------------------------------------- | ------------ |
| Y: Strong signs of nonlinear trend in ETV | 6            |
| Y?: Might be nonlinear trend in ETV       | 10           |
| ?: Uncertain ETV results                  | 22           |
| N?: Most likely no ETV                    | 4            |
| Skip: insufficient TESS data for ETV      | 9            |

See methods section on how these candidates are selected and how ETV analysis was carried out.

### Candidates with strong signs of nonlinear trend in ETV

Those with `has_etv = Y` in candidates table.

#### TIC 353894978

EB to SB Period ~ 1:3 . O-C period: possibly ~141 days, ~184 days

<img src="results_assets/oc_plot_tic353894978.png" width=429 height=323 alt="O-C plot TIC 353894978">


#### TIC 30034081

EB to SB Period ~ 2.8:1 . O-C period: no clear one.

<img src="results_assets/oc_plot_tic30034081.png" width=429 height=323 alt="O-C plot TIC 30034081">


#### TIC 154453805

EB : SB Period ~ 1:3.07 . O-C period: possibly 1000+ days or ~215 days.

<img src="results_assets/oc_plot_tic154453805.png" width=429 height=323 alt="O-C plot TIC 154453805">


#### TIC 290035575

EB : SB Period ~= 9:10 . O-C period: ~7.6 days

It is worth noting that the orbital period from eclipses is ~0.95 day, SB Period ~1.05 day, and O-C period is ~7.6 days.

<img src="results_assets/oc_plot_tic290035575.png" width=429 height=323 alt="O-C plot TIC 290035575">

The periodic ETV variation is visible in a phase plot of the O-C:

<img src="results_assets/oc_plot_phase_tic290035575.png" width=429 height=323 alt="O-C phase plot TIC 290035575">


#### TIC 410354930

EB : SB Period ~= 1.7 : 1 . O-C Period: unclear.

<img src="results_assets/oc_plot_tic410354930.png" width=429 height=323 alt="O-C plot TIC 410354930">


#### TIC 167304040

EB : SB Period ~= 7.1 : 1 . O-C Period: unclear.

<img src="results_assets/oc_plot_tic167304040.png" width=429 height=323 alt="O-C plot TIC 167304040">


### Samples of other candidates

#### TIC 38383256

`has_etv = Y?`. Anti-correlated timings: The trend of primary is opposite of the one of secondary. Possibly long term non-linear trend as well.

<img src="results_assets/oc_plot_tic38383256.png" width=429 height=323 alt="oc_plot_tic38383256.png">

<img src="results_assets/oc_plot_zoom1_tic38383256.png" width=200 height=323 alt="oc_plot_zoom1_tic38383256.png"> <img src="results_assets/oc_plot_zoom2_tic38383256.png" width=200 height=323 alt="oc_plot_zoom2_tic38383256.png"> <img src="results_assets/oc_plot_zoom3_tic38383256.png" width=200 height=323 alt="oc_plot_zoom3_tic38383256.png">


#### TIC 36883123

`has_etv = Y?`. Anti-correlated timings: The trend of primary is opposite of the one of secondary.
No clear long-term non-linear trend (a small possibility of a ~14 d period). Only 2 sectors of data.

Other targets with anti-correlated timings, but have no clear long-term non-linear trend (possibly due to lack of data)
TIC 57297550 (possibly O-C period of ~21 d), TIC 440051305, TIC 11546041

<img src="results_assets/oc_plot_tic36883123.png" width=429 height=323 alt="oc_plot_tic36883123.png">

<img src="results_assets/oc_plot_zoom1_tic36883123.png" width=300 height=323 alt="oc_plot_zoom1_tic36883123.png"> <img src="results_assets/oc_plot_zoom2_tic36883123.png" width=300 height=323 alt="oc_plot_zoom2_tic36883123.png">


#### TIC 52368472

`has_etv = Y?`. Unsure if the trend is non-linear.

<img src="results_assets/oc_plot_tic52368472.png" width=429 height=323 alt="oc_plot_tic52368472.png">

<img src="results_assets/oc_plot_zoom1_tic52368472.png" width=200 height=323 alt="oc_plot_zoom1_tic52368472.png"> <img src="results_assets/oc_plot_zoom2_tic52368472.png" width=200 height=323 alt="oc_plot_zoom2_tic52368472.png"> <img src="results_assets/oc_plot_zoom3_tic52368472.png" width=200 height=323 alt="oc_plot_zoom3_tic52368472.png">


#### TIC 231922417

`has_etv = Y?`. Apparent non-linear trend. Huge error in O-C, especially in the 30 minute cadence data in the early sectors.

<img src="results_assets/oc_plot_tic231922417.png" width=429 height=323 alt="oc_plot_tic231922417.png">

<img src="results_assets/oc_plot_zoom1_tic231922417.png" width=300 height=323 alt="oc_plot_zoom1_tic231922417.png"> <img src="results_assets/oc_plot_zoom2_tic231922417.png" width=300 height=323 alt="oc_plot_zoom2_tic231922417.png">


#### TIC 273874851

`has_etv = Y?`. Possibly long-term (1000d+) non-linear trend.

<img src="results_assets/oc_plot_tic273874851.png" width=429 height=323 alt="oc_plot_tic273874851.png">

<img src="results_assets/oc_plot_zoom1_tic273874851.png" width=200 height=323 alt="oc_plot_zoom3_tic273874851.png"> <img src="results_assets/oc_plot_zoom2_tic273874851.png" width=200 height=323 alt="oc_plot_zoom3_tic273874851.png"> <img src="results_assets/oc_plot_zoom3_tic273874851.png" width=200 height=323 alt="oc_plot_zoom3_tic273874851.png">


---

## Possible Future Work


---

## Methods

### Candidate Selection

- Cross matched 3485 TICs [*] known to be an EB in VSX, with Gaia DR3 NSS SB type.
- 515 TICs are SBs per Gaia DR3.
- 65 TICs: period from eclipses in VSX is different from the period from spectroscopy in Gaia DR3 NSS or aliases.
    - the periods are considered the same if the ratio of VSX Period to SB Period are in [0.99, 1.01], [0.495, 0.505], or [1.98, 2.02].
- 51 TICs: period from eclipses are indeed different from the period from spectroscopy after inspecting data in TESS.
    - 14 are deemed false positives. Most of them are because the period in VSX are not accurate.

- [*] The initial 3485 TICs come from TICs tagged as eclipsing binary by the volunteers of [Zooniverse Planet Hunters TESS](https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/talk/tags/eclipsingbinary) that are also known EB in VSX. The choice is primarily a matter of expedience. The dataset can be seen as a sample of known EB in VSX which also has some TESS data in 2-minute cadence.


### Eclipse Timing Variations (ETV) Analysis

- ETV analysis has been carried out for the 51 TICs based on their TESS data.
    - The TESS data used are based on 2-minute cadence data provided by SPOC, supplemented with 10-minute / 30-minute cadence data from TESS-SPOC or QLP in sectors where 2-minute cadence data is not available.
- Eclipse midpoint ephemeris is identified by fitting a geometric `cosh` Gaussian model to individual eclipses with [`emcee`](https://emcee.readthedocs.io/), using a custom version of [`https://github.com/noraeisner/ETVs`](https://github.com/noraeisner/ETVs) .
    - The initial eclipse midpoint and period are determined using a variety of means,including consulting the values from VSX and TESS TCEs, manual fitting, and MCMC-based fitting of a geometric `cosh` Gaussian model with period as an additional free parameter.
    - In most cases, the model fitting is done on the central portion of eclipses, excluding ingress/egress, as fitting on ingress/egress tends to be less accurate, and is not essential for the purpose of determining eclipse midpoints.
    - For cases the underlying TESS data are 30 minute cadence and eclipse duration is short, the modeled midpoints often have large errors (or are altogether excluded) because of limited number of data points per eclipse. For example, an eclipse with 4 hour duration would imply at most 8 data points in 30 minute cadence data.


---

## TODOs

