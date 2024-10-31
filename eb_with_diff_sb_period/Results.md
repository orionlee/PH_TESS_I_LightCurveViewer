<!-- Numbering (does not work on github)
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

-->

<!--
regex for images
^!(.+)$
<img src="results_assets/$1" width=429 height=323 alt="$1">
-->

# Identifying Multi Star System Candidates using VSX, Gaia DR3 NSS and TESS

## 1. Results

- Identified 51 candidates (from over 3000 known EBs), of which the orbital period from eclipses (VSX / TESS) is different from the period from spectroscopy (SB in Gaia DR3 NSS)
- ETV analysis was attempted for the candidates. The breakdown:

| Disposition                               | Num. of TICs |
| ----------------------------------------- | ------------ |
| `Y`: Strong signs of nonlinear trend in ETV | 6            |
| `Y?`: Might be nonlinear trend in ETV       | 10           |
| `?`: Uncertain ETV results                  | 22           |
| `N?`: Most likely no ETV                    | 4            |
| `Skip`: insufficient TESS data for ETV      | 9            |

<!--
    For the ETV table below, use the link to colab version of the notebook, so that for the notebook links of the individual targets, colab version is used.
    It's done so because github's rendering for individual target notebooks is not good.
-->
- See the table in Eclipse Timing Variation (ETV) Analysis section of this [notebook](https://colab.research.google.com/github/orionlee/PH_TESS_I_LightCurveViewer/blob/main/eb_with_diff_sb_period/Dashboard_EB_with_SB.ipynb).
    - Some of them, including the most noteworthy ones, are presented in the next section.
- See methods section on how these candidates are selected and how ETV analysis was carried out.


### 1.1 Candidates with strong signs of nonlinear trend in ETV

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


### 1.2 Candidates with some signs of nonlinear trend in ETV

Those with `has_etv = Y?` in candidates table.

#### TIC 38383256

 Anti-correlated timings: The trend of primary is opposite of the one of secondary. Possibly long term non-linear trend as well.

<img src="results_assets/oc_plot_tic38383256.png" width=429 height=323 alt="oc_plot_tic38383256.png">

<img src="results_assets/oc_plot_zoom1_tic38383256.png" width=200 height=323 alt="oc_plot_zoom1_tic38383256.png"> <img src="results_assets/oc_plot_zoom2_tic38383256.png" width=200 height=323 alt="oc_plot_zoom2_tic38383256.png"> <img src="results_assets/oc_plot_zoom3_tic38383256.png" width=200 height=323 alt="oc_plot_zoom3_tic38383256.png">


#### TIC 36883123

Anti-correlated timings: The trend of primary is opposite of the one of secondary.
No clear long-term non-linear trend (a long-shot possibility of a ~14 d O-C period). Only 2 sectors of data.

<img src="results_assets/oc_plot_tic36883123.png" width=429 height=323 alt="oc_plot_tic36883123.png">

<img src="results_assets/oc_plot_zoom1_tic36883123.png" width=300 height=323 alt="oc_plot_zoom1_tic36883123.png"> <img src="results_assets/oc_plot_zoom2_tic36883123.png" width=300 height=323 alt="oc_plot_zoom2_tic36883123.png">

Other targets with similar anti-correlated timings, but have no clear long-term non-linear trend (possibly due to lack of data):

| Target        | EB Period (d) | SB Period (d) | (Long-shot) O-C Period (d) |
| ------------- | ------------- | ------------- | -------------------------- |
| TIC 57297550  | 0.38          | 292.9         | ~21                        |
| TIC 440051305 | 0.33          | 63.3          | ~48                        |
| TIC 11546041  | 0.41          | 1.45          | ~6                         |


#### TIC 52368472

 Unsure if the trend is non-linear (a long-shot O-C period of ~211 d)

<img src="results_assets/oc_plot_tic52368472.png" width=429 height=323 alt="oc_plot_tic52368472.png">

<img src="results_assets/oc_plot_zoom1_tic52368472.png" width=200 height=323 alt="oc_plot_zoom1_tic52368472.png"> <img src="results_assets/oc_plot_zoom2_tic52368472.png" width=200 height=323 alt="oc_plot_zoom2_tic52368472.png"> <img src="results_assets/oc_plot_zoom3_tic52368472.png" width=200 height=323 alt="oc_plot_zoom3_tic52368472.png">


#### TIC 231922417

Apparent non-linear trend. Huge error in O-C, especially in the 30 minute cadence data in the early sectors.

<img src="results_assets/oc_plot_tic231922417.png" width=429 height=323 alt="oc_plot_tic231922417.png">

<img src="results_assets/oc_plot_zoom1_tic231922417.png" width=300 height=323 alt="oc_plot_zoom1_tic231922417.png"> <img src="results_assets/oc_plot_zoom2_tic231922417.png" width=300 height=323 alt="oc_plot_zoom2_tic231922417.png">


#### TIC 273874851

Possibly long-term (100d+ or 1000d+) non-linear trend.

<img src="results_assets/oc_plot_tic273874851.png" width=429 height=323 alt="oc_plot_tic273874851.png">

<img src="results_assets/oc_plot_zoom1_tic273874851.png" width=200 height=323 alt="oc_plot_zoom3_tic273874851.png"> <img src="results_assets/oc_plot_zoom2_tic273874851.png" width=200 height=323 alt="oc_plot_zoom3_tic273874851.png"> <img src="results_assets/oc_plot_zoom3_tic273874851.png" width=200 height=323 alt="oc_plot_zoom3_tic273874851.png"> <img src="results_assets/oc_plot_zoom4_tic273874851.png" width=200 height=323 alt="oc_plot_zoom4_tic273874851.png">


### 1.3 Samples of other candidates

The remaining candidates do not appear to show strong non linear trend, but most of them cannot be discounted (from having ETVs) outright either. Some of them are presented here.

TODO:

---

## 2. Discussion

### 2.1 Open Issues

- Validity of ETVs, especially those with apparent non-linear trends (`has_etv = Y or Y?`)
    - For those with possible O-C period (or rough O-C period range) identified, none of them match the correspond SB period.
- Any values for the candidates wit no apparent ETV but EB Period appears to be different from SB Period


### 2.2 Possible Future Work

- Go in depth of some interesting candidates.
- Go in breath with a more systematic and complete survey, e.g., all known EBs in VSX and in TESS continuous viewing zone (CVZ) with different SB periods. (The CVZ constraint would ensure the identified candidates have ample of TESS data for ETV analysis.)

---

## 3. Methods

### 3.1 Candidate Selection

- Cross matched 3485 TICs [*] known to be EBs in VSX, with Gaia DR3 NSS SB type.
- 515 TICs are SBs per Gaia DR3.
- 65 TICs: period from eclipses in VSX is different from the period from spectroscopy in Gaia DR3 NSS or aliases.
    - the periods are considered the same or are aliases, if the ratio of VSX Period to SB Period are in [0.99, 1.01], [0.495, 0.505], or [1.98, 2.02].
- 51 TICs: period from eclipses are indeed different from the period from spectroscopy after inspecting data in TESS.
    - The remaining 14 TICs are deemed false positives, i.e., no significant difference between EB Period and SB Period. Most of them are because the periods in VSX are not accurate.

- [*] The initial 3485 TICs come from TICs tagged as eclipsing binary by the volunteers of [Zooniverse Planet Hunters TESS](https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/talk/tags/eclipsingbinary) that are also known EB in VSX. The choice is primarily a matter of expedience. The dataset can be seen as a sample of known EBs in VSX which also have some TESS data in 2-minute cadence.


### 3.2 Eclipse Timing Variations (ETV) Analysis

- ETV analysis has been carried out for the 51 TICs based on their TESS data.
    - The TESS data used are based on 2-minute cadence data provided by SPOC, supplemented with 10-minute / 30-minute cadence data from TESS-SPOC or QLP in sectors where 2-minute cadence data is not available.
- Eclipse midpoint ephemeris is identified by fitting a geometric `cosh` Gaussian model to individual eclipses with [`emcee`](https://emcee.readthedocs.io/), using a custom version of [`https://github.com/noraeisner/ETVs`](https://github.com/noraeisner/ETVs) .
    - The initial eclipse midpoint and period are determined using a variety of means,including consulting the values from VSX and TESS TCEs, manual fitting, and MCMC-based fitting of a geometric `cosh` Gaussian model with period as an additional free parameter.
    - In most cases, the model fitting is done on the central portion of eclipses, excluding ingress/egress, as fitting on ingress/egress tends to be less accurate, and is not essential for the purpose of determining eclipse midpoints.
    - For cases the underlying TESS data are 30 minute cadence and eclipse duration is short, the modeled midpoints often have large errors (or are altogether excluded) because of limited number of data points per eclipse. For example, an eclipse with 4 hour duration would imply at most 8 data points in 30 minute cadence data.

