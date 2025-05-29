# ETV Analysis Notebooks

- de-facto template: ETV_extraction_coshgauss_342645261.ipynb
- Files with prefix `ETV_extraction_`: Analysis for the candidates from EB - SB period difference.
- Files with prefix `Misc_`: Miscellaneous targets
- Files with prefix that signifies a paper, e.g, `2013ApJ_`, `2019AJ_`, etc.: reproducing results of the correspond papers.


## Notebooks with special code paths

- `Misc_extraction_coshgauss_352830705_B.ipynb`: group multiple eclipses together (roughly by TESS sectors) for fitting due to low SNR
- `Misc_ETV_extraction_coshgauss_156846286.ipynb`: incorporate ASAS-SN data; grouping multiple eclipses together due to sparse data
- `Misc_extraction_coshgauss_122994468.ipynb`: codes to mask out a second set of eclipses
- `Misc_ETV_extraction_coshgauss_378275980.ipynb`:
  - both primary and secondary have varying depths over time. Custom indiviudal fit needed (as the standard one assumes fixed shape. (the varying depth make shape changes as well)
  - The period for primary is different from the period in secondary.
  - ETV in both per-sector and per-cycle basis.
  - Custom log_prior so that emcee will not sample amplitude (`alpha1`) that is way too large (on secondary)
- `Misc_ETV_extraction_coshgauss_220410224.ipynb`: a) handle varying eclipse depth, b) include sparse ASAS-SN data in ETV

