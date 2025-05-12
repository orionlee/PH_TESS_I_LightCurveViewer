from types import SimpleNamespace


def _tce_reference_name_with_year(publication_year):
    return f"TESS Threshold Crossing Event (TCE), {publication_year}, (online data)"


# commonly used references
BIBS = SimpleNamespace(
    K2_N="Howell, S. B.; et al., 2014, The K2 Mission: Characterization and Early Results",
    K2_B="2014PASP..126..398H",
    TESS_N="Ricker, G. R.; et al., 2014, Transiting Exoplanet Survey Satellite (TESS)",
    TESS_B="2014SPIE.9143E..20R",
    TESS_SPOC_N="Caldwell, D. A.; et al., 2020, TESS Science Processing Operations Center FFI Target List Products",
    TESS_SPOC_B="2020RNAAS...4..201C",
    QLP_N="Huang, C. X.; et al., 2020, Photometry of 10 Million Stars from the First Two Years of TESS Full Frame Images: Part I",
    QLP_B="2020RNAAS...4..204H",
    GSFC_ELEANOR_LITE_N="Powell, B. P.; et al., 2022, The NASA GSFC TESS Full Frame Image Light Curve Data Set",
    GSFC_ELEANOR_LITE_B="2022RNAAS...6..111P",
    TCE_N=_tce_reference_name_with_year,
    # links to TCE is case specific
    TIC_N="Stassun, K. G.; et al., 2019, The Revised TESS Input Catalog and Candidate Target List",  # the paper describing TIC v8, the subsequent paper for v8.1/8.2 focuses mainly on the changes and is not as helpful
    TIC_B="2019AJ....158..138S",
    ASAS_SN_N="Kochanek, C. S.; et al., 2017, The All-Sky Automated Survey for Supernovae (ASAS-SN) Light Curve Server v1.0",
    ASAS_SN_B="2017PASP..129j4502K",
    ASAS_SNv2a_N="Hart, K.; et al, 2023, ASAS-SN Sky Patrol V2.0",
    ASAS_SNv2a_B="2023arXiv230403791H",
    ASAS_SNv2b_N="Shappee, B. J.; et al, 2014, The Man behind the Curtain: X-Rays Drive the UV through NIR Variability in the 2013 Active Galactic Nucleus Outburst in NGC 2617",
    ASAS_SNv2b_B="2014ApJ...788...48S",
    ASAS3_N="Pojmanski, G., 2002, The All Sky Automated Survey",
    ASAS3_B="2002AcA....52..397P",
    ZTF_N="Masci, F. J.; et al., 2019, The Zwicky Transient Facility: Data Processing, Products, and Archive",
    ZTF_B="2019PASP..131a8003M",
    GAIA_DR3_N="Gaia collaboration; et al., 2022, Gaia Data Release 3 (Gaia DR3) Part 1 Main source",
    GAIA_DR3_B="2022yCat.1355....0G",
    GAIA_DR3_ASTROPHY_N="Creevey, O. L.; et al., 2023, Gaia Data Release 3: Astrophysical parameters inference system (Apsis) I -- methods and content overview",
    GAIA_DR3_ASTROPHY_B="2023A&A...674A..26C",
    GAIA_DR3_NSS_N="Gaia collaboration; et al., 2022, Gaia Data Release 3 (Gaia DR3) Part 3 Non-single stars",
    GAIA_DR3_NSS_B="2022yCat.1357....0G",
    GAIA_DR3_VAR_N="Gaia collaboration; et al., 2022, Gaia Data Release 3 (Gaia DR3) Part 4 Variability",
    GAIA_DR3_VAR_B="2022yCat.1358....0G",
    TESSEB_N="Prša, A.; et al., 2022, TESS Eclipsing Binary Stars. I. Short-cadence Observations of 4584 Eclipsing Binaries in Sectors 1-26",
    TESSEB_B="2022ApJS..258...16P",
    TESSEB_LIVE_N="TESS Eclipsing Binary Catalogue (online data)",
    # links to live TESS EB is case specific
    # PHt II paper discussed methodlogy, with (indirect) mentions of eclipsingbinary tagging
    PHT_II_N="Eisner, N. L.; et al., 2021, Planet Hunters TESS II: findings from the first two years of TESS",
    PHT_II_B="2021MNRAS.501.4669E",
)
