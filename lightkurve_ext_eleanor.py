import traceback
import warnings

from astropy.units import Quantity
from astropy.wcs import WCS

import numpy as np

from lightkurve import TessTargetPixelFile, LightkurveWarning
from lightkurve.targetpixelfile import HduToMetaMapping
from lightkurve.io.eleanor import read_eleanor_lightcurve


class EleanorTargetPixelFile(TessTargetPixelFile):
    """Read Eleanor-generated TargetPixelFile."""

    def __init__(self, path, quality_bitmask="default", **kwargs):
        super().__init__(path, quality_bitmask=quality_bitmask, **kwargs)

        tic_id = self.get_header().get("TIC_ID")
        if self.targetid is None:
            self.targetid = tic_id
        # mimic SPOC-generated TPF headers for compatibilities
        self.hdu[0].header["TICID"] = tic_id
        self.hdu[0].header["TESSMAG"] = self.hdu[0].header.get("TMAG")
        # meta needs to be updated after hdu[0].header is changed (ugly)
        self.meta = HduToMetaMapping(self.hdu[0])

    def hdu(self, value, keys=("TPF", "QUALITY")):
        super().hdu(value, keys)

    # TODO:
    # - override  _parse_aperture_mask(self, aperture_mask):
    #   to support names of eleanor aperture
    # - mimic TICID and TESSMAG header
    # - handle tpf.background_mask

    @property
    def cadenceno(self):
        """Return the cadence number for all good-quality cadences."""
        return self.hdu[1].data["FFIINDEX"][self.quality_mask]

    @property
    def flux(self) -> Quantity:
        """Returns the flux for all good-quality cadences."""
        return Quantity(self.hdu[1].data["TPF"][self.quality_mask], unit="electron/s")

    @property
    def flux_err(self) -> Quantity:
        """Returns the flux uncertainty for all good-quality cadences."""
        return Quantity(self.hdu[1].data["TPF_ERR"][self.quality_mask], unit="electron/s")

    @property
    def flux_bkg(self) -> Quantity:
        """Returns the background flux for all good-quality cadences."""
        return Quantity(np.full_like(self.flux, np.nan), unit="electron/s")

    @property
    def flux_bkg_err(self) -> Quantity:
        return Quantity(np.full_like(self.flux, np.nan), unit="electron/s")

    @property
    def column(self):
        """CCD pixel column number ('1CRV5P' header keyword)."""
        return self.get_keyword("1CRV5P", hdu=0, default=0)

    @property
    def row(self):
        """CCD pixel row number ('2CRV5P' header keyword)."""
        return self.get_keyword("2CRV5P", hdu=0, default=0)

    @property
    def pos_corr1(self):
        """Returns the column position correction. N/A for Eleanor."""
        return np.full_like(self.quality_mask, 0.0)

    @property
    def pos_corr2(self):
        """Returns the row position correction. N/A for Eleanor."""
        return np.full_like(self.quality_mask, 0.0)

    @property
    def wcs(self) -> WCS:
        """Returns an `astropy.wcs.WCS` object with the World Coordinate System
        solution for the target pixel file.

        Returns
        -------
        w : `astropy.wcs.WCS` object
            WCS solution
        """
        # the same as the SPOC branch in TessTargetPixelFile.wcs(),
        # except the hdu index is different

        # For standard (Ames-pipeline-produced) TPF files, we use the WCS
        # keywords provided in the first extension (the data table extension).
        # Specifically, we use the WCS keywords for the 5th data column (FLUX).
        wcs_keywords = {
            "1CTYP5": "CTYPE1",
            "2CTYP5": "CTYPE2",
            "1CRPX5": "CRPIX1",
            "2CRPX5": "CRPIX2",
            "1CRVL5": "CRVAL1",
            "2CRVL5": "CRVAL2",
            "1CUNI5": "CUNIT1",
            "2CUNI5": "CUNIT2",
            "1CDLT5": "CDELT1",
            "2CDLT5": "CDELT2",
            "11PC5": "PC1_1",
            "12PC5": "PC1_2",
            "21PC5": "PC2_1",
            "22PC5": "PC2_2",
            "NAXIS1": "NAXIS1",
            "NAXIS2": "NAXIS2",
        }
        mywcs = {}
        for oldkey, newkey in wcs_keywords.items():
            if self.hdu[0].header.get(oldkey, None) is not None:
                mywcs[newkey] = self.hdu[0].header[oldkey]
        return WCS(mywcs)

    @property
    def pipeline_mask(self):
        """Returns the optimal aperture mask used by the pipeline.

        If the aperture extension is missing from the file, a mask
        composed of all `True` values will be returned.
        """
        # Both Kepler and TESS flag the pixels in the optimal aperture using
        # bit number 2 in the aperture mask extension, e.g. see Section 6 of
        # the TESS Data Products documentation (EXP-TESS-ARC-ICD-TM-0014.pdf).

        try:
            aperture_col = self.hdu[0].header.get("APERTURE")
            return self.get_mask_of_aperture(aperture_col)
        except Exception as err:
            # should not happen with Eleanor tpf
            err_str = f"{type(err).__name__}: {err}\n" + "".join(traceback.format_exc())
            warnings.warn(
                f"Unexpected error in returning pipeline_mask. Return an empty one. The error: {err_str}",
                LightkurveWarning,
            )
            return np.ones(self.hdu[1].data["TPF"][0].shape, dtype=bool)

    def get_bkg_lightcurve(self, aperture_mask=None):
        if not (aperture_mask is None or aperture_mask == "pipeline"):
            return NotImplementedError("Eleanor TargetPixelFile supports background lightcurve from pipeline only.")
        # user read_eleanor_lightcurve() so that
        # I can pass self.hdu instead of self.path
        # (self.path might not be valid if the tpf is inited with a HDUList object)
        lc = read_eleanor_lightcurve(self.hdu).select_flux("flux_bkg")
        lc.flux_err = lc.flux_err.value * lc.flux.unit  # workaround lightkurve bug
        return lc["time", "flux", "flux_err", "centroid_col", "centroid_row", "cadenceno", "quality"]

    #
    # new methods specific to Eleanor
    #

    def aperture_names(self):
        """Return the list of aperture names available."""
        return [c.name for c in self.hdu[2].data.columns]

    def get_mask_of_aperture(self, name):
        return self.hdu[2].data[name] > 0
