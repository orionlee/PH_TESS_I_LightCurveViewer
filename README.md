# TESS Target (TIC) Vetting

  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/orionlee/PH_TESS_I_LightCurveViewer/blob/main/TIC_Vetting.ipynb)

[TIC_Vetting.ipynb](TIC_Vetting.ipynb) - a Jupyter notebook that can be used a template for vetting / inspecting a particular TIC.

- A [visual walkthrough](TIC_Vetting_Walkthrough.md)
- it downloads the lightcurves of a given TIC, and let you inspect it in a variety of ways, including:
    - show zoomed-in plots
    - visualize transit period
    - check against background spikes
- it can also download the Target Pixel File, so that one can inspect it in per-pixel fashion.
    - useful for checking against contamination from nearby stars.


---

## Others: Interactive Subject Lightcurve Viewer

  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/orionlee/PH_TESS_I_LightCurveViewer/blob/main/PH_TESS_I_LightCurveViewer.ipynb)
  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/orionlee/PH_TESS_I_LightCurveViewer/main?filepath=PH_TESS_I_LightCurveViewer.ipynb)

An interactive (but more limited) one that uses the data available on the Planet Hunters TESS subject page, without relying on MAST Portal.

Instructions:

1. Goto the subject page of a star ([example](https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/talk/subjects/36971891)) on Planet Hunters TESS
2. Click the `2` button below the lightcurve. The data in json format will be shown.
3. Copy the data and paste it to the `Lightcurve data` textbox in the notebook.
4. The lightcurve will be plotted. You can then explore it using the controls to zoom in, show the curve in moving average, etc.

---

## Others: TCE Finder

  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/orionlee/PH_TESS_I_LightCurveViewer/blob/main/TCE_Finder.ipynb)

[TCE_Finder.ipynb](TCE_Finder.ipynb) - a Jupyter notebook that can be used a web-based tool to look up TCEs. Useful when TCEs are not yet available on exo.mast.

---

## Learn More

- [About Planet Hunters TESS](https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/about/research)
- [Light curve](https://en.wikipedia.org/wiki/Light_curve)
- [Exoplanets detection methods](https://en.wikipedia.org/wiki/Methods_of_detecting_exoplanets)

