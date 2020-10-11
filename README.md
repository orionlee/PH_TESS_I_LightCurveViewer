# Planet Hunters TESS Subject Lightcurve Viewers
Jupyter Notebooks that lets you interactively visualize the a star's lightcurve provided by Planet Hunters TESS:
    https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess

## TIC Lightcurve Viewer 
  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/orionlee/PH_TESS_I_LightCurveViewer/blob/master/TIC_Lightcurve_Viewer.ipynb)

A notebook that given a TIC, produces lightcurves, with many other options.
- It is essentially a template that you're likely need to modify to suit the specific needs.
- In simple cases, you  just need to replace the TIC, but you do need to know some Python / Jupyter notebooks.


## Interactive Subject Lightcurve Viewer 
  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/orionlee/PH_TESS_I_LightCurveViewer/master?filepath=PH_TESS_I_LightCurveViewer.ipynb)

An interactive (but more limited) one that uses the data available on the Planet Hunters TESS subject page, without relying on MAST Portal. 

### Instructions
1. Goto the subject page of a star ([example](https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/talk/subjects/36971891)) on Planet Hunters TESS
2. Click the `2` button below the lightcurve. The data in json format will be shown.
3. Copy the data and paste it to the `Lightcurve data` textbox in the notebook.
4. The lightcurve will be plotted. You can then explore it using the controls to zoom in, show the curve in moving average, etc.

---
## Learn More
- [About Planet Hunters TESS](https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess/about/research)
- [Light curve](https://en.wikipedia.org/wiki/Light_curve)
- [Exoplanets detection methods](https://en.wikipedia.org/wiki/Methods_of_detecting_exoplanets)

