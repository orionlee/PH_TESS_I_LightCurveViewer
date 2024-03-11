# TESS TCEs Webapp for Google Cloud Run

- Contains additional files needed for Google Cloud Run deployment


## Misc Notes

- requires entrypoint to be `main.py` with `app` attribute (of the flask app)
- triggered warning that memory usage exceed the default of 512Mb, increased
  the memory limit to 1Gb.  (Unclear on the root cause of the high memory usage.)
