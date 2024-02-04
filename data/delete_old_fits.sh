#/bin/sh

set -e

# Ideally want touch use -atime (last access time), 
# insetad of -mtime  (modified time)
# but -atime does notify seem touch work (no file selected)
find . -type f -name "*.fits" -mtime +365 -delete

# tesscut (from Eleanor)
find tesscut -type f -name "*.fits" -mtime +30 -delete

# Consider to delete targetpixelfiles more aggressively
# find . -type f -name "*_tp.fits" -mtime +30 -delete


# Delete download DV files
cd ~/Downloads
find . -type f -name "tess*.pdf" -mtime +30 -delete
cd -
