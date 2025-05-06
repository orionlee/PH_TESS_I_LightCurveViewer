#/bin/sh

set -e

echo "Disk usage before cleanup"
df -h .
du -h -s mastDownload tesscut .

# Ideally want touch use -atime (last access time), 
# insetad of -mtime  (modified time)
# but -atime does notify seem touch work (no file selected)
find . -type f -name "*.fits" -mtime +120 -delete

# tesscut (from Eleanor, raw TessCut,  etc)
find tesscut -type f -name "*.fits" -mtime +14 -delete

# delete targetpixelfiles more aggressively
find . -type f -name "*_tp.fits" -mtime +14 -delete

# aggressively deletes the tesscut fits and the processed pickle dumps there
#   (they're huge)
cd /c/dev/_juypter/transit-APP/example_output/
find . -type f -name "*.fits" -mtime +2 -delete 
find . -type f -name "*.pickle" -mtime +2 -delete
cd - 

# aggressively deletes these copies (they're nost liked copied from .)
cd /c/dev/_juypter/TESSPositionalProbability/Lightcurves/
find . -type f -name "*.fits" -mtime +2 -delete 
cd -

# Delete download DV files
cd ~/Downloads
find . -type f -name "tess*.pdf" -mtime +2 -delete
cd -

echo "Disk usage after cleanup"
df -h .
du -h -s mastDownload tesscut .
