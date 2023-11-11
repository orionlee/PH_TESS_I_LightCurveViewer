#/bin/sh

set -e

# Ideally want touch use -atime (last access time), 
# insetad of -mtime  (modified time)
# but -atime does notify seem touch work (no file selected)
find . -type f -name "*.fits" -mtime +365 -delete

# Delete download DV files
cd ~/Downloads
find . -type f -name "tess*.pdf" -mtime +30 -delete
cd -
