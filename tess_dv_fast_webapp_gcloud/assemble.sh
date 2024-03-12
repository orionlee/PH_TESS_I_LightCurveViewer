#!/bin/sh

# Assemble TESS TCEs webapp that can be source-deployed to Google Cloud Run
#

base=`dirname $0`
dest=$1

set -e

mkdir -p $dest

mkdir -p $dest/data/tess_dv_fast
# --update --archive
cp --update --archive  $base/../data/tess_dv_fast/tess_tcestats.db  $dest/data/tess_dv_fast

cp --update --archive  $base/*.*  $dest
cp --update --archive  $base/../tess_dv_fast.py $base/../tess_dv_fast_webapp.py  $dest

ls -l $dest/data/tess_dv_fast/tess_tcestats.db

echo Source assembled. You can do the following for actual deployment:
echo cd $dest
# for local verification
echo python main.py
echo gcloud run deploy --source .
