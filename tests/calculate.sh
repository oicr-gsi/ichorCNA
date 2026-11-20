#!/bin/bash
set -o errexit
set -o nounset
set -o pipefail

#1. list output files to detect new or missing files
ls -1

#2. check on ID.cna.seq, from which other files are derived . : cat ID.cna.seq | cut -f 1-8 | md5sum
cat *.cna.seg | cut -f 1-8 | md5sum


#3. check metrics json, ignore loglik
cat *_metrics.json | tr , '\n' | grep -v loglik | md5sum

#4. check size of RData binary file
find . -iname "*.RData" -size +0 -printf '%p found and non-zero file size\n'

#5. calculate what report files are archived
tar --list --file *_plots.tar.gz
