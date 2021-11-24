#!/bin/bash
set -o errexit
set -o nounset
set -o pipefail

# list output files to detect new or missing files
ls -1

# 2. check on ID.cna.seq, from which other files are derived . : cat ID.cna.seq | cut -f 1-8 | md5sum
cat *.cna.seg | cut -f 1-8 | md5sum


#3. check on params.txt, from which ploidy and purity values are taken : cat id.params.txt | cut -f1-5 | md5sum
cat *.params.txt | cut -f 1-5 | md5sum


# check size of RData binary file
find . -iname "*.RData" -size +0 -printf '%p found and non-zero file size\n'

# calculate what report files are archived
tar --list --file *_plots.tar.gz
