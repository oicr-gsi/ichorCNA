#!/bin/bash
set -o errexit
set -o nounset
set -o pipefail

# list output files to detect new or missing files
ls -1

# count lines for .seg, .seg.txt, correctedDepth.txt and params.txt due to nondeterministic numeric fields
find . \( -iname "*.seg" -o -iname "*.seg.txt" -o -iname "*.correctedDepth.txt" -o -iname "*params.txt" \) -exec wc -l {} \;

# check size of RData binary file
find . -iname "*.RData" -size +0 -printf '%p found and non-zero file size\n'

# calculate what report files are archived
tar --list --file *_plots.tar.gz
