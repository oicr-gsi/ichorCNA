#!/bin/bash
set -o errexit
set -o nounset
set -o pipefail

# list output files to detect new or missing files
ls -1

# params.txt output should be identical
md5sum *params.txt

# only count lines for .seg, .seg.txt and correctedDepth.txt due to nondeterministic numeric fields
find . \( -iname "*.seg" -o -iname "*.seg.txt" -o -iname "*.correctedDepth.txt" \) -exec wc -l {} \;

# check size of RData binary file
du -h *.RData

# calculate what report files are archived
tar --list --file *_plots.tar.gz
