#!/bin/sh

INFILE=$1
OUTFILE=$2
isPbPb=$3

root -l -b<<EOF
.x loop.C+("${INFILE}","${OUTFILE}", true, ${isPbPb}, 0, -1, false, true, true)
.q
EOF
