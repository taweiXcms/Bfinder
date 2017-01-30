#!/bin/sh

CONFIGFILE=$1
INFILE=$2
OUTFILE=$3
startEntries=$4
endEntries=$5
#DESTINATION=$6

#nohup root -l -b -q  'loophltpp.C+("${INFILE}","${OUTFILE}", false, ${startEntries}, ${endEntries}, false, true)' >& out_${startEntries}\_${endEntries} &
#mv ${OUTFILE} ${DESTINATION}/${OUTFILE}

root -l -b<<EOF
.x loop.C+("${INFILE}","${OUTFILE}", true, true, ${startEntries}, ${endEntries}, false, true)
.q
EOF
