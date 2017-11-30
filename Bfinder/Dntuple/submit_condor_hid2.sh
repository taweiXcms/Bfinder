#!/bin/bash
###This is a modification of the default condor submission script for doing submission from the new hidisk2
###Initialize your grid certificate before submitting jobs
###TA-WEI WANG, 06/08/2017 created

DATASET=$1
DESTINATION=$2
MAXFILES=$3
LOGDIR=$4

PROXYFILE=$(ls /tmp/ -lt | grep $USER | grep -m 1 x509 | awk '{print $NF}')
CONFIGFILE="loop.C"
TRANSFERFILE="loop.C,loop.h,Dntuple.h,format.h"
OUTFILE="ntuple"

rm mylistfinal.txt
ls $DATASET  | awk '{print "" $0}' | grep .root >> mylistfinal.txt

SRM_PREFIX="/mnt/hadoop/"
SRM_PATH=${DESTINATION#${SRM_PREFIX}}

if [ ! -d $DESTINATION ]
then
    gfal-mkdir -p gsiftp://se01.cmsaf.mit.edu:2811/${SRM_PATH}
fi

if [ ! -d $LOGDIR ]
then
    mkdir -p $LOGDIR
fi

counter=0
for i in `cat mylistfinal.txt`
do
    if [ $counter -ge $MAXFILES ]
    then
	break
    fi
    # ifexist=`ls $DESTINATION/ntuple_$i`
    # if [ -z $ifexist ]
    # then
    if [ ! -f ${DESTINATION}/${OUTFILE}_$i ] && [ -f ${DATASET}/$i ]
    then
        echo -e "\033[38;5;242mSubmitting a job for output\033[0m ${DESTINATION}/${OUTFILE}_$i"
	infn=`echo $i | awk -F "." '{print $1}'`
	INFILE="${DATASET}/$i"
	
# make the condor file
	cat > subfile <<EOF

Universe = vanilla
Initialdir = .
Executable = exec_condor_hid2.sh
+AccountingGroup = "group_cmshi.$(whoami)"
Arguments =  $CONFIGFILE $DESTINATION ntuple_${infn}.root $INFILE $PROXYFILE
GetEnv       = True
Input = /dev/null

# log files
Output       = $LOGDIR/log-${infn}.out
Error        = $LOGDIR/log-${infn}.err
Log          = $LOGDIR/log-${infn}.log

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
requirements = GLIDEIN_Site == "MIT_CampusFactory" && BOSCOGroup == "bosco_cmshi" && HAS_CVMFS_cms_cern_ch && BOSCOCluster == "ce03.cmsaf.mit.edu"
transfer_input_files = $TRANSFERFILE,/tmp/$PROXYFILE

Queue
EOF

############################ Submit ###############################

#cat subfile
condor_submit subfile -name submit.mit.edu
# condor_submit subfile -pool submit.mit.edu:9615 -name submit.mit.edu -spool
mv subfile $LOGDIR/log-${infn}.subfile
counter=$((counter+1))
    fi
done

echo -e "Submitted \033[1;36m$counter\033[0m jobs to Condor."
