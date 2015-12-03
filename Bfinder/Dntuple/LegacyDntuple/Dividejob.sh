#!/bin/sh

#./Dividejob.sh root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/Run2015E/PromptReco/MinimumBias/HiForest_pp_run262163PromptReco_7M.root ntD_pp_run262163PromptReco_MBPD_7MEvents 500000 16 

INFILE=$1
OUTFILE=$2
EvtPerJob=$3
Njobs=$4
echo 'EvtPerJob: ' $EvtPerJob

for ((count=1; count <= $Njobs; count++))
do
 echo "Job $count"
 startevt=$(((count-1)*EvtPerJob))
 endevt=$((count*EvtPerJob))
 echo "from $startevt to $endevt"
 nohup ./exec_test.sh loop.C $INFILE $OUTFILE $startevt $endevt >& out_${OUTFILE}_${startevt}_${endevt} &
 sleep 5 
done
