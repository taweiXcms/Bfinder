#!/bin/bash

if [[ $0 != *.sh ]]
then
    echo -e "\e[31;1merror:\e[0m use \e[32;1m./script.sh\e[0m instead of \e[32;1msource script.sh\e[0m"
    return 1
fi

MAXFILENO=5

INPUTDIR="/mnt/hadoop/cms/store/user/wangj/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/crab_DfinderMC_pp_20171129_dPt1tkPt0p5dDecay2p5_D0/171130_013258/0000/"
OUTPUTPRIDIR="/mnt/hadoop/cms/store/user/jwang/DntupleRpAgfal"
OUTPUTSUBDIR="Dntuple_20171129_DfinderMC_PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8_20171129_dPt1tkPt0p5dDecay2p5_D0"
WORKDIR="/work/$USER/Dntuple"

OUTPUTDIR="${OUTPUTPRIDIR}/${OUTPUTSUBDIR}"
LOGDIR="logs/log_${OUTPUTSUBDIR}"

movetosubmit=${1:-0}
runjobs=${2:-0}

##

if [ ! -d $WORKDIR ]
then
    mkdir -p $WORKDIR
fi

##

if [ "$movetosubmit" -eq 1 ]
then
    if [[ $(hostname) == "submit-hi2.mit.edu" || $(hostname) == "submit.mit.edu" ]]
    then
        cp loop.C $WORKDIR/
        cp loop.h $WORKDIR/
        cp Dntuple.h $WORKDIR/
        cp format.h $WORKDIR/
        cp $0 $WORKDIR/
        cp exec_condor_hid2.sh $WORKDIR/
        cp submit_condor_hid2.sh $WORKDIR/
    else
        echo -e "\e[31;1merror:\e[0m compile macros on \e[32;1msubmit-hiX.mit.edu\e[0m or \e[32;1msubmit.mit.edu\e[0m."
    fi
fi

if [ "$runjobs" -eq 1 ]
then
    if [[ $(hostname) == "submit.mit.edu" ]]
    then
        ./submit_condor_hid2.sh $INPUTDIR $OUTPUTDIR $MAXFILENO $LOGDIR
    else
        echo -e "\e[31;1merror:\e[0m submit jobs on \e[32;1msubmit.mit.edu\e[0m."
    fi
fi