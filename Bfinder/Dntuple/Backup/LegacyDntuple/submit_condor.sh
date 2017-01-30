###Condor submitting template for plain root jobs. Run on all the files separately in a given folder, the DATASET folder below
###TA-WEI WANG, 02/20/2014 created
###Your plain root file needs to be modified , see example loop.C
###Please compile the .C before submission to make sure your code is working.
###Checking condor jobs status: condor_q <username> 

###Plain root .C to be run
CONFIGFILE="loop.C"

###All the header/related files needed
TRANSFERFILE="loop.C,loop.h"

###Folder location within which files are to be run
#DATASET=/mnt/hadoop/cms/store/user/twang/HI_Dfinder/DfinderData_HIMinBiasUPC_HIRun2011-14Mar2014-v2_20150912/*
#DATASET=/mnt/hadoop/cms/store/user/twang/HI_DfinderNtuple/DinderMC_Pyquen_D0tokaonpion_D0pt1p0_Pthat0_TuneZ2_Unquenched_2760GeV_20150912/*
#DATASET=/mnt/hadoop/cms/store/user/twang/HI_DfinderNtuple/DinderMC_Pyquen_D0tokaonpion_D0pt1p0_Pthat15_TuneZ2_Unquenched_2760GeV_20150912/*
#DATASET=/mnt/hadoop/cms/store/user/twang/HI_DfinderNtuple/DinderMC_Pyquen_D0tokaonpion_D0pt1p0_Pthat30_TuneZ2_Unquenched_2760GeV_20150912/*
#DATASET=/mnt/hadoop/cms/store/user/twang/HI_DfinderNtuple/DinderMC_Pyquen_D0tokaonpion_D0pt1p0_Pthat50_TuneZ2_Unquenched_2760GeV_20150912/*
#DATASET=/mnt/hadoop/cms/store/user/twang/HI_Dfinder/DinderMC_richard-HydjetMB5020_750_75X_mcRun2_20150919/*
#DATASET=/mnt/hadoop/cms/store/user/twang/HI_Dfinder/DinderMC_richard-HydjetMB5020_750_75X_mcRun2_centrality30_100_20150927/*
DATASET=/mnt/hadoop/cms/store/user/twang/HI_DfinderNtuple_HLTemulation/DinderMC_Pyquen_D0tokaonpion_D0pt15p0_Pthat15_TuneZ2_Unquenched_5020GeV_GENSIM_75x_v2_20151010/*

###Output file location
#DESTINATION=/mnt/hadoop/cms/store/user/jwang/Dmeson/ntD_20150924_DfinderData_HIMinBiasUPC_HIRun2011-14Mar2014-v2_20150912_EvtBase_skim
#DESTINATION=/mnt/hadoop/cms/store/user/jwang/Dmeson/ntD_20150924_DinderMC_Pyquen_D0tokaonpion_D0pt1p0_Pthat0_TuneZ2_Unquenched_2760GeV_20150912_EvtBase_skim
#DESTINATION=/mnt/hadoop/cms/store/user/jwang/Dmeson/ntD_20150924_DinderMC_Pyquen_D0tokaonpion_D0pt1p0_Pthat15_TuneZ2_Unquenched_2760GeV_20150912_EvtBase_skim
#DESTINATION=/mnt/hadoop/cms/store/user/jwang/Dmeson/ntD_20150924_DinderMC_Pyquen_D0tokaonpion_D0pt1p0_Pthat30_TuneZ2_Unquenched_2760GeV_20150912_EvtBase_skim
#DESTINATION=/mnt/hadoop/cms/store/user/jwang/Dmeson/ntD_20150924_DinderMC_Pyquen_D0tokaonpion_D0pt1p0_Pthat50_TuneZ2_Unquenched_2760GeV_20150912_EvtBase_skim
#DESTINATION=/mnt/hadoop/cms/store/user/jwang/Dmeson/ntD_20150924_DinderMC_richard-HydjetMB5020_750_75X_mcRun2_20150919_EvtBase
#DESTINATION=/mnt/hadoop/cms/store/user/jwang/Dmeson/ntD_20151001_DinderMC_richard-HydjetMB5020_750_75X_mcRun2_centrality30_100_20150927_EvtBase_skim
DESTINATION=/mnt/hadoop/cms/store/user/jwang/Dmeson/ntD_20151012_HLTemulation_DinderMC_Pyquen_D0tokaonpion_D0pt15p0_Pthat15_TuneZ2_Unquenched_5020GeV_GENSIM_75x_v2_20151010_EvtBase

###Output file name
#OUTFILE="ntD_20150924_data_20150912_evtbase_skim"
#OUTFILE="ntD_20150924_mc_pthat0_20150912_evtbase_skim"
#OUTFILE="ntD_20150924_mc_pthat15_20150912_evtbase_skim"
#OUTFILE="ntD_20150924_mc_pthat30_20150912_evtbase_skim"
#OUTFILE="ntD_20150924_mc_pthat50_20150912_evtbase_skim"
#OUTFILE="ntD_20150924_DinderMC_richard-HydjetMB5020_750_75X_mcRun2_20150919_EvtBase"
#OUTFILE="ntD_20151001_DinderMC_richard-HydjetMB5020_750_75X_mcRun2_centrality30_100_20150927_EvtBase_skim"
OUTFILE="ntD_20151012_HLTemulation_DinderMC_Pyquen_D0tokaonpion_D0pt15p0_Pthat15_TuneZ2_Unquenched_5020GeV_GENSIM_75x_v2_20151010_EvtBase"

###Maximum number of files to be run
MAXFILES=6000

###Log file location and it's name
#LOGDIR=/export/d00/scratch/jwang/hadooplogs/ntD_20150924_data_20150912_EvtBase_skim
#LOGDIR=/export/d00/scratch/jwang/hadooplogs/ntD_20150924_mc_pthat0_20150912_EvtBase_skim
#LOGDIR=/export/d00/scratch/jwang/hadooplogs/ntD_20150924_mc_pthat15_20150912_EvtBase_skim
#LOGDIR=/export/d00/scratch/jwang/hadooplogs/ntD_20150924_mc_pthat30_20150912_EvtBase_skim
#LOGDIR=/export/d00/scratch/jwang/hadooplogs/ntD_20150924_mc_pthat50_20150912_EvtBase_skim
#LOGDIR=/export/d00/scratch/jwang/hadooplogs/ntD_20150924_DinderMC_richard-HydjetMB5020_750_75X_mcRun2_20150919_EvtBase
#LOGDIR=/export/d00/scratch/jwang/hadooplogs/ntD_20151001_DinderMC_richard-HydjetMB5020_750_75X_mcRun2_centrality30_100_20150927_EvtBase_skim
LOGDIR=/export/d00/scratch/jwang/hadooplogs/ntD_20151012_HLTemulation_DinderMC_Pyquen_D0tokaonpion_D0pt15p0_Pthat15_TuneZ2_Unquenched_5020GeV_GENSIM_75x_v2_20151010_EvtBase

#LOGNAME=log_ntD_20150924_data_20150912_EvtBase_skim
#LOGNAME=log_ntD_20150924_mc_pthat0_20150912_EvtBase_skim
#LOGNAME=log_ntD_20150924_mc_pthat15_20150912_EvtBase_skim
#LOGNAME=log_ntD_20150924_mc_pthat30_20150912_EvtBase_skim
#LOGNAME=log_ntD_20150924_mc_pthat50_20150912_EvtBase_skim
#LOGNAME=log_ntD_20150924_DinderMC_richard-HydjetMB5020_750_75X_mcRun2_20150919_EvtBase
#LOGNAME=log_ntD_20151001_DinderMC_richard-HydjetMB5020_750_75X_mcRun2_centrality30_100_20150927_EvtBase_skim
LOGNAME=log_ntD_20151012_HLTemulation_DinderMC_Pyquen_D0tokaonpion_D0pt15p0_Pthat15_TuneZ2_Unquenched_5020GeV_GENSIM_75x_v2_20151010_EvtBase

########################## Create subfile ###############################
dateTime=$(date +%Y%m%d%H%M)
fileCounter=0
INFILE=""
mkdir -p $DESTINATION
mkdir -p $LOGDIR

for file in $DATASET
do
    if [ $fileCounter -ge $MAXFILES ]
    then
	break
    fi

    INFILE="$file"
    fileCounter=$((fileCounter+1))

# make the condor file
cat > subfile <<EOF

Universe = vanilla
Initialdir = .
Executable = exec_condor.sh
+AccountingGroup = "group_cmshi.$(whoami)"
Arguments =  $CONFIGFILE $DESTINATION ${OUTFILE}_${fileCounter}.root $INFILE
GetEnv       = True
Input = /dev/null

# log files
Output       = $LOGDIR/$LOGNAME-$dateTime-${fileCounter}.out
Error        = $LOGDIR/$LOGNAME-$dateTime-${fileCounter}.err
Log          = $LOGDIR/$LOGNAME-$dateTime-${fileCounter}.log

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = $TRANSFERFILE

Queue
EOF

############################ Submit ###############################

#cat subfile
condor_submit subfile
mv subfile $LOGDIR/$LOGNAME-$dateTime-$fileCounter.subfile
done
echo "Submitted $fileCounter jobs to Condor."