###Condor submitting template for plain root jobs. Run on all the files separately in a given folder, the DATASET folder below
###TA-WEI WANG, 02/20/2014 created
###Your plain root file needs to be modified , see example loop.C
###Please compile the .C before submission to make sure your code is working.
###Checking condor jobs status: condor_q <username> 

###Plain root .C to be run
CONFIGFILE="loop.C"

###All the header/related files needed
TRANSFERFILE="loop.C,loop.h,Dntuple.h,format.h"

###Folder location within which files are to be run
DATASET=/mnt/hadoop/cms/store/user/wangj/PAMinimumBias1/crab_DfinderData_PAMinimumBias1_pPb_20170407_PARun2016B_PromptReco_v1_Dstar_dPt0tkPt0p5/170407_213459/0000/

###Output file location
#DESTINATION=/mnt/hadoop/cms/store/user/jwang/Dmeson/ntD_20151012_HLTemulation_DinderMC_Pyquen_D0tokaonpion_D0pt15p0_Pthat15_TuneZ2_Unquenched_5020GeV_GENSIM_75x_v2_20151010_EvtBase
DESTINATION=/mnt/hadoop/cms/store/user/jwang/DntupleRpA/ntD_EvtBase_20170410_DfinderData_20170407_pPb_MinimumBias_PARun2016B_PromptReco_v1_AOD_MinimumBias1_dPt2tkPt0p5_Dstar_skim

###Maximum number of files to be run
MAXFILES=1000000000

###Log file location and it's name
LOGDIR=/export/d00/scratch/jwang/hadooplogs/ntD_EvtBase_20170410_DfinderData_20170407_pPb_MinimumBias_PARun2016B_PromptReco_v1_AOD_MinimumBias1_dPt2tkPt0p5_Dstar_skim

########################## Create subfile ###############################
rm mylistfinal.txt
ls $DATASET  | awk '{print "" $0}' | grep .root >> mylistfinal.txt

if [ ! -d $DESTINATION ]
then
    mkdir -p $DESTINATION
fi
if [ ! -d $LOGDIR ]
then
    mkdir -p $LOGDIR
fi

dateTime=$(date +%Y%m%d%H%M)
INFILE=""
fileCounter=0
for i in `cat mylistfinal.txt`
do
    if [ $fileCounter -ge $MAXFILES ]
    then
	break
    fi
    ifexist=`ls $DESTINATION/ntuple_$i`
    if [ -z $ifexist ]
    then
	infn=`echo $i | awk -F "." '{print $1}'`
	INFILE="$DATASET$i"
	
# make the condor file
	cat > subfile <<EOF

Universe = vanilla
Initialdir = .
Executable = exec_condor.sh
+AccountingGroup = "group_cmshi.$(whoami)"
Arguments =  $CONFIGFILE $DESTINATION ntuple_${infn}.root $INFILE
GetEnv       = True
Input = /dev/null

# log files
Output       = $LOGDIR/log-${infn}.out
Error        = $LOGDIR/log-${infn}.err
Log          = $LOGDIR/log-${infn}.log

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = $TRANSFERFILE

Queue
EOF
	
############################ Submit ###############################
	
#cat subfile
	condor_submit subfile
	mv subfile $LOGDIR/log-${infn}.subfile
	fileCounter=$((fileCounter+1))
    fi
done
echo "Submitted $fileCounter jobs to Condor."
