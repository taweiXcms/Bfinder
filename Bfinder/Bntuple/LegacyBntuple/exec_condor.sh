CONFIGFILE=$1
DESTINATION=$2
OUTFILE=$3
INFILE=$4

source /osg/app/cmssoft/cms/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
cd /cvmfs/cms.cern.ch/slc5_amd64_gcc462/cms/cmssw-patch/CMSSW_5_3_2_patch4/src
eval `scramv1 runtime -sh` 
cd -

#root -l -b -q  $CONFIGFILE++\(\"${INFILE}\",\"${OUTFILE}\"\)
root -l -b -q  $CONFIGFILE\(\"${INFILE}\",\"${OUTFILE}\"\)
mv ${OUTFILE} ${DESTINATION}/${OUTFILE}