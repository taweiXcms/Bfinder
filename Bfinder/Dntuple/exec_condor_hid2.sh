CONFIGFILE=$1
DESTINATION=$2
OUTFILE=$3
INFILE=$4
export X509_USER_PROXY=${PWD}/$5

source /osg/app/cmssoft/cms/cmsset_default.sh
#export SCRAM_ARCH=slc5_amd64_gcc462
#cd /cvmfs/cms.cern.ch/slc5_amd64_gcc462/cms/cmssw-patch/CMSSW_5_3_2_patch4/src
#export SCRAM_ARCH=slc6_amd64_gcc472
#cd /cvmfs/cms.cern.ch/slc6_amd64_gcc472/cms/cmssw/CMSSW_5_3_24/src
#export SCRAM_ARCH=slc6_amd64_gcc491
#cd /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_5_0/src
export SCRAM_ARCH=slc6_amd64_gcc491
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_5_8/src

eval `scramv1 runtime -sh` 
cd -

#root -l -b -q  $CONFIGFILE++\(\"${INFILE}\",\"${OUTFILE}\"\)
root -l -b -q  $CONFIGFILE\(\"${INFILE}\",\"${OUTFILE}\"\)
if [[ $? -eq 0 ]]; then
	#mv ${OUTFILE} ${DESTINATION}/${OUTFILE}
    gfal-copy file://${PWD}/${OUTFILE}  srm://se01.cmsaf.mit.edu:8443/srm/v2/server?SFN=${DESTINATION}/${OUTFILE}
fi
