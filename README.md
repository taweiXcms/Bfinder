To setup Bfinder
=====

Branch for CMSSW_10XX Recommended using `CMSSW_10_3_1`

Check forest version in: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForestSetup#Setup_for_10_3_1_2018_PbPb_data

```
cmsrel CMSSW_10_3_1
cd CMSSW_10_3_1/src
cmsenv
git cms-merge-topic -u CmsHI:forest_CMSSW_10_3_1
# Switch to the branch HEAD
git remote add cmshi git@github.com:CmsHI/cmssw.git
git fetch cmshi --no-tags # don't fetch tags unless you have 20 mins to burn
git checkout -b forest_CMSSW_10_3_1 remotes/cmshi/forest_CMSSW_10_3_1
cd HeavyIonsAnalysis/JetAnalysis/python/jets
./makeJetSequences.sh
cd ../../../..
scram b -j4
```

To add D/Bfinder to forest:
=====

```
cd $CMSSW_BASE/src
cmsenv
git clone --branch CMSSW_10XX-F20190115_MC https://github.com/boundino/Bfinder.git --depth 1
source Bfinder/test/DnBfinder_to_Forest_103X.sh
scram b -j4
# mkdir -p bfinder && cp HeavyIonsAnalysis/JetAnalysis/test/runForestAOD_pponAA_DATA_103X_onlyDfinder.py bfinder/runForestAOD_pponAA_DATA_103X_onlyBfinder.py # Bfinder data
mkdir -p bfinder && cp HeavyIonsAnalysis/JetAnalysis/test/runForestAOD_pponAA_DATA_103X_onlyBfinder.py bfinder/runForestAOD_pponAA_MIX_103X_onlyBfinder.py # Bfinder MC
cd bfinder/
cmsRun runForestAOD_pponAA_MIX_103X_onlyBfinder.py
```
