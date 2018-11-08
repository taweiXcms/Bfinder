To setup Bfinder
=====

Branch for CMSSW_10XX Recommended using `CMSSW_10_3_1`

Check forest version in: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForestSetup#Setup%20for%2010_3_1%20(%202018%20_PbPb%20da

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
scram build -j4
```

To add Dfinder to forest:
=====

```
cd $CMSSW_BASE/src
cmsenv
git clone -b Dfinder_10XX https://github.com/taweiXcms/Bfinder.git
scram b -j4
source Bfinder/test/Dfinder_to_Forest.sh

# config to run Forest+Dfinder: HeavyIonsAnalysis/JetAnalysis/test/runForestAOD_pponAA_DATA_103X_wDfinder.py
# config to run Dfinder only: HeavyIonsAnalysis/JetAnalysis/test/runForestAOD_pponAA_DATA_103X_onlyDfinder.py

```
