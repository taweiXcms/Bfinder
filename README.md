D meson reconstruction is added to this version in CMSSW_5_3_20

for CMSSW5 (CMSSW_5_3_20)
git cms-merge-topic -u CmsHI:forest_CMSSW_5_3_20
git clone -b Dfinder https://github.com/taweiXcms/Bfinder.git

for CMSSW7 (CMSSW_7_5_5_patch4)
check forest setup version in: https://github.com/CmsHI/cmssw/tree/forest_CMSSW_7_5_5_patch4/HeavyIonsAnalysis

git cms-merge-topic -u CmsHI:forest_$CMSSW_VERSION
git clone -b Dfinder https://github.com/taweiXcms/Bfinder.git

