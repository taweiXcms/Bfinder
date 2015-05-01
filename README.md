This is a special temparary version for CMSSW741
Reminder: git checkout temp_v1_CMSSW_7_4_1

Before compile:
mv HeavyIonsAnalysis ../


Update log:
Remove DataFormats, RecoHI, HeavyIonsAnalysis/EventAnalysis, EventAnalysis/src/HiEvtAnalyzer.cc, EventAnalysis/python/hievtanalyzer_data_cfi.py 
comment removeAllPATObjectsBut from Bfinder_cfg.py

Bfinder_PbPb_cfg.py
Added several ad-hoc fixs.
