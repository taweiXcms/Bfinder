D meson reconstruction is added to this version in CMSSW_5_3_20
git checkout Dfinder

note:
D vertex probability: 0.05 (in B we use 0.01)

=================
Few actions now needed, please do:
mv DataFormats  HeavyIonsAnalysis RecoHI ../
Before scram b (compile).

Comments==========
In:
DataFormats/PatCandidates/src/Muon.cc
DataFormats/MuonReco/src/MuonCocktails.cc
Comment out the error log which seems to be delaying the program

Add:
HeavyIonsAnalysis/Configuration/python/collisionEventSelection_cff.py
For HI event selection
