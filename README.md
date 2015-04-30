test branch

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
