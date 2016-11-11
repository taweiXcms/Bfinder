CMSSW_8_0_22
=====

0.check forest setup version in: https://github.com/CmsHI/cmssw/tree/forest_CMSSW_7_5_5_patch4/HeavyIonsAnalysis
`git cms-merge-topic -u CmsHI:forest_$CMSSW_VERSION`
`git clone -b test80forpPb https://github.com/taweiXcms/Bfinder.git`


To add B/D finder, paste the following block:
-----

```
AddCaloMuon = False
runOnMC = False
HIFormat = False
UseGenPlusSim = False
VtxLabel = "offlinePrimaryVerticesWithBS"
TrkLabel = "generalTracks"
from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X
finderMaker_75X(process, AddCaloMuon, runOnMC, HIFormat, UseGenPlusSim, VtxLabel, TrkLabel)
process.Dfinder.alphaCut = cms.vdouble(0.2, 0.2, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0)
process.Dfinder.svpvDistanceCut_lowptD = cms.vdouble(2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
process.Dfinder.dPtCut = cms.vdouble(2.0, 2.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0)
process.Bfinder.tkPtCut = cms.double(0.5)
process.Bfinder.bPtCut = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
```

And add to your cms.Path()
-----

```
process.DfinderSequence+process.BfinderSequence
```
