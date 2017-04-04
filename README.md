To setup Bfinder
=====

Recommended using `CMSSW_7_5_8_patch3` for 2015 5TeV pp and PbPb data and `CMSSW_8_0_22` for 2016 pPb data

Check forest version in: https://github.com/CmsHI/cmssw/tree/forest_CMSSW_7_5_8_patch3/HeavyIonsAnalysis

```
git cms-merge-topic -u CmsHI:forest_$CMSSW_VERSION
git clone -b Dfinder https://github.com/taweiXcms/Bfinder.git
```

To add B/D finder to forest, paste the following block:
=====

pp
-----

```python
#################
### D/B finder
#################
AddCaloMuon = False
runOnMC = False
HIFormat = False
UseGenPlusSim = False
VtxLabel = "offlinePrimaryVerticesWithBS"
TrkLabel = "generalTracks"
from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X
finderMaker_75X(process, AddCaloMuon, runOnMC, HIFormat, UseGenPlusSim, VtxLabel, TrkLabel)
```

pPb
-----

```python
#################
### D/B finder
#################
AddCaloMuon = False
runOnMC = False
HIFormat = False
UseGenPlusSim = False
VtxLabel = "offlinePrimaryVerticesWithBS"
TrkLabel = "generalTracks"
from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X
finderMaker_75X(process, AddCaloMuon, runOnMC, HIFormat, UseGenPlusSim, VtxLabel, TrkLabel)
### MVA label changed in pPb data CMSSW8XX
process.Bfinder.MVAMapLabel = cms.InputTag(TrkLabel,"MVAValues")
process.Dfinder.MVAMapLabel = cms.InputTag(TrkLabel,"MVAValues")
```

PbPb
-----

```python
#################
### D/B finder
#################
AddCaloMuon = False
runOnMC = False
HIFormat = False
UseGenPlusSim = False
from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X
finderMaker_75X(process, AddCaloMuon, runOnMC, HIFormat, UseGenPlusSim)
```

And add to your cms.Path()
-----

```python
process.DfinderSequence+process.BfinderSequence
```

And comment the following line in
-----
`Bfinder/finderMaker/python/finderMaker_75X_cff.py`

```python
process.patTrigger.collections.remove("hltL3MuonCandidates")
```

If running on MC
-----

Remember to set
```python
runOnMC = True
```

Note:
-----

By default it will only run B to Jpsi Pi and D to K Pi channel.

Find more customizations in `Bfinder/finderMaker/python/finderMaker_75X_cff.py`

and add lines like

```python
process.Dfinder.Dchannel = cms.vint32(1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0)
```

to customize your own selection values and channels
