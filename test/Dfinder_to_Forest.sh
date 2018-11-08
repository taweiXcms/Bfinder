#!/bin/bash

PATHTOTEST=$CMSSW_BASE/src/HeavyIonsAnalysis/JetAnalysis/test

cp ${PATHTOTEST}/runForestAOD_pponAA_DATA_103X.py ${PATHTOTEST}/runForestAOD_pponAA_DATA_103X_wDfinder.py

echo '
#################### D/B finder ################# 
AddCaloMuon = False 
runOnMC = False 
HIFormat = False 
UseGenPlusSim = False 
VtxLabel = "offlinePrimaryVerticesWithBS" 
TrkLabel = "generalTracks" 
from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X 
finderMaker_75X(process, AddCaloMuon, runOnMC, HIFormat, UseGenPlusSim, VtxLabel, TrkLabel)
process.Dfinder.MVAMapLabel = cms.InputTag(TrkLabel,"MVAValues")
process.Dfinder.makeDntuple = cms.bool(True)
process.Dfinder.tkPtCut = cms.double(1.) # before fit
process.Dfinder.dPtCut = cms.vdouble(2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0) # before fit
process.Dfinder.dCutSeparating_PtVal = cms.vdouble(5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.)
process.Dfinder.tktkRes_svpvDistanceCut_lowptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0.)
process.Dfinder.tktkRes_svpvDistanceCut_highptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0.)
process.Dfinder.svpvDistanceCut_lowptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0., 0., 0., 2.5, 2.5)
process.Dfinder.svpvDistanceCut_highptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0., 0., 0., 2.5, 2.5)
process.Dfinder.Dchannel = cms.vint32(1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1)

process.dfinder = cms.Path(process.DfinderSequence)

' >> ${PATHTOTEST}/runForestAOD_pponAA_DATA_103X_wDfinder.py

cp ${PATHTOTEST}/runForestAOD_pponAA_DATA_103X_wDfinder.py ${PATHTOTEST}/runForestAOD_pponAA_DATA_103X_onlyDfinder.py

echo '
process.ana_step = cms.Path(
    process.HiForest +
    process.hltanalysis +
    process.hltobject +
    process.centralityBin +
    process.hiEvtAnalyzer 
    )
' >> ${PATHTOTEST}/runForestAOD_pponAA_DATA_103X_onlyDfinder.py