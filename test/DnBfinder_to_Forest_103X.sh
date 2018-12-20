#!/bin/bash

PATHTOTEST=$CMSSW_BASE/src/HeavyIonsAnalysis/JetAnalysis/test
FOREST=runForestAOD_pponAA_DATA_103X

##
cp ${PATHTOTEST}/${FOREST}.py ${PATHTOTEST}/${FOREST}_wDfinder.py

echo '
#################### D/B finder ################# 
AddCaloMuon = False 
runOnMC = False 
HIFormat = False 
UseGenPlusSim = False 
VtxLabel = "offlinePrimaryVerticesRecovery" 
TrkLabel = "generalTracks" 
from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X 
finderMaker_75X(process, AddCaloMuon, runOnMC, HIFormat, UseGenPlusSim, VtxLabel, TrkLabel)
process.Dfinder.MVAMapLabel = cms.InputTag(TrkLabel,"MVAValues")
process.Dfinder.makeDntuple = cms.bool(True)
process.Dfinder.tkPtCut = cms.double(1.) # before fit
process.Dfinder.dPtCut = cms.vdouble(2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0) # before fit
process.Dfinder.VtxChiProbCut = cms.vdouble(0.05, 0.05, 0.0, 0.0, 0.0, 0.0, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.05)
process.Dfinder.dCutSeparating_PtVal = cms.vdouble(5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.)
process.Dfinder.tktkRes_svpvDistanceCut_lowptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0.)
process.Dfinder.tktkRes_svpvDistanceCut_highptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0.)
process.Dfinder.svpvDistanceCut_lowptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0., 0., 0., 2.5, 2.5)
process.Dfinder.svpvDistanceCut_highptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0., 0., 0., 2.5, 2.5)
process.Dfinder.Dchannel = cms.vint32(1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1)

process.dfinder = cms.Path(process.DfinderSequence)

' >> ${PATHTOTEST}/${FOREST}_wDfinder.py

##
cp ${PATHTOTEST}/${FOREST}.py ${PATHTOTEST}/${FOREST}_wBfinder.py

echo '
#################### D/B finder #################
AddCaloMuon = False
runOnMC = False
HIFormat = False
UseGenPlusSim = False
VtxLabel = "offlinePrimaryVerticesRecovery"
TrkLabel = "generalTracks"
from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X
finderMaker_75X(process, AddCaloMuon, runOnMC, HIFormat, UseGenPlusSim, VtxLabel, TrkLabel)

process.Bfinder.MVAMapLabel = cms.InputTag(TrkLabel,"MVAValues")
process.Bfinder.makeBntuple = cms.bool(True)
process.Bfinder.tkPtCut = cms.double(1.0)#before fit
process.Bfinder.jpsiPtCut = cms.double(0.0)#before fit
process.Bfinder.bPtCut = cms.vdouble(5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0)#before fit
process.Bfinder.Bchannel = cms.vint32(1, 0, 0, 1, 1, 1, 1)
process.Bfinder.VtxChiProbCut = cms.vdouble(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
process.Bfinder.svpvDistanceCut = cms.vdouble(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
process.Bfinder.doTkPreCut = cms.bool(True)
process.Bfinder.MuonTriggerMatchingPath = cms.vstring(
    "HLT_HIL1DoubleMu0_v1",
    "HLT_HIL1DoubleMu0_part1_v1",
    "HLT_HIL1DoubleMu0_part2_v1",
    "HLT_HIL1DoubleMu0_part3_v1",
    "HLT_HIL1DoubleMu0_2HF_v1",
    "HLT_HIL1DoubleMu0_2HF0_v1",
    "HLT_HIL1DoubleMu10_v1",
    "HLT_HIL2DoubleMu0_NHitQ_v2",
    "HLT_HIL2DoubleMu0_NHitQ_2HF_v1",
    "HLT_HIL2DoubleMu0_NHitQ_2HF0_v1",
    "HLT_HIL1DoubleMu0_2HF_Cent30100_v1",
    "HLT_HIL1DoubleMu0_2HF0_Cent30100_v1",
    "HLT_HIL2DoubleMu0_2HF_Cent30100_NHitQ_v1",
    "HLT_HIL1DoubleMu0_Cent30_v1",
    "HLT_HIL2DoubleMu0_2HF0_Cent30100_NHitQ_v1",
    "HLT_HIL2DoubleMu0_Cent30_NHitQ_v1",
    "HLT_HIL2DoubleMu0_Cent30_OS_NHitQ_v1",
    "HLT_HIL3DoubleMu0_Cent30_v1",
    "HLT_HIL3DoubleMu0_Cent30_OS_m2p5to4p5_v1",
    "HLT_HIL3DoubleMu0_Cent30_OS_m7to14_v1",
    "HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1",
    "HLT_HIL3DoubleMu0_OS_m7to14_v1")
process.Bfinder.MuonTriggerMatchingFilter = cms.vstring("hltHIDoubleMu0L1Filtered",
                                                        "hltHIDoubleMu0MinBiasL1Filtered",
                                                        "hltHIDoubleMu0HFTower0Filtered",
                                                        "hltHIDoubleMu10L1Filtered",
                                                        "hltHIL2DoubleMu0NHitQFiltered",
                                                        "hltHIL2DoubleMu0NHitQ2HFFiltered",
                                                        "hltHIL2DoubleMu0NHitQ2HF0Filtered",
                                                        "hltHIDoubleMu0MinBiasCent30to100L1Filtered",
                                                        "hltHIDoubleMu0HFTower0Cent30to100L1Filtered",
                                                        "hltHIL2DoubleMu02HFcent30100NHitQFiltered",
                                                        "hltHIDoubleMu0MinBiasCent30L1Filtered",
                                                        "hltHIL2DoubleMu02HF0cent30100NHitQFiltered",
                                                        "hltHIL2DoubleMu0cent30NHitQFiltered",
                                                        "hltHIL2DoubleMu0cent30OSNHitQFiltered",
                                                        "hltHIDimuonOpenCentrality30L3Filter",
                                                        "hltHIDimuonOpenCentrality30OSm2p5to4p5L3Filter",
                                                        "hltHIDimuonOpenCentrality30OSm7to14L3Filter",
                                                        "hltHIDimuonOpenOSm2p5to4p5L3Filter",
                                                        "hltHIDimuonOpenOSm7to14L3Filter")
process.p = cms.Path(process.BfinderSequence)

' >> ${PATHTOTEST}/${FOREST}_wBfinder.py

#
cp ${PATHTOTEST}/${FOREST}_wDfinder.py ${PATHTOTEST}/${FOREST}_onlyDfinder.py
cp ${PATHTOTEST}/${FOREST}_wBfinder.py ${PATHTOTEST}/${FOREST}_onlyBfinder.py

for ifile in ${PATHTOTEST}/${FOREST}_onlyBfinder.py ${PATHTOTEST}/${FOREST}_onlyDfinder.py
do
    sed -i "/D\/B finder/i \\
process.ana_step = cms.Path( \\
    process.offlinePrimaryVerticesRecovery + \\
    process.HiForest + \\
    process.hltanalysis + \\
    process.hltobject + \\
    process.centralityBin + \\
    process.hiEvtAnalyzer  \\
    ) \\
" $ifile
done

##