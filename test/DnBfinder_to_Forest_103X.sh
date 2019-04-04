#!/bin/bash

PATHTOTEST=$CMSSW_BASE/src/HeavyIonsAnalysis/JetAnalysis/test
FORESTS=(runForestAOD_pponAA_DATA_103X runForestAOD_pponAA_MIX_103X)
RUNONMC=(False True)
# DIFFPATH=("process.hltanalysisReco *" "process.hltanalysis * process.runAnalyzer *")
INFILES=(
    "file:/eos/cms/store/group/phys_heavyions/wangj/AOD/HIDoubleMuon_PromptReco-v1/CF143D2D-4992-8040-9717-F6ADA30B914C.root"
    "file:/afs/cern.ch/work/w/wangj/public/Hydjet_Pythia8_Psi2SToJpsiPiPi_prompt_Pthat30_TuneCP5_5020GeV_Drum5Ev8/MC_20181231_Psipt0p0_103X_upgrade2018_realistic_HI_v7_RECO/step2_reco_121.root"
)

cc=0
for FOREST in ${FORESTS[@]}
do
##
    cp ${PATHTOTEST}/${FOREST}.py ${PATHTOTEST}/${FOREST}_wDfinder.py

    echo '
#################### D/B finder ################# 
AddCaloMuon = False 
runOnMC = '${RUNONMC[cc]}' ## !!
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
runOnMC = '${RUNONMC[cc]}' ## !!
HIFormat = False
UseGenPlusSim = False
VtxLabel = "offlinePrimaryVerticesRecovery"
TrkLabel = "generalTracks"
from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X
finderMaker_75X(process, AddCaloMuon, runOnMC, HIFormat, UseGenPlusSim, VtxLabel, TrkLabel)

process.Bfinder.MVAMapLabel = cms.InputTag(TrkLabel,"MVAValues")
process.Bfinder.makeBntuple = cms.bool(True)
process.Bfinder.tkPtCut = cms.double(0.7)#before fit
process.Bfinder.jpsiPtCut = cms.double(0.0)#before fit
process.Bfinder.bPtCut = cms.vdouble(5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0)#before fit
process.Bfinder.Bchannel = cms.vint32(1, 0, 0, 1, 1, 1, 1)
process.Bfinder.VtxChiProbCut = cms.vdouble(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.10)
process.Bfinder.svpvDistanceCut = cms.vdouble(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0)
process.Bfinder.doTkPreCut = cms.bool(True)
process.Bfinder.MuonTriggerMatchingPath = cms.vstring(
    "HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1",
    "HLT_HIL1DoubleMuOpen_v1",
    "HLT_HIL1DoubleMu10_v1")
process.Bfinder.MuonTriggerMatchingFilter = cms.vstring(
    "hltL2fDoubleMuOpenL2DR3p5PreFiltered0",
    "hltL1fL1sL1DoubleMuOpenL1Filtered0",
    "hltL1fL1sL1DoubleMu10L1Filtered0")
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
    process.runAnalyzer + \\
    process.hltanalysis + \\
    # process.hltobject + \\
    process.centralityBin + \\
    process.hiEvtAnalyzer  \\
    ) \\
" $ifile
    done

#
    for ifile in ${PATHTOTEST}/${FOREST}_onlyBfinder.py ${PATHTOTEST}/${FOREST}_onlyDfinder.py ${PATHTOTEST}/${FOREST}_wBfinder.py ${PATHTOTEST}/${FOREST}_wDfinder.py
    do
        echo '
###############################
import FWCore.ParameterSet.VarParsing as VarParsing
ivars = VarParsing.VarParsing('"'"'analysis'"'"')

ivars.maxEvents = -1
ivars.outputFile='"'"'HiForestAOD.root'"'"'
ivars.inputFiles='"'${INFILES[cc]}'"'
ivars.parseArguments() # get and parse the command line arguments

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(ivars.inputFiles)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(ivars.maxEvents)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(ivars.outputFile))
' >> $ifile
    done
##
    cc=$((cc+1))
done