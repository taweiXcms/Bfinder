import FWCore.ParameterSet.Config as cms
process = cms.Process('HiForest')
import FWCore.ParameterSet.VarParsing as VarParsing
ivars = VarParsing.VarParsing('analysis')
ivars.inputFiles='file:/afs/cern.ch/work/g/ginnocen/HF2016TriggerDevelop_New/testofDfinder80X/CMSSW_8_0_9/src/Bfinder/test/example/step3_MinBias_pPb_RAW2DIGI_L1Reco_RECO_1.root'
ivars.outputFile='finder_pp.root'
ivars.parseArguments()# get and parse the command line arguments

### Custom options
########## MUST CUSTOMIZE THE FOLLOWING THREE ##########
### pp B/Dfinder recommended setting, choose only one from them or set all to false and made your own setting
ppBdefault = 0
ppDHFdefault = 0
ppDMBdefault = 1
ppD0DstarV2 = 0

### Run on MC?
runOnMC = True

### Use AOD event filter
RunOnAOD = True
########## MUST CUSTOMIZE THE FOLLOWING THREE ##########

### More custom options
### Add Calo muons
AddCaloMuon = False

### Switching between "hiGenParticles"(pPb MC) and "genParticles" (pp MC)
HIFormat = False

### Include SIM tracks for matching?
UseGenPlusSim = False

### Add centrality filter
CentralityFilter = False

### Vertex/Track label 
VtxLabel = "offlinePrimaryVerticesWithBS"
TrkLabel = "generalTracks"

### Set maxEvents
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))

### output module
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('keep *',
    )
)
#process.e = cms.EndPath(process.out)

### Set output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(ivars.outputFile)
)

### PoolSource will be ignored when running crab
process.source = cms.Source("PoolSource",
    skipEvents=cms.untracked.uint32(0),
	fileNames = cms.untracked.vstring(ivars.inputFiles)
)

### Using JSON file
#if not runOnMC:
#    #import PhysicsTools.PythonAnalysis.LumiList as LumiList
#    import FWCore.PythonUtilities.LumiList as LumiList
#    process.source.lumisToProcess = LumiList.LumiList(filename =
#    '/net/hisrv0001/home/tawei/HeavyFlavor_20131030/Bfinder/CMSSW_5_3_20/src/Bfinder/JSON/Cert_181530-183126_HI7TeV_25Oct2012ReReco_Collisions11_JSON_MuonPhys_HF_manualPatch.txt'
#    ).getVLuminosityBlockRange()

### General setups, Geometry/GlobalTag/BField
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

### All relevent GlobalTags
globalTag = ""
#MC
if runOnMC:
    #globalTag = 'START53_V7F::All'#Summer12_DR53X
    #globalTag = 'STARTHI53_V26::All'
    #globalTag = 'START52_V5::All'
    #globalTag = 'START52_V7::All'
    #globalTag = 'START53_V17::All'
    #globalTag = 'STARTHI53_LV1::All'##PbPb
    #globalTag = 'START53_V27::All'##pPb
    #globalTag = 'MCHI1_74_V4::All'##PbPb for 7_4_0_pre8
    #globalTag = 'MCHI1_74_V6::All'##PbPb for 7_4_2
    #globalTag = '75X_mcRun2_HeavyIon_v1'##PbPb for 7_5_0
    #globalTag = '75X_mcRun2_HeavyIon_v4'##PbPb for 7_5_3_patch1
    #globalTag = '75X_mcRun2_asymptotic_v5'##pp for 7_5_3_patch1
    #globalTag = 'auto:run2_mc'
    globalTag = '80X_mcRun2_asymptotic_v12'
	
#Data
else:
    #globalTag = 'FT_53_V6_AN2::All'#for 2012AB
    #globalTag = 'FT_53_V10_AN2::All'#for 2012C
    #globalTag = 'FT_P_V42_AN2::All'#for 2012D
    #globalTag = 'GR_R_53_LV6::All'##PbPb
    #globalTag = 'GR_P_V43D::All'##pp
    #globalTag = 'GR_P_V43F::All'##pPb: /PAMuon/HIRun2013-28Sep2013-v1/RECO
    #globalTag = 'GR_P_V43D::All'##pPb: /PAMuon/HIRun2013-PromptReco-v1/RECO
    #globalTag = 'GR_R_74_V8A::All'##CMSSW_7_4_0_pre8 PbPb
    #globalTag = 'GR_R_74_V12A::All'##CMSSW_7_4_2 PbPb
    #globalTag = '75X_dataRun2_v2'##CMSSW_7_5_0 PbPb
    globalTag = 'auto:run2_data'

process.GlobalTag.globaltag = cms.string(globalTag)
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag, '')

#### HI infomation
from GeneratorInterface.HiGenCommon.HeavyIon_cff import *
process.load('GeneratorInterface.HiGenCommon.HeavyIon_cff')

#### SetUp Evt Info (centrality)
#for pp data create centrality object and bin
process.load("RecoHI.HiCentralityAlgos.pACentrality_cfi")
process.pACentrality.producePixelTracks = cms.bool(False)
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("pACentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
if not RunOnAOD:
	process.centrality_path = cms.Path(process.siPixelRecHits*process.pACentrality*process.centralityBin)

### Run the hiEvtAnalyzer sequence
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
### pp centrality
process.hiEvtAnalyzer.CentralitySrc = cms.InputTag("pACentrality")
process.hiEvtAnalyzer.Vertex = cms.InputTag("offlinePrimaryVertices")
process.hiEvtAnalyzer.doEvtPlane = cms.bool(False)
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
process.evtAna = cms.Path(process.hiEvtAnalyzer)
if RunOnAOD:
	process.hiEvtAnalyzer.doCentrality = cms.bool(False) 
if runOnMC:
	process.hiEvtAnalyzer.doMC = cms.bool(True)
	#process.evtAna = cms.Path(process.heavyIon*process.hiEvtAnalyzer)

### Run HLT info sequence
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
from HeavyIonsAnalysis.EventAnalysis.dummybranches_cff import addHLTdummybranchesForPP
addHLTdummybranchesForPP(process)
process.hltanalysis.OfflinePrimaryVertices0 = cms.InputTag(VtxLabel)
#if HIFormat:
    #process.hltanalysis.mctruth = cms.InputTag("hiGenParticles")# Will cause segmentation violation
    #process.hltanalysis.HLTProcessName = cms.string("HISIGNAL")
    #process.hltanalysis.hltresults = cms.InputTag("TriggerResults","","HISIGNAL")
    #process.hltanalysis.l1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap::HISIGNAL")
process.hltAna = cms.Path(process.hltanalysis)

### Set basic filter
process.load('HeavyIonsAnalysis.JetAnalysis.EventSelection_cff')
process.pHBHENoiseFilterResultProducer = cms.Path( process.HBHENoiseFilterResultProducer )
process.HBHENoiseFilterResult = cms.Path(process.fHBHENoiseFilterResult)
process.HBHENoiseFilterResultRun1 = cms.Path(process.fHBHENoiseFilterResultRun1)
process.HBHENoiseFilterResultRun2Loose = cms.Path(process.fHBHENoiseFilterResultRun2Loose)
process.HBHENoiseFilterResultRun2Tight = cms.Path(process.fHBHENoiseFilterResultRun2Tight)
process.HBHEIsoNoiseFilterResult = cms.Path(process.fHBHEIsoNoiseFilterResult)

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True), # otherwise it won't filter the events
)

process.NoScraping = cms.EDFilter("FilterOutScraping",
 applyfilter = cms.untracked.bool(True),
 debugOn = cms.untracked.bool(False),
 numtrack = cms.untracked.uint32(10),
 thresh = cms.untracked.double(0.25)
)

process.pPAprimaryVertexFilter = cms.Path(process.PAprimaryVertexFilter)
process.pBeamScrapingFilter=cms.Path(process.NoScraping)

process.load("HeavyIonsAnalysis.VertexAnalysis.PAPileUpVertexFilter_cff")

process.pVertexFilterCutG = cms.Path(process.pileupVertexFilterCutG)
process.pVertexFilterCutGloose = cms.Path(process.pileupVertexFilterCutGloose)
process.pVertexFilterCutGtight = cms.Path(process.pileupVertexFilterCutGtight)
process.pVertexFilterCutGplus = cms.Path(process.pileupVertexFilterCutGplus)
process.pVertexFilterCutE = cms.Path(process.pileupVertexFilterCutE)
process.pVertexFilterCutEandG = cms.Path(process.pileupVertexFilterCutEandG)

process.pAna = cms.EndPath(process.skimanalysis)

### finder building block
from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X
finderMaker_75X(process, AddCaloMuon, runOnMC, HIFormat, UseGenPlusSim, VtxLabel, TrkLabel)
process.p = cms.Path(process.finderSequence)

process.Bfinder.Bchannel = cms.vint32(
    0,#RECONSTRUCTION: J/psi + K
    0,#RECONSTRUCTION: J/psi + Pi
    0,#RECONSTRUCTION: J/psi + Ks
    0,#RECONSTRUCTION: J/psi + K* (K+, Pi-)
    0,#RECONSTRUCTION: J/psi + K* (K-, Pi+)
    0,#RECONSTRUCTION: J/psi + phi
    0,#RECONSTRUCTION: J/psi + pi pi <= psi', X(3872), Bs->J/psi f0
)
process.Dfinder.Dchannel = cms.vint32(
    0,#RECONSTRUCTION: K+pi- : D0bar
    0,#RECONSTRUCTION: K-pi+ : D0
    0,#RECONSTRUCTION: K-pi+pi+ : D+
    0,#RECONSTRUCTION: K+pi-pi- : D-
    0,#RECONSTRUCTION: K-pi-pi+pi+ : D0
    0,#RECONSTRUCTION: K+pi+pi-pi- : D0bar
    0,#RECONSTRUCTION: K+K-(Phi)pi+ : Ds+
    0,#RECONSTRUCTION: K+K-(Phi)pi- : Ds-
    0,#RECONSTRUCTION: D0(K-pi+)pi+ : D+*
    0,#RECONSTRUCTION: D0bar(K+pi-)pi- : D-*
    0,#RECONSTRUCTION: D0(K-pi-pi+pi+)pi+ : D+*
    0,#RECONSTRUCTION: D0bar(K+pi+pi-pi-)pi- : D-*
)
## pp Bfinder setting on DoubleMu
if ppBdefault and not ppDMBdefault and not ppDHFdefault and not ppD0DstarV2:
    process.Bfinder.tkPtCut = cms.double(0.5)#before fit
    process.Bfinder.jpsiPtCut = cms.double(0.0)#before fit
    process.Bfinder.bPtCut = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)#before fit
    process.Bfinder.Bchannel = cms.vint32(1, 0, 0, 0, 0, 0, 0)
    process.Bfinder.doTkPreCut = cms.bool(False)
    process.Bfinder.MuonTriggerMatchingPath = cms.vstring("HLT_HIL1DoubleMu0_v1",
                                                          "HLT_HIL1DoubleMu10_v1",
                                                          "HLT_HIL2DoubleMu0_NHitQ_v1",
                                                          "HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1",
                                                          "HLT_HIL3DoubleMu0_OS_m7to14_v1")
    process.Bfinder.MuonTriggerMatchingFilter = cms.vstring("hltHIDoubleMu0L1Filtered",
                                                            "hltHIDoubleMu10MinBiasL1Filtered",
                                                            "hltHIL2DoubleMu0NHitQFiltered",
                                                            "hltHIDimuonOpenOSm2p5to4p5L3Filter",
                                                            "hltHIDimuonOpenOSm7to14L3Filter")
    process.Bfinder.makeBntuple = cms.bool(False)
    process.p = cms.Path(process.BfinderSequence)
## pp Dfinder setting on HeavyFlavor
if ppDHFdefault and not ppDMBdefault and not ppBdefault and not ppD0DstarV2:
    process.Dfinder.tkPtCut = cms.double(1.)#before fit
    process.Dfinder.dPtCut = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)#before fit
    process.Dfinder.dCutSeparating_PtVal = cms.vdouble(5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.)
    process.Dfinder.tktkRes_svpvDistanceCut_lowptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5)
    process.Dfinder.tktkRes_svpvDistanceCut_highptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5)
    process.Dfinder.svpvDistanceCut_lowptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0.)
    process.Dfinder.svpvDistanceCut_highptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0.)
    process.Dfinder.Dchannel = cms.vint32(1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
    process.p = cms.Path(process.DfinderSequence)
##pp Dfinder setting on MB
if ppDMBdefault and not ppDHFdefault and not ppBdefault and not ppD0DstarV2:
    process.Dfinder.tkPtCut = cms.double(1.)#before fit
    process.Dfinder.dPtCut = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)#before fit
    process.Dfinder.dCutSeparating_PtVal = cms.vdouble(5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.)
    process.Dfinder.tktkRes_svpvDistanceCut_lowptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5)
    process.Dfinder.tktkRes_svpvDistanceCut_highptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5)
    process.Dfinder.svpvDistanceCut_lowptD = cms.vdouble(4.0, 4.0, 2.5, 2.5, 2.5, 2.5, 4.0, 4.0, 0., 0., 0., 0.)
    process.Dfinder.svpvDistanceCut_highptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0.)
    process.Dfinder.Dchannel = cms.vint32(1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0)
    process.p = cms.Path(process.DfinderSequence)
## default cut version 2
if ppD0DstarV2 and not ppDHFdefault and not ppDMBdefault and not ppBdefault:
    process.Dfinder.tkPtCut = cms.double(0.5)#before fit
    process.Dfinder.dPtCut = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)#before fit
    process.Dfinder.dCutSeparating_PtVal = cms.vdouble(8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8.)
    process.Dfinder.tktkRes_svpvDistanceCut_lowptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 3.0, 3.0, 3.0, 3.0 )
    process.Dfinder.tktkRes_svpvDistanceCut_highptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 1.5, 1.5, 1.5, 1.5)
    process.Dfinder.svpvDistanceCut_lowptD = cms.vdouble(3.0, 3.0, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0.)
    process.Dfinder.svpvDistanceCut_highptD = cms.vdouble(1.5, 1.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0.)
    process.Dfinder.Dchannel = cms.vint32(1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
    process.p = cms.Path(process.DfinderSequence)

### Add centrality filter
if CentralityFilter:
    process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")
    #process.cenfilterClone = process.centralityFilter.clone(selectedBins = [0,1,2,3,4])
    process.cenfilterClone = process.centralityFilter.clone(selectedBins = range(59,201))
    process.filter = cms.Sequence(process.cenfilterClone)
    for path in process.paths:
       getattr(process,path)._seq = process.filter * getattr(process,path)._seq
       
