import FWCore.ParameterSet.Config as cms
process = cms.Process('HiForest')
import FWCore.ParameterSet.VarParsing as VarParsing
ivars = VarParsing.VarParsing('analysis')
#ivars.inputFiles='file:/data/twang/MC_samples/MinBias_TuneCUETP8M1_5p02TeV-pythia8/MinBias_TuneCUETP8M1_5p02TeV_pythia8_pp502Fall15_MCRUN2_71_V1_v1_AOD_CMSSW_7_5_4_20151113/step3_RAW2DIGI_L1Reco_RECO_993_1_q1f.root'
#ivars.inputFiles='file:/data/twang/Data_samples/Run2015E/DoubleMu/AOD/PromptReco-v1/000/262/235/00000/0E9E6AA6-F394-E511-B74B-02163E01474F.root'#AOD DoubleMu
#ivars.inputFiles='file:/data/twang/Data_samples/Run2015E/DoubleMu/RECO/PromptReco-v1/000/262/163/00000/14A3BF17-D591-E511-868F-02163E014117.root'#RECO DoubleMu
#ivars.inputFiles='file:/data/twang/Data_samples/Run2015E/HeavyFlavor/AOD/PromptReco-v1/000/262/273/00000/06090E4E-0C97-E511-856C-02163E0142D2.root'#AOD HeavyFlavor
#ivars.inputFiles='file:/data/twang/Data_samples/Run2015E/MinimumBias1/AOD/PromptReco-v1/000/262/274/00000/0E443E25-3C9A-E511-9263-02163E013626.root'#AOD MB1
ivars.outputFile='finder_pp.root'
ivars.parseArguments()# get and parse the command line arguments

### Custom options
########## MUST CUSTOMIZE THE FOLLOWING THREE ##########
### pp B/Dfinder recommended setting, choose only one from them or set all to false and made your own setting
ppBdefault = True
ppDHFdefault = False
ppDMBdefault = False

### Run on MC?
runOnMC = False

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
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

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
    globalTag = 'auto:run2_mc'
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
	process.evtAna = cms.Path(process.heavyIon*process.hiEvtAnalyzer)

### Set basic filter
process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True), # otherwise it won't filter the events
)
process.pPAprimaryVertexFilter = cms.Path(process.PAprimaryVertexFilter)

process.NoScraping = cms.EDFilter("FilterOutScraping",
 applyfilter = cms.untracked.bool(True),
 debugOn = cms.untracked.bool(False),
 numtrack = cms.untracked.uint32(10),
 thresh = cms.untracked.double(0.25)
)
process.pBeamScrapingFilter=cms.Path(process.NoScraping)

# Common offline event selection
process.load('HeavyIonsAnalysis.Configuration.hfCoincFilter_cff')

process.PAcollisionEventSelection = cms.Sequence(process.hfCoincFilter *
                                         process.PAprimaryVertexFilter *
                                         process.NoScraping 
                                         )
if not RunOnAOD:
	process.phfCoincFilter = cms.Path(process.hfCoincFilter )
	process.phfCoincFilter3 = cms.Path(process.hfCoincFilter3 )
	process.PAcollisionEventSelection = cms.Path(process.PAcollisionEventSelection)

### Run HLT info sequence
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
process.hltanalysis.dummyBranches = cms.untracked.vstring()
process.hltanalysis.OfflinePrimaryVertices0 = cms.InputTag(VtxLabel)
#if HIFormat:
    #process.hltanalysis.mctruth = cms.InputTag("hiGenParticles")# Will cause segmentation violation
    #process.hltanalysis.HLTProcessName = cms.string("HISIGNAL")
    #process.hltanalysis.hltresults = cms.InputTag("TriggerResults","","HISIGNAL")
    #process.hltanalysis.l1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap::HISIGNAL")
process.hltAna = cms.Path(process.hltanalysis)

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
if ppBdefault and not ppDMBdefault and not ppDHFdefault:
    process.Bfinder.tkPtCut = cms.double(0.5)#before fit
    process.Bfinder.jpsiPtCut = cms.double(0.0)#before fit
    process.Bfinder.bPtCut = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)#before fit
    process.Bfinder.Bchannel = cms.vint32(1, 0, 0, 1, 1, 1, 1)
    process.p = cms.Path(process.BfinderSequence)
## pp Dfinder setting on HeavyFlavor
if ppDHFdefault and not ppDMBdefault and not ppBdefault:
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
if ppDMBdefault and not ppDHFdefault and not ppBdefault:
    process.Dfinder.tkPtCut = cms.double(0.5)#before fit
    process.Dfinder.dPtCut = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)#before fit
    process.Dfinder.dCutSeparating_PtVal = cms.vdouble(5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.)
    process.Dfinder.tktkRes_svpvDistanceCut_lowptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5)
    process.Dfinder.tktkRes_svpvDistanceCut_highptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5)
    process.Dfinder.svpvDistanceCut_lowptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0.)
    process.Dfinder.svpvDistanceCut_highptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0.)
    process.Dfinder.Dchannel = cms.vint32(1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0)
    process.p = cms.Path(process.DfinderSequence)

process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.pHBHENoiseFilterResultProducer = cms.Path( process.HBHENoiseFilterResultProducer )

process.pAna = cms.EndPath(process.skimanalysis)

### Add centrality filter
if CentralityFilter:
    process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")
    #process.cenfilterClone = process.centralityFilter.clone(selectedBins = [0,1,2,3,4])
    process.cenfilterClone = process.centralityFilter.clone(selectedBins = range(59,201))
    process.filter = cms.Sequence(process.cenfilterClone)
    for path in process.paths:
       getattr(process,path)._seq = process.filter * getattr(process,path)._seq
