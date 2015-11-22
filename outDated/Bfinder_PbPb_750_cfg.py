import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
ivars = VarParsing.VarParsing('analysis')
#ivars.inputFiles='file:/mnt/hadoop/cms/store/user/twang/HIDiMuon/RECO_HIDiMuon_L2DoubleMu3Skim_v10_JpsiFilter_v1_CMSSW740pre8_20150428/3c3450dda05abb66de621932774972fa/hiRecoData_RAW2DIGI_L1Reco_RECO_filter_975_1_PTa.root'
#ivars.inputFiles='file:/mnt/hadoop/cms/store/user/twang/Pyquen_CMSSW742_Unquenched_PbPb_2760GeV_GEN_SIM_PU_BuKp_20150526_100kevt/Pyquen_CMSSW742_Unquenched_PbPb_2760GeV_step3_BuKp_20150526_100kevt/27ff3fcdfd0b68d12bfbb80768287940/step3_RAW2DIGI_L1Reco_RECO_PU_90_1_Ole.root'
ivars.inputFiles='file:/mnt/hadoop/cms/store/user/richard/MBHydjet5020/Hydjet_Quenched_MinBias_5020GeV/HydjetMB5020_750_75X_mcRun2_HeavyIon_v1_RealisticHICollisions2011_STARTHI50_mc_RECOSIM_v3/150729_144407/0000/step3_98.root'

ivars.outputFile='Bfinder_PbPb_all.root'
# get and parse the command line arguments
ivars.parseArguments()

### Add Calo muons
AddCaloMuon = False

### Run on MC?
runOnMC = True
#runOnMC = False

### Switching between "hiGenParticles"(pPb MC) and "genParticles" (pp MC)
#HIFormat = True
HIFormat = False

### Include SIM tracks for matching?
UseGenPlusSim = False

### Add centrality filter
CentralityFilter = False

process = cms.Process("demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
### Set TransientTrackBuilder 
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
### Set Geometry/GlobalTag/BField
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

### output module
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('keep *',
    )
)

### Set output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(ivars.outputFile)
)

### Set maxEvents
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#### HI infomation
from GeneratorInterface.HiGenCommon.HeavyIon_cff import *
process.load('GeneratorInterface.HiGenCommon.HeavyIon_cff')

#### SetUp Evt Info (centrality)
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v1', '')

process.GlobalTag.toGet.extend([
 cms.PSet(record = cms.string("HeavyIonRcd"),
tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5_v740x01_mc"),
connect = cms.string("frontier://FrontierProd/CMS_COND_31X_PHYSICSTOOLS"),
label = cms.untracked.string("HFtowersHydjetDrum5")
 ),
])

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("HydjetDrum5")

process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')

### Set global tag
if runOnMC:
    #process.GlobalTag.globaltag = cms.string( 'START53_V7F::All' )  #Summer12_DR53X
    #process.GlobalTag.globaltag = cms.string( 'STARTHI53_V26::All' ) 
    #process.GlobalTag.globaltag = cms.string( 'START52_V5::All' ) 
    #process.GlobalTag.globaltag = cms.string( 'START52_V7::All' )
    #process.GlobalTag.globaltag = cms.string( 'START53_V17::All' )
    #process.GlobalTag.globaltag = cms.string( 'STARTHI53_LV1::All' ) ##PbPb
    #process.GlobalTag.globaltag = cms.string( 'START53_V27::All' ) ##pPb
    #process.GlobalTag.globaltag = cms.string( 'MCHI1_74_V4::All' ) ##PbPb for 7_4_0_pre8
    #process.GlobalTag.globaltag = cms.string( 'MCHI1_74_V6::All' ) ##PbPb for 7_4_2
    process.GlobalTag.globaltag = cms.string( '75X_mcRun2_HeavyIon_v1' ) ##PbPb for 7_5_0
else:
    #process.GlobalTag.globaltag = cms.string( 'FT_53_V6_AN2::All' ) #for 2012AB
    #process.GlobalTag.globaltag = cms.string( 'FT_53_V10_AN2::All' )#for 2012C
    #process.GlobalTag.globaltag = cms.string( 'FT_P_V42_AN2::All' ) #for 2012D
    #process.GlobalTag.globaltag = cms.string( 'GR_R_53_LV6::All' ) ##PbPb
    #process.GlobalTag.globaltag = cms.string( 'GR_P_V43D::All' ) ##pp
    #process.GlobalTag.globaltag = cms.string( 'GR_P_V43F::All' ) ##pPb: /PAMuon/HIRun2013-28Sep2013-v1/RECO
    #process.GlobalTag.globaltag = cms.string( 'GR_P_V43D::All' ) ##pPb: /PAMuon/HIRun2013-PromptReco-v1/RECO
    #process.GlobalTag.globaltag = cms.string( 'GR_R_74_V8A::All' ) ##CMSSW_7_4_0_pre8 PbPb
    #process.GlobalTag.globaltag = cms.string( 'GR_R_74_V12A::All' ) ##CMSSW_7_4_2 PbPb
    process.GlobalTag.globaltag = cms.string( '75X_dataRun2_v2' ) ##CMSSW_7_5_0 PbPb

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

### Set basic filter
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
	#vertexCollection = cms.InputTag('offlinePrimaryVertices'),
	#vertexCollection = cms.InputTag('offlinePrimaryVerticesWithBS'),
	vertexCollection = cms.InputTag('hiSelectedVertex'),
	minimumNDOF = cms.uint32(4) ,
	maxAbsZ = cms.double(24),	
	maxd0 = cms.double(2)	
)

process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
	debugOn = cms.untracked.bool(False),
	numtrack = cms.untracked.uint32(10),
	thresh = cms.untracked.double(0.25)
)

# Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

#process.filter = cms.Sequence(process.primaryVertexFilter+process.noscraping)
#process.filter = cms.Sequence(process.noscraping)
#process.filter = cms.Sequence(process.PAcollisionEventSelection)
process.filter = cms.Sequence(process.collisionEventSelection)

## Add centrality filter
if CentralityFilter:
	process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")
	#process.cenfilter = process.centralityFilter.clone(selectedBins = [0,1,2,3,4])
	process.cenfilter = process.centralityFilter.clone(selectedBins = range(59,201))
	process.filter = cms.Sequence(process.centralityBin*process.cenfilter*process.collisionEventSelection)

##Producing Gen list with SIM particles
process.genParticlePlusGEANT = cms.EDProducer("GenPlusSimParticleProducer",
        src           = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
        setStatus     = cms.int32(8),             # set status = 8 for GEANT GPs
        filter        = cms.vstring("pt > 0.0"),  # just for testing (optional)
	    genParticles   = cms.InputTag("genParticles") # original genParticle list
)

### Setup Pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")
###### Needed in CMSSW7
process.particleFlowPtrs.src = "particleFlowTmp"
process.pfPileUpIsoPFBRECO.Vertices = cms.InputTag("hiSelectedVertex")
process.pfPileUpPFBRECO.Vertices = cms.InputTag("hiSelectedVertex")
###### Needed in CMSSW7

if HIFormat:
	process.muonMatch.matched = cms.InputTag("hiGenParticles")
	process.genParticlePlusGEANT.genParticles = cms.InputTag("hiGenParticles")

##Using GEN plus SIM list for matching
if UseGenPlusSim:
	process.muonMatch.matched = cms.InputTag("genParticlePlusGEANT")

## TrackCand
from PhysicsTools.PatAlgos.tools.trackTools import *
if runOnMC:
    makeTrackCandidates(process,              # patAODTrackCands
        label='TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
        #tracks=cms.InputTag('generalTracks'), # input track collection
        tracks=cms.InputTag('hiGeneralTracks'), # input track collection
    	particleType='pi+',                   # particle type (for assigning a mass)
        preselection='pt > 0.3',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
        selection='pt > 0.3',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
    	isolation={},                         # Isolations to use ('source':deltaR; set to {} for None)
       	isoDeposits=[],
        mcAs='muon'                           # Replicate MC match as the one used for Muons
    );                                        # you can specify more than one collection for this
    ### MC+mcAs+Match/pat_label options
    #process.patTrackCandsMCMatch.matched = cms.InputTag("hiGenParticles")
    process.patTrackCandsMCMatch.resolveByMatchQuality = cms.bool(True)
    process.patTrackCandsMCMatch.resolveAmbiguities = cms.bool(True)
    process.patTrackCandsMCMatch.checkCharge = cms.bool(True)
    process.patTrackCandsMCMatch.maxDPtRel = cms.double(0.5)
    process.patTrackCandsMCMatch.maxDeltaR = cms.double(0.7)
    process.patTrackCandsMCMatch.mcPdgId = cms.vint32(111, 211, 311, 321)
    process.patTrackCandsMCMatch.mcStatus = cms.vint32(1)
    l1cands = getattr(process, 'patTrackCands')
    l1cands.addGenMatch = True

else :
    makeTrackCandidates(process,              # patAODTrackCands
        label='TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
        #tracks=cms.InputTag('generalTracks'), # input track collection
        tracks=cms.InputTag('hiGeneralTracks'), # input track collection
        particleType='pi+',                   # particle type (for assigning a mass)
        preselection='pt > 0.3',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
        selection='pt > 0.3',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
        isolation={},                         # Isolations to use ('source':deltaR; set to {} for None)
        isoDeposits=[],
        mcAs=None                             # Replicate MC match as the one used for Muons
    );                                        # you can specify more than one collection for this
    l1cands = getattr(process, 'patTrackCands')
    l1cands.addGenMatch = False
process.TrackCandSequence = cms.Sequence(process.patAODTrackCandsUnfiltered*process.patAODTrackCands*process.patTrackCandsMCMatch*process.patTrackCands*process.selectedPatTrackCands)

## patMuonsWithTrigger
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
process.patMuonsWithTriggerSequence = cms.Sequence(process.pfParticleSelectionForIsoSequence*process.muonPFIsolationPATSequence*process.patMuonsWithTriggerSequence)
if runOnMC:
	addMCinfo(process)
	process.muonMatch.resolveByMatchQuality = True
changeTriggerProcessName(process, "HLT")
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
###Criterias from Hyunchul's 
process.muonL1Info.maxDeltaR = 0.3
process.muonL1Info.fallbackToME1 = True
process.muonMatchHLTL1.maxDeltaR = 0.3
process.muonMatchHLTL1.fallbackToME1 = True
process.muonMatchHLTL2.maxDeltaR = 0.3
process.muonMatchHLTL2.maxDPtRel = 10.0
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0

# Merge muons, calomuons in a single collection for T&P
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
process.mergedMuons = cms.EDProducer("CaloMuonMerger",                                                                                                                                                  
    muons     = cms.InputTag("muons"),
    mergeCaloMuons = cms.bool(True),  ### NEEDED TO RUN ON AOD
    caloMuons = cms.InputTag("calomuons"),
    minCaloCompatibility = cms.double(0.6),
    mergeTracks = cms.bool(False),
    #tracks = cms.InputTag("generalTracks"),
    tracks = cms.InputTag("hiGeneralTracks"),
)
if AddCaloMuon:
    #changeRecoMuonInput(process, "mergedMuons")#Add calo muon to the collection
    #process.patMuons.muonSource = cms.InputTag("mergedMuons")#Need to use the same collection as they are internally entengled
    #process.patMuons.embedCaloMETMuonCorrs = cms.bool(False)
    #process.patMuons.embedTcMETMuonCorrs   = cms.bool(False)

    #Or we change the muonMatch source of our patMuonsWithoutTrigger
    process.patMuonsWithoutTrigger.muonSource = cms.InputTag("mergedMuons")
    process.patMuonsWithoutTriggerMatch = PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi.muonMatch.clone( src = cms.InputTag("mergedMuons"))
    if runOnMC:
        process.patMuonsWithTriggerSequence.replace(process.patMuonsWithoutTrigger, process.patMuonsWithoutTriggerMatch + process.patMuonsWithoutTrigger)
        process.patMuonsWithoutTrigger.genParticleMatch = 'patMuonsWithoutTriggerMatch'
    process.patMuonsWithTriggerSequence = cms.Sequence(process.mergedMuons*process.patMuonsWithTriggerSequence)

### Set Bfinder option
process.demo = cms.EDAnalyzer('Bfinder',
	Bchannel 		= cms.vint32(
		0,#RECONSTRUCTION: J/psi + K
		0,#RECONSTRUCTION: J/psi + Pi
		0,#RECONSTRUCTION: J/psi + Ks 
		0,#RECONSTRUCTION: J/psi + K* (K+, Pi-)
		0,#RECONSTRUCTION: J/psi + K* (K-, Pi+)
		0,#RECONSTRUCTION: J/psi + phi
		0,#RECONSTRUCTION: J/psi + pi pi <= psi', X(3872), Bs->J/psi f0
	),
    #MuonTriggerMatchingPath = cms.vstring("HLT_PAMu3_v*"),
    MuonTriggerMatchingPath = cms.vstring("HLT_HIL2DoubleMu3_v*"),
    #MuonTriggerMatchingPath = cms.vstring("HLT_PAMu3_v*", "HLT_PAMu7_v*", "HLT_PAMu12_v*"),
	HLTLabel        = cms.InputTag('TriggerResults::HLT'),
    GenLabel        = cms.InputTag('genParticles'),
	MuonLabel       = cms.InputTag('patMuonsWithTrigger'),
	TrackLabel      = cms.InputTag('patTrackCands'),
    PUInfoLabel     = cms.InputTag("addPileupInfo"),
    BSLabel     = cms.InputTag("offlineBeamSpot"),
    #PVLabel     = cms.InputTag("offlinePrimaryVerticesWithBS"),
    PVLabel     = cms.InputTag("hiSelectedVertex"),
    tkPtCut = cms.double(1.0),#before fit
    tkEtaCut = cms.double(999.0),#before fit
    jpsiPtCut = cms.double(3.0),#before fit
    bPtCut = cms.double(5.0),#before fit
    bEtaCut = cms.double(2.4),#before fit, not used currently
	VtxChiProbCut = cms.double(0.01),
    RunOnMC = cms.bool(False),
    doTkPreCut = cms.bool(False),
    doMuPreCut = cms.bool(True)
)
### Set Dfinder option
process.Dfinder = cms.EDAnalyzer('Dfinder',
	Dchannel 		= cms.vint32(
		0,#RECONSTRUCTION: K+pi-
		0,#RECONSTRUCTION: K-pi+
		0,#RECONSTRUCTION: K-pi+pi+
		0,#RECONSTRUCTION: K+pi-pi-
		0,#RECONSTRUCTION: K-pi-pi+pi+
		0,#RECONSTRUCTION: K+pi+pi-pi-
		0,#RECONSTRUCTION: K+K-(Phi)pi+
		0,#RECONSTRUCTION: K+K-(Phi)pi-
	),
	HLTLabel        = cms.InputTag('TriggerResults::HLT'),
    GenLabel        = cms.InputTag('genParticles'),
	TrackLabel      = cms.InputTag('patTrackCands'),
    PUInfoLabel     = cms.InputTag("addPileupInfo"),
    BSLabel     = cms.InputTag("offlineBeamSpot"),
    #PVLabel     = cms.InputTag("offlinePrimaryVerticesWithBS"),
    PVLabel     = cms.InputTag("hiSelectedVertex"),
    tkPtCut = cms.double(1.),#before fit
    tkEtaCut = cms.double(1.1),#before fit
    dPtCut = cms.double(3.0),#before fit
    dEtaCut = cms.double(1.5),#before fit, not used currently
	VtxChiProbCut = cms.double(0.05),
    svpvDistanceCut = cms.double(2.0),
    MaxDocaCut = cms.double(999.),
    alphaCut = cms.double(999.),
    RunOnMC = cms.bool(False),
    doTkPreCut = cms.bool(True),
)
if runOnMC:
    process.demo.RunOnMC = cms.bool(True)
    process.Dfinder.RunOnMC = cms.bool(True)
if HIFormat:
	process.demo.GenLabel = cms.InputTag('hiGenParticles')
	process.Dfinder.GenLabel = cms.InputTag('hiGenParticles')
if UseGenPlusSim:
	process.demo.GenLabel = cms.InputTag('genParticlePlusGEANT')
	process.Dfinder.GenLabel = cms.InputTag('genParticlePlusGEANT')

### SetUp HLT info
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
process.hltanalysis.dummyBranches = cms.untracked.vstring()
process.hltanalysis.OfflinePrimaryVertices0 = cms.InputTag("hiSelectedVertex")
#if HIFormat:
	#process.hltanalysis.mctruth = cms.InputTag("hiGenParticles")# Will cause segmentation violation
	#process.hltanalysis.HLTProcessName = cms.string("HISIGNAL")
	#process.hltanalysis.hltresults = cms.InputTag("TriggerResults","","HISIGNAL")
	#process.hltanalysis.l1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap::HISIGNAL")
process.hltAna = cms.Path(process.filter*process.hltanalysis)

### Run the hiEvtAnalyzer
process.evtAna = cms.Path(process.filter*process.hiEvtAnalyzer)
if runOnMC:
    process.hiEvtAnalyzer.doMC = cms.bool(True)
    process.evtAna = cms.Path(process.filter*process.heavyIon*process.hiEvtAnalyzer)

if runOnMC and UseGenPlusSim:
	process.patMuonsWithTriggerSequence *= process.genParticlePlusGEANT

process.p = cms.Path(	
    process.filter*process.patMuonsWithTriggerSequence*process.TrackCandSequence*process.demo*process.Dfinder
)

#process.e = cms.EndPath(process.out)
process.schedule = cms.Schedule(
	process.p
	,process.hltAna
	,process.evtAna
#    ,process.e
)
