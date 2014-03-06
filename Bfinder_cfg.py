# 2014Feb05   twang   fix HI MC HLT process name
import FWCore.ParameterSet.Config as cms

### Run on MC?
runOnMC = True

### HI label?
HIFormat = True
#HIFormat = False

### Include SIM tracks for matching?
UseGenPlusSim = False

process = cms.Process("demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

### Set TransientTrackBuilder 
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

### Set Geometry/GlobalTag/BField
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

### keep only Pat:: part 
from PhysicsTools.PatAlgos.patEventContent_cff import *
### output module
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *',
    )
)

### Set maxEvents
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

### Set global tag
if runOnMC:
    #process.GlobalTag.globaltag = cms.string( 'START53_V7F::All' )  #Summer12_DR53X
    #process.GlobalTag.globaltag = cms.string( 'STARTHI53_V26::All' ) 
    #process.GlobalTag.globaltag = cms.string( 'START52_V5::All' ) 
    #process.GlobalTag.globaltag = cms.string( 'START52_V7::All' )
    #process.GlobalTag.globaltag = cms.string( 'START53_V17::All' )
    process.GlobalTag.globaltag = cms.string( 'START53_V27::All' )
else:
    #process.GlobalTag.globaltag = cms.string( 'FT_53_V6_AN2::All' ) #for 2012AB
    #process.GlobalTag.globaltag = cms.string( 'FT_53_V10_AN2::All' )#for 2012C
    #process.GlobalTag.globaltag = cms.string( 'FT_P_V42_AN2::All' ) #for 2012D
    process.GlobalTag.globaltag = cms.string( 'GR_P_V43D::All' )

### PoolSource will be ignored when running crab
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
#'root://eoscms//eos/cms/store/mc/Summer12/BuToJPsiK_K2MuPtEtaEtaFilter_8TeV-pythia6-evtgen/AODSIM/PU_S7_START52_V9-v1/0000/00D5E95C-7798-E111-B5B5-00237DA13CA2.root'
#'file:/raid2/w/twang/JpsiKp/JpsiKp/PyquenMix_embedHIJING_Bp2JpsiKp_5TeV_61_1_hye.root'
#'root://eoscms//eos/cms/store/user/twang/HIBmeson_20131220/test_20140106/JpsiKp/PyquenMix_embedHIJING_Bp2JpsiKp_5TeV_102_2_VpA.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBdKs/HIJINGemb_BdJpsiKs_TuneZ2star_5TeV_cff_GEN_SIM_Bd_JpsiKs_mumu.root'
'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBdKs/step4/HIJINGemb_BdJpsiKs_TuneZ2star_5TeV_cff_step4_RAW2DIGI_L1Reco_RECO_Bd_JpsiKs_mumu.root'
   )
)
#process.load("_eos_cms_store_user_twang_HIBmeson_20131220_test_20140106_JpsiKp_cff")

### Set basic filter
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
	#vertexCollection = cms.InputTag('offlinePrimaryVertices'),
	vertexCollection = cms.InputTag('offlinePrimaryVerticesWithBS'),
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

#process.filter = cms.Sequence(process.primaryVertexFilter+process.noscraping)
process.filter = cms.Sequence(process.noscraping)

##Producing Gen list with SIM particles
process.genParticlePlusGEANT = cms.EDProducer("GenPlusSimParticleProducer",
        src           = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
        setStatus     = cms.int32(8),             # set status = 8 for GEANT GPs
        filter        = cms.vstring("pt > 0.0"),  # just for testing (optional)
	    genParticles   = cms.InputTag("genParticles") # original genParticle list
)

### Setup Pat
### Ref: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATMCMatching
process.load("PhysicsTools.PatAlgos.patSequences_cff")
if HIFormat:
	process.muonMatch.matched = cms.InputTag("hiGenParticles")
	process.genParticlePlusGEANT.genParticles = cms.InputTag("hiGenParticles")

##Using GEN plus SIM list for matching
if UseGenPlusSim:
	process.muonMatch.matched = cms.InputTag("genParticlePlusGEANT")

#process.allLayer1Jets.addJetCorrFactors = False
from PhysicsTools.PatAlgos.tools.trackTools import *
#process.load( 'PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cff' )
if runOnMC:
    makeTrackCandidates(process,              # patAODTrackCands
        label='TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
        tracks=cms.InputTag('generalTracks'), # input track collection
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
        tracks=cms.InputTag('generalTracks'), # input track collection
        particleType='pi+',                   # particle type (for assigning a mass)
        preselection='pt > 0.3',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
        selection='pt > 0.3',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
        isolation={},                         # Isolations to use ('source':deltaR; set to {} for None)
        isoDeposits=[],
        mcAs=None                             # Replicate MC match as the one used for Muons
    );                                        # you can specify more than one collection for this
    l1cands = getattr(process, 'patTrackCands')
    l1cands.addGenMatch = False
from PhysicsTools.PatAlgos.tools.coreTools import *
removeAllPATObjectsBut(process, ['Muons'])
#removeSpecificPATObjects(process, ['Jets'])

if not runOnMC :
	removeMCMatching(process, ['All'] )

### Set Bfinder option
process.demo = cms.EDAnalyzer('Bfinder',
	Bchannel 		= cms.vint32(
		0,#RECONSTRUCTION: J/psi + K
		0,#RECONSTRUCTION: J/psi + Pi
		1,#RECONSTRUCTION: J/psi + Ks 
		0,#RECONSTRUCTION: J/psi + K* (K+, Pi-)
		0,#RECONSTRUCTION: J/psi + K* (K-, Pi+)
		0,#RECONSTRUCTION: J/psi + phi
		0,),#RECONSTRUCTION: J/psi + pi pi <= psi', X(3872), Bs->J/psi f0
	HLTLabel        = cms.InputTag('TriggerResults::HLT'),
    GenLabel        = cms.InputTag('genParticles'),
	MuonLabel       = cms.InputTag('selectedPatMuons'),         #selectedPatMuons
	TrackLabel      = cms.InputTag('selectedPatTrackCands'),    #selectedPat
    PUInfoLabel     = cms.InputTag("addPileupInfo"),
    BSLabel     = cms.InputTag("offlineBeamSpot"),
    PVLabel     = cms.InputTag("offlinePrimaryVerticesWithBS")
)
if HIFormat:
	process.demo.GenLabel = cms.InputTag('hiGenParticles')

if UseGenPlusSim:
	process.demo.GenLabel = cms.InputTag('genParticlePlusGEANT')

### SetUp HLT info
process.load('Bfinder.HiHLTAlgos.hltanalysis_cff')
process.hltanalysis.dummyBranches = cms.untracked.vstring()
#if HIFormat:
	#process.hltanalysis.HLTProcessName = cms.string("HISIGNAL")
	#process.hltanalysis.hltresults = cms.InputTag("TriggerResults","","HISIGNAL")
	#process.hltanalysis.l1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap::HISIGNAL")
process.hltAna = cms.Path(process.hltanalysis)

### Set output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('Bfinder_all.root')
)

if runOnMC:
	process.patDefaultSequence *= process.genParticlePlusGEANT

process.p = cms.Path(	
#	process.filter*process.genParticlePlusGEANT*process.patDefaultSequence*process.demo
	process.filter*process.patDefaultSequence*process.demo
)

process.schedule = cms.Schedule(
	process.p
	,process.hltAna
)
