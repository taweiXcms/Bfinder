import FWCore.ParameterSet.Config as cms

### Run on MC?
runOnMC = True

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
    process.GlobalTag.globaltag = cms.string( 'START52_V5::All' ) 
else:
    #process.GlobalTag.globaltag = cms.string( 'FT_53_V6_AN2::All' ) #for 2012AB
    #process.GlobalTag.globaltag = cms.string( 'FT_53_V10_AN2::All' )#for 2012C
    process.GlobalTag.globaltag = cms.string( 'FT_P_V42_AN2::All' ) #for 2012D


### PoolSource will be ignored when running crab
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
       #'root://eoscms//eos/cms/store/user/pchen/dumpData/MuOnia_Run2012A-13Jul2012-v1_AOD/fetch_1_1_jBm.root'
       #'file:///afs/cern.ch/user/p/pchen/Upsilon2S_RECO_9_1_BHo.root'
#'/store/data/Run2012A/MuOnia/AOD/13Jul2012-v1/00000/B8AA89D2-2ECF-E111-8084-003048FFCB84.root'
'root://eoscms//eos/cms/store/mc/Summer12/BuToJPsiK_K2MuPtEtaEtaFilter_8TeV-pythia6-evtgen/AODSIM/PU_S7_START52_V9-v1/0000/00D5E95C-7798-E111-B5B5-00237DA13CA2.root'
   )
)

### Set basic filter
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
vertexCollection = cms.InputTag('offlinePrimaryVertices'),
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

process.filter = cms.Sequence(process.primaryVertexFilter+process.noscraping)

### Setup Pat
    ### Ref: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATMCMatching
process.load("PhysicsTools.PatAlgos.patSequences_cff")
#process.allLayer1Jets.addJetCorrFactors = False
from PhysicsTools.PatAlgos.tools.trackTools import *
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
    process.patTrackCandsMCMatch.resolveByMatchQuality = cms.bool(True)
    process.patTrackCandsMCMatch.resolveAmbiguities = cms.bool(True)
    process.patTrackCandsMCMatch.checkCharge = cms.bool(True)
    process.patTrackCandsMCMatch.maxDPtRel = cms.double(0.5)
    process.patTrackCandsMCMatch.maxDeltaR = cms.double(0.7)
    process.patTrackCandsMCMatch.mcPdgId = cms.vint32(111, 211, 311, 321)
#    process.patTrackCandsMCMatch.mcPdgId = cms.vint32(all)
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
	HLTLabel        = cms.InputTag('TriggerResults::HLT'),
    GenLabel        = cms.InputTag('genParticles'),
	MuonLabel       = cms.InputTag('selectedPatMuons'),         #selectedPatMuons
	TrackLabel      = cms.InputTag('selectedPatTrackCands'),    #selectedPat
    PUInfoLabel     = cms.InputTag("addPileupInfo"),
	NtupleType      = cms.untracked.string('jpsi')               #['all','upsilon','jpsi','no_c'], which mass constraint
)

### SetUp HLT info
process.load('Bfinder.HiHLTAlgos.hltanalysis_cff')
process.hltanalysis.dummyBranches = cms.untracked.vstring()
process.hltAna = cms.Path(process.hltanalysis)

### Set output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('Bfinder_all.root')
)

process.p = cms.Path(	
	process.filter*process.patDefaultSequence*process.demo
)

process.schedule = cms.Schedule(
	process.p
	,process.hltAna
)
