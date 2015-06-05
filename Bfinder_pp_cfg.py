import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
ivars = VarParsing.VarParsing('analysis')
#ivars.inputFiles='root://eoscms//eos/cms/store/user/twang/HIBmeson_20131220/test_20140106/JpsiKp/PyquenMix_embedHIJING_Bp2JpsiKp_5TeV_102_2_VpA.root'
#ivars.inputFiles='file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBu/step4/HIJINGemb_BdJpsiKstar_TuneZ2star_5TeV_cff_step4_RAW2DIGI_L1Reco_RECO.root'
#ivars.inputFiles='file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/localRun/RunMore/PyquenMix_embedHIJING_Bp2JpsiKp_Bpt5_5TeV_boostedMC.root'
#ivars.inputFiles='file:/mnt/hadoop/cms/store/user/twang/Hijing_PPb502_MinimumBias/PyquenMix_STARTHI53_V27_HIJINGembed_pPb_step4_RAW2DIGI_L1Reco_RECO_Bpt5_BuJpsiK_20140225/5dc89fb1319c58a400229c5d020a3799/HIJINGemb_BuJpsiK_TuneZ2star_5TeV_cff_step4_RAW2DIGI_L1Reco_RECO_92_1_Yd0.root'
#ivars.inputFiles='file:/mnt/hadoop/cms/store/himc/HiWinter13/PYTHIA6_inclBtoPsiMuMu_5TeV02/GEN-SIM-RECO/pa_STARTHI53_V27-v1/20000/F66E8E9B-AD56-E311-9214-848F69FD3D0D.root'
#ivars.inputFiles='file:/mnt/hadoop/cms/store/user/tawei/Data_samples/HIRun2013/PAMuon/RECO/PromptReco-v1/000/210/676/00000/38291323-E567-E211-8DD9-5404A63886C5.root'
#ivars.inputFiles='file:/mnt/hadoop/cms/store/user/tawei/MC_samples/hckim-HIJINGemb_inclBtoPsiMuMu_5TeV_boost_FEVTDEBUGHLT_v7All/HIJINGemb_inclBtoPsiMuMu_5TeV_boost_RECO_STARTHI53_V27_evt600_3446_1_4Qm.root'
#ivars.inputFiles='file:/net/hisrv0001/home/tawei/HeavyFlavor_20131030/Gen_HiBmeson_20150121/PbPbtest/CMSSW_7_4_0_pre8/src/GenPP/step3_RAW2DIGI_L1Reco_RECO_VALIDATION.root'
#ivars.inputFiles='file:/mnt/hadoop/cms/store/user/twang/HIDiMuon/RECO_HIDiMuon_L2DoubleMu3Skim_v10_JpsiFilter_v1_CMSSW740pre8_20150428/3c3450dda05abb66de621932774972fa/hiRecoData_RAW2DIGI_L1Reco_RECO_filter_97_1_SUD.root'
ivars.inputFiles='file:/net/hisrv0001/home/tawei/HeavyFlavor_20131030/Bfinder/TempFor_CMSSW_7_4_1_20150430/test/CMSSW_7_4_1/src/data/store/relval/CMSSW_7_4_1/RelValBuJpsiK_13/GEN-SIM-RECO/MCRUN2_74_V9_gensim_740pre7-v1/00000/382A5F31-12EC-E411-879A-0025905A48EC.root'
ivars.outputFile='Bfinder_pp_all.root'
# get and parse the command line arguments
ivars.parseArguments()

### Add Calo muons
#AddCaloMuon = True
AddCaloMuon = False

### Run on MC?
runOnMC = True
#runOnMC = False

### Switching between "hiGenParticles"(pPb MC) and "genParticles" (pp MC)
#HIFormat = True
HIFormat = False

### Include SIM tracks for matching?
UseGenPlusSim = False

### Using pat muon with trigger or not
UsepatMuonsWithTrigger = True

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
    outputCommands = cms.untracked.vstring('drop *',
    )
)

### Set maxEvents
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

### Set global tag
if runOnMC:
    #process.GlobalTag.globaltag = cms.string( 'START53_V7F::All' )  #Summer12_DR53X
    #process.GlobalTag.globaltag = cms.string( 'STARTHI53_V26::All' ) 
    #process.GlobalTag.globaltag = cms.string( 'START52_V5::All' ) 
    #process.GlobalTag.globaltag = cms.string( 'START52_V7::All' )
    #process.GlobalTag.globaltag = cms.string( 'START53_V17::All' )
#    process.GlobalTag.globaltag = cms.string( 'START53_V27::All' ) ##pPb
    process.GlobalTag.globaltag = cms.string( 'GR_R_74_V8A::All' ) ##PbPb 740pre8 data
else:
    #process.GlobalTag.globaltag = cms.string( 'FT_53_V6_AN2::All' ) #for 2012AB
    #process.GlobalTag.globaltag = cms.string( 'FT_53_V10_AN2::All' )#for 2012C
#    process.GlobalTag.globaltag = cms.string( 'FT_P_V42_AN2::All' ) #for 2012D
#    process.GlobalTag.globaltag = cms.string( 'GR_P_V43D::All' ) ##pp
    process.GlobalTag.globaltag = cms.string( 'GR_R_74_V8A::All' ) ##PbPb 740pre8 data
#    process.GlobalTag.globaltag = cms.string( 'GR_P_V43F::All' ) ##pPb: /PAMuon/HIRun2013-28Sep2013-v1/RECO
#    process.GlobalTag.globaltag = cms.string( 'GR_P_V43D::All' ) ##pPb: /PAMuon/HIRun2013-PromptReco-v1/RECO

### PoolSource will be ignored when running crab
process.source = cms.Source("PoolSource",
    skipEvents=cms.untracked.uint32(0),
	fileNames = cms.untracked.vstring(ivars.inputFiles)
)

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

# Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

#process.filter = cms.Sequence(process.primaryVertexFilter+process.noscraping)
#process.filter = cms.Sequence(process.noscraping)
process.filter = cms.Sequence(process.PAcollisionEventSelection)

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
### keep only Pat:: part 
#from PhysicsTools.PatAlgos.patEventContent_cff import *
if HIFormat:
	process.muonMatch.matched = cms.InputTag("hiGenParticles")
	process.genParticlePlusGEANT.genParticles = cms.InputTag("hiGenParticles")

##Using GEN plus SIM list for matching
if UseGenPlusSim:
	process.muonMatch.matched = cms.InputTag("genParticlePlusGEANT")

#process.allLayer1Jets.addJetCorrFactors = False
#from PhysicsTools.PatAlgos.tools.trackTools import *######MODIFIED
from Bfinder.tempTools.trackTools import *######MODIFIED
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
from Bfinder.tempTools.tempTools import *######MODIFIED
removeSpecificPATObjects(process, ['Photons', 'Electrons', 'Taus', 'Jets', 'METs', 'Muons'])######MODIFIED
#removeAllPATObjectsBut(process, ['Muons'])
#removeSpecificPATObjects(process, ['Jets'])

if not runOnMC :
	removeMCMatching(process, ['All'] )

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
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
    tracks = cms.InputTag("generalTracks"),
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
    process.patDefaultSequence = cms.Sequence(process.mergedMuons*process.patDefaultSequence)

### Set Bfinder option
process.demo = cms.EDAnalyzer('Bfinder',
	Bchannel 		= cms.vint32(
		1,#RECONSTRUCTION: J/psi + K
		1,#RECONSTRUCTION: J/psi + Pi
		1,#RECONSTRUCTION: J/psi + Ks 
		1,#RECONSTRUCTION: J/psi + K* (K+, Pi-)
		1,#RECONSTRUCTION: J/psi + K* (K-, Pi+)
		1,#RECONSTRUCTION: J/psi + phi
		1,),#RECONSTRUCTION: J/psi + pi pi <= psi', X(3872), Bs->J/psi f0
#    MuonTriggerMatchingPath = cms.vstring("HLT_PAMu3_v1"),
    MuonTriggerMatchingPath = cms.vstring("HLT_PAMu3_v*"),
#    MuonTriggerMatchingPath = cms.vstring("HLT_PAMu3_v*", "HLT_PAMu7_v*", "HLT_PAMu12_v*"),
	HLTLabel        = cms.InputTag('TriggerResults::HLT'),
    GenLabel        = cms.InputTag('genParticles'),
	MuonLabel       = cms.InputTag('selectedPatMuons'),         #selectedPatMuons
	TrackLabel      = cms.InputTag('selectedPatTrackCands'),    #selectedPat
    PUInfoLabel     = cms.InputTag("addPileupInfo"),
    BSLabel     = cms.InputTag("offlineBeamSpot"),
    PVLabel     = cms.InputTag("offlinePrimaryVerticesWithBS"),
    tkPtCut = cms.double(0.4),
    jpsiPtCut = cms.double(0.0),
    bPtCut = cms.double(0.0),
    RunOnMC = cms.bool(False),
    doTkPreCut = cms.bool(True),
    doMuPreCut = cms.bool(True)
)
if HIFormat:
	process.demo.GenLabel = cms.InputTag('hiGenParticles')
if UseGenPlusSim:
	process.demo.GenLabel = cms.InputTag('genParticlePlusGEANT')
if UsepatMuonsWithTrigger:
	process.demo.MuonLabel = cms.InputTag('patMuonsWithTrigger')	

### SetUp HLT info
#process.load('Bfinder.HiHLTAlgos.hltanalysis_cff')
process.load('Bfinder.EventAnalysis.hltanalysis_cff')
process.hltanalysis.dummyBranches = cms.untracked.vstring()
#if HIFormat:
	#process.hltanalysis.mctruth = cms.InputTag("hiGenParticles")
	#process.hltanalysis.HLTProcessName = cms.string("HISIGNAL")
	#process.hltanalysis.hltresults = cms.InputTag("TriggerResults","","HISIGNAL")
	#process.hltanalysis.l1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap::HISIGNAL")
process.hltAna = cms.Path(process.filter*process.hltanalysis)

### Set output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(ivars.outputFile)
)

if runOnMC and UseGenPlusSim:
	process.patDefaultSequence *= process.genParticlePlusGEANT
if UsepatMuonsWithTrigger:
	process.patDefaultSequence *= process.patMuonsWithTriggerSequence

process.p = cms.Path(	
#    process.patDefaultSequence*process.demo
    process.filter*process.patDefaultSequence*process.demo
)
#process.e = cms.EndPath(process.out)
process.schedule = cms.Schedule(
	process.p
	,process.hltAna
#    ,process.e
)
