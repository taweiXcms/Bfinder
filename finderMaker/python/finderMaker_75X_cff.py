import FWCore.ParameterSet.Config as cms

def finderMaker_75X(process, AddCaloMuon = False, runOnMC = True, HIFormat = False, UseGenPlusSim = False, VtxLabel = "hiSelectedVertex", TrkLabel = "hiGeneralTracks"):
	### Set TransientTrackBuilder 
	process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
	
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
	process.pfPileUpIsoPFBRECO.Vertices = cms.InputTag(VtxLabel)
	process.pfPileUpPFBRECO.Vertices = cms.InputTag(VtxLabel)
	###### Needed in CMSSW7
	
	if HIFormat:
		process.muonMatch.matched = cms.InputTag("hiGenParticles")
		process.genParticlePlusGEANT.genParticles = cms.InputTag("hiGenParticles")
	
	##Using GEN plus SIM list for matching
	if UseGenPlusSim:
		process.muonMatch.matched = cms.InputTag("genParticlePlusGEANT")
	
	## TrackCand
	from PhysicsTools.PatAlgos.tools.trackTools import makeTrackCandidates
	if runOnMC:
	    makeTrackCandidates(process,              # patAODTrackCands
	        label='TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
	        tracks=cms.InputTag(TrkLabel), # input track collection
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
	        tracks=cms.InputTag(TrkLabel), # input track collection
	        particleType='pi+',                   # particle type (for assigning a mass)
	        preselection='pt > 0.3',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
	        selection='pt > 0.3',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
	        isolation={},                         # Isolations to use ('source':deltaR; set to {} for None)
	        isoDeposits=[],
	        mcAs=None                             # Replicate MC match as the one used for Muons
	    );                                        # you can specify more than one collection for this
	    l1cands = getattr(process, 'patTrackCands')
	    l1cands.addGenMatch = False
	if runOnMC:
		process.TrackCandSequence = cms.Sequence(process.patAODTrackCandsUnfiltered*process.patAODTrackCands*process.patTrackCandsMCMatch*process.patTrackCands*process.selectedPatTrackCands)
	else:
		process.TrackCandSequence = cms.Sequence(process.patAODTrackCandsUnfiltered*process.patAODTrackCands*process.patTrackCands*process.selectedPatTrackCands)
	
	## patMuonsWithTrigger
	process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
	from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, changeTriggerProcessName, switchOffAmbiguityResolution
	#process.patMuonsWithTriggerSequence = cms.Sequence(process.pfParticleSelectionForIsoSequence*process.muonPFIsolationPATSequence*process.patMuonsWithTriggerSequence)
	process.patMuonsWithTriggerSequence = cms.Sequence(process.patMuonsWithTriggerSequence)
	process.patMuonsWithoutTrigger.isoDeposits = cms.PSet()
	process.patMuonsWithoutTrigger.isolationValues = cms.PSet()
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
	    tracks = cms.InputTag(TrkLabel),
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
	process.Bfinder = cms.EDAnalyzer('Bfinder',
		Bchannel 		= cms.vint32(
			1,#RECONSTRUCTION: J/psi + K
			0,#RECONSTRUCTION: J/psi + Pi
			0,#RECONSTRUCTION: J/psi + Ks 
			0,#RECONSTRUCTION: J/psi + K* (K+, Pi-)
			0,#RECONSTRUCTION: J/psi + K* (K-, Pi+)
			0,#RECONSTRUCTION: J/psi + phi
			0,#RECONSTRUCTION: J/psi + pi pi <= psi', X(3872), Bs->J/psi f0
		),
        detailMode = cms.bool(True),
        dropUnusedTracks = cms.bool(True),
	    #MuonTriggerMatchingPath = cms.vstring("HLT_PAMu3_v*"),
	    MuonTriggerMatchingPath = cms.vstring("HLT_HIL2DoubleMu3_v*"),
	    #MuonTriggerMatchingPath = cms.vstring("HLT_PAMu3_v*", "HLT_PAMu7_v*", "HLT_PAMu12_v*"),
		HLTLabel = cms.InputTag('TriggerResults::HLT'),
	    GenLabel = cms.InputTag('genParticles'),
		MuonLabel = cms.InputTag('patMuonsWithTrigger'),
		TrackLabel = cms.InputTag('patTrackCands'),
        MVAMapLabel = cms.string(TrkLabel),
	    PUInfoLabel = cms.InputTag("addPileupInfo"),
	    BSLabel = cms.InputTag("offlineBeamSpot"),
	    PVLabel = cms.InputTag(VtxLabel),
	    tkPtCut = cms.double(1.0),#before fit
	    tkEtaCut = cms.double(999.0),#before fit
	    jpsiPtCut = cms.double(3.0),#before fit
	    bPtCut = cms.double(5.0),#before fit
	    bEtaCut = cms.double(2.4),#before fit, not used currently
		VtxChiProbCut = cms.double(0.01),
	    svpvDistanceCut = cms.double(0.0),
	    MaxDocaCut = cms.double(999.),
	    alphaCut = cms.double(999.),
	    RunOnMC = cms.bool(False),
	    doTkPreCut = cms.bool(True),
	    doMuPreCut = cms.bool(True),
	    makeBntuple = cms.bool(True),
	    doBntupleSkim = cms.bool(False),
	)
	### Set Dfinder option
	process.Dfinder = cms.EDAnalyzer('Dfinder',
		Dchannel 		= cms.vint32(
	        1,#RECONSTRUCTION: K+pi- : D0bar
	        1,#RECONSTRUCTION: K-pi+ : D0
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
		),
        detailMode = cms.bool(False),
        dropUnusedTracks = cms.bool(True),
		HLTLabel = cms.InputTag('TriggerResults::HLT'),
	    GenLabel = cms.InputTag('genParticles'),
		TrackLabel = cms.InputTag('patTrackCands'),
		MVAMapLabel = cms.string(TrkLabel),
	    PUInfoLabel = cms.InputTag("addPileupInfo"),
	    BSLabel = cms.InputTag("offlineBeamSpot"),
	    PVLabel = cms.InputTag(VtxLabel),
	    tkPtCut = cms.double(1.),#before fit
	    tkEtaCut = cms.double(2.0),#before fit
	    dPtCut = cms.double(8.0),#before fit
	    dEtaCut = cms.double(1.5),#before fit, not used currently
		VtxChiProbCut = cms.double(0.05),
	    svpvDistanceCut = cms.double(0.0),
	    MaxDocaCut = cms.double(999.),
	    alphaCut = cms.double(999.),
	    RunOnMC = cms.bool(False),
	    doTkPreCut = cms.bool(True),
	    makeDntuple = cms.bool(True),
	    doDntupleSkim = cms.bool(False),
	)
	if runOnMC:
	    process.Bfinder.RunOnMC = cms.bool(True)
	    process.Dfinder.RunOnMC = cms.bool(True)
	if HIFormat:
		process.Bfinder.GenLabel = cms.InputTag('hiGenParticles')
		process.Dfinder.GenLabel = cms.InputTag('hiGenParticles')
	if UseGenPlusSim:
		process.Bfinder.GenLabel = cms.InputTag('genParticlePlusGEANT')
		process.Dfinder.GenLabel = cms.InputTag('genParticlePlusGEANT')
	
	if runOnMC and UseGenPlusSim:
		process.patMuonsWithTriggerSequence *= process.genParticlePlusGEANT
	
	process.BfinderSequence = cms.Sequence(process.patMuonsWithTriggerSequence*process.TrackCandSequence*process.Bfinder)
	process.DfinderSequence = cms.Sequence(process.TrackCandSequence*process.Dfinder)
	process.finderSequence = cms.Sequence(process.patMuonsWithTriggerSequence*process.TrackCandSequence*process.Bfinder*process.Dfinder)
