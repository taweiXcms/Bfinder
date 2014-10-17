import FWCore.ParameterSet.Config as cms

process = cms.Process('EvtAna')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
			    fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/hidata/HIRun2013/PAMinBiasUPC/RECO/28Sep2013-v1/10000/001397FC-462D-E311-A034-782BCB3BCADD.root'),
)

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_P_V43D::All'
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')

process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')

process.load('HeavyIonsAnalysis.VertexAnalysis.PAPileUpVertexFilter_cff')

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltZeroBiasSingleTrack = process.hltHighLevel.clone()
process.hltZeroBiasSingleTrack.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]

process.load('RecoHI.HiCentralityAlgos.pACentrality_cfi')

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
	centralityVariable = cms.string("HFtowersTrunc"), #or HFtowersPlusTrunc
	nonDefaultGlauberModel = cms.string(""),
	centralitySrc = cms.InputTag("pACentrality"),
	pPbRunFlip = cms.untracked.uint32(211313)
	)

process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string("test.root"))

process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.hiEvtAnalyzer.Centrality = cms.InputTag("pACentrality")
process.hiEvtAnalyzer.Vertex = cms.InputTag("offlinePrimaryVertices")

process.p = cms.Path(process.hltZeroBiasSingleTrack * process.PAcollisionEventSelection * process.pileupVertexFilterCutGplus * process.siPixelRecHits * process.pACentrality * process.hiEvtAnalyzer)
