import FWCore.ParameterSet.Config as cms

process = cms.Process('EvtAna')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/hidata/HIRun2011/HIMinBiasUPC/RECO/14Mar2014-v2/00000/0011854A-F0AD-E311-885B-FA163E32A814.root', 'root://xrootd.unl.edu//store/hidata/HIRun2011/HIMinBiasUPC/RECO/14Mar2014-v2/00000/0018A8E7-F9AF-E311-ADAB-FA163E565820.root', 'root://xrootd.unl.edu//store/hidata/HIRun2011/HIMinBiasUPC/RECO/14Mar2014-v2/00000/002456EE-FCAF-E311-87B3-FA163E632CDA.root'),
)

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_53_LV6::All'
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')

process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltMinBiasHFOrBSC = process.hltHighLevel.clone()
process.hltMinBiasHFOrBSC.HLTPaths = ["HLT_HIMinBiasHfOrBSC_v1"]

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string("test.root"))

process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')

process.p = cms.Path(process.hltMinBiasHFOrBSC * process.collisionEventSelection * process.hiEvtAnalyzer)
