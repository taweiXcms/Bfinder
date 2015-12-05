import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
ivars = VarParsing.VarParsing('analysis')
#ivars.inputFiles='file:/mnt/hadoop/cms/store/user/tawei/Data_samples/HIRun2011/HIDiMuon/RAW/v1/000/183/013/02A69E73-C01F-E111-A008-00237DDC5AF6.root'
#ivars.inputFiles='file:/net/hisrv0001/home/tawei/HeavyFlavor_20131030/Skim/skim20141215/CMSSW_5_3_20/src/test/NoHLTfilter.root'
ivars.inputFiles='file:/net/hisrv0001/home/tawei/HeavyFlavor_20131030/Skim/skim20141215/CMSSW_5_3_20/src/test/HLTfilter.root'
ivars.outputFile='NoHLTfilter_all.root'
ivars.parseArguments()
process = cms.Process("demo")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.TFileService = cms.Service("TFileService",fileName = cms.string(ivars.outputFile))
process.source = cms.Source("PoolSource",
    skipEvents=cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(ivars.inputFiles)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
### Set Geometry/GlobalTag/BField
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

process.GlobalTag.globaltag = cms.string( 'GR_R_53_LV6::All' ) ##PbPb
process.load('Bfinder.EventAnalysis.hltanalysis_cff')
process.hltanalysis.dummyBranches = cms.untracked.vstring()
#if HIFormat:
    #process.hltanalysis.HLTProcessName = cms.string("HISIGNAL")
    #process.hltanalysis.hltresults = cms.InputTag("TriggerResults","","HISIGNAL")
    #process.hltanalysis.l1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap::HISIGNAL")
process.hltAna = cms.Path(process.hltanalysis)
process.schedule = cms.Schedule(process.hltAna)
