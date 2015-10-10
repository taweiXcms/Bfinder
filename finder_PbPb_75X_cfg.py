import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
ivars = VarParsing.VarParsing('analysis')
#ivars.inputFiles='file:/mnt/hadoop/cms/store/user/twang/HIDiMuon/RECO_HIDiMuon_L2DoubleMu3Skim_v10_JpsiFilter_v1_CMSSW740pre8_20150428/3c3450dda05abb66de621932774972fa/hiRecoData_RAW2DIGI_L1Reco_RECO_filter_975_1_PTa.root'
#ivars.inputFiles='file:/mnt/hadoop/cms/store/user/twang/Pyquen_CMSSW742_Unquenched_PbPb_2760GeV_GEN_SIM_PU_BuKp_20150526_100kevt/Pyquen_CMSSW742_Unquenched_PbPb_2760GeV_step3_BuKp_20150526_100kevt/27ff3fcdfd0b68d12bfbb80768287940/step3_RAW2DIGI_L1Reco_RECO_PU_90_1_Ole.root'
ivars.inputFiles='file:/mnt/hadoop/cms/store/user/richard/MBHydjet5020/Hydjet_Quenched_MinBias_5020GeV/HydjetMB5020_750_75X_mcRun2_HeavyIon_v1_RealisticHICollisions2011_STARTHI50_mc_RECOSIM_v3/150729_144407/0000/step3_98.root'

ivars.outputFile='Bfinder_PbPb_all.root'
# get and parse the command line arguments
ivars.parseArguments()

process = cms.Process("demo")

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

### Set maxEvents
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

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
process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

### All relevent GlobalTags
#MC
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
#Data
#process.GlobalTag.globaltag = cms.string( 'FT_53_V6_AN2::All' ) #for 2012AB
#process.GlobalTag.globaltag = cms.string( 'FT_53_V10_AN2::All' )#for 2012C
#process.GlobalTag.globaltag = cms.string( 'FT_P_V42_AN2::All' ) #for 2012D
#process.GlobalTag.globaltag = cms.string( 'GR_R_53_LV6::All' ) ##PbPb
#process.GlobalTag.globaltag = cms.string( 'GR_P_V43D::All' ) ##pp
#process.GlobalTag.globaltag = cms.string( 'GR_P_V43F::All' ) ##pPb: /PAMuon/HIRun2013-28Sep2013-v1/RECO
#process.GlobalTag.globaltag = cms.string( 'GR_P_V43D::All' ) ##pPb: /PAMuon/HIRun2013-PromptReco-v1/RECO
#process.GlobalTag.globaltag = cms.string( 'GR_R_74_V8A::All' ) ##CMSSW_7_4_0_pre8 PbPb
#process.GlobalTag.globaltag = cms.string( 'GR_R_74_V12A::All' ) ##CMSSW_7_4_2 PbPb
#process.GlobalTag.globaltag = cms.string( '75X_dataRun2_v2' ) ##CMSSW_7_5_0 PbPb
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v1', '')

#### HI infomation
from GeneratorInterface.HiGenCommon.HeavyIon_cff import *
process.load('GeneratorInterface.HiGenCommon.HeavyIon_cff')

#### SetUp Evt Info (centrality)
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

### Custom options
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
CentralityFilter = True

### Add centrality filter
if CentralityFilter:
    process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")
    #process.cenfilter = process.centralityFilter.clone(selectedBins = [0,1,2,3,4])
    process.cenfilter = process.centralityFilter.clone(selectedBins = range(59,201))
    process.filter = cms.Sequence(process.centralityBin*process.cenfilter*process.collisionEventSelection)

### Run the hiEvtAnalyzer sequence
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.evtAna = cms.Path(process.filter*process.hiEvtAnalyzer)

if runOnMC:
	process.hiEvtAnalyzer.doMC = cms.bool(True)
	process.evtAna = cms.Path(process.filter*process.heavyIon*process.hiEvtAnalyzer)

### Run HLT info sequence
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
process.hltanalysis.dummyBranches = cms.untracked.vstring()
process.hltanalysis.OfflinePrimaryVertices0 = cms.InputTag("hiSelectedVertex")
#if HIFormat:
    #process.hltanalysis.mctruth = cms.InputTag("hiGenParticles")# Will cause segmentation violation
    #process.hltanalysis.HLTProcessName = cms.string("HISIGNAL")
    #process.hltanalysis.hltresults = cms.InputTag("TriggerResults","","HISIGNAL")
    #process.hltanalysis.l1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap::HISIGNAL")
process.hltAna = cms.Path(process.filter*process.hltanalysis)

### finder building block
from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X
finderMaker_75X(process, AddCaloMuon, runOnMC, HIFormat, UseGenPlusSim)
process.p = cms.Path(process.filter*process.finderSequence)

process.schedule = cms.Schedule(
	process.p
	,process.hltAna
	,process.evtAna
#    ,process.e
)
