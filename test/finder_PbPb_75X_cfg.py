import FWCore.ParameterSet.Config as cms
process = cms.Process('HiForest')
import FWCore.ParameterSet.VarParsing as VarParsing
ivars = VarParsing.VarParsing('analysis')

#ivars.inputFiles='file:/mnt/hadoop/cms/store/user/twang/HIDiMuon/RECO_HIDiMuon_L2DoubleMu3Skim_v10_JpsiFilter_v1_CMSSW740pre8_20150428/3c3450dda05abb66de621932774972fa/hiRecoData_RAW2DIGI_L1Reco_RECO_filter_975_1_PTa.root'
#ivars.inputFiles='file:/mnt/hadoop/cms/store/user/twang/Pyquen_CMSSW742_Unquenched_PbPb_2760GeV_GEN_SIM_PU_BuKp_20150526_100kevt/Pyquen_CMSSW742_Unquenched_PbPb_2760GeV_step3_BuKp_20150526_100kevt/27ff3fcdfd0b68d12bfbb80768287940/step3_RAW2DIGI_L1Reco_RECO_PU_90_1_Ole.root'
#ivars.inputFiles='file:/mnt/hadoop/cms/store/user/richard/MBHydjet5020/Hydjet_Quenched_MinBias_5020GeV/HydjetMB5020_750_75X_mcRun2_HeavyIon_v1_RealisticHICollisions2011_STARTHI50_mc_RECOSIM_v3/150729_144407/0000/step3_98.root'
#ivars.inputFiles='file:/data/twang/MC_samples/Hydjet_Quenched_MinBias_5020GeV_750/Hydjet_Quenched_MinBias_5020GeV_750_HiFall15_step3_20151110/step3_RAW2DIGI_L1Reco_RECO_2054_1_2fm.root'
#ivars.inputFiles='file:/data/twang/MC_samples/Pythia8_BuToJpsiK_TuneCUEP8M1_5020GeV_BPHMod_filter_GEN_SIM_PU_20151105/Pythia8_BuToJpsiK_TuneCUEP8M1_5020GeV_BPHMod_filter_step3_20151105/step3_RAW2DIGI_L1Reco_RECO_756_1_ecf.root'
#ivars.inputFiles='file:/data/twang/MC_samples/Pythia8_BdToJpsiKstar_TuneCUEP8M1_5020GeV_BPHMod_filter_GEN_SIM_PU_20151105/Pythia8_BdToJpsiKstar_TuneCUEP8M1_5020GeV_BPHMod_filter_step3_20151105/step3_RAW2DIGI_L1Reco_RECO_698_2_gsq.root'
#ivars.inputFiles='file:/data/twang/MC_samples/Pythia8_BsToJpsiPhi_TuneCUEP8M1_5020GeV_BPHMod_filter_GEN_SIM_PU_20151105/Pythia8_BsToJpsiPhi_TuneCUEP8M1_5020GeV_BPHMod_filter_step3_PU_20151105/step3_RAW2DIGI_L1Reco_RECO_721_1_JPr.root'
#ivars.inputFiles='file:/data/twang/MC_samples/Pythia8_5020GeV_DstarD0kpi_755patch3_GEN_SIM_PU_20151120/Pythia8_5020GeV_DstarD0kpi_755patch3_step3_20151120/step3_RAW2DIGI_L1Reco_RECO_623_1_QKI.root'
#ivars.inputFiles='file:/data/twang/MC_samples/Pythia8_5020GeV_DstarD0kpipipi_755patch3_GEN_SIM_PU_20151120/Pythia8_5020GeV_DstarD0kpipipi_755patch3_step3_20151120/step3_RAW2DIGI_L1Reco_RECO_614_1_jV8.root'
#ivars.inputFiles='file:/data/twang/Data_samples/HIRun2015/HIOniaL1DoubleMu0/AOD/PromptReco-v1/000/262/735/00000/B0EC6FA1-4E99-E511-B663-02163E013910.root'#HIOniaL1DoubleMu0
#ivars.inputFiles='file:/data/twang/Data_samples/HIRun2015/HIHardProbes/AOD/PromptReco-v1/000/262/735/00000/FC577170-6C99-E511-9848-02163E014120.root'#HIHardProbes
#ivars.inputFiles='file:/data/twang/Data_samples/HIRun2015/HIHardProbes/RECO/D0Meson-PromptReco-v1/000/262/735/00000/6E423E98-5C99-E511-B72B-02163E0138EE.root'#HIHardProbes
#ivars.inputFiles='file:/data/twang/Data_samples/HIRun2015/HIHardProbes/RECO/D0Meson-PromptReco-v1/000/262/735/00000/E67242E4-5E99-E511-947B-02163E0127B4.root'#HIHardProbes
#ivars.inputFiles='file:/data/twang/Data_samples/HIRun2015/HIMinimumBias1/AOD/PromptReco-v1/000/262/726/00000/EE7F4A63-4599-E511-9CE5-02163E013850.root'#HIMinimumBias1
ivars.outputFile='finder_PbPb.root'
ivars.parseArguments()# get and parse the command line arguments

### Custom options
########## MUST CUSTOMIZE THE FOLLOWING THREE ##########
### PbPb B/Dfinder recommended setting, choose only one from them or set all to false and made your own setting
PbPbBdefault = True
PbPbDHPdefault = False
PbPbDMBdefault = False

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
VtxLabel = "hiSelectedVertex"
#TrkLabel = "hiGeneralTracks"

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
	globalTag = 'auto:run2_mc_HIon'
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

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag, '')

#### HI infomation
from GeneratorInterface.HiGenCommon.HeavyIon_cff import *
process.load('GeneratorInterface.HiGenCommon.HeavyIon_cff')

#### SetUp Evt Info (centrality)
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
if not runOnMC:
	process.centralityBin.Centrality = cms.InputTag("hiCentrality")
	process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centrality_path = cms.Path(process.centralityBin)

### Run the hiEvtAnalyzer sequence
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.evtAna = cms.Path(process.hiEvtAnalyzer)

if runOnMC:
	process.hiEvtAnalyzer.doMC = cms.bool(True)
	process.hiEvtAnalyzer.doHiMC = cms.bool(True)
	#process.evtAna = cms.Path(process.heavyIon*process.hiEvtAnalyzer)

### Set basic filter
# Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.phfCoincFilter = cms.Path(process.hfCoincFilter )
process.phfCoincFilter3 = cms.Path(process.hfCoincFilter3 )
process.pcollisionEventSelection = cms.Path(process.collisionEventSelection)
if RunOnAOD:
	process.pcollisionEventSelection = cms.Path(process.collisionEventSelectionAOD)

### Run HLT info sequence
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
from HeavyIonsAnalysis.EventAnalysis.dummybranches_cff import addHLTdummybranches
addHLTdummybranches(process)
process.hltanalysis.OfflinePrimaryVertices0 = cms.InputTag(VtxLabel)
#if HIFormat:
    #process.hltanalysis.mctruth = cms.InputTag("hiGenParticles")# Will cause segmentation violation
    #process.hltanalysis.HLTProcessName = cms.string("HISIGNAL")
    #process.hltanalysis.hltresults = cms.InputTag("TriggerResults","","HISIGNAL")
    #process.hltanalysis.l1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap::HISIGNAL")
process.hltAna = cms.Path(process.hltanalysis)

### finder building block
from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X
finderMaker_75X(process, AddCaloMuon, runOnMC, HIFormat, UseGenPlusSim)
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
## PbPb Bfinder setting on DoubleMu
if PbPbBdefault and not PbPbDMBdefault and not PbPbDHPdefault:
    process.Bfinder.tkPtCut = cms.double(0.8)#before fit
    process.Bfinder.jpsiPtCut = cms.double(0.0)#before fit
    process.Bfinder.bPtCut = cms.vdouble(5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0)#before fit
    process.Bfinder.Bchannel = cms.vint32(1, 0, 0, 1, 1, 1, 1)
    process.p = cms.Path(process.BfinderSequence)
## PbPb Dfinder setting on HardProbe
if PbPbDHPdefault and not PbPbDMBdefault and not PbPbBdefault:
    process.Dfinder.tkPtCut = cms.double(2.5)#before fit
    process.Dfinder.dPtCut = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)#before fit
    process.Dfinder.dCutSeparating_PtVal = cms.vdouble(5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.)
    process.Dfinder.tktkRes_svpvDistanceCut_lowptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5)
    process.Dfinder.tktkRes_svpvDistanceCut_highptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5)
    process.Dfinder.svpvDistanceCut_lowptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0.)
    process.Dfinder.svpvDistanceCut_highptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0.)
    process.Dfinder.Dchannel = cms.vint32(1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
    process.p = cms.Path(process.DfinderSequence)
## PbPb Dfinder setting on MB
if PbPbDMBdefault and not PbPbDHPdefault and not PbPbBdefault:
    process.Dfinder.tkPtCut = cms.double(0.7)#before fit
    process.Dfinder.dPtCut = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)#before fit
    process.Dfinder.dCutSeparating_PtVal = cms.vdouble(5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.)
    process.Dfinder.tktkRes_svpvDistanceCut_lowptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5)
    process.Dfinder.tktkRes_svpvDistanceCut_highptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5)
    process.Dfinder.svpvDistanceCut_lowptD = cms.vdouble(4.0, 4.0, 2.5, 2.5, 2.5, 2.5, 4.0, 4.0, 0., 0., 0., 0.)
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
