#file path will be:
#outLFNDirBase/inputDataset/requestName/time_tag/...
from CRABClient.UserUtilities import config
config = config()
config.General.transferOutputs = True
config.General.requestName = 'BfinderMC_Bu_20151124'
#config.General.requestName = 'BfinderMC_Bd_20151124'
#config.General.requestName = 'BfinderMC_Bs_20151124'
config.General.workArea = 'crab_projects'
config.JobType.psetName = 'finder_PbPb_75X_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['rssLimit']
config.JobType.pyCfgParams = ['noprint']
config.JobType.outputFiles = ['finder_PbPb.root']
config.Data.inputDataset = '/Pythia8_BuToJpsiK_TuneCUEP8M1_5020GeV_BPHMod_filter_GEN_SIM_PU_20151105/twang-Pythia8_BuToJpsiK_TuneCUEP8M1_5020GeV_BPHMod_filter_step3_20151105-61b88e8365cfdee4d0bfdfd20dfa5ba1/USER'
#config.Data.inputDataset = '/Pythia8_BdToJpsiKstar_TuneCUEP8M1_5020GeV_BPHMod_filter_GEN_SIM_PU_20151105/twang-Pythia8_BdToJpsiKstar_TuneCUEP8M1_5020GeV_BPHMod_filter_step3_20151105-61b88e8365cfdee4d0bfdfd20dfa5ba1/USER'
#config.Data.inputDataset = '/Pythia8_BsToJpsiPhi_TuneCUEP8M1_5020GeV_BPHMod_filter_GEN_SIM_PU_20151105/twang-Pythia8_BsToJpsiPhi_TuneCUEP8M1_5020GeV_BPHMod_filter_step3_PU_20151105-61b88e8365cfdee4d0bfdfd20dfa5ba1/USER'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 500
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.outLFNDirBase = '/store/user/twang/BfinderRun2'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/HeavyFlavorRun2/BfinderRun2'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/5TeV/Cert_262081-262273_5TeV_PromptReco_Collisions15_25ns_JSON_MuonPhys_v2.txt'
#config.Site.storageSite = 'T2_US_MIT'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_US_MIT']
