#file path will be:
#outLFNDirBase/inputDataset/requestName/time_tag/...
from CRABClient.UserUtilities import config
config = config()
config.General.transferOutputs = True
config.General.transferLogs = True

#config.General.requestName = 'BfinderMC_PbPb_Pythia8_BuToJpsiK_Bpt0_Pthat0_TuneCUEP8M1_20160816_bPt5jpsiPt0tkPt0p8_Bp'
#config.General.requestName = 'BfinderMC_PbPb_Pythia8_BuToJpsiK_Bpt0_Pthat5_TuneCUEP8M1_20160816_bPt5jpsiPt0tkPt0p8_Bp'
#config.General.requestName = 'BfinderMC_PbPb_Pythia8_BuToJpsiK_Bpt0_Pthat15_TuneCUEP8M1_20160816_bPt5jpsiPt0tkPt0p8_Bp'
#config.General.requestName = 'BfinderMC_PbPb_Pythia8_BuToJpsiK_Bpt0_Pthat30_TuneCUEP8M1_20160816_bPt5jpsiPt0tkPt0p8_Bp'
#config.General.requestName = 'BfinderMC_PbPb_Pythia8_BuToJpsiK_Bpt0_Pthat50_TuneCUEP8M1_20160816_bPt5jpsiPt0tkPt0p8_Bp'
#config.General.requestName = 'BfinderMC_PbPb_Pythia8_BJpsiMM_ptJpsi_00_03_Hydjet_MB_20160816_bPt5jpsiPt0tkPt0p8_Bp'
#config.General.requestName = 'BfinderMC_PbPb_Pythia8_BJpsiMM_ptJpsi_03_06_Hydjet_MB_20160816_bPt5jpsiPt0tkPt0p8_Bp'
#config.General.requestName = 'BfinderMC_PbPb_Pythia8_BJpsiMM_ptJpsi_06_09_Hydjet_MB_20160816_bPt5jpsiPt0tkPt0p8_Bp'
#config.General.requestName = 'BfinderMC_PbPb_Pythia8_BJpsiMM_ptJpsi_09_12_Hydjet_MB_20160816_bPt5jpsiPt0tkPt0p8_Bp'
#config.General.requestName = 'BfinderMC_PbPb_Pythia8_BJpsiMM_ptJpsi_12_15_Hydjet_MB_20160816_bPt5jpsiPt0tkPt0p8_Bp'
#config.General.requestName = 'BfinderMC_PbPb_Pythia8_BJpsiMM_ptJpsi_15_30_Hydjet_MB_20160816_bPt5jpsiPt0tkPt0p8_Bp'
#config.General.requestName = 'BfinderMC_PbPb_Pythia8_BJpsiMM_ptJpsi_30_inf_Hydjet_MB_20160816_bPt5jpsiPt0tkPt0p8_Bp'

config.General.workArea = 'crab_projects'
config.JobType.psetName = 'finder_PbPb_75X_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['rssLimit']
config.JobType.pyCfgParams = ['noprint']
config.JobType.outputFiles = ['finder_PbPb.root']

#config.Data.inputDataset = '/Pythia8_BuToJpsiK_Bpt0_Pthat0_TuneCUEP8M1/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_BuToJpsiK_Bpt0_Pthat5_TuneCUEP8M1/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_BuToJpsiK_Bpt0_Pthat15_TuneCUEP8M1/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_BuToJpsiK_Bpt0_Pthat30_TuneCUEP8M1/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_BuToJpsiK_Bpt0_Pthat50_TuneCUEP8M1/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_BJpsiMM_ptJpsi_00_03_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_BJpsiMM_ptJpsi_03_06_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_BJpsiMM_ptJpsi_06_09_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_BJpsiMM_ptJpsi_09_12_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_BJpsiMM_ptJpsi_12_15_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_BJpsiMM_ptJpsi_15_30_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_BJpsiMM_ptJpsi_30_inf_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'

config.Data.totalUnits = -1
config.Data.unitsPerJob = 2000
#config.Data.totalUnits = 18000
#config.Data.unitsPerJob = 9000
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/HeavyFlavourRun2/BfinderRun2/MC_official'
#config.Data.allowNonValidInputDataset = True
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/HI/Cert_262548-263757_PromptReco_HICollisions15_JSON_MuonPhys_v2.txt'
#config.Site.storageSite = 'T2_US_MIT'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_US_MIT']
