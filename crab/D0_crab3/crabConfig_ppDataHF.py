#file path will be:
#outLFNDirBase/inputDataset/requestName/time_tag/...
from CRABClient.UserUtilities import config
config = config()
config.General.transferOutputs = True
#config.General.requestName = 'DfinderData_pp_20160326_dPt0tkPt0p5_D0Dstar'
#config.General.requestName = 'DfinderData_pp_20160327_dPt0tkPt0p5_D0Dstar'
config.General.requestName = 'DfinderData_pp_20160328_dPt0tkPt0p5_D0Dstar'
config.General.workArea = 'crab_projects'
config.JobType.psetName = 'finder_pp_75X_cfg.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.inputFiles = ['rssLimit']
config.JobType.pyCfgParams = ['noprint']
config.JobType.outputFiles = ['finder_pp.root']
config.Data.inputDataset = '/HeavyFlavor/Run2015E-PromptReco-v1/AOD'

config.Data.totalUnits = -1
config.Data.unitsPerJob = 1
#config.Data.inputDBS = 'phys03'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.splitting = 'LumiBased'
#config.Data.outLFNDirBase = '/store/user/twang/BfinderRun2'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/HeavyFlavor'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/5TeV/Cert_262081-262273_5TeV_PromptReco_Collisions15_25ns_JSON_MuonPhys_v2.txt'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/5TeV/Cert_262081-262328_5TeV_PromptReco_Collisions15_25ns_JSON.txt'
#config.Site.storageSite = 'T2_US_MIT'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_US_MIT']
