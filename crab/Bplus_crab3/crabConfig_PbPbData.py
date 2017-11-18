#file path will be:
#outLFNDirBase/inputDataset/requestName/time_tag/...
from CRABClient.UserUtilities import config
config = config()
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'BfinderData_PbPb_20160816_bPt5jpsiPt0tkPt0p8_Bp'
config.General.workArea = 'crab_projects'

config.JobType.psetName = 'finder_PbPb_75X_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['rssLimit']
config.JobType.pyCfgParams = ['noprint']
config.JobType.outputFiles = ['finder_PbPb.root']
config.JobType.maxMemoryMB = 4000

config.Data.inputDataset = '/HIOniaL1DoubleMu0/HIRun2015-PromptReco-v1/AOD'
#config.Data.inputDataset = '/HIOniaL1DoubleMu0B/HIRun2015-PromptReco-v1/AOD'
#config.Data.inputDataset = '/HIOniaL1DoubleMu0C/HIRun2015-PromptReco-v1/AOD'
#config.Data.inputDataset = '/HIOniaL1DoubleMu0D/HIRun2015-PromptReco-v1/AOD'

config.Data.totalUnits = -1
config.Data.unitsPerJob = 7000
#config.Data.unitsPerJob = 5000
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/HeavyFlavourRun2/BfinderRun2/'
#config.Data.allowNonValidInputDataset = True
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/HI/Cert_262548-263757_PromptReco_HICollisions15_JSON_MuonPhys_v2.txt'
#config.Site.storageSite = 'T2_US_MIT'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_US_MIT']
