#file path will be:
#outLFNDirBase/inputDataset/requestName/time_tag/...
from CRABClient.UserUtilities import config
config = config()
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'BfinderData_pp_20160606_bPt0jpsiPt0tkPt0p5_Bp'
config.General.workArea = 'crab_projects'

config.JobType.psetName = 'finder_pp_75X_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['noprint']
config.JobType.outputFiles = ['finder_pp.root']

config.Data.inputDataset = '/DoubleMu/Run2015E-PromptReco-v1/AOD'

config.Data.totalUnits = -1
config.Data.unitsPerJob = 30000
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/HeavyFlavourRun2/BfinderRun2/'
#config.Data.allowNonValidInputDataset = True
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/5TeV/Cert_262081-262328_5TeV_PromptReco_Collisions15_25ns_JSON_MuonPhys.txt'
#config.Site.storageSite = 'T2_US_MIT'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_US_MIT']
