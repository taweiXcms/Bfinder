#

from CRABClient.UserUtilities import config
config = config()
config.General.transferOutputs = True
config.General.requestName = 'DfinderData_PAMinimumBias5_pPb_20170328_PARun2016B_PromptReco_v1_D0_dPt0tkPt0p2'
config.General.workArea = 'crab_projects'
config.JobType.psetName = 'finder_pPb_80X_cfg.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.inputFiles = ['rssLimit']
config.JobType.pyCfgParams = ['noprint']
config.JobType.outputFiles = ['finder_pPb.root']
config.Data.inputDataset = '/PAMinimumBias5/PARun2016B-PromptReco-v1/AOD'

config.Data.totalUnits = -1
config.Data.unitsPerJob = 40000
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'LumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/wangj/DfinderData2016'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/5TeV/Cert_262081-262328_5TeV_PromptReco_Collisions15_25ns_JSON.txt'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/HI/Cert_285090-285383_HI5TeV_PromptReco_pPb_Collisions16_JSON_noL1T.txt'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_US_MIT']
#config.Site.ignoreGlobalBlacklist = True
