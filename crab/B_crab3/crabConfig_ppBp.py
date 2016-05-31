#file path will be:
#outLFNDirBase/inputDataset/requestName/time_tag/...
from CRABClient.UserUtilities import config
config = config()
config.General.transferOutputs = True
config.General.requestName = 'BfinderMC_pp_20160420_bPt0jpsiPt0tkPt0p5_BpB0BsX'
#config.General.workArea = 'crab_projects_ppBupt5'
#config.General.workArea = 'crab_projects_ppBupt15'
#config.General.workArea = 'crab_projects_ppBupt30'
config.General.workArea = 'crab_projects_ppBupt50'
config.JobType.psetName = 'finder_pp_75X_cfg.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.inputFiles = ['rssLimit']
config.JobType.pyCfgParams = ['noprint']
config.JobType.outputFiles = ['finder_pp.root']

#config.Data.inputDataset = '/pp_Pythia8_BuToJpsiK_Bpt5p0_Pthat5_TuneCUEP8M1_5020GeV_GEN_SIM_20151212_v2/twang-pp_Pythia8_BuToJpsiK_Bpt5p0_Pthat5_TuneCUEP8M1_5020GeV_step3_20151212-e70ee3b406076082d57525c487e745ee/USER'
#config.Data.inputDataset = '/pp_Pythia8_BuToJpsiK_Bpt15p0_Pthat15_TuneCUEP8M1_5020GeV_GEN_SIM_20151212_v2/twang-pp_Pythia8_BuToJpsiK_Bpt15p0_Pthat15_TuneCUEP8M1_5020GeV_step3_20151212-e70ee3b406076082d57525c487e745ee/USER'
#config.Data.inputDataset = '/pp_Pythia8_BuToJpsiK_Bpt30p0_Pthat30_TuneCUEP8M1_5020GeV_GEN_SIM_20151212_v2/twang-pp_Pythia8_BuToJpsiK_Bpt30p0_Pthat30_TuneCUEP8M1_5020GeV_step3_20151212-e70ee3b406076082d57525c487e745ee/USER'
config.Data.inputDataset = '/pp_Pythia8_BuToJpsiK_Bpt50p0_Pthat50_TuneCUEP8M1_5020GeV_GEN_SIM_20151212_v2/twang-pp_Pythia8_BuToJpsiK_Bpt50p0_Pthat50_TuneCUEP8M1_5020GeV_step3_20151212-e70ee3b406076082d57525c487e745ee/USER'

config.Data.totalUnits = -1
config.Data.unitsPerJob = 1000
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/HeavyFlavourRun2/BfinderRun2/MC'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/5TeV/Cert_262081-262273_5TeV_PromptReco_Collisions15_25ns_JSON_MuonPhys_v2.txt'
#config.Site.storageSite = 'T2_US_MIT'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_US_MIT']
