#file path will be:
#outLFNDirBase/inputDataset/requestName/time_tag/...
from CRABClient.UserUtilities import config
config = config()
config.General.transferOutputs = True
config.General.requestName = 'DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0'

#config.General.workArea = 'crab_projects_PbPbD0pt0'
#config.General.workArea = 'crab_projects_PbPbD0pt5'
#config.General.workArea = 'crab_projects_PbPbD0pt10'
#config.General.workArea = 'crab_projects_PbPbD0pt15'
#config.General.workArea = 'crab_projects_PbPbD0pt30'
#config.General.workArea = 'crab_projects_PbPbD0pt50'
#config.General.workArea = 'crab_projects_PbPbD0pt80'
#config.General.workArea = 'crab_projects_PbPbD0pt120'
#config.General.workArea = 'crab_projects_PbPbD0pt170'

#config.General.workArea = 'crab_projects_PbPbNPD0pt0'
#config.General.workArea = 'crab_projects_PbPbNPD0pt5'
#config.General.workArea = 'crab_projects_PbPbNPD0pt10'
#config.General.workArea = 'crab_projects_PbPbNPD0pt15'
#config.General.workArea = 'crab_projects_PbPbNPD0pt30'
#config.General.workArea = 'crab_projects_PbPbNPD0pt50'
#config.General.workArea = 'crab_projects_PbPbNPD0pt80'
#config.General.workArea = 'crab_projects_PbPbNPD0pt120'
#config.General.workArea = 'crab_projects_PbPbNPD0pt170'

config.JobType.psetName = 'finder_PbPb_75X_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['rssLimit']
config.JobType.pyCfgParams = ['noprint']
config.JobType.outputFiles = ['finder_PbPb.root']

#config.Data.inputDataset = '/Pythia8_prompt_D0pt0p0_Pthat0_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_prompt_D0pt0p0_Pthat5_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_prompt_D0pt0p0_Pthat10_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_prompt_D0pt0p0_Pthat15_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_prompt_D0pt0p0_Pthat30_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_prompt_D0pt0p0_Pthat50_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_prompt_D0pt0p0_Pthat80_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_prompt_D0pt0p0_Pthat120_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_prompt_D0pt0p0_Pthat170_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'

#
#config.Data.inputDataset = '/Pythia8_nonprompt_D0pt0p0_Pthat0_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_nonprompt_D0pt0p0_Pthat5_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_nonprompt_D0pt0p0_Pthat10_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_nonprompt_D0pt0p0_Pthat15_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_nonprompt_D0pt0p0_Pthat30_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_nonprompt_D0pt0p0_Pthat50_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_nonprompt_D0pt0p0_Pthat80_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_nonprompt_D0pt0p0_Pthat120_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'
#config.Data.inputDataset = '/Pythia8_nonprompt_D0pt0p0_Pthat170_Hydjet_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1/AODSIM'

config.Data.totalUnits = -1
config.Data.unitsPerJob = 200
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/MC_official'
#config.Data.allowNonValidInputDataset = True
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/5TeV/Cert_262081-262273_5TeV_PromptReco_Collisions15_25ns_JSON_MuonPhys_v2.txt'
#config.Site.storageSite = 'T2_US_MIT'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_US_MIT']
