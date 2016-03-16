#file path will be:
#outLFNDirBase/inputDataset/requestName/time_tag/...
from CRABClient.UserUtilities import config
config = config()
config.General.transferOutputs = True
config.General.requestName = 'DfinderMC_pp_20160302_dPt0tkPt0p5_D0Dstar'
config.General.workArea = 'crab_projects_ppD0pt0'
#config.General.workArea = 'crab_projects_ppD0pt5'
#config.General.workArea = 'crab_projects_ppD0pt10'
#config.General.workArea = 'crab_projects_ppD0pt15'
#config.General.workArea = 'crab_projects_ppD0pt30'
#config.General.workArea = 'crab_projects_ppD0pt50'
#config.General.workArea = 'crab_projects_ppD0pt80'
#config.General.workArea = 'crab_projects_ppD0pt120'
#config.General.workArea = 'crab_projects_ppD0pt170'
config.JobType.psetName = 'finder_pp_75X_cfg.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.inputFiles = ['rssLimit']
config.JobType.pyCfgParams = ['noprint']
config.JobType.outputFiles = ['finder_pp.root']

config.Data.inputDataset = '/pp_Pythia8_prompt_D0pt0p0_Pthat0_TuneCUETP8M1_5020GeV_evtgen130_GEN_SIM_20160229/twang-pp_Pythia8_prompt_D0pt0p0_Pthat0_TuneCUETP8M1_5020GeV_evtgen130_step3_20160229-5f793ffd6c0245cf30cde3c9fc53e2e3/USER'
#config.Data.inputDataset = '/pp_Pythia8_prompt_D0pt0p0_Pthat5_TuneCUETP8M1_5020GeV_evtgen130_GEN_SIM_20160229/twang-pp_Pythia8_prompt_D0pt0p0_Pthat5_TuneCUETP8M1_5020GeV_evtgen130_step3_20160229-5f793ffd6c0245cf30cde3c9fc53e2e3/USER'
#config.Data.inputDataset = '/pp_Pythia8_prompt_D0pt0p0_Pthat10_TuneCUETP8M1_5020GeV_evtgen130_GEN_SIM_20160229/twang-pp_Pythia8_prompt_D0pt0p0_Pthat10_TuneCUETP8M1_5020GeV_evtgen130_step3_20160229-5f793ffd6c0245cf30cde3c9fc53e2e3/USER'
#config.Data.inputDataset = '/pp_Pythia8_prompt_D0pt0p0_Pthat15_TuneCUETP8M1_5020GeV_evtgen130_GEN_SIM_20160229/twang-pp_Pythia8_prompt_D0pt0p0_Pthat15_TuneCUETP8M1_5020GeV_evtgen130_step3_20160229-5f793ffd6c0245cf30cde3c9fc53e2e3/USER'
#config.Data.inputDataset = '/pp_Pythia8_prompt_D0pt0p0_Pthat30_TuneCUETP8M1_5020GeV_evtgen130_GEN_SIM_20160229/twang-pp_Pythia8_prompt_D0pt0p0_Pthat30_TuneCUETP8M1_5020GeV_evtgen130_step3_20160229-5f793ffd6c0245cf30cde3c9fc53e2e3/USER'
#config.Data.inputDataset = '/pp_Pythia8_prompt_D0pt0p0_Pthat50_TuneCUETP8M1_5020GeV_evtgen130_GEN_SIM_20160229/twang-pp_Pythia8_prompt_D0pt0p0_Pthat50_TuneCUETP8M1_5020GeV_evtgen130_step3_20160229-5f793ffd6c0245cf30cde3c9fc53e2e3/USER'
#config.Data.inputDataset = '/pp_Pythia8_prompt_D0pt0p0_Pthat80_TuneCUETP8M1_5020GeV_evtgen130_GEN_SIM_20160229/twang-pp_Pythia8_prompt_D0pt0p0_Pthat80_TuneCUETP8M1_5020GeV_evtgen130_step3_20160229-5f793ffd6c0245cf30cde3c9fc53e2e3/USER'
#config.Data.inputDataset = '/pp_Pythia8_prompt_D0pt0p0_Pthat120_TuneCUETP8M1_5020GeV_evtgen130_GEN_SIM_20160229/twang-pp_Pythia8_prompt_D0pt0p0_Pthat120_TuneCUETP8M1_5020GeV_evtgen130_step3_20160229-5f793ffd6c0245cf30cde3c9fc53e2e3/USER'
#config.Data.inputDataset = '/pp_Pythia8_prompt_D0pt0p0_Pthat170_TuneCUETP8M1_5020GeV_evtgen130_GEN_SIM_20160229/twang-pp_Pythia8_prompt_D0pt0p0_Pthat170_TuneCUETP8M1_5020GeV_evtgen130_step3_20160229-5f793ffd6c0245cf30cde3c9fc53e2e3/USER'

config.Data.totalUnits = -1
config.Data.unitsPerJob = 2000
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/MC'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/5TeV/Cert_262081-262273_5TeV_PromptReco_Collisions15_25ns_JSON_MuonPhys_v2.txt'
#config.Site.storageSite = 'T2_US_MIT'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_US_MIT']
