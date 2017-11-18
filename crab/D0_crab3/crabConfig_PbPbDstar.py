#file path will be:
#outLFNDirBase/inputDataset/requestName/time_tag/...
from CRABClient.UserUtilities import config
config = config()
config.General.transferOutputs = True
config.General.requestName = 'DfinderMC_PbPb_20160215_dPt0tkPt2p5_D0Dstar3p5p'
config.General.workArea = 'crab_projects'
config.JobType.psetName = 'finder_PbPb_75X_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['rssLimit']
config.JobType.pyCfgParams = ['noprint']
config.JobType.outputFiles = ['finder_PbPb.root']

config.Data.inputDataset = '/Pythia8D0kpi_Dstarpt5p0_Pthat5_TuneCUETP8M1_5020GeV_GEN_SIM_PU_20151212/twang-Pythia8D0kpi_Dstarpt5p0_Pthat5_TuneCUETP8M1_5020GeV_step3_20151221-bafe009c9accba23246002568a520694/USER'
#config.Data.inputDataset = '/Pythia8D0kpi_Dstarpt10p0_Pthat10_TuneCUETP8M1_5020GeV_GEN_SIM_PU_20151212/twang-Pythia8D0kpi_Dstarpt10p0_Pthat10_TuneCUETP8M1_5020GeV_step3_20151221-bafe009c9accba23246002568a520694/USER'
#config.Data.inputDataset = '/Pythia8D0kpi_Dstarpt15p0_Pthat15_TuneCUETP8M1_5020GeV_GEN_SIM_PU_20151212/twang-Pythia8D0kpi_Dstarpt15p0_Pthat15_TuneCUETP8M1_5020GeV_step3_20151221-bafe009c9accba23246002568a520694/USER'
#config.Data.inputDataset = '/Pythia8D0kpi_Dstarpt30p0_Pthat30_TuneCUETP8M1_5020GeV_GEN_SIM_PU_20151212/twang-Pythia8D0kpi_Dstarpt30p0_Pthat30_TuneCUETP8M1_5020GeV_step3_20151221-bafe009c9accba23246002568a520694/USER'
#config.Data.inputDataset = '/Pythia8D0kpi_Dstarpt50p0_Pthat50_TuneCUETP8M1_5020GeV_GEN_SIM_PU_20151212/twang-Pythia8D0kpi_Dstarpt50p0_Pthat50_TuneCUETP8M1_5020GeV_step3_20151221-bafe009c9accba23246002568a520694/USER'
#config.Data.inputDataset = '/Pythia8D0kpipipi_Dstarpt5p0_Pthat5_TuneCUETP8M1_5020GeV_GEN_SIM_PU_20151212/twang-Pythia8D0kpipipi_Dstarpt5p0_Pthat5_TuneCUETP8M1_5020GeV_step3_20151221-bafe009c9accba23246002568a520694/USER'
#config.Data.inputDataset = '/Pythia8D0kpipipi_Dstarpt10p0_Pthat10_TuneCUETP8M1_5020GeV_GEN_SIM_PU_20151212/twang-Pythia8D0kpipipi_Dstarpt10p0_Pthat10_TuneCUETP8M1_5020GeV_step3_20151221-bafe009c9accba23246002568a520694/USER'
#config.Data.inputDataset = '/Pythia8D0kpipipi_Dstarpt15p0_Pthat15_TuneCUETP8M1_5020GeV_GEN_SIM_PU_20151212/twang-Pythia8D0kpipipi_Dstarpt15p0_Pthat15_TuneCUETP8M1_5020GeV_step3_20151221-bafe009c9accba23246002568a520694/USER'
#config.Data.inputDataset = '/Pythia8D0kpipipi_Dstarpt30p0_Pthat30_TuneCUETP8M1_5020GeV_GEN_SIM_PU_20151212/twang-Pythia8D0kpipipi_Dstarpt30p0_Pthat30_TuneCUETP8M1_5020GeV_step3_20151221-bafe009c9accba23246002568a520694/USER'
#config.Data.inputDataset = '/Pythia8D0kpipipi_Dstarpt50p0_Pthat50_TuneCUETP8M1_5020GeV_GEN_SIM_PU_20151212/twang-Pythia8D0kpipipi_Dstarpt50p0_Pthat50_TuneCUETP8M1_5020GeV_step3_20151221-bafe009c9accba23246002568a520694/USER'


config.Data.totalUnits = -1
config.Data.unitsPerJob = 1000
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/HeavyFlavourRun2/DfinderRun2/MC'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/5TeV/Cert_262081-262273_5TeV_PromptReco_Collisions15_25ns_JSON_MuonPhys_v2.txt'
#config.Site.storageSite = 'T2_US_MIT'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_US_MIT']
