#include "format.h"
#include "loop.h"

int loop(TString infile="root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/Run2015E/HIMinimumBias2/Merged/HIForestExpress_run262620.root",
         TString outfile="./ntD_HIForestExpress_run262620", bool REAL=true, bool isPbPb=true, int startEntries=0, int endEntries=-1, bool skim=false, bool gskim=true)
{
  infile="/afs/cern.ch/user/t/twang/work/MITHIG/HeavyFlavor/Bfinder/DfinderDev_20150813/Dev_20151202/CMSSW_7_5_5_patch4/src/test3/finder_PbPb.root";
  outfile="test_noskim";
  REAL=false;
//skim = true;

  cout<<endl;
  if(REAL) cout<<"--- Processing - REAL DATA"<<endl;
  else cout<<"--- Processing - MC"<<endl;
  
  TFile* f = TFile::Open(infile);
  TTree* root = (TTree*)f->Get("Dfinder/root");
  TTree* hltroot = (TTree*)f->Get("hltanalysis/HltTree");
  TTree* skimroot = (TTree*)f->Get("skimanalysis/HltTree");
  TTree* hiroot;
  hiroot = (TTree*)f->Get("hiEvtAnalyzer/HiTree");

  DntupleBranches     *Dntuple = new DntupleBranches;
  EvtInfoBranches     *EvtInfo = new EvtInfoBranches;
  VtxInfoBranches     *VtxInfo = new VtxInfoBranches;
  TrackInfoBranches   *TrackInfo = new TrackInfoBranches;
  DInfoBranches       *DInfo = new DInfoBranches;
  GenInfoBranches     *GenInfo = new GenInfoBranches;

  setHltTreeBranch(hltroot);
  if(isPbPb) setHiTreeBranch(hiroot);

  EvtInfo->setbranchadd(root);
  VtxInfo->setbranchadd(root);
  TrackInfo->setbranchadd(root);
  DInfo->setbranchadd(root);
  GenInfo->setbranchadd(root);

  Long64_t nentries = root->GetEntries();
  if( endEntries > nentries )
      endEntries = nentries;
  if( endEntries == -1 )  endEntries = nentries;
  TFile *outf = new TFile(Form("%s_Evtfrom%dto%d.root", outfile.Data(), startEntries, endEntries),"recreate");

  int isDchannel[6];
  isDchannel[0] = 1; //k+pi-
  isDchannel[1] = 1; //k-pi+
  isDchannel[2] = 0; //k-pi+pi+
  isDchannel[3] = 0; //k+pi-pi-
  isDchannel[4] = 0; //k-pi-pi+pi+
  isDchannel[5] = 0; //k+pi+pi-pi-

  cout<<"--- Building trees"<<endl;
  TTree* ntD1 = new TTree("ntDkpi","");       Dntuple->buildDBranch(ntD1);
  TTree* ntD2 = new TTree("ntDkpipi","");     Dntuple->buildDBranch(ntD2);
  TTree* ntD3 = new TTree("ntDkpipipi","");   Dntuple->buildDBranch(ntD3);
  TTree* ntGen = new TTree("ntGen","");       Dntuple->buildGenBranch(ntGen);
  TTree* ntHlt = hltroot->CloneTree(0);
  ntHlt->SetName("ntHlt");
  TTree* ntSkim = skimroot->CloneTree(0);
  ntSkim->SetName("ntSkim");
  TTree* ntHi;
  ntHi = hiroot->CloneTree(0);
  cout<<"--- Building trees finished"<<endl;

  int flagEvt=0, offsetHltTree=0;
  cout<<"--- Check the number of events for three trees"<<endl;
  cout<<root->GetEntries()<<" "<<hltroot->GetEntries();
  if(isPbPb) cout<<" "<<hiroot->GetEntries();
  cout<<endl;
  cout<<"--- Processing events"<<endl;
  for(int i=startEntries;i<endEntries;i++)
    {
      root->GetEntry(i);
      hltroot->GetEntry(i);
	  skimroot->GetEntry(i);
      if(isPbPb) hiroot->GetEntry(i);
      if(i%100000==0) cout<<setw(7)<<i<<" / "<<nentries<<endl;
      if((int)Df_HLT_Event!=EvtInfo->EvtNo||Df_HLT_Run!=EvtInfo->RunNo||Df_HLT_LumiBlock!=EvtInfo->LumiNo)
	{
	  if(!isPbPb||(isPbPb&&(Df_HiTree_Evt!=EvtInfo->EvtNo||Df_HiTree_Run!=EvtInfo->RunNo||Df_HiTree_Lumi!=EvtInfo->LumiNo)))
	    {
	      cout<<"Error: not matched "<<i<<" | ";
	      cout<<Df_HLT_Event<<","<<EvtInfo->EvtNo<<"   "<<Df_HLT_Run<<","<<EvtInfo->RunNo<<"   "<<Df_HLT_LumiBlock<<","<<EvtInfo->LumiNo<<" | "<<Df_HiTree_Evt<<","<<EvtInfo->EvtNo<<"   "<<Df_HiTree_Run<<","<<EvtInfo->RunNo<<"   "<<Df_HiTree_Lumi<<","<<EvtInfo->LumiNo<<endl;
	      continue;
	    }
	}

      ntHlt->Fill();
	  ntSkim->Fill();
      if(isPbPb) ntHi->Fill();

      Dntuple->makeDNtuple(isDchannel, REAL, skim, EvtInfo, VtxInfo, TrackInfo, DInfo, GenInfo, ntD1, ntD2, ntD3);
      if(!REAL) Dntuple->fillDGenTree(ntGen, GenInfo, gskim);
    }
  outf->Write();
  cout<<"--- Writing finished"<<endl;
  outf->Close();

  return 1;
}
