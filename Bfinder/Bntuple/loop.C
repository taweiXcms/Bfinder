#include "format.h"
#include "../interface/Bntuple.h"
#include "loop.h"

Bool_t iseos = false;
int loop(TString infile="/data/twang/BfinderRun2/DoubleMu/BfinderData_pp_20151130/finder_pp_merged.root", TString outfile="/data/wangj/Data2015/Bntuple/ntB_DoubleMu_pp_20151130.root", bool REAL=true, bool isPbPb=false, int startEntries=0, int endEntries=-1,  bool skim=false, bool gskim=true, bool testMatching=true)
{
  infile = "/data/twang/BfinderRun2/DoubleMu/BfinderData_pp_20151130/finder_pp_99_1_07c.root";
  //infile = "/afs/cern.ch/user/t/twang/work/MITHIG/HeavyFlavor/Bfinder/DfinderDev_20150813/Dev_20151130/CMSSW_7_5_5_patch4/src/test2/finder_PbPb_oldSize.root";
  //infile = "/afs/cern.ch/user/t/twang/work/MITHIG/HeavyFlavor/Bfinder/DfinderDev_20150813/Dev_20151130/CMSSW_7_5_5_patch4/src/test2/finder_PbPb.root";
  outfile = "test.root";
  REAL = true;
  iseos = false;

  cout<<endl;
  if(REAL) cout<<"--- Processing - REAL DATA"<<endl;
  else cout<<"--- Processing - MC"<<endl;

  TString ifname;
  if(iseos) ifname = Form("root://eoscms.cern.ch//eos/cms%s",infile.Data());
  else ifname = infile;
  TFile* f = TFile::Open(ifname);
  TTree* root = (TTree*)f->Get("Bfinder/root");
  TTree* hltroot = (TTree*)f->Get("hltanalysis/HltTree");
  TTree* hiroot  = (TTree*)f->Get("hiEvtAnalyzer/HiTree");

  setHltBranch(hltroot);
  if(isPbPb) setHiTreeBranch(hiroot);

  BntupleBranches     *Bntuple = new BntupleBranches;
  EvtInfoBranches     *EvtInfo = new EvtInfoBranches;
  VtxInfoBranches     *VtxInfo = new VtxInfoBranches;
  MuonInfoBranches    *MuonInfo = new MuonInfoBranches;
  TrackInfoBranches   *TrackInfo = new TrackInfoBranches;
  BInfoBranches       *BInfo = new BInfoBranches;
  GenInfoBranches     *GenInfo = new GenInfoBranches;
  EvtInfo->setbranchadd(root);
  VtxInfo->setbranchadd(root);
  MuonInfo->setbranchadd(root);
  TrackInfo->setbranchadd(root);
  BInfo->setbranchadd(root);
  GenInfo->setbranchadd(root);

  Long64_t nentries = root->GetEntries();
  if(endEntries>nentries || endEntries == -1) endEntries = nentries;
  TFile *outf = new TFile(Form("%s_Evtfrom%dto%d.root", outfile.Data(), startEntries, endEntries),"recreate");

  int ifchannel[7];
  ifchannel[0] = 1; //jpsi+Kp
  ifchannel[1] = 1; //jpsi+pi
  ifchannel[2] = 1; //jpsi+Ks(pi+,pi-)
  ifchannel[3] = 1; //jpsi+K*(K+,pi-)
  ifchannel[4] = 1; //jpsi+K*(K-,pi+)
  ifchannel[5] = 1; //jpsi+phi(K+,K-)
  ifchannel[6] = 1; //jpsi+pi pi <= psi', X(3872), Bs->J/psi f0
  
  cout<<"--- Building trees"<<endl;
  TTree* nt0 = new TTree("ntKp","");     Bntuple->buildBranch(nt0);
  TTree* nt1 = new TTree("ntpi","");     Bntuple->buildBranch(nt1);
  TTree* nt2 = new TTree("ntKs","");     Bntuple->buildBranch(nt2);
  TTree* nt3 = new TTree("ntKstar","");  Bntuple->buildBranch(nt3);
  TTree* nt5 = new TTree("ntphi","");    Bntuple->buildBranch(nt5);
  TTree* nt6 = new TTree("ntmix","");    Bntuple->buildBranch(nt6);
  TTree* ntGen = new TTree("ntGen","");  Bntuple->buildGenBranch(ntGen);
  TTree* ntHlt = hltroot->CloneTree(0);
  TTree* ntHi = hiroot->CloneTree(0);
  cout<<"--- Building trees finished"<<endl;

  cout<<"--- Check the number of events for two trees"<<endl;
  cout<<root->GetEntries()<<" "<<hltroot->GetEntries();
  if(isPbPb) cout<<" "<<hiroot->GetEntries();
  cout<<endl;
  cout<<"--- Processing events"<<endl;

  for(int i=startEntries;i<endEntries;i++)
  {
    root->GetEntry(i);
    hltroot->GetEntry(i);
    if(isPbPb) hiroot->GetEntry(i);
    if(i%100000==0) cout<<setw(8)<<i<<" / "<<(endEntries-startEntries)<<endl;
    if(testMatching)
	{
	  if(((int)Bf_HLT_Event!=EvtInfo->EvtNo||(int)Bf_HLT_Run!=EvtInfo->RunNo||(int)Bf_HLT_LumiBlock!=EvtInfo->LumiNo) || (isPbPb&&((int)Bf_HiTree_Evt!=EvtInfo->EvtNo||(int)Bf_HiTree_Run!=EvtInfo->RunNo||(int)Bf_HiTree_Lumi!=EvtInfo->LumiNo)))
	  {
	    cout<<"Error: not matched "<<i<<" | ";
	    cout<<"EvtNo("<<Bf_HLT_Event<<","<<EvtInfo->EvtNo<<") RunNo("<<Bf_HLT_Run<<","<<EvtInfo->RunNo<<") LumiNo("<<Bf_HLT_LumiBlock<<","<<EvtInfo->LumiNo<<") | EvtNo("<<Bf_HiTree_Evt<<","<<EvtInfo->EvtNo<<") RunNo("<<Bf_HiTree_Run<<","<<EvtInfo->RunNo<<") LumiNo("<<Bf_HiTree_Lumi<<","<<EvtInfo->LumiNo<<")"<<endl;
	    continue;
	  }
	}
    ntHlt->Fill();
    if(isPbPb) ntHi->Fill();
    Bntuple->makeNtuple(ifchannel, REAL, skim, EvtInfo, VtxInfo, MuonInfo, TrackInfo, BInfo, GenInfo, nt0, nt1, nt2, nt3, nt5, nt6);
    if(!REAL) Bntuple->fillGenTree(ntGen, GenInfo, gskim);
  }//entries loop
  outf->Write();
  outf->Close();
  return 1;
}
