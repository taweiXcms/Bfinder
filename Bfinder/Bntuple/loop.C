using namespace std;

#include "format.h"
#include "Bntuple.h"
#include "loop.h"

Bool_t istest = false;
int loop(TString infile="/data/twang/BfinderRun2/DoubleMu/BfinderData_pp_20151130/finder_pp_merged.root", 
         TString outfile="test.root", 
         Bool_t REAL=false, Bool_t isPbPb=false, Int_t startEntries=0, Int_t endEntries=-1,  Bool_t skim=true, Bool_t gskim=true, Bool_t checkMatching=true, Bool_t iseos=false, Bool_t SkimHLTtree=true)
{
  if(istest)
    {
      infile="/store/group/phys_heavyions/HeavyFlavourRun2/BfinderRun2/MC_official/Pythia8_BuToJpsiK_Bpt0_Pthat15_TuneCUEP8M1/crab_BfinderMC_pp_Pythia8_BuToJpsiK_Bpt0_Pthat15_TuneCUEP8M1_20160816_bPt0jpsiPt0tkPt0p5_Bp/160816_214911/0000/finder_pp_5.root";
      outfile="test.root";
      REAL=false;
      isPbPb=false;
      skim=false;
      checkMatching=true;
      iseos=true;
    }

  cout<<endl;
  if(REAL) cout<<"--- Processing - REAL DATA";
  else cout<<"--- Processing - MC";
  if(isPbPb) cout<<" - PbPb";
  else cout<<" - pp";
  cout<<endl;

  TString ifname;
  if(iseos) ifname = Form("root://eoscms.cern.ch//eos/cms%s",infile.Data());
  else ifname = infile;
  TFile* f = TFile::Open(ifname);
  TTree* root = (TTree*)f->Get("Bfinder/root");
  TTree* hltroot = (TTree*)f->Get("hltanalysis/HltTree");
  TTree* skimroot = (TTree*)f->Get("skimanalysis/HltTree");
  TTree* hiroot = (TTree*)f->Get("hiEvtAnalyzer/HiTree");

  BntupleBranches     *Bntuple = new BntupleBranches;
  EvtInfoBranches     *EvtInfo = new EvtInfoBranches;
  VtxInfoBranches     *VtxInfo = new VtxInfoBranches;
  MuonInfoBranches    *MuonInfo = new MuonInfoBranches;
  TrackInfoBranches   *TrackInfo = new TrackInfoBranches;
  BInfoBranches       *BInfo = new BInfoBranches;
  GenInfoBranches     *GenInfo = new GenInfoBranches;

  if(SkimHLTtree) SetHlttreestatus(hltroot, isPbPb);
  setHltBranch(hltroot);
  setHiTreeBranch(hiroot);

  EvtInfo->setbranchadd(root);
  VtxInfo->setbranchadd(root);
  MuonInfo->setbranchadd(root);
  TrackInfo->setbranchadd(root);
  BInfo->setbranchadd(root);
  GenInfo->setbranchadd(root);

  Long64_t nentries = root->GetEntries();
  if(endEntries>nentries || endEntries == -1) endEntries = nentries;
  TFile *outf = TFile::Open(Form("%s", outfile.Data()),"recreate");

  Int_t ifchannel[8];
  ifchannel[0] = 1; //jpsi+Kp
  ifchannel[1] = 1; //jpsi+pi
  ifchannel[2] = 1; //jpsi+Ks(pi+,pi-)
  ifchannel[3] = 1; //jpsi+K*(K+,pi-)
  ifchannel[4] = 1; //jpsi+K*(K-,pi+)
  ifchannel[5] = 1; //jpsi+phi(K+,K-)
  ifchannel[6] = 1; //jpsi+pi pi <= psi', X(3872), Bs->J/psi f0
  ifchannel[7] = 1; //inclusive jpsi
  
  cout<<"--- Building trees"<<endl;
  TTree* nt0 = new TTree("ntKp","");      Bntuple->buildBranch(nt0);
  TTree* nt1 = new TTree("ntpi","");      Bntuple->buildBranch(nt1);
  TTree* nt2 = new TTree("ntKs","");      Bntuple->buildBranch(nt2);
  TTree* nt3 = new TTree("ntKstar","");   Bntuple->buildBranch(nt3);
  TTree* nt5 = new TTree("ntphi","");     Bntuple->buildBranch(nt5);
  TTree* nt6 = new TTree("ntmix","");     Bntuple->buildBranch(nt6);
  TTree* nt7 = new TTree("ntJpsi","");    Bntuple->buildBranch(nt7,true);
  TTree* ntGen = new TTree("ntGen","");   Bntuple->buildGenBranch(ntGen);
  TTree* ntHlt = hltroot->CloneTree(0);
  ntHlt->SetName("ntHlt");
  TTree* ntSkim = skimroot->CloneTree(0);
  ntSkim->SetName("ntSkim");
  TTree* ntHi = hiroot->CloneTree(0);
  ntHi->SetName("ntHi");
  cout<<"--- Building trees finished"<<endl;

  cout<<"--- Check the number of events for four trees"<<endl;
  cout<<root->GetEntries()<<" "<<hltroot->GetEntries()<<" "<<hiroot->GetEntries();
  cout<<" "<<skimroot->GetEntries()<<endl;
  cout<<endl;

  cout<<"--- Processing events"<<endl;
  for(Int_t i=startEntries;i<endEntries;i++)
    {
      root->GetEntry(i);
      hltroot->GetEntry(i);
      skimroot->GetEntry(i);
      hiroot->GetEntry(i);
      
      if(i%100000==0) cout<<setw(7)<<i<<" / "<<endEntries<<endl;
      if(checkMatching)
	{
          if(((int)Bf_HLT_Event!=EvtInfo->EvtNo||(int)Bf_HLT_Run!=EvtInfo->RunNo||(int)Bf_HLT_LumiBlock!=EvtInfo->LumiNo) || 
             ((int)Bf_HiTree_Evt!=EvtInfo->EvtNo||(int)Bf_HiTree_Run!=EvtInfo->RunNo||(int)Bf_HiTree_Lumi!=EvtInfo->LumiNo))
            {
              cout<<"Error: not matched "<<i<<" | (Hlt,Bfr,Hi) | ";
              cout<<"EvtNo("<<Bf_HLT_Event<<","<<EvtInfo->EvtNo<<","<<Bf_HiTree_Evt<<") ";
              cout<<"RunNo("<<Bf_HLT_Run<<","<<EvtInfo->RunNo<<","<<Bf_HiTree_Run<<") ";
              cout<<"LumiNo("<<Bf_HLT_LumiBlock<<","<<EvtInfo->LumiNo<<","<<Bf_HiTree_Lumi<<")"<<endl;
              continue;
            }
	}
      ntHlt->Fill();
      ntSkim->Fill();
      ntHi->Fill();
      Bntuple->makeNtuple(ifchannel, REAL, skim, EvtInfo, VtxInfo, MuonInfo, TrackInfo, BInfo, GenInfo, nt0, nt1, nt2, nt3, nt5, nt6, nt7);
      if(!REAL) Bntuple->fillGenTree(ntGen, GenInfo, gskim);
    }
  outf->Write();
  cout<<"--- Writing finished"<<endl;
  outf->Close();

  cout<<"--- In/Output files"<<endl;
  cout<<ifname<<endl;
  cout<<outfile<<endl;
  cout<<endl;

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc==3)
    {
      loop(argv[1], argv[2]);
    }
  else
    {
      std::cout << "Usage: mergeForest <input_collection> <output_file>" << std::endl;
      return 1;
    }
  return 0;
}


