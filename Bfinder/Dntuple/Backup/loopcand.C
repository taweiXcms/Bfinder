using namespace std;

#include "loop.h"
#include "format.h"
#include "DntupleCand.h"

Bool_t iseos = false;
Bool_t istest = false;
int loopcand(TString infile="root://eoscms//eos/cms//store/group/phys_heavyions/jisun/PbPb2015_HeavyFlavor/PbPb_2015_HIHardProbes_AOD_tkpt0p8_D0pt0p8_D3d2_Prob0p05_alpha0p3_Dstarpt6_1210/finder_962.root",
         TString outfile="testCand.root", Bool_t REAL=true, Bool_t isPbPb=false, Int_t startEntries=0, Int_t endEntries=-1, Bool_t skim=false, Bool_t gskim=true, Bool_t checkMatching=true)
{
  if(istest)
    {
      //this testing sample doesn't have skimtree and hitree
      infile="/data/twang/DfinderRun2/Pythia8_5020GeV_DstarD0kpipipi_755patch3_GEN_SIM_PU_20151120/crab_DfinderMC_Dstar5p_tkPt2_20151126/151127_005816/0000/finder_PbPb_253.root";
      outfile="test.root";
      REAL=false;
      isPbPb=false;
      iseos=false;
      checkMatching=true;
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
  TTree* root = (TTree*)f->Get("Dfinder/root");
  TTree* hltroot = (TTree*)f->Get("hltanalysis/HltTree");
  TTree* skimroot = (TTree*)f->Get("skimanalysis/HltTree");
  TTree* hiroot = (TTree*)f->Get("hiEvtAnalyzer/HiTree");

  DntupleBranchesCand *DntupleCand = new DntupleBranchesCand;
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
  if(endEntries>nentries || endEntries == -1) endEntries = nentries;
  TFile *outf = new TFile(Form("%s", outfile.Data()),"recreate");

  int isDchannel[12];
  isDchannel[0] = 1; //k+pi-
  isDchannel[1] = 1; //k-pi+
  isDchannel[2] = 1; //k-pi+pi+
  isDchannel[3] = 1; //k+pi-pi-
  isDchannel[4] = 1; //k-pi-pi+pi+
  isDchannel[5] = 1; //k+pi+pi-pi-
  isDchannel[6] = 1; 
  isDchannel[7] = 1; 
  isDchannel[8] = 1; 
  isDchannel[9] = 1; 
  isDchannel[10] = 1; 
  isDchannel[11] = 1;

  cout<<"--- Building trees"<<endl;
  TTree* ntD1 = new TTree("ntDkpi","");         DntupleCand->buildDBranch(ntD1);
  TTree* ntD2 = new TTree("ntDkpipi","");       DntupleCand->buildDBranch(ntD2);
  TTree* ntD3 = new TTree("ntDkpipipi","");     DntupleCand->buildDBranch(ntD3);
  TTree* ntD4 = new TTree("ntDPhikkpi","");     DntupleCand->buildDBranch(ntD4);
  TTree* ntD5 = new TTree("ntDD0kpipi","");     DntupleCand->buildDBranch(ntD5);
  TTree* ntD6 = new TTree("ntDD0kpipipipi",""); DntupleCand->buildDBranch(ntD6);
  TTree* ntGen = new TTree("ntGen","");         DntupleCand->buildGenBranch(ntGen);
  /*
  TTree* ntHlt = hltroot->CloneTree(0);
  ntHlt->SetName("ntHlt");
  TTree* ntSkim = skimroot->CloneTree(0);
  ntSkim->SetName("ntSkim");
  TTree* ntHi = hiroot->CloneTree(0);
  ntHi->SetName("ntHi");
  */
  cout<<"--- Building trees finished"<<endl;

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
      if(i%100000==0) cout<<setw(7)<<i<<" / "<<endEntries<<endl;
      if(checkMatching)
        {
          if(((int)Df_HLT_Event!=EvtInfo->EvtNo||(int)Df_HLT_Run!=EvtInfo->RunNo||(int)Df_HLT_LumiBlock!=EvtInfo->LumiNo) || 
             (isPbPb&&((int)Df_HiTree_Evt!=EvtInfo->EvtNo||(int)Df_HiTree_Run!=EvtInfo->RunNo||(int)Df_HiTree_Lumi!=EvtInfo->LumiNo)))
            {
              cout<<"Error: not matched "<<i<<" | (Hlt,Dfr,Hi) | ";
              cout<<"EvtNo("<<Df_HLT_Event<<","<<EvtInfo->EvtNo<<","<<Df_HiTree_Evt<<") ";
              cout<<"RunNo("<<Df_HLT_Run<<","<<EvtInfo->RunNo<<","<<Df_HiTree_Run<<") ";
              cout<<"LumiNo("<<Df_HLT_LumiBlock<<","<<EvtInfo->LumiNo<<","<<Df_HiTree_Lumi<<")"<<endl;
              continue;
            }
        }
      /*
      ntHlt->Fill();
      ntSkim->Fill();
      if(isPbPb) ntHi->Fill();
      */
      DntupleCand->makeDNtuple(isDchannel, REAL, skim, EvtInfo, VtxInfo, TrackInfo, DInfo, GenInfo, ntD1, ntD2, ntD3, ntD4, ntD5, ntD6);
      if(!REAL)
        {
          DntupleCand->fillDGenTree(ntGen, GenInfo, gskim);
        }
    }
  outf->Write();
  cout<<"--- Writing finished"<<endl;
  outf->Close();

  cout<<"--- In/Output files"<<endl;
  cout<<ifname<<endl;
  cout<<outfile<<endl;
  cout<<endl;

  return 1;
}

int main(int argc, char *argv[])
{
  if((argc != 3) && (argc != 4))
  {
    std::cout << "Usage: mergeForest <input_collection> <output_file>" << std::endl;
    return 1;
  }
  
  if(argc == 3)
    loopcand(argv[1], argv[2]);
  //else if (argc == 4)
  //  loop(argv[1], argv[2], argv[3]);
  return 0;
}
