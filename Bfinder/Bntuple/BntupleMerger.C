#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "iostream"
#include "iomanip"
using namespace std;

string inFileNames[] = {
//  "ntuple_finder_PbPb_989_1_vSW.root",
//  "ntuple_finder_PbPb_995_1_5Vd.root",
	"HIOniaL1DoubleMu0/Bntuple20160610_crab_BfinderData_PbPb_HIOniaL1DoubleMu0_20160607_bPt5jpsiPt0tkPt0p8_Bp.root",
	"HIOniaL1DoubleMu0B/Bntuple20160610_crab_BfinderData_PbPb_HIOniaL1DoubleMu0B_20160607_bPt5jpsiPt0tkPt0p8_Bp.root",
	"HIOniaL1DoubleMu0C/Bntuple20160610_crab_BfinderData_PbPb_HIOniaL1DoubleMu0C_20160607_bPt5jpsiPt0tkPt0p8_Bp.root",
	"HIOniaL1DoubleMu0D/Bntuple20160610_crab_BfinderData_PbPb_HIOniaL1DoubleMu0D_20160607_bPt5jpsiPt0tkPt0p8_Bp.root"
};
int nFiles = sizeof(inFileNames)/sizeof(string);

void addtochain(TChain* root){
  for (int i=0; i<nFiles; i++) {
    root->Add(inFileNames[i].c_str());
  }	    
}

int BntupleMerger()
{
  TChain* ntKp= new TChain("ntKp");
  TChain* ntGen = new TChain("ntGen");
  TChain* ntHlt = new TChain("ntHlt");
  TChain* ntHi = new TChain("ntHi");
  TChain* ntSkim = new TChain("ntSkim");
  addtochain(ntKp);
  addtochain(ntGen);
  addtochain(ntHlt);
  addtochain(ntHi);
  addtochain(ntSkim);

  cout<<" -- Check evt no. for three trees"<<endl;
  cout<<"    "<<ntKp->GetEntries()<<", "<<ntGen->GetEntries()<<", "<<ntHlt->GetEntries()<<", "<<ntHi->GetEntries()<<", "<<ntSkim->GetEntries()<<endl;
  if(ntKp->GetEntries()!=ntHlt->GetEntries())
    {
      cout<<"    Error: Event numbers are different in three trees."<<endl;
      return 0;
    }

  int      Bsize;
  int      RunNo;
  int      EvtNo;
  int      LumiNo;
  ntKp->SetBranchAddress("Bsize",&Bsize);
  ntKp->SetBranchAddress("RunNo",&RunNo);
  ntKp->SetBranchAddress("EvtNo",&EvtNo);
  ntKp->SetBranchAddress("LumiNo",&LumiNo);
  Bool_t skimevents=false;

  TString ofname="test.root";
  TFile* outf = TFile::Open(ofname,"recreate");
  TTree* ntKp_new = ntKp->CloneTree(0);
  TTree* ntGen_new = ntGen->CloneTree(0);
  TTree* ntHlt_new = ntHlt->CloneTree(0);
  TTree* ntHi_new = ntHi->CloneTree(0);
  TTree* ntSkim_new = ntSkim->CloneTree(0);

  Int_t fCurrent = -1;
  map< pair<int, int>, int> eList;
  map< pair<int, int>, int>::iterator it;
  int nDuplicate = 0;

  Long64_t nentries = ntKp->GetEntries();
  cout<<" -- Event reading"<<endl;
  for(Long64_t i=0;i<nentries;i++)
    {
      if(i%100000==0) cout<<setiosflags(ios::left)<<"    "<<setw(8)<<i<<" / "<<nentries<<endl;

      Long64_t centry = ntKp->LoadTree(i);
      if (centry < 0) break;
      if (ntKp->GetTreeNumber() != fCurrent) {
        fCurrent = ntKp->GetTreeNumber();
        cout << "[INFO] Changed to Tree number: " << inFileNames[fCurrent] << endl;
      }

      ntKp->GetEntry(i);
      ntGen->GetEntry(i);
      ntHlt->GetEntry(i);
      ntHi->GetEntry(i);
      ntSkim->GetEntry(i);

      //if ((HLTriggers&(ULong64_t)pow(2,0))==(ULong64_t)pow(2,0)) {
      bool isDuplicate = false;
      pair<int, int> element = make_pair(EvtNo, RunNo);
      pair<pair<int, int>, int> element2 = make_pair(element, i);
      pair<map<pair<int, int>, int>::iterator, bool> result = eList.insert(element2);
      if (result.second == 0) {
        cout << " Duplicated event in File " << inFileNames[fCurrent] << " : " << RunNo << " " << EvtNo << endl;
	    isDuplicate = true;
        nDuplicate += 1;
      }
      //}
      //cout << " Duplicated event: " << runNb << " " << eventNb << endl;

	  bool flag = false;
	  if(Bsize > 0){
        flag = true;
      }
      if(!isDuplicate && (!skimevents || flag))
	  {
	    ntKp_new->Fill();
	    ntGen_new->Fill();
	    ntHlt_new->Fill(); 
	    ntHi_new->Fill(); 
	    ntSkim_new->Fill(); 
	  } 
    }
  outf->Write();
  cout<<"# of duplicate: "<<nDuplicate<<endl;
  cout<<" -- Writing new trees done"<<endl;
  outf->Close();

  return 1;  
}
