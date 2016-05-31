using namespace std;
#include "skimTree_Dzero.h"

int skimTree_Dzero(
//                   TString ifname="root://eoscms//eos/cms//store/group/phys_heavyions/wangj/DntupleData/ntD_EvtBase_20160425_HighPtJet80_DfinderData_pp_20160329_dPt0tkPt0p5_D0Dstar/ntuple_69.root",
                   TString ifname="root://eoscms//eos/cms//store/group/phys_heavyions/wangj/DntupleData/ntD_EvtBase_20160405_HIHardProbes_DfinderData_PbPb_20160402_dPt0tkPt2p5_D0Dstar3p5p_FINALJSON/ntuple_987.root",
                   TString ofname="test.root",
                   Int_t isPbPb=1,
                   Bool_t skimbranch=true,
                   Bool_t skimevents=true)
{
  cout<<endl;
  cout<<" -- Checking if input and output files are same"<<endl;
  if(ifname==ofname)
    {
      cout<<"    Error: Input file will be overwritten."<<endl;
      return 0;
    }
  TFile* inf = TFile::Open(ifname);
  TTree* ntDkpi = (TTree*)inf->Get("ntDkpi");
  TTree* ntGen = (TTree*)inf->Get("ntGen");
  TTree* ntHlt = (TTree*)inf->Get("ntHlt");
  TTree* ntHi = (TTree*)inf->Get("ntHi");
  TTree* ntSkim = (TTree*)inf->Get("ntSkim");

  Int_t Dsize_ntDkpi;      ntDkpi->SetBranchAddress("Dsize",&Dsize_ntDkpi);

  if(skimbranch)
    {
      SelectRecoBranches(ntDkpi);
      SelectHltBranches(ntHlt,isPbPb);
    }
  TFile* outf = TFile::Open(ofname,"recreate");
  TTree* ntDkpi_new = ntDkpi->CloneTree(0);
  TTree* ntGen_new = ntGen->CloneTree(0);
  TTree* ntHlt_new = ntHlt->CloneTree(0);
  TTree* ntHi_new = ntHi->CloneTree(0);
  TTree* ntSkim_new = ntSkim->CloneTree(0);

  cout<<" -- Check evt no. for three trees"<<endl;
  cout<<"    "<<ntDkpi->GetEntries()<<", "<<ntGen->GetEntries()<<", "<<ntHlt->GetEntries()<<", "<<ntHi->GetEntries()<<", "<<ntSkim->GetEntries()<<endl;
  if(ntDkpi->GetEntries()!=ntHlt->GetEntries())
    {
      cout<<"    Error: Event numbers are different in three trees."<<endl;
      return 0;
    }
  Long64_t nentries = ntDkpi->GetEntries();

  cout<<" -- Event reading"<<endl;
  for(Long64_t i=0;i<nentries;i++)
    {
      if(i%100000==0) cout<<setiosflags(ios::left)<<"    "<<setw(8)<<i<<" / "<<nentries<<endl;
      ntDkpi->GetEntry(i);
      ntGen->GetEntry(i);
      ntHlt->GetEntry(i);
      ntHi->GetEntry(i);
      ntSkim->GetEntry(i);
      if(!skimevents||Dsize_ntDkpi>0)
	{
	  ntDkpi_new->Fill();
	  ntGen_new->Fill();
	  ntHlt_new->Fill(); 
	  ntHi_new->Fill(); 
	  ntSkim_new->Fill(); 
	} 
    }
  outf->Write();
  cout<<" -- Writing new trees done"<<endl;
  outf->Close();

  return 1;  
}

int main(int argc, char *argv[])
{
  if(argc==4)
    {
      skimTree_Dzero(argv[1], argv[2], atoi(argv[3]));
      return 1;
    }
  else
    {
      std::cout<<"Error: invalid parameter number."<<std::endl;
      return 0;
    }
}
