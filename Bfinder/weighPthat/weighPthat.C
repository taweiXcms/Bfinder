using namespace std;
#include "uti.h"
#define MAX_GEN      6000

int weighPurePthat(TString ifname = "",
                   TString ofname = "")
{
  Bool_t isInsidebin(Float_t xpthat, Float_t maxBgenpt, Int_t i);
  cout<<endl;
  cout<<" -- Checking if input and output files are same"<<endl;
  if(ifname==ofname)
    {
      cout<<"    Error: Input file will be overwritten."<<endl;
      return 0;
    }
  cout<<" -- Opening unweighed sample"<<endl;
  TFile* inf = TFile::Open(ifname);
  TTree* ntGen = (TTree*)inf->Get("ntGen");
  TTree* ntHi = (TTree*)inf->Get("ntHi");
  Int_t Gsize; ntGen->SetBranchAddress("Gsize",&Gsize);
  Float_t Gpt[MAX_GEN]; ntGen->SetBranchAddress("Gpt",Gpt);
  Float_t GisSignal[MAX_GEN]; ntGen->SetBranchAddress("GisSignal",GisSignal);
  Float_t pthat; ntHi->SetBranchAddress("pthat",&pthat);

  Float_t weight[nBins],nweight[nBins];
  for(Int_t j=0;j<nBins;j++)
    {
      weight[j]=0;
      nweight[j]=0;
    }
  cout<<" -- Checking event number"<<endl;
  if(ntGen->GetEntries()!=ntHi->GetEntries())
    {
      cout<<"    Error: Gen tree and Hi tree have different event number."<<endl;
      return 0;
    }
  Int_t nentries = ntGen->GetEntries();
  cout<<" -- Calculating weights"<<endl;
  for(Int_t i=0;i<nentries;i++)
    {
      ntGen->GetEntry(i);
      ntHi->GetEntry(i);
      if(i%100000==0) cout<<"    Processing event "<<setiosflags(ios::left)<<setw(7)<<i<<" / "<<nentries<<endl;
      Float_t maxpt=0;
      for(Int_t k=0;k<Gsize;k++)
        {
          if((GisSignal[k]==genSignal[0]||GisSignal[k]==genSignal[1])&&Gpt[k]>maxpt) maxpt=Gpt[k];
        }
      maxpt;
      for(Int_t j=0;j<nBins;j++)
        {
          if(isInsidebin(pthat,maxpt,j)) nweight[j]++;
        }
    }
  cout<<" -- Weight results"<<endl;
  for(Int_t j=0;j<nBins;j++)
    {
      if(nweight[j]==0)
        {
          cout<<"    Error: Weight fails."<<endl;
          return 0;
        }
      weight[j] = (crosssec[j]-crosssec[j+1])/nweight[j];
      cout<<"    Pthat"<<setiosflags(ios::left)<<setw(3)<<pthatBin[j]<<": "<<weight[j]<<endl;
    }

  cout<<" -- Building weight branch"<<endl;
  TFile* otf = TFile::Open(ofname,"update");
  TTree* ntHinew = (TTree*)otf->Get("ntHi");
  Float_t pthatweight,maxBgenpt;
  TBranch* newBr_pthatweight = ntHinew->Branch("pthatweight", &pthatweight, "pthatweight/F");
  TBranch* newBr_maxBgenpt = ntHinew->Branch("maxBgenpt", &maxBgenpt, "maxBgenpt/F");
  cout<<" -- Filling weight branch"<<endl;
  for(Int_t i=0;i<nentries;i++)
    {
      ntGen->GetEntry(i);
      ntHi->GetEntry(i);
      if(i%100000==0) cout<<"    Processing event "<<setiosflags(ios::left)<<setw(7)<<i<<" / "<<nentries<<endl;
      pthatweight=0;
      Float_t maxpt=0;
      for(Int_t k=0;k<Gsize;k++)
        {
          if((GisSignal[k]==genSignal[0]||GisSignal[k]==genSignal[1])&&Gpt[k]>maxpt) maxpt=Gpt[k];
        }
      maxBgenpt = maxpt;
      for(Int_t j=0;j<nBins;j++)
        {
          if(isInsidebin(pthat,maxBgenpt,j))
            {
              pthatweight = weight[j];
            }
        }
      newBr_pthatweight->Fill();
      newBr_maxBgenpt->Fill();
    }
  ntHinew->Write("", TObject::kOverwrite);

  cout<<" -- End"<<endl;
  cout<<endl;
  return 1;
}

Bool_t isInsidebin(Float_t xpthat, Float_t maxBgenpt, Int_t i)
{
  if(i>=nBins)
    {
      cout<<"    Error: invalid input"<<endl;
      return false;
    }
  if(i<(nBins-1)
	&& (xpthat>pthatBin[i]&&maxBgenpt>pthatBin[i])
    &&!(xpthat>pthatBin[i+1]&&maxBgenpt>pthatBin[i+1])) return true;
  else if(i==(nBins-1)&&xpthat>pthatBin[i]&&maxBgenpt>pthatBin[i]) return true;
  else return false;
}

int main(int argc, char *argv[])
{
  if(argc==3)
    {
      weighPurePthat(argv[1], argv[2]);
      return 1;
    }
  else
    {
      std::cout<<"Invalid parameter number"<<std::endl;
      return 0;
    }
}
