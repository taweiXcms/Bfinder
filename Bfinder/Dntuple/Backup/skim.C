#include <iostream>
#include <vector>
#include <algorithm>
#include <TTree.h>
#include <TFile.h>


// Get entry from the forest
void GetEntry(std::vector<TTree*> forest, int j)
{
  for (unsigned int i=0; i<forest.size(); i++)
  {
    forest[i]->GetEntry(j);
  }

}

// Fill the output for each tree in the forest
void FillOutput(std::vector<TTree*> cloneForest)
{
  for (unsigned int i=0; i<cloneForest.size(); i++)
  {
    cloneForest[i]->Fill();
  }
}

// Clone a tree
void AddCloneTree(std::vector<TTree*> &cloneForest,TFile *outf, TTree* t, const char *treeName)
{
  // Make directory
  outf->cd();
  //outf->mkdir(dirName);
  //outf->cd(dirName);

  // Add a clone tree to the clone forest
  TTree *tClone = t->CloneTree(0);
  tClone->SetMaxTreeSize(40000000000);
  tClone->SetName(treeName);

  cloneForest.push_back(tClone);
  std::cout <<"size"<<" "<<cloneForest.size();
}

// main routine
void skim(char *infname,char *outfname)
{
   std::vector<TTree*> cloneForest;
   std::vector<TTree*> forest;

   TFile *inf = new TFile(infname);
   
   // take the relevant trees from the file
   TTree *ntDkpi = (TTree*)inf->Get("ntDkpi");
   TTree *ntHlt = (TTree*)inf->Get("ntHlt");
   //TTree *ntHi = (TTree*)inf->Get("ntHi");
   TTree *ntSkim = (TTree*)inf->Get("ntSkim");
   TTree *ntGen = (TTree*)inf->Get("ntGen");
      
   forest.push_back(ntDkpi);
   //forest.push_back(ntHlt);
   forest.push_back(ntSkim);
   forest.push_back(ntGen);
   
   // define an output file
   TFile *outf = new TFile(outfname,"recreate");
   TTree *ntHltclone=new TTree("ntHltclone","");
   
   Int_t HLT_HIDmesonHITrackingGlobal_Dpt20_v1;
   Int_t HLT_HIDmesonHITrackingGlobal_Dpt40_v1;
   Int_t HLT_HIDmesonHITrackingGlobal_Dpt60_v1;
   Int_t L1_MinimumBiasHF1_AND;
   Int_t L1_MinimumBiasHF1_AND_Prescl;
   Int_t L1_MinimumBiasHF2_AND;
   Int_t L1_MinimumBiasHF2_AND_Prescl;
   Int_t L1_SingleS1Jet28_BptxAND;
   Int_t L1_SingleS1Jet28_BptxAND_Prescl;
   Int_t L1_SingleJet44_BptxAND;
   Int_t L1_SingleJet44_BptxAND_Prescl;
   
   ntHltclone->Branch("HLT_HIDmesonHITrackingGlobal_Dpt20_v1",&HLT_HIDmesonHITrackingGlobal_Dpt20_v1,"HLT_HIDmesonHITrackingGlobal_Dpt20_v1/I");
   ntHltclone->Branch("HLT_HIDmesonHITrackingGlobal_Dpt40_v1",&HLT_HIDmesonHITrackingGlobal_Dpt40_v1,"HLT_HIDmesonHITrackingGlobal_Dpt40_v1/I");
   ntHltclone->Branch("HLT_HIDmesonHITrackingGlobal_Dpt60_v1",&HLT_HIDmesonHITrackingGlobal_Dpt60_v1,"HLT_HIDmesonHITrackingGlobal_Dpt60_v1/I");
   ntHltclone->Branch("L1_MinimumBiasHF1_AND",&L1_MinimumBiasHF1_AND,"L1_MinimumBiasHF1_AND/I");
   ntHltclone->Branch("L1_MinimumBiasHF2_AND",&L1_MinimumBiasHF2_AND,"L1_MinimumBiasHF2_AND/I");
   ntHltclone->Branch("L1_SingleS1Jet28_BptxAND",&L1_SingleS1Jet28_BptxAND,"L1_SingleS1Jet28_BptxAND/I");
   ntHltclone->Branch("L1_SingleJet44_BptxAND",&L1_SingleJet44_BptxAND,"L1_SingleJet44_BptxAND/I");
   ntHltclone->Branch("L1_MinimumBiasHF1_AND_Prescl",&L1_MinimumBiasHF1_AND_Prescl,"L1_MinimumBiasHF1_AND_Prescl/I");
   ntHltclone->Branch("L1_MinimumBiasHF2_AND_Prescl",&L1_MinimumBiasHF2_AND_Prescl,"L1_MinimumBiasHF2_AND_Prescl/I");
   ntHltclone->Branch("L1_SingleS1Jet28_BptxAND_Prescl",&L1_SingleS1Jet28_BptxAND_Prescl,"L1_SingleS1Jet28_BptxAND_Prescl/I");
   ntHltclone->Branch("L1_SingleJet44_BptxAND_Prescl",&L1_SingleJet44_BptxAND_Prescl,"L1_SingleJet44_BptxAND_Prescl/I");

   AddCloneTree(cloneForest,outf,ntDkpi,"ntDkpi");
   //AddCloneTree(cloneForest,outf,ntHlt,"ntHlt");
   AddCloneTree(cloneForest,outf,ntSkim,"ntSkim");
   AddCloneTree(cloneForest,outf,ntGen,"ntGen");
   
   // You only need the branches which can help you decide if you want to keep the event
   int Dsize; 
     
   Int_t myHLT_HIDmesonHITrackingGlobal_Dpt20_v1;
   Int_t myHLT_HIDmesonHITrackingGlobal_Dpt40_v1;
   Int_t myHLT_HIDmesonHITrackingGlobal_Dpt60_v1;
   Int_t myL1_MinimumBiasHF1_AND;
   Int_t myL1_MinimumBiasHF1_AND_Prescl;
   Int_t myL1_MinimumBiasHF2_AND;
   Int_t myL1_MinimumBiasHF2_AND_Prescl;
   Int_t myL1_SingleS1Jet28_BptxAND;
   Int_t myL1_SingleS1Jet28_BptxAND_Prescl;
   Int_t myL1_SingleJet44_BptxAND;
   Int_t myL1_SingleJet44_BptxAND_Prescl;

   ntDkpi->SetBranchAddress("Dsize",&Dsize);
   ntHlt->SetBranchAddress("HLT_HIDmesonHITrackingGlobal_Dpt20_v1",&myHLT_HIDmesonHITrackingGlobal_Dpt20_v1);
   ntHlt->SetBranchAddress("HLT_HIDmesonHITrackingGlobal_Dpt40_v1",&myHLT_HIDmesonHITrackingGlobal_Dpt40_v1);
   ntHlt->SetBranchAddress("HLT_HIDmesonHITrackingGlobal_Dpt60_v1",&myHLT_HIDmesonHITrackingGlobal_Dpt60_v1);
   ntHlt->SetBranchAddress("L1_MinimumBiasHF1_AND",&myL1_MinimumBiasHF1_AND);
   ntHlt->SetBranchAddress("L1_MinimumBiasHF1_AND_Prescl",&myL1_MinimumBiasHF1_AND_Prescl);
   ntHlt->SetBranchAddress("L1_MinimumBiasHF2_AND",&myL1_MinimumBiasHF2_AND);
   ntHlt->SetBranchAddress("L1_MinimumBiasHF2_AND_Prescl",&myL1_MinimumBiasHF2_AND_Prescl);
   ntHlt->SetBranchAddress("L1_SingleS1Jet28_BptxAND",&myL1_SingleS1Jet28_BptxAND);
   ntHlt->SetBranchAddress("L1_SingleS1Jet28_BptxAND_Prescl",&myL1_SingleS1Jet28_BptxAND_Prescl);
   ntHlt->SetBranchAddress("L1_SingleJet44_BptxAND",&myL1_SingleJet44_BptxAND);
   ntHlt->SetBranchAddress("L1_SingleJet44_BptxAND_Prescl",&myL1_SingleJet44_BptxAND_Prescl);

   // main loop
   for (Long64_t i=0;i<ntDkpi->GetEntries();i++)
   {
         if (i%10000==0) std::cout <<ntDkpi->GetEntries()<<"/"<<i<<std::endl;
	 ntDkpi->GetEntry(i);
	 ntHlt->GetEntry(i);
	 
	 HLT_HIDmesonHITrackingGlobal_Dpt20_v1=myHLT_HIDmesonHITrackingGlobal_Dpt20_v1;
	 HLT_HIDmesonHITrackingGlobal_Dpt40_v1=myHLT_HIDmesonHITrackingGlobal_Dpt40_v1;
	 HLT_HIDmesonHITrackingGlobal_Dpt60_v1=myHLT_HIDmesonHITrackingGlobal_Dpt60_v1;
	 L1_MinimumBiasHF1_AND=myL1_MinimumBiasHF1_AND;
	 L1_MinimumBiasHF1_AND_Prescl=myL1_MinimumBiasHF1_AND_Prescl;
	 L1_MinimumBiasHF2_AND_Prescl=myL1_MinimumBiasHF2_AND;
	 L1_MinimumBiasHF2_AND_Prescl=myL1_MinimumBiasHF2_AND_Prescl;
	 L1_SingleS1Jet28_BptxAND=myL1_SingleS1Jet28_BptxAND;
	 L1_SingleS1Jet28_BptxAND_Prescl=myL1_SingleS1Jet28_BptxAND_Prescl;
	 L1_SingleJet44_BptxAND=myL1_SingleJet44_BptxAND;
	 L1_SingleJet44_BptxAND_Prescl=myL1_SingleJet44_BptxAND_Prescl;
	
	 ntHltclone->Fill();
	 	 
	 // if pass the selection -> accept this event
	 //if (Dsize>0) {
     GetEntry(forest,i);
	 FillOutput(cloneForest);
	 //} 
   }
   
   for (unsigned int i=0;i<cloneForest.size();i++)
   {
      cloneForest[i]->AutoSave();
   }
   ntHltclone->SetName("ntHlt");
   ntHltclone->AutoSave();
   
   outf->Close();
}


int main(int argc, char *argv[])
{
  if((argc != 3) && (argc != 4))
  {
    std::cout << "Usage: mergeForest <input_collection> <output_file>" << std::endl;
    return 1;
  }
  
  if(argc == 3)
    skim(argv[1], argv[2]);
  //else if (argc == 4)
  //  loop(argv[1], argv[2], argv[3]);
  return 0;
}



