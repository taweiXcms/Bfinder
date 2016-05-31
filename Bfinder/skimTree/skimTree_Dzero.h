#include "uti.h"

void SelectRecoBranches(TTree* nt)
{
  nt->SetBranchStatus("DtktkRes*",0);
  nt->SetBranchStatus("Dtrk3*",0);
  nt->SetBranchStatus("Dtrk4*",0);
  nt->SetBranchStatus("DRestrk*",0);
}

void SelectHltBranches(TTree* nthlt, Int_t isPbPb)
{
  nthlt->SetBranchStatus("*",0);
  if(isPbPb==1)
    {
      //nthlt->SetBranchStatus("HLT*",1);
      nthlt->SetBranchStatus("HLT_HIL1MinimumBias*",1);
      nthlt->SetBranchStatus("HLT_HIDmeson*",1);
      nthlt->SetBranchStatus("HLT_HIPuAK4CaloJet*",1);
      nthlt->SetBranchStatus("L1_MinimumBiasHF1_AND",1);
      nthlt->SetBranchStatus("L1_MinimumBiasHF2_AND",1);
      nthlt->SetBranchStatus("L1_SingleS1Jet16_BptxAND*", 1);
      nthlt->SetBranchStatus("L1_SingleS1Jet28_BptxAND*", 1);
      nthlt->SetBranchStatus("L1_SingleS1Jet32_BptxAND*", 1);
      nthlt->SetBranchStatus("L1_SingleJet44_BptxAND*",   1);
      nthlt->SetBranchStatus("L1_SingleS1Jet52_BptxAND*", 1);
      nthlt->SetBranchStatus("L1_SingleS1Jet56_BptxAND*", 1);
      nthlt->SetBranchStatus("L1_SingleS1Jet64_BptxAND*", 1);
    }
  else
    {
      nthlt->SetBranchStatus("HLT_L1MinimumBias*",1);
      nthlt->SetBranchStatus("HLT_Dmeson*",1);
      nthlt->SetBranchStatus("HLT_AK4CaloJet40_Eta5p1_v1",1);
      nthlt->SetBranchStatus("HLT_AK4CaloJet60_Eta5p1_v1",1);
      nthlt->SetBranchStatus("HLT_AK4CaloJet80_Eta5p1_v1",1);
      nthlt->SetBranchStatus("L1_SingleJet28_BptxAND",1);
      nthlt->SetBranchStatus("L1_SingleJet40_BptxAND",1);
      nthlt->SetBranchStatus("L1_SingleJet48_BptxAND",1);
    }
}
