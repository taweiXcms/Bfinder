#include "TTree.h"

//HltInfo
int           Df_HLT_Run;
ULong64_t     Df_HLT_Event;
int           Df_HLT_LumiBlock;
void setHltTreeBranch(TTree* hltroot)
{
  hltroot->SetBranchAddress("Run",&Df_HLT_Run);
  hltroot->SetBranchAddress("Event",&Df_HLT_Event);
  hltroot->SetBranchAddress("LumiBlock",&Df_HLT_LumiBlock);
}

//hiEvtInfo
unsigned int       Df_HiTree_Run;
unsigned long long Df_HiTree_Evt;
unsigned int       Df_HiTree_Lumi;
void setHiTreeBranch(TTree* hitreeroot)
{
  hitreeroot->SetBranchAddress("run",&Df_HiTree_Run);
  hitreeroot->SetBranchAddress("evt",&Df_HiTree_Evt);
  hitreeroot->SetBranchAddress("lumi",&Df_HiTree_Lumi);
}

void SetHlttreestatus(TTree * hltroot, bool isPbPb)
{
    hltroot->SetBranchStatus("*",0);
    hltroot->SetBranchStatus("Event",1);
    hltroot->SetBranchStatus("LumiBlock",1);
    hltroot->SetBranchStatus("Run",1);

    if( isPbPb )
    {   
        hltroot->SetBranchStatus("HLT_HIL1MinimumBiasHF1AND_*",1);
        hltroot->SetBranchStatus("HLT_HIL1MinimumBiasHF2AND_*",1);
        hltroot->SetBranchStatus("HLT_HIL1Centralityext*",1);
        hltroot->SetBranchStatus("HLT_HIDmesonHITrackingGlobal_Dpt*",1);
        hltroot->SetBranchStatus("HLT_HIPuAK4CaloJet*_Eta5p1_v*",1);
        hltroot->SetBranchStatus("HLT_HIPuAK4CaloDJet*_Eta2p1_v*",1);
        hltroot->SetBranchStatus("L1_MinimumBiasHF*_AND*",1);
        hltroot->SetBranchStatus("L1_Centrality_ext*",1);
        hltroot->SetBranchStatus("L1_SingleS1Jet*",1);
        hltroot->SetBranchStatus("L1_SingleJet*",1);
        hltroot->SetBranchStatus("L1_ZeroBias*",1);
    }   
    else
    {   
        hltroot->SetBranchStatus("HLT_L1MinimumBiasHF1OR_part*",1);
        hltroot->SetBranchStatus("HLT_DmesonPPTrackingGlobal_Dpt*",1);
        hltroot->SetBranchStatus("HLT_AK4*Jet*_*_v*",1);
        hltroot->SetBranchStatus("L1_MinimumBiasHF1_OR*",1);
        hltroot->SetBranchStatus("L1_SingleJet*",1);
        hltroot->SetBranchStatus("L1_ZeroBias*",1);
    }   
}
