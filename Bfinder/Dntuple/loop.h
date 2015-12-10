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
