//HltInfo
int           Df_HLT_Run;
ULong64_t       Df_HLT_Event;
int           Df_HLT_LumiBlock;
void setHltTreeBranch(TTree* hltroot)
{
  hltroot->SetBranchAddress("Run",&Df_HLT_Run);
  hltroot->SetBranchAddress("Event",&Df_HLT_Event);
  hltroot->SetBranchAddress("LumiBlock",&Df_HLT_LumiBlock);
}

//hiEvtInfo
int           Df_HiTree_Run;
int           Df_HiTree_Evt;
int           Df_HiTree_Lumi;
void setHiTreeBranch(TTree* hitreeroot)
{
  hitreeroot->SetBranchAddress("run",&Df_HiTree_Run);
  hitreeroot->SetBranchAddress("evt",&Df_HiTree_Evt);
  hitreeroot->SetBranchAddress("lumi",&Df_HiTree_Lumi);
}
