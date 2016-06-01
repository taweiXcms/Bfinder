//HltInfo
int           Bf_HLT_Run;
ULong64_t     Bf_HLT_Event;
int           Bf_HLT_LumiBlock;
void setHltBranch(TTree* hltroot)
{
  hltroot->SetBranchAddress("Run",&Bf_HLT_Run);
  hltroot->SetBranchAddress("Event",&Bf_HLT_Event);
  hltroot->SetBranchAddress("LumiBlock",&Bf_HLT_LumiBlock);
}

//hiEvtInfo
unsigned int       Bf_HiTree_Run;
unsigned long long Bf_HiTree_Evt;
unsigned int       Bf_HiTree_Lumi;
void setHiTreeBranch(TTree* hitreeroot)
{
  hitreeroot->SetBranchAddress("run",&Bf_HiTree_Run);
  hitreeroot->SetBranchAddress("evt",&Bf_HiTree_Evt);
  hitreeroot->SetBranchAddress("lumi",&Bf_HiTree_Lumi);
}
void SetHlttreestatus(TTree * hltroot, bool isPbPb)
{
    hltroot->SetBranchStatus("*",0);
    hltroot->SetBranchStatus("Event",1);
    hltroot->SetBranchStatus("LumiBlock",1);
    hltroot->SetBranchStatus("Run",1);

    if( isPbPb )
    {
        hltroot->SetBranchStatus("*DoubleMu*",1);
    }
    else
    {
        hltroot->SetBranchStatus("*DoubleMu*",1);
    }
}
