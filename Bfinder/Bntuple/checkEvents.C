#include <iostream>
#include <map>
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;

int checkEvents(void) {

  string inFileNames[] = { 
    "/afs/cern.ch/user/a/anstahll/work/PromptReco/OniaTree_A.root",
    "/afs/cern.ch/user/a/anstahll/work/PromptReco/OniaTree_B.root",
    "/afs/cern.ch/user/a/anstahll/work/PromptReco/OniaTree_C.root",
    "/afs/cern.ch/user/a/anstahll/work/PromptReco/OniaTree_D.root"
  };  
  int nFiles = sizeof(inFileNames)/sizeof(string);

  TChain *chain = new TChain("hionia/myTree");
  for (int i=0; i<nFiles; i++) {
    chain->Add(inFileNames[i].c_str());
  }

  int eventNb,runNb;
  ULong64_t HLTriggers;
  map< pair<int, int>, int> eList;
  map< pair<int, int>, int>::iterator it;

  Int_t fCurrent = -1;
  chain->SetBranchAddress("eventNb",&eventNb);
  chain->SetBranchAddress("runNb",&runNb);
  chain->SetBranchAddress("HLTriggers",&HLTriggers);

  int nentries = chain->GetEntries();
  for (int ev=0; ev<nentries; ev++) {
    if (ev%1000000 == 0) cout << "Event: " << ev << " / " << nentries << endl;
   
    Long64_t centry = chain->LoadTree(ev);
    if (centry < 0) break;
    if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      cout << "[INFO] Changed to Tree number: " << inFileNames[fCurrent] << endl; 
    }
    chain->GetEntry(ev);
    //   if ((HLTriggers&(ULong64_t)pow(2,0))==(ULong64_t)pow(2,0)) {
      pair<int, int> element = make_pair(eventNb, runNb);
      pair<pair<int, int>, int> element2 = make_pair(element, ev);
      pair<map<pair<int, int>, int>::iterator, bool> result = eList.insert(element2);
      if (result.second == 0) cout << " Duplicated event in File " << inFileNames[fCurrent] << " : " << runNb << " " << eventNb << endl;
      // }
//    cout << " Duplicated event: " << runNb << " " << eventNb << endl;
  }

  cout << "Total map size: " << eList.size() << " | Total tree event size: " << chain->GetEntries()<< endl;

  delete chain;
  return 0;
}
