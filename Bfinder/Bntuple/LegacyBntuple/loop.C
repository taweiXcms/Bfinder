#include "loop.h"

Bool_t iseos = true;
int loop(TString infile="/store/group/phys_heavyions/velicanu/forest/Run2015E/ExpressPhysics/Merged/HiForestExpress_baobab.root",
	 TString outfile="/data/wangj/Data2015/Bntuple/ntB_20151130_HiForestExpress_baobab.root", Bool_t REAL=true, Bool_t isPbPb=false, Int_t startEntries=0, Bool_t skim=false, Bool_t gskim=true, Bool_t checkMatching=false)
{
  void fillTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, Int_t j, Int_t typesize, Double_t track_mass1, Double_t track_mass2, Bool_t REAL);
  bool signalGen(Int_t Btype, Int_t j);

  cout<<endl;
  if(REAL) cout<<"--- Processing - REAL DATA"<<endl;
  else cout<<"--- Processing - MC"<<endl;

  TString ifname;
  if(iseos) ifname = Form("root://eoscms.cern.ch//eos/cms%s",infile.Data());
  else ifname = infile;
  TFile* f = TFile::Open(ifname);
  TTree* root = (TTree*)f->Get("Bfinder/root");
  TTree* hltroot = (TTree*)f->Get("hltanalysis/HltTree");
  TTree* hiroot  = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
  TFile* outf = new TFile(outfile,"recreate");
  setBBranch(root);
  setHltBranch(hltroot);
  if(isPbPb) setHiTreeBranch(hiroot);

  int ifchannel[7];
  ifchannel[0] = 1; //jpsi+Kp
  ifchannel[1] = 1; //jpsi+pi
  ifchannel[2] = 1; //jpsi+Ks(pi+,pi-)
  ifchannel[3] = 1; //jpsi+K*(K+,pi-)
  ifchannel[4] = 1; //jpsi+K*(K-,pi+)
  ifchannel[5] = 1; //jpsi+phi(K+,K-)
  ifchannel[6] = 1; //jpsi+pi pi <= psi', X(3872), Bs->J/psi f0
  
  cout<<"--- Building trees"<<endl;
  TTree* nt0 = new TTree("ntKp","");     buildBranch(nt0);
  TTree* nt1 = new TTree("ntpi","");     buildBranch(nt1);
  TTree* nt2 = new TTree("ntKs","");     buildBranch(nt2);
  TTree* nt3 = new TTree("ntKstar","");  buildBranch(nt3);
  TTree* nt5 = new TTree("ntphi","");    buildBranch(nt5);
  TTree* nt6 = new TTree("ntmix","");    buildBranch(nt6);
  TTree* ntGen = new TTree("ntGen","");  buildGenBranch(ntGen);
  TTree* ntHlt = hltroot->CloneTree(0);
  TTree* ntHi = hiroot->CloneTree(0);
  cout<<"--- Building trees finished"<<endl;

  Long64_t nentries = root->GetEntries();
  TVector3* bP = new TVector3;
  TVector3* bVtx = new TVector3;
  TLorentzVector* b4P = new TLorentzVector;
  TLorentzVector* b4Pout = new TLorentzVector;
  TLorentzVector* bGen = new TLorentzVector;
  cout<<"--- Check the number of events for two trees"<<endl;
  cout<<root->GetEntries()<<" "<<hltroot->GetEntries();
  if(isPbPb) cout<<" "<<hiroot->GetEntries();
  cout<<endl;
  cout<<"--- Processing events"<<endl;
  //nentries=1000;
  for(Int_t i=startEntries;i<nentries;i++)
    {
      root->GetEntry(i);
      hltroot->GetEntry(i);
      if(isPbPb) hiroot->GetEntry(i);
      if(i%100000==0) cout<<setw(8)<<i<<" / "<<nentries<<endl;
      if(checkMatching)
	{
	  if((Int_t)Bf_HLT_Event!=EvtInfo_EvtNo||Bf_HLT_Run!=EvtInfo_RunNo||Bf_HLT_LumiBlock!=EvtInfo_LumiNo || (isPbPb&&(Bf_HiTree_Evt!=EvtInfo_EvtNo||Bf_HiTree_Run!=EvtInfo_RunNo||Bf_HiTree_Lumi!=EvtInfo_LumiNo)))
	    {
	      cout<<"Error: not matched "<<i<<" | ";
	      cout<<"EvtNo("<<Bf_HLT_Event<<","<<EvtInfo_EvtNo<<") RunNo("<<Bf_HLT_Run<<","<<EvtInfo_RunNo<<") LumiNo("<<Bf_HLT_LumiBlock<<","<<EvtInfo_LumiNo<<") | EvtNo("<<Bf_HiTree_Evt<<","<<EvtInfo_EvtNo<<") RunNo("<<Bf_HiTree_Run<<","<<EvtInfo_RunNo<<") LumiNo("<<Bf_HiTree_Lumi<<","<<EvtInfo_LumiNo<<")"<<endl;
	      continue;
	    }
	}
      Int_t Btypesize[7]={0,0,0,0,0,0,0};
      Int_t ptflag=-1,ptMatchedflag=-1,probflag=-1,probMatchedflag=-1,tktkflag=-1,tktkMatchedflag=-1;
      Double_t pttem=0,ptMatchedtem=0,probtem=0,probMatchedtem=0,tktktem=0,tktkMatchedtem=0;
      for(Int_t t=0;t<7;t++)
	{
	  Int_t tidx = t-1;
	  if(t!=4)
	    {
	      tidx = t;
	      Bsize = 0;
	      ptflag = -1;
	      pttem = 0;
	      ptMatchedflag = -1;
	      ptMatchedtem = 0;
	      probflag = -1;
	      probtem = 0;
	      probMatchedflag = -1;
	      probMatchedtem = 0;
	      tktkflag = -1;
	      tktktem = 1000000.;
	      tktkMatchedflag = -1;
	      tktkMatchedtem = 1000000.;
	    }
	  if(ifchannel[t]==1)
	    {
	      for(int j=0;j<BInfo_size;j++)
		{
		  if(BInfo_type[j]==(t+1))
		    {
		      fillTree(bP,bVtx,b4P,j,Btypesize[tidx],tk1mass[t],tk2mass[t],REAL);
		      if(BInfo_pt[j]>pttem)
			{
			  ptflag = Btypesize[tidx];
			  pttem = BInfo_pt[j];
			}
		      if(TMath::Prob(BInfo_vtxchi2[j],BInfo_vtxdof[j])>probtem)
			{
			  probflag = Btypesize[tidx];
			  probtem = TMath::Prob(BInfo_vtxchi2[j],BInfo_vtxdof[j]);
			}
		      if(BInfo_type[j]>2&&BInfo_type[j]<7)
			{
			  if(TMath::Abs(BInfo_tktk_mass[j]-midmass[t])<tktktem)
			    {
			      tktkflag = Btypesize[tidx];
			      tktktem = TMath::Abs(BInfo_tktk_mass[j]-midmass[t]);
			    }
			}
		      if((!REAL&&(Bgen[Btypesize[tidx]]==23333||Bgen[Btypesize[tidx]]==41000))||REAL)//////////////////////////////
			{
			  if(BInfo_pt[j]>ptMatchedtem)
			    {
			      ptMatchedflag = Btypesize[tidx];
			      ptMatchedtem = BInfo_pt[j];
			    }
			  if(TMath::Prob(BInfo_vtxchi2[j],BInfo_vtxdof[j])>probMatchedtem)
			    {
			      probMatchedflag = Btypesize[tidx];
			      probMatchedtem = TMath::Prob(BInfo_vtxchi2[j],BInfo_vtxdof[j]);
			    }
			  if(BInfo_type[j]>2&&BInfo_type[j]<7)
			    {
			      if(TMath::Abs(BInfo_tktk_mass[j]-midmass[t])<tktkMatchedtem)
				{
				  tktkMatchedflag = Btypesize[tidx];
				  tktkMatchedtem = TMath::Abs(BInfo_tktk_mass[j]-midmass[t]);
				}
			    }
			}
		      Btypesize[tidx]++;
		    }
		}
	      if(t!=3)
		{
		  if(ptflag>=0) Bmaxpt[ptflag] = true;
		  if(probflag>=0) Bmaxprob[probflag] = true;
		  if(tktkflag>=0) Bbesttktkmass[tktkflag] = true;
		  if(ptMatchedflag>=0) BmaxptMatched[ptMatchedflag] = true;
		  if(probMatchedflag>=0) BmaxprobMatched[probMatchedflag] = true;
		  if(tktkMatchedflag>=0) BbesttktkmassMatched[tktkMatchedflag] = true;
		}
	      if(t==0)      nt0->Fill();
	      else if(t==1) nt1->Fill();
	      else if(t==2) nt2->Fill();
	      else if(t==4) nt3->Fill();
	      else if(t==5) nt5->Fill();
	      else if(t==6) nt6->Fill();
	    }
	}

      ntHlt->Fill();
      if(isPbPb) ntHi->Fill();

      if(!REAL)
	{
	  Int_t gt=0,sigtype=0;
	  Int_t gsize=0;
	  Gsize = 0;
	  for(int j=0;j<GenInfo_size;j++)
	    {
	      if((TMath::Abs(GenInfo_pdgId[j])!=BPLUS_PDGID&&TMath::Abs(GenInfo_pdgId[j])!=BZERO_PDGID&&TMath::Abs(GenInfo_pdgId[j])!=BSUBS_PDGID) && gskim) continue;
	      Gsize = gsize+1;
	      Gpt[gsize] = GenInfo_pt[j];
	      Geta[gsize] = GenInfo_eta[j];
	      Gphi[gsize] = GenInfo_phi[j];
	      GpdgId[gsize] = GenInfo_pdgId[j];
	      bGen->SetPtEtaPhiM(GenInfo_pt[j],GenInfo_eta[j],GenInfo_phi[j],GenInfo_mass[j]);
	      Gy[gsize] = bGen->Rapidity();
	      sigtype=0;
	      for(gt=1;gt<8;gt++)
		{
		  if(signalGen(gt,j))
		    {
		      sigtype=gt;
		      break;
		    }
		}
	      GisSignal[gsize] = sigtype;
	      Gmu1pt[gsize] = -1;
	      Gmu1eta[gsize] = -20;
	      Gmu1phi[gsize] = -20;
	      Gmu1p[gsize] = -1;
	      Gmu2pt[gsize] = -1;
	      Gmu2eta[gsize] = -20;
	      Gmu2phi[gsize] = -20;
	      Gmu2p[gsize] = -1;
	      Gtk1pt[gsize] = -1;
	      Gtk1eta[gsize] = -20;
	      Gtk1phi[gsize] = -20;
	      Gtk2pt[gsize] = -1;
	      Gtk2eta[gsize] = -20;
	      Gtk2phi[gsize] = -20;
	      if(sigtype!=0)
		{
		  Gmu1pt[gsize] = GenInfo_pt[GenInfo_da1[GenInfo_da1[j]]];
		  Gmu1eta[gsize] = GenInfo_eta[GenInfo_da1[GenInfo_da1[j]]];
		  Gmu1phi[gsize] = GenInfo_phi[GenInfo_da1[GenInfo_da1[j]]];
		  Gmu1p[gsize] = Gmu1pt[gsize]*cosh(Gmu1eta[gsize]);
		  Gmu2pt[gsize] = GenInfo_pt[GenInfo_da2[GenInfo_da1[j]]];
		  Gmu2eta[gsize] = GenInfo_eta[GenInfo_da2[GenInfo_da1[j]]];
		  Gmu2phi[gsize] = GenInfo_phi[GenInfo_da2[GenInfo_da1[j]]];
		  Gmu2p[gsize] = Gmu2pt[gsize]*cosh(Gmu2eta[gsize]);
		  if(sigtype==1||sigtype==2)
		    {
		      Gtk1pt[gsize] = GenInfo_pt[GenInfo_da2[j]];
		      Gtk1eta[gsize] = GenInfo_eta[GenInfo_da2[j]];
		      Gtk1phi[gsize] = GenInfo_phi[GenInfo_da2[j]];
		    }
		  else
		    {
		      Gtk1pt[gsize] = GenInfo_pt[GenInfo_da1[GenInfo_da2[j]]];
		      Gtk1eta[gsize] = GenInfo_eta[GenInfo_da1[GenInfo_da2[j]]];
		      Gtk1phi[gsize] = GenInfo_phi[GenInfo_da1[GenInfo_da2[j]]];
		      Gtk2pt[gsize] = GenInfo_pt[GenInfo_da2[GenInfo_da2[j]]];
		      Gtk2eta[gsize] = GenInfo_eta[GenInfo_da2[GenInfo_da2[j]]];
		      Gtk2phi[gsize] = GenInfo_phi[GenInfo_da2[GenInfo_da2[j]]];
		    }
		}
	      gsize++;
	    }
	  ntGen->Fill();
	}
    }

  outf->Write();
  outf->Close();
  return 1;
}


void fillTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, Int_t j, Int_t typesize, Double_t track_mass1, Double_t track_mass2, Bool_t REAL)
{

  //Event Info
  RunNo = EvtInfo_RunNo;
  EvtNo = EvtInfo_EvtNo;
  LumiNo = EvtInfo_LumiNo;
  Bsize = typesize+1;
  PVx = EvtInfo_PVx;
  PVy = EvtInfo_PVy;
  PVz = EvtInfo_PVz;
  PVxE = EvtInfo_PVxE;
  PVyE = EvtInfo_PVyE;
  PVzE = EvtInfo_PVzE;
  PVnchi2 = EvtInfo_PVnchi2;
  PVchi2 = EvtInfo_PVchi2;
  BSx = EvtInfo_BSx;
  BSy = EvtInfo_BSy;
  BSz = EvtInfo_BSz;
  BSxErr = EvtInfo_BSxErr;
  BSyErr = EvtInfo_BSyErr;
  BSzErr = EvtInfo_BSzErr;
  BSdxdz = EvtInfo_BSdxdz;
  BSdydz = EvtInfo_BSdydz;
  BSdxdzErr = EvtInfo_BSdxdzErr;
  BSdydzErr = EvtInfo_BSdydzErr;
  BSWidthX = EvtInfo_BSWidthX;
  BSWidthXErr = EvtInfo_BSWidthXErr;
  BSWidthY = EvtInfo_BSWidthY;
  BSWidthYErr = EvtInfo_BSWidthYErr;
  bP->SetXYZ(BInfo_px[j],BInfo_py[j],BInfo_pz[j]*0);
  bVtx->SetXYZ(BInfo_vtxX[j]-EvtInfo_PVx,
	       BInfo_vtxY[j]-EvtInfo_PVy,
	       BInfo_vtxZ[j]*0-EvtInfo_PVz*0);
  b4P->SetXYZM(BInfo_px[j],BInfo_py[j],BInfo_pz[j],BInfo_mass[j]);

  Bindex[typesize] = typesize;
  Btype[typesize] = BInfo_type[j];
  Bmass[typesize] = BInfo_mass[j];
  Bpt[typesize] = BInfo_pt[j];
  Beta[typesize] = BInfo_eta[j];
  Bphi[typesize] = BInfo_phi[j];
  By[typesize] = b4P->Rapidity();
  BvtxX[typesize] = BInfo_vtxX[j] - EvtInfo_PVx;
  BvtxY[typesize] = BInfo_vtxY[j] - EvtInfo_PVy;
  Bd0[typesize] = TMath::Sqrt((BInfo_vtxX[j]-EvtInfo_PVx)*(BInfo_vtxX[j]-EvtInfo_PVx)+(BInfo_vtxY[j]-EvtInfo_PVy)*(BInfo_vtxY[j]-EvtInfo_PVy));
  Bd0Err[typesize] = TMath::Sqrt(BInfo_vtxXErr[j]*BInfo_vtxXErr[j]+BInfo_vtxYErr[j]*BInfo_vtxYErr[j]);
  Bdxyz[typesize] = TMath::Sqrt((BInfo_vtxX[j]-EvtInfo_PVx)*(BInfo_vtxX[j]-EvtInfo_PVx)+(BInfo_vtxY[j]-EvtInfo_PVy)*(BInfo_vtxY[j]-EvtInfo_PVy)+(BInfo_vtxZ[j]-EvtInfo_PVz)*(BInfo_vtxZ[j]-EvtInfo_PVz));
  BdxyzErr[typesize] = TMath::Sqrt(BInfo_vtxXErr[j]*BInfo_vtxXErr[j]+BInfo_vtxYErr[j]*BInfo_vtxYErr[j]+BInfo_vtxZErr[j]*BInfo_vtxZErr[j]);
  Bchi2ndf[typesize] = BInfo_vtxchi2[j]/BInfo_vtxdof[j];
  Bchi2cl[typesize] = TMath::Prob(BInfo_vtxchi2[j],BInfo_vtxdof[j]);
  Bdtheta[typesize] = bP->Angle(*bVtx);
  Blxy[typesize] = ((BInfo_vtxX[j]-EvtInfo_PVx)*BInfo_px[j] + (BInfo_vtxY[j]-EvtInfo_PVy)*BInfo_py[j])/BInfo_pt[j];
  Double_t r2lxyBS = (BInfo_vtxX[j]-EvtInfo_BSx+(BInfo_vtxZ[j]-EvtInfo_BSz)*EvtInfo_BSdxdz) * (BInfo_vtxX[j]-EvtInfo_BSx+(BInfo_vtxZ[j]-EvtInfo_BSz)*EvtInfo_BSdxdz)
    + (BInfo_vtxY[j]-EvtInfo_BSy+(BInfo_vtxZ[j]-EvtInfo_BSz)*EvtInfo_BSdydz) * (BInfo_vtxY[j]-EvtInfo_BSy+(BInfo_vtxZ[j]-EvtInfo_BSz)*EvtInfo_BSdydz);
  Double_t xlxyBS = BInfo_vtxX[j]-EvtInfo_BSx + (BInfo_vtxZ[j]-EvtInfo_BSz)*EvtInfo_BSdxdz;
  Double_t ylxyBS = BInfo_vtxY[j]-EvtInfo_BSy + (BInfo_vtxZ[j]-EvtInfo_BSz)*EvtInfo_BSdydz;
  BlxyBS[typesize] = TMath::Sqrt(r2lxyBS);
  //BlxyBSErr[typesize] = 0;
  BlxyBSErr[typesize] = (1./r2lxyBS) * ((xlxyBS*xlxyBS)*BInfo_vtxXErr[j] + (2*xlxyBS*ylxyBS)*BInfo_vtxYXErr[j] + (ylxyBS*ylxyBS)*BInfo_vtxYErr[j]);
  Bmaxpt[typesize] = false;
  Bmaxprob[typesize] = false;
  Bbesttktkmass[typesize] = false;
  BmaxptMatched[typesize] = false;
  BmaxprobMatched[typesize] = false;
  BbesttktkmassMatched[typesize] = false;

  Double_t mu1px,mu1py,mu1pz,mu1E;
  mu1px = MuonInfo_pt[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]*cos(MuonInfo_phi[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]);
  mu1py = MuonInfo_pt[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]*sin(MuonInfo_phi[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]);
  mu1pz = MuonInfo_pt[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]*sinh(MuonInfo_eta[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]);
  b4P->SetXYZM(mu1px,mu1py,mu1pz,MUON_MASS);
  mu1E = b4P->E();
  Bmu1pt[typesize] = b4P->Pt();
  Bmu1p[typesize] = b4P->P();
  Bmu1eta[typesize] = b4P->Eta();
  Bmu1phi[typesize] = b4P->Phi();
  Bmu1y[typesize] = b4P->Rapidity();
  Double_t mu2px,mu2py,mu2pz,mu2E;
  mu2px = MuonInfo_pt[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]*cos(MuonInfo_phi[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]);
  mu2py = MuonInfo_pt[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]*sin(MuonInfo_phi[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]);
  mu2pz = MuonInfo_pt[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]*sinh(MuonInfo_eta[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]);
  b4P->SetXYZM(mu2px,mu2py,mu2pz,MUON_MASS);
  mu2E = b4P->E();
  Bmu2pt[typesize] = b4P->Pt();
  Bmu2p[typesize] = b4P->P();
  Bmu2eta[typesize] = b4P->Eta();
  Bmu2phi[typesize] = b4P->Phi();
  Bmu2y[typesize] = b4P->Rapidity();

  Bmu1dzPV[typesize] = MuonInfo_dzPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2dzPV[typesize] = MuonInfo_dzPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1dxyPV[typesize] = MuonInfo_dxyPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2dxyPV[typesize] = MuonInfo_dxyPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1normchi2[typesize] = MuonInfo_normchi2[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2normchi2[typesize] = MuonInfo_normchi2[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1Chi2ndf[typesize] = MuonInfo_i_chi2[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]/MuonInfo_i_ndf[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2Chi2ndf[typesize] = MuonInfo_i_chi2[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]/MuonInfo_i_ndf[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1muqual[typesize] = MuonInfo_muqual[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2muqual[typesize] = MuonInfo_muqual[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1TrackerMuArbitrated[typesize] = MuonInfo_TrackerMuonArbitrated[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2TrackerMuArbitrated[typesize] = MuonInfo_TrackerMuonArbitrated[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1isTrackerMuon[typesize] = MuonInfo_isTrackerMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2isTrackerMuon[typesize] = MuonInfo_isTrackerMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1isGlobalMuon[typesize] = MuonInfo_isGlobalMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2isGlobalMuon[typesize] = MuonInfo_isGlobalMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1TMOneStationTight[typesize] = MuonInfo_TMOneStationTight[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2TMOneStationTight[typesize] = MuonInfo_TMOneStationTight[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1InPixelLayer[typesize] = MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2InPixelLayer[typesize] = MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1InStripLayer[typesize] = MuonInfo_i_nStripLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2InStripLayer[typesize] = MuonInfo_i_nStripLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1InTrackerLayer[typesize] = MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]] + MuonInfo_i_nStripLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2InTrackerLayer[typesize] = MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]] + MuonInfo_i_nStripLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  b4P->SetPxPyPzE(mu1px+mu2px,
		  mu1py+mu2py,
		  mu1pz+mu2pz,
		  mu1E+mu2E);
  Bmumumass[typesize] = b4P->Mag();
  Bmumueta[typesize] = b4P->Eta();
  Bmumuphi[typesize] = b4P->Phi();
  Bmumuy[typesize] = b4P->Rapidity();
  Bmumupt[typesize] = b4P->Pt();
  Bujmass[typesize] = BInfo_uj_mass[BInfo_rfuj_index[j]];
  BujvProb[typesize] = TMath::Prob(BInfo_uj_vtxchi2[BInfo_rfuj_index[j]],BInfo_uj_vtxdof[BInfo_rfuj_index[j]]);
  b4P->SetXYZM(BInfo_uj_px[BInfo_rfuj_index[j]],
	       BInfo_uj_py[BInfo_rfuj_index[j]],
	       BInfo_uj_pz[BInfo_rfuj_index[j]],
	       BInfo_uj_mass[BInfo_rfuj_index[j]]);
  Bujpt[typesize] = b4P->Pt();
  Bujeta[typesize] = b4P->PseudoRapidity();
  Bujphi[typesize] = b4P->Phi();
  Bujy[typesize] = b4P->Rapidity();
  Bujlxy[typesize] = ((BInfo_uj_vtxX[BInfo_rfuj_index[j]]-EvtInfo_PVx)*BInfo_uj_px[BInfo_rfuj_index[j]] + (BInfo_uj_vtxY[BInfo_rfuj_index[j]]-EvtInfo_PVy)*BInfo_uj_py[BInfo_rfuj_index[j]])/Bujpt[typesize];
  /*
  if(MuonInfo_muqual[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]&16) mu1TrackerMuArbitrated[typesize] = 1;
  else mu1TrackerMuArbitrated[typesize] = 0;
  if(MuonInfo_muqual[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]&4096) mu1TMOneStationTight[typesize] = 1;
  else mu1TMOneStationTight[typesize] = 0;
  if(MuonInfo_muqual[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]&16) mu2TrackerMuArbitrated[typesize] = 1;
  else mu2TrackerMuArbitrated[typesize] = 0;
  if(MuonInfo_muqual[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]&4096) mu2TMOneStationTight[typesize] = 1;
  else mu2TMOneStationTight[typesize] = 0;
  */

  Double_t tk1px,tk1py,tk1pz,tk1E;
  Double_t tk2px,tk2py,tk2pz,tk2E;
  Btrk1Idx[typesize] = BInfo_rftk1_index[j];
  Btrk2Idx[typesize] = BInfo_rftk2_index[j];
  if(BInfo_type[j]==1 || BInfo_type[j]==2)
    {
      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk1_index[j]],TrackInfo_eta[BInfo_rftk1_index[j]],TrackInfo_phi[BInfo_rftk1_index[j]],track_mass1);
      Btrk1Pt[typesize] = TrackInfo_pt[BInfo_rftk1_index[j]];
      Btrk1Eta[typesize] = TrackInfo_eta[BInfo_rftk1_index[j]];
      Btrk1Phi[typesize] = TrackInfo_phi[BInfo_rftk1_index[j]];
      Btrk1PtErr[typesize] = TrackInfo_ptErr[BInfo_rftk1_index[j]];
      Btrk1EtaErr[typesize] = TrackInfo_etaErr[BInfo_rftk1_index[j]];
      Btrk1PhiErr[typesize] = TrackInfo_phiErr[BInfo_rftk1_index[j]];
      Btrk1Y[typesize] = b4P->Rapidity();
      Btrk1Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk1_index[j]];
      Btrk1D0Err[typesize] = TrackInfo_d0error[BInfo_rftk1_index[j]];
      Btrk1PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk1_index[j]];
      Btrk1StripHit[typesize] = TrackInfo_striphit[BInfo_rftk1_index[j]];
      Btrk1nPixelLayer[typesize] = TrackInfo_nPixelLayer[BInfo_rftk1_index[j]];
      Btrk1nStripLayer[typesize] = TrackInfo_nStripLayer[BInfo_rftk1_index[j]];
      Btrk1Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk1_index[j]]/TrackInfo_ndf[BInfo_rftk1_index[j]];
      Btrk1MVAVal[typesize] = TrackInfo_trkMVAVal[BInfo_rftk1_index[j]];
      Btrk1Algo[typesize] = TrackInfo_trkAlgo[BInfo_rftk1_index[j]];
      Btrk1highPurity[typesize] = TrackInfo_highPurity[BInfo_rftk1_index[j]];
      Btrk1Quality[typesize] = TrackInfo_trackQuality[BInfo_rftk1_index[j]];
      Btrk2Pt[typesize] = -1;
      Btrk2Eta[typesize] = -20;
      Btrk2Phi[typesize] = -20;
      Btrk2PtErr[typesize] = 0;
      Btrk2EtaErr[typesize] = 0;
      Btrk2PhiErr[typesize] = 0;
      Btrk2Y[typesize] = -1;
      Btrk2Dxy[typesize] = -1;
      Btrk2D0Err[typesize] = -1;
      Btrk2PixelHit[typesize] = -1;
      Btrk2StripHit[typesize] = -1;
      Btrk2nPixelLayer[typesize] = -1;
      Btrk2nStripLayer[typesize] = -1;
      Btrk2Chi2ndf[typesize] = -1;
      Btrk2MVAVal[typesize] = -100;
      Btrk2Algo[typesize] = 0;
      Btrk2highPurity[typesize] = false;
      Btrk2Quality[typesize] = 0;
      Btktkmass[typesize] = -1;
      BtktkvProb[typesize] = -1;
      Btktkpt[typesize] = -1;
      Btktketa[typesize] = -20;
      Btktkphi[typesize] = -20;
      Btktky[typesize] = -1;
      Bdoubletmass[typesize] = -1;
      Bdoubletpt[typesize] = -1;
      Bdoubleteta[typesize] = -20;
      Bdoubletphi[typesize] = -20;
      Bdoublety[typesize] = -1;
    }  
  else if(BInfo_type[j]==5)
    {
      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk2_index[j]],TrackInfo_eta[BInfo_rftk2_index[j]],TrackInfo_phi[BInfo_rftk2_index[j]],track_mass1);
      Btrk1Pt[typesize] = TrackInfo_pt[BInfo_rftk2_index[j]];
      Btrk1Eta[typesize] = TrackInfo_eta[BInfo_rftk2_index[j]];
      Btrk1Phi[typesize] = TrackInfo_phi[BInfo_rftk2_index[j]];
      Btrk1PtErr[typesize] = TrackInfo_ptErr[BInfo_rftk2_index[j]];
      Btrk1EtaErr[typesize] = TrackInfo_etaErr[BInfo_rftk2_index[j]];
      Btrk1PhiErr[typesize] = TrackInfo_phiErr[BInfo_rftk2_index[j]];
      Btrk1Y[typesize] = b4P->Rapidity();
      Btrk1Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk2_index[j]];
      Btrk1D0Err[typesize] = TrackInfo_d0error[BInfo_rftk2_index[j]];
      Btrk1PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk2_index[j]];
      Btrk1StripHit[typesize] = TrackInfo_striphit[BInfo_rftk2_index[j]];
      Btrk1nPixelLayer[typesize] = TrackInfo_nPixelLayer[BInfo_rftk2_index[j]];
      Btrk1nStripLayer[typesize] = TrackInfo_nStripLayer[BInfo_rftk2_index[j]];
      Btrk1Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk2_index[j]]/TrackInfo_ndf[BInfo_rftk2_index[j]];
      Btrk1MVAVal[typesize] = TrackInfo_trkMVAVal[BInfo_rftk2_index[j]];
      Btrk1Algo[typesize] = TrackInfo_trkAlgo[BInfo_rftk2_index[j]];
      Btrk1highPurity[typesize] = TrackInfo_highPurity[BInfo_rftk2_index[j]];
      Btrk1Quality[typesize] = TrackInfo_trackQuality[BInfo_rftk2_index[j]];
      tk1px = b4P->Px();
      tk1py = b4P->Py();
      tk1pz = b4P->Pz();
      tk1E = b4P->E();
      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk1_index[j]],TrackInfo_eta[BInfo_rftk1_index[j]],TrackInfo_phi[BInfo_rftk1_index[j]],track_mass2);
      Btrk2Pt[typesize] = TrackInfo_pt[BInfo_rftk1_index[j]];
      Btrk2Eta[typesize] = TrackInfo_eta[BInfo_rftk1_index[j]];
      Btrk2Phi[typesize] = TrackInfo_phi[BInfo_rftk1_index[j]];
      Btrk2PtErr[typesize] = TrackInfo_ptErr[BInfo_rftk1_index[j]];
      Btrk2EtaErr[typesize] = TrackInfo_etaErr[BInfo_rftk1_index[j]];
      Btrk2PhiErr[typesize] = TrackInfo_phiErr[BInfo_rftk1_index[j]];
      Btrk2Y[typesize] = b4P->Rapidity();
      Btrk2Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk1_index[j]];
      Btrk2D0Err[typesize] = TrackInfo_d0error[BInfo_rftk1_index[j]];
      Btrk2PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk1_index[j]];
      Btrk2StripHit[typesize] = TrackInfo_striphit[BInfo_rftk1_index[j]];
      Btrk2nPixelLayer[typesize] = TrackInfo_nPixelLayer[BInfo_rftk1_index[j]];
      Btrk2nStripLayer[typesize] = TrackInfo_nStripLayer[BInfo_rftk1_index[j]];
      Btrk2Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk1_index[j]]/TrackInfo_ndf[BInfo_rftk1_index[j]];
      Btrk2MVAVal[typesize] = TrackInfo_trkMVAVal[BInfo_rftk1_index[j]];
      Btrk2Algo[typesize] = TrackInfo_trkAlgo[BInfo_rftk1_index[j]];
      Btrk2highPurity[typesize] = TrackInfo_highPurity[BInfo_rftk1_index[j]];
      Btrk2Quality[typesize] = TrackInfo_trackQuality[BInfo_rftk1_index[j]];
      tk2px = b4P->Px();
      tk2py = b4P->Py();
      tk2pz = b4P->Pz();
      tk2E = b4P->E();

      b4P->SetPxPyPzE(tk1px+tk2px,
		      tk1py+tk2py,
		      tk1pz+tk2pz,
		      tk1E+tk2E);
      Btktkmass[typesize] = b4P->Mag();
      Btktketa[typesize] = b4P->Eta();
      Btktkphi[typesize] = b4P->Phi();
      Btktky[typesize] = b4P->Rapidity();
      Btktkpt[typesize] = b4P->Pt();
      BtktkvProb[typesize] = TMath::Prob(BInfo_tktk_vtxchi2[j],BInfo_tktk_vtxdof[j]);
      Bdoubletmass[typesize] = BInfo_tktk_mass[j];
      b4P->SetXYZM(BInfo_tktk_px[j],BInfo_tktk_py[j],BInfo_tktk_pz[j],BInfo_tktk_mass[j]);
      Bdoubletpt[typesize] = b4P->Pt();
      Bdoubleteta[typesize] = b4P->PseudoRapidity();
      Bdoubletphi[typesize] = b4P->Phi();
      Bdoublety[typesize] = b4P->Rapidity();

      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk1_index[j]],TrackInfo_eta[BInfo_rftk1_index[j]],TrackInfo_phi[BInfo_rftk1_index[j]],KAON_MASS);
      double tk1EK = b4P->E();
      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk2_index[j]],TrackInfo_eta[BInfo_rftk2_index[j]],TrackInfo_phi[BInfo_rftk2_index[j]],KAON_MASS);
      double tk2EK = b4P->E();
      b4P->SetPxPyPzE(tk1px+tk2px,
		      tk1py+tk2py,
		      tk1pz+tk2pz,
		      tk1EK+tk2EK);
      BtktkmassKK[typesize] = b4P->Mag();
    }
  else
    {
      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk1_index[j]],TrackInfo_eta[BInfo_rftk1_index[j]],TrackInfo_phi[BInfo_rftk1_index[j]],track_mass1);
      Btrk1Pt[typesize] = TrackInfo_pt[BInfo_rftk1_index[j]];
      Btrk1Eta[typesize] = TrackInfo_eta[BInfo_rftk1_index[j]];
      Btrk1Phi[typesize] = TrackInfo_phi[BInfo_rftk1_index[j]];
      Btrk1PtErr[typesize] = TrackInfo_ptErr[BInfo_rftk1_index[j]];
      Btrk1EtaErr[typesize] = TrackInfo_etaErr[BInfo_rftk1_index[j]];
      Btrk1PhiErr[typesize] = TrackInfo_phiErr[BInfo_rftk1_index[j]];
      Btrk1Y[typesize] = b4P->Rapidity();
      Btrk1Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk1_index[j]];
      Btrk1D0Err[typesize] = TrackInfo_d0error[BInfo_rftk1_index[j]];
      Btrk1PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk1_index[j]];
      Btrk1StripHit[typesize] = TrackInfo_striphit[BInfo_rftk1_index[j]];
      Btrk1nPixelLayer[typesize] = TrackInfo_nPixelLayer[BInfo_rftk1_index[j]];
      Btrk1nStripLayer[typesize] = TrackInfo_nStripLayer[BInfo_rftk1_index[j]];
      Btrk1Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk1_index[j]]/TrackInfo_ndf[BInfo_rftk1_index[j]];
      Btrk1MVAVal[typesize] = TrackInfo_trkMVAVal[BInfo_rftk1_index[j]];
      Btrk1Algo[typesize] = TrackInfo_trkAlgo[BInfo_rftk1_index[j]];
      Btrk1highPurity[typesize] = TrackInfo_highPurity[BInfo_rftk1_index[j]];
      Btrk1Quality[typesize] = TrackInfo_trackQuality[BInfo_rftk1_index[j]];
      tk1px = b4P->Px();
      tk1py = b4P->Py();
      tk1pz = b4P->Pz();
      tk1E = b4P->E();
      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk2_index[j]],TrackInfo_eta[BInfo_rftk2_index[j]],TrackInfo_phi[BInfo_rftk2_index[j]],track_mass2);
      Btrk2Pt[typesize] = TrackInfo_pt[BInfo_rftk2_index[j]];
      Btrk2Eta[typesize] = TrackInfo_eta[BInfo_rftk2_index[j]];
      Btrk2Phi[typesize] = TrackInfo_phi[BInfo_rftk2_index[j]];
      Btrk2PtErr[typesize] = TrackInfo_ptErr[BInfo_rftk2_index[j]];
      Btrk2EtaErr[typesize] = TrackInfo_etaErr[BInfo_rftk2_index[j]];
      Btrk2PhiErr[typesize] = TrackInfo_phiErr[BInfo_rftk2_index[j]];
      Btrk2Y[typesize] = b4P->Rapidity();
      Btrk2Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk2_index[j]];
      Btrk2D0Err[typesize] = TrackInfo_d0error[BInfo_rftk2_index[j]];
      Btrk2PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk2_index[j]];
      Btrk2StripHit[typesize] = TrackInfo_striphit[BInfo_rftk2_index[j]];
      Btrk2nPixelLayer[typesize] = TrackInfo_nPixelLayer[BInfo_rftk2_index[j]];
      Btrk2nStripLayer[typesize] = TrackInfo_nStripLayer[BInfo_rftk2_index[j]];
      Btrk2Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk2_index[j]]/TrackInfo_ndf[BInfo_rftk2_index[j]];
      Btrk2MVAVal[typesize] = TrackInfo_trkMVAVal[BInfo_rftk2_index[j]];
      Btrk2Algo[typesize] = TrackInfo_trkAlgo[BInfo_rftk2_index[j]];
      Btrk2highPurity[typesize] = TrackInfo_highPurity[BInfo_rftk2_index[j]];
      Btrk2Quality[typesize] = TrackInfo_trackQuality[BInfo_rftk2_index[j]];
      tk2px = b4P->Px();
      tk2py = b4P->Py();
      tk2pz = b4P->Pz();
      tk2E = b4P->E();

      b4P->SetPxPyPzE(tk1px+tk2px,
		      tk1py+tk2py,
		      tk1pz+tk2pz,
		      tk1E+tk2E);
      Btktkmass[typesize] = b4P->Mag();
      Btktketa[typesize] = b4P->Eta();
      Btktkphi[typesize] = b4P->Phi();
      Btktky[typesize] = b4P->Rapidity();
      Btktkpt[typesize] = b4P->Pt();
      BtktkvProb[typesize] = TMath::Prob(BInfo_tktk_vtxchi2[j],BInfo_tktk_vtxdof[j]);
      Bdoubletmass[typesize] = BInfo_tktk_mass[j];
      b4P->SetXYZM(BInfo_tktk_px[j],BInfo_tktk_py[j],BInfo_tktk_pz[j],BInfo_tktk_mass[j]);
      Bdoubletpt[typesize] = b4P->Pt();
      Bdoubleteta[typesize] = b4P->PseudoRapidity();
      Bdoubletphi[typesize] = b4P->Phi();
      Bdoublety[typesize] = b4P->Rapidity();

      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk1_index[j]],TrackInfo_eta[BInfo_rftk1_index[j]],TrackInfo_phi[BInfo_rftk1_index[j]],KAON_MASS);
      double tk1EK = b4P->E();
      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk2_index[j]],TrackInfo_eta[BInfo_rftk2_index[j]],TrackInfo_phi[BInfo_rftk2_index[j]],KAON_MASS);
      double tk2EK = b4P->E();
      b4P->SetPxPyPzE(tk1px+tk2px,
		      tk1py+tk2py,
		      tk1pz+tk2pz,
		      tk1EK+tk2EK);
      BtktkmassKK[typesize] = b4P->Mag();
    }

  //gen info judgement
  if(!REAL)
    {
      Bgen[typesize] = 0;
      BgenIndex[typesize] = -1;
      Bgenpt[typesize] = -1;
      Bgeneta[typesize] = -20;
      Bgenphi[typesize] = -20;
      Bgeny[typesize] = -1;
      int mGenIdxTk1=-1;
      int mGenIdxTk2=-1;
      int bGenIdxTk1=-1;
      int bGenIdxTk2=-1;
      int bGenIdxMu1=-1;
      int bGenIdxMu2=-1;
      int ujGenIdxMu1=-1;
      int ujGenIdxMu2=-1;
      
      Double_t BId,MId,tk1Id,tk2Id;
      //tk1:positive, tk2:negtive
      if(BInfo_type[j]==1)
	{
	  BId = 521;//B+-
	  MId = -1;
	  tk1Id = 321;//K+-
	  tk2Id = -1;
	}
      if(BInfo_type[j]==2)
	{
	  BId = 521;//B+-
	  MId = -1;
	  tk1Id = 211;//pi+-
	  tk2Id = -1;
	}
      if(BInfo_type[j]==3)
	{
	  BId = 511;//B0
	  MId = 310;//Ks
	  tk1Id = 211;//pi+
	  tk2Id = 211;//pi-
	}
      if(BInfo_type[j]==4)
	{
	  BId = 511;//B0
	  MId = 313;//K*0
	  tk1Id = 321;//K+
	  tk2Id = 211;//pi-
	}
      if(BInfo_type[j]==5)
	{
	  BId = 511;//B0
	  MId = 313;//K*0
	  tk1Id = 211;//pi+
	  tk2Id = 321;//K-
	}
      if(BInfo_type[j]==6)
	{
	  BId = 531;//Bs
	  MId = 333;//phi
	  tk1Id = 321;//K+
	  tk2Id = 321;//K-
	}

      int twoTks,kStar,flagkstar=0;
      if(BInfo_type[j]==1 || BInfo_type[j]==2) twoTks=0;
      else twoTks=1;
      if(BInfo_type[j]==4 || BInfo_type[j]==5) kStar=1;
      else kStar=0;
      int nonprompt=0,prompt=0;

      // tk1
      if(TrackInfo_geninfo_index[BInfo_rftk1_index[j]]>-1)
	{
	  int level =0;
	  if(abs(GenInfo_pdgId[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]])==tk1Id)
	    {
	      level = 1;
	      if(GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]>-1)
		{
		  if(!twoTks)//one trk channel
		    {
		      mGenIdxTk1=0;
		      if(abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]])==BId)
			{
			  level = 3;
			  bGenIdxTk1=GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]];
			}		  
		    }
		  else//two trk channel
		    {
		      if(abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]])==MId)
			{
			  level = 2;
			  if(GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]]>-1)
			    {
			      if(abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]]])==BId)
				{
				  level = 3;
				  bGenIdxTk1=GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]];
				}
			    }
			  mGenIdxTk1=GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]];
			}
		    }
		}
	    }
	  Bgen[typesize]=level;
	}
      
      //tk2
      if(!twoTks)//one trk channel
	{
	  Bgen[typesize]+=30;
	  mGenIdxTk2=0;
	  bGenIdxTk2=0;
	}
      else//two trk channel
	{
	  if(TrackInfo_geninfo_index[BInfo_rftk2_index[j]]>-1)
	    {
	      int level =0;
	      if(abs(GenInfo_pdgId[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]])==tk2Id)
		{
		  level = 1;
		  if(GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]>-1)
		    {
		      if(abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]])==MId)
			{
			  level = 2;
			  if(GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]]>-1)
			    {
			      if(abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]]])==BId)
				{
				  level = 3;
				  bGenIdxTk2 = GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]];
				}
			    }
			  mGenIdxTk2 = GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]];
			}
		    }
		}
	      Bgen[typesize]+=(level*10);
	    }
	}

      
      //mu1
      //cout<<MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]<<endl;
      if(MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]>-1)
	{
	  int level =0;
	  if(abs(GenInfo_pdgId[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]])==13)
	    {
	      level=1;
	      if(GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]>-1)
		{
		  if(GenInfo_pdgId[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]]==443)
		    {
		      ujGenIdxMu1 = GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]];
		      level=2;
		      if(GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]]>-1)
			{
			  if(abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]]])==BId)
			    {
			      //nonprompt=1;
			      level = 3;
			      bGenIdxMu1=GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]];
			      flagkstar++;///////////////////////////////////////////////=1
			    }
			}
		      else 
			{
			  prompt=1;
			}
		    } 
		}
	    }
	  Bgen[typesize]+=(level*100);
	}
      
      //mu2
      if(MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]>-1)
	{  
	  int level =0;
	  if(abs(GenInfo_pdgId[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]])==13)
	    {
	      level = 1;
	      if(GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]>-1)
		{
		  if(GenInfo_pdgId[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]]==443)
		    {
		      ujGenIdxMu2 = GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]];
		      level = 2;
		      if(GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]]>-1)
			{
			  if(abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]]])==BId)
			    {
			      level = 3;
			      bGenIdxMu2=GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]];
			      flagkstar++;///////////////////////////////////////////////////=2
			    }
			}
		    }
		}
	    }
	  Bgen[typesize]+=(level*1000);
	}

      int level=0;
      if(mGenIdxTk1!=-1 && mGenIdxTk2!=-1)
	{
	  if(!twoTks) level=1;
	  else
	    {
	      if(mGenIdxTk1==mGenIdxTk2) level=1;
	    }
	}
      if(bGenIdxMu1!=-1 && bGenIdxMu1==bGenIdxMu2 && bGenIdxMu1==bGenIdxTk1)
	{
	  if(!twoTks)
	    {
	      level=2;
	      BgenIndex[typesize] = bGenIdxMu1;
	    }
	  else if(bGenIdxMu1==bGenIdxTk2)
	    {
	      level=2;
	      BgenIndex[typesize] = bGenIdxMu1;
	    }
	}
      Bgen[typesize]+=(level*10000);

      //kstar#############################################################################
      if(kStar)
	{
	  //tk1
	  if(TrackInfo_geninfo_index[BInfo_rftk1_index[j]]>-1)
	    {
	      if(abs(GenInfo_pdgId[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]])==tk2Id)
		{
		  if(GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]>-1)
		    {
		      if(abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]])==MId)
			{
			  if(GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]]>-1)
			    {
			      if(abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]]])==BId)
				{
				  flagkstar++;//////////////////////////////////////////////=3
				  bGenIdxTk1=GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]];
				}
			    }
			  mGenIdxTk1=GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]];
			}
		    }
		}
	    }
	  
	  //tk2
	  if(TrackInfo_geninfo_index[BInfo_rftk2_index[j]]>-1)
	    {
	      if(abs(GenInfo_pdgId[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]])==tk1Id)
		{
		  if(GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]>-1)
		    {
		      if(abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]])==MId)
			{
			  if(GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]]>-1)
			    {
			      if(abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]]])==BId)
				{
				  flagkstar++;////////////////////////////////////////////////////=4
				  bGenIdxTk2 = GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]];
				}
			    }
			  mGenIdxTk2 = GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]];
			}
		    }
		}
	    }
	  if(flagkstar==4)
	    {
	      if((bGenIdxMu1!=-1) 
		 && (bGenIdxMu1==bGenIdxMu2)
		 && (bGenIdxMu1==bGenIdxTk1)
		 && (bGenIdxMu1==bGenIdxTk2)
		 )
		{
		  Bgen[typesize]=41000;
		}
	    }
	}//kstar End#############################################################################

      int tgenIndex=BgenIndex[typesize];
      if(Bgen[typesize]==23333 || Bgen[typesize]==41000)
	{
	  Bgenpt[typesize] = GenInfo_pt[tgenIndex];
	  Bgeneta[typesize] = GenInfo_eta[tgenIndex];
	  Bgenphi[typesize] = GenInfo_phi[tgenIndex];
	  b4P->SetXYZM(GenInfo_pt[tgenIndex]*cos(GenInfo_phi[tgenIndex]),
		       GenInfo_pt[tgenIndex]*sin(GenInfo_phi[tgenIndex]),
		       GenInfo_pt[tgenIndex]*sinh(GenInfo_eta[tgenIndex]),
		       GenInfo_mass[tgenIndex]);
	  Bgeny[typesize] = b4P->Rapidity();
	}
    }
}

bool signalGen(Int_t Btype, Int_t j)
{
  Double_t BId,MId,tk1Id,tk2Id;
  int twoTks;
  //tk1:positive, tk2:negtive
  if(Btype==1)
    {
      BId = 521;//B+-
      MId = -1;
      tk1Id = 321;//K+-
      tk2Id = -1;
      twoTks = 0;
    }
  if(Btype==2)
    {
      BId = 521;//B+-
      MId = -1;
      tk1Id = 211;//pi+-
      tk2Id = -1;
      twoTks = 0;
    }
  if(Btype==3)
    {
      BId = 511;//B0
      MId = 310;//Ks
      tk1Id = 211;//pi+
      tk2Id = -211;//pi-
      twoTks = 1;
    }
  if(Btype==4)
    {
      BId = 511;//B0
      MId = 313;//K*0
      tk1Id = 321;//K+
      tk2Id = -211;//pi-
      twoTks = 1;
    }
  if(Btype==5)
    {
      BId = 511;//B0
      MId = 313;//K*0
      tk1Id = -321;//pi+
      tk2Id = 211;//K-
      twoTks = 1;
    }
  if(Btype==6)
    {
      BId = 531;//Bs
      MId = 333;//phi
      tk1Id = 321;//K+
      tk2Id = -321;//K-
      twoTks = 1;
    }

  int flag=0;
  if (abs(GenInfo_pdgId[j])==BId&&GenInfo_nDa[j]==2&&GenInfo_da1[j]!=-1&&GenInfo_da2[j]!=-1)
    {
      if (abs(GenInfo_pdgId[GenInfo_da1[j]])==443)//jpsi
	{
	  if(GenInfo_da1[GenInfo_da1[j]]!=-1&&GenInfo_da2[GenInfo_da1[j]]!=-1)
	    {
	      if(abs(GenInfo_pdgId[GenInfo_da1[GenInfo_da1[j]]])==13&&abs(GenInfo_pdgId[GenInfo_da2[GenInfo_da1[j]]])==13)
		{
		  if(!twoTks)
		    {
		      if(abs(GenInfo_pdgId[GenInfo_da2[j]])==tk1Id) flag++;
		    }
		  else
		    {
		      if (abs(GenInfo_pdgId[GenInfo_da2[j]])==MId) 
			{
			  if(GenInfo_da1[GenInfo_da2[j]]!=-1 && GenInfo_da2[GenInfo_da2[j]]!=-1)
			    {
			      if(GenInfo_pdgId[GenInfo_da1[GenInfo_da2[j]]]==tk1Id && GenInfo_pdgId[GenInfo_da2[GenInfo_da2[j]]]==tk2Id) flag++;
			    }
			}
		    }
		}
	    }
	}
    }
  return flag;
}
