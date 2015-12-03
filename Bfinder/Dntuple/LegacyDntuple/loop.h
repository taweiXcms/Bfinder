#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TString.h>
#include <TNtuple.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <iostream>
#include <iomanip>
#include <cmath>

#define MUON_MASS   0.10565837
#define PION_MASS   0.13957018
#define KAON_MASS   0.493677
#define KSHORT_MASS 0.497614
#define KSTAR_MASS  0.89594
#define PHI_MASS    1.019455
#define JPSI_MASS   3.096916

Int_t PION_PDGID = 211;
Int_t KAON_PDGID = 321;
Int_t DZERO_PDGID = 421;
Int_t DPLUS_PDGID = 411;
Int_t DSUBS_PDGID = 431;

#define MAX_XB 16384
#define MAX_MUON 512
#define MAX_TRACK 4096
#define MAX_GEN 8192
#define MAX_BX 128
#define MAX_Vertices 4096
//#define N_TRIGGER_BOOKINGS 5842

//EvtInfo
Int_t      RunNo;
Int_t      EvtNo;
Int_t      Bsize;
Int_t      Dsize;
Double_t   PVx;
Double_t   PVy;
Double_t   PVz;
Double_t   PVxE;
Double_t   PVyE;
Double_t   PVzE;
Double_t   PVnchi2;
Double_t   PVchi2;
Double_t   BSx;
Double_t   BSy;
Double_t   BSz;
Double_t   BSxErr;
Double_t   BSyErr;
Double_t   BSzErr;
Double_t   BSdxdz;
Double_t   BSdydz;
Double_t   BSdxdzErr;
Double_t   BSdydzErr;
Double_t   BSWidthX;
Double_t   BSWidthXErr;
Double_t   BSWidthY;
Double_t   BSWidthYErr;
//DInfo
Int_t      Dindex[MAX_XB];
Int_t      Dtype[MAX_XB];
Double_t   Dmass[MAX_XB];
Double_t   Dpt[MAX_XB];
Double_t   Deta[MAX_XB];
Double_t   Dphi[MAX_XB];
Double_t   Dy[MAX_XB];
Double_t   DvtxX[MAX_XB];
Double_t   DvtxY[MAX_XB];
Double_t   Dd0[MAX_XB];
Double_t   Dd0Err[MAX_XB];
Double_t   Ddxyz[MAX_XB];
Double_t   DdxyzErr[MAX_XB];
Double_t   Dchi2ndf[MAX_XB];
Double_t   Dchi2cl[MAX_XB];
Double_t   Ddtheta[MAX_XB];
Double_t   Dlxy[MAX_XB];
Double_t   Dalpha[MAX_XB];
Double_t   DsvpvDistance[MAX_XB];
Double_t   DsvpvDisErr[MAX_XB];
Double_t   DsvpvDistance_2D[MAX_XB];
Double_t   DsvpvDisErr_2D[MAX_XB];
Double_t   DlxyBS[MAX_XB];
Double_t   DlxyBSErr[MAX_XB];
Double_t   DMaxDoca[MAX_XB];
Bool_t     Dmaxpt[MAX_XB];
Bool_t     Dmaxprob[MAX_XB];
Bool_t     DmaxptMatched[MAX_XB];
Bool_t     DmaxprobMatched[MAX_XB];

//DInfo.trkInfo
Double_t   Dtrk1Idx[MAX_XB];
Double_t   Dtrk2Idx[MAX_XB];
Double_t   Dtrk3Idx[MAX_XB];
Double_t   Dtrk4Idx[MAX_XB];
Double_t   Dtrk1Pt[MAX_XB];
Double_t   Dtrk2Pt[MAX_XB];
Double_t   Dtrk3Pt[MAX_XB];
Double_t   Dtrk4Pt[MAX_XB];
Double_t   Dtrk1Eta[MAX_XB];
Double_t   Dtrk2Eta[MAX_XB];
Double_t   Dtrk3Eta[MAX_XB];
Double_t   Dtrk4Eta[MAX_XB];
Double_t   Dtrk1Phi[MAX_XB];
Double_t   Dtrk2Phi[MAX_XB];
Double_t   Dtrk3Phi[MAX_XB];
Double_t   Dtrk4Phi[MAX_XB];
Double_t   Dtrk1PtErr[MAX_XB];
Double_t   Dtrk2PtErr[MAX_XB];
Double_t   Dtrk3PtErr[MAX_XB];
Double_t   Dtrk4PtErr[MAX_XB];
Double_t   Dtrk1EtaErr[MAX_XB];
Double_t   Dtrk2EtaErr[MAX_XB];
Double_t   Dtrk3EtaErr[MAX_XB];
Double_t   Dtrk4EtaErr[MAX_XB];
Double_t   Dtrk1PhiErr[MAX_XB];
Double_t   Dtrk2PhiErr[MAX_XB];
Double_t   Dtrk3PhiErr[MAX_XB];
Double_t   Dtrk4PhiErr[MAX_XB];
Double_t   Dtrk1Y[MAX_XB];
Double_t   Dtrk2Y[MAX_XB];
Double_t   Dtrk3Y[MAX_XB];
Double_t   Dtrk4Y[MAX_XB];
Double_t   Dtrk1Dxy[MAX_XB];
Double_t   Dtrk2Dxy[MAX_XB];
Double_t   Dtrk3Dxy[MAX_XB];
Double_t   Dtrk4Dxy[MAX_XB];
Double_t   Dtrk1D0Err[MAX_XB];
Double_t   Dtrk2D0Err[MAX_XB];
Double_t   Dtrk3D0Err[MAX_XB];
Double_t   Dtrk4D0Err[MAX_XB];
Double_t   Dtrk1PixelHit[MAX_XB];
Double_t   Dtrk2PixelHit[MAX_XB];
Double_t   Dtrk3PixelHit[MAX_XB];
Double_t   Dtrk4PixelHit[MAX_XB];
Double_t   Dtrk1StripHit[MAX_XB];
Double_t   Dtrk2StripHit[MAX_XB];
Double_t   Dtrk3StripHit[MAX_XB];
Double_t   Dtrk4StripHit[MAX_XB];
Double_t   Dtrk1nStripLayer[MAX_XB];
Double_t   Dtrk2nStripLayer[MAX_XB];
Double_t   Dtrk3nStripLayer[MAX_XB];
Double_t   Dtrk4nStripLayer[MAX_XB];
Double_t   Dtrk1nPixelLayer[MAX_XB];
Double_t   Dtrk2nPixelLayer[MAX_XB];
Double_t   Dtrk3nPixelLayer[MAX_XB];
Double_t   Dtrk4nPixelLayer[MAX_XB];
Double_t   Dtrk1Chi2ndf[MAX_XB];
Double_t   Dtrk2Chi2ndf[MAX_XB];
Double_t   Dtrk3Chi2ndf[MAX_XB];
Double_t   Dtrk4Chi2ndf[MAX_XB];
Double_t   Dtrk1MassHypo[MAX_XB];
Double_t   Dtrk2MassHypo[MAX_XB];
Double_t   Dtrk3MassHypo[MAX_XB];
Double_t   Dtrk4MassHypo[MAX_XB];
Double_t   Dtrkminpt[MAX_XB];
Double_t   Dtrkmaxpt[MAX_XB];
Int_t      Dtrkminptindex[MAX_XB];
Int_t      Dtrkmaxptindex[MAX_XB];
Double_t   Dtrk1MVAVal[MAX_XB];
Double_t   Dtrk2MVAVal[MAX_XB];
Double_t   Dtrk3MVAVal[MAX_XB];
Double_t   Dtrk4MVAVal[MAX_XB];
Int_t      Dtrk1Algo[MAX_XB];
Int_t      Dtrk2Algo[MAX_XB];
Int_t      Dtrk3Algo[MAX_XB];
Int_t      Dtrk4Algo[MAX_XB];
Bool_t     Dtrk1highPurity[MAX_XB];
Bool_t     Dtrk2highPurity[MAX_XB];
Bool_t     Dtrk3highPurity[MAX_XB];
Bool_t     Dtrk4highPurity[MAX_XB];
Int_t      Dtrk1Quality[MAX_XB];
Int_t      Dtrk2Quality[MAX_XB];
Int_t      Dtrk3Quality[MAX_XB];
Int_t      Dtrk4Quality[MAX_XB];
//DInfo.tktkResInfo
Double_t   DtktkResmass[MAX_XB];
Double_t   DtktkRespt[MAX_XB];
Double_t   DtktkReseta[MAX_XB];
Double_t   DtktkResphi[MAX_XB];
//DInfo.genInfo
Double_t   Dgen[MAX_XB];
Int_t      DgennDa[MAX_XB];
Int_t      DgenIndex[MAX_XB];
Double_t   Dgenpt[MAX_XB];
Double_t   Dgeneta[MAX_XB];
Double_t   Dgenphi[MAX_XB];
Double_t   Dgeny[MAX_XB];

void buildDBranch(TTree* dnt)
{
  //EvtInfo
  dnt->Branch("RunNo",&RunNo);
  dnt->Branch("EvtNo",&EvtNo);
  dnt->Branch("Dsize",&Dsize);
  dnt->Branch("PVx",&PVx);
  dnt->Branch("PVy",&PVy);
  dnt->Branch("PVz",&PVz);
  dnt->Branch("PVxE",&PVxE);
  dnt->Branch("PVyE",&PVyE);
  dnt->Branch("PVzE",&PVzE);
  dnt->Branch("PVnchi2",&PVnchi2);
  dnt->Branch("BSx",&BSx);
  dnt->Branch("BSy",&BSy);
  dnt->Branch("BSz",&BSz);
  dnt->Branch("BSxErr",&BSxErr);
  dnt->Branch("BSyErr",&BSyErr);
  dnt->Branch("BSzErr",&BSzErr);
  dnt->Branch("BSdxdz",&BSdxdz);
  dnt->Branch("BSdydz",&BSdydz);
  dnt->Branch("BSdxdzErr",&BSdxdzErr);
  dnt->Branch("BSdydzErr",&BSdydzErr);
  dnt->Branch("BSWidthX",&BSWidthX);
  dnt->Branch("BSWidthXErr",&BSWidthXErr);
  dnt->Branch("BSWidthY",&BSWidthY);
  dnt->Branch("BSWidthYErr",&BSWidthYErr);
  //DInfo
  dnt->Branch("Dindex",Dindex,"Dindex[Dsize]/I");
  dnt->Branch("Dtype",Dtype,"Dtype[Dsize]/I");
  dnt->Branch("Dmass",Dmass,"Dmass[Dsize]/D");
  dnt->Branch("Dpt",Dpt,"Dpt[Dsize]/D");
  dnt->Branch("Deta",Deta,"Deta[Dsize]/D");
  dnt->Branch("Dphi",Dphi,"Dphi[Dsize]/D");
  dnt->Branch("Dy",Dy,"Dy[Dsize]/D");
  dnt->Branch("DvtxX",DvtxX,"DvtxX[Dsize]/D");
  dnt->Branch("DvtxY",DvtxY,"DvtxY[Dsize]/D");
  dnt->Branch("Dd0",Dd0,"Dd0[Dsize]/D");
  dnt->Branch("Dd0Err",Dd0Err,"Dd0Err[Dsize]/D");
  dnt->Branch("Ddxyz",Ddxyz,"Ddxyz[Dsize]/D");
  dnt->Branch("DdxyzErr",DdxyzErr,"DdxyzErr[Dsize]/D");
  dnt->Branch("Dchi2ndf",Dchi2ndf,"Dchi2ndf[Dsize]/D");
  dnt->Branch("Dchi2cl",Dchi2cl,"Dchi2cl[Dsize]/D");
  dnt->Branch("Ddtheta",Ddtheta,"Ddtheta[Dsize]/D");
  dnt->Branch("Dlxy",Dlxy,"Dlxy[Dsize]/D");
  dnt->Branch("Dalpha",Dalpha,"Dalpha[Dsize]/D");
  dnt->Branch("DsvpvDistance",DsvpvDistance,"DsvpvDistance[Dsize]/D");
  dnt->Branch("DsvpvDisErr",DsvpvDisErr,"DsvpvDisErr[Dsize]/D");
  dnt->Branch("DsvpvDistance_2D",DsvpvDistance_2D,"DsvpvDistance_2D[Dsize]/D");
  dnt->Branch("DsvpvDisErr_2D",DsvpvDisErr_2D,"DsvpvDisErr_2D[Dsize]/D");
  dnt->Branch("DlxyBS",DlxyBS,"DlxyBS[Dsize]/D");
  dnt->Branch("DlxyBSErr",DlxyBSErr,"DlxyBSErr[Dsize]/D");
  dnt->Branch("DMaxDoca",DMaxDoca,"DMaxDoca[Dsize]/D");
  dnt->Branch("Dmaxpt",Dmaxpt,"Dmaxpt[Dsize]/O");
  dnt->Branch("Dmaxprob",Dmaxprob,"Dmaxprob[Dsize]/O");
  dnt->Branch("DmaxptMatched",DmaxptMatched,"DmaxptMatched[Dsize]/O");
  dnt->Branch("DmaxprobMatched",DmaxprobMatched,"DmaxprobMatched[Dsize]/O");
  //DInfo.trkInfo
  dnt->Branch("Dtrk1Idx",Dtrk1Idx,"Dtrk1Idx[Dsize]/I");
  dnt->Branch("Dtrk2Idx",Dtrk2Idx,"Dtrk2Idx[Dsize]/I");
  dnt->Branch("Dtrk3Idx",Dtrk3Idx,"Dtrk3Idx[Dsize]/I");
  dnt->Branch("Dtrk4Idx",Dtrk4Idx,"Dtrk4Idx[Dsize]/I");
  dnt->Branch("Dtrk1Pt",Dtrk1Pt,"Dtrk1Pt[Dsize]/D");
  dnt->Branch("Dtrk2Pt",Dtrk2Pt,"Dtrk2Pt[Dsize]/D");
  dnt->Branch("Dtrk3Pt",Dtrk3Pt,"Dtrk3Pt[Dsize]/D");
  dnt->Branch("Dtrk4Pt",Dtrk4Pt,"Dtrk4Pt[Dsize]/D");
  dnt->Branch("Dtrk1Eta",Dtrk1Eta,"Dtrk1Eta[Dsize]/D");
  dnt->Branch("Dtrk2Eta",Dtrk2Eta,"Dtrk2Eta[Dsize]/D");
  dnt->Branch("Dtrk3Eta",Dtrk3Eta,"Dtrk3Eta[Dsize]/D");
  dnt->Branch("Dtrk4Eta",Dtrk4Eta,"Dtrk4Eta[Dsize]/D");
  dnt->Branch("Dtrk1Phi",Dtrk1Phi,"Dtrk1Phi[Dsize]/D");
  dnt->Branch("Dtrk2Phi",Dtrk2Phi,"Dtrk2Phi[Dsize]/D");
  dnt->Branch("Dtrk3Phi",Dtrk3Phi,"Dtrk3Phi[Dsize]/D");
  dnt->Branch("Dtrk4Phi",Dtrk4Phi,"Dtrk4Phi[Dsize]/D");
  dnt->Branch("Dtrk1PtErr",Dtrk1PtErr,"Dtrk1PtErr[Dsize]/D");
  dnt->Branch("Dtrk2PtErr",Dtrk2PtErr,"Dtrk2PtErr[Dsize]/D");
  dnt->Branch("Dtrk3PtErr",Dtrk3PtErr,"Dtrk3PtErr[Dsize]/D");
  dnt->Branch("Dtrk4PtErr",Dtrk4PtErr,"Dtrk4PtErr[Dsize]/D");
  dnt->Branch("Dtrk1EtaErr",Dtrk1EtaErr,"Dtrk1EtaErr[Dsize]/D");
  dnt->Branch("Dtrk2EtaErr",Dtrk2EtaErr,"Dtrk2EtaErr[Dsize]/D");
  dnt->Branch("Dtrk3EtaErr",Dtrk3EtaErr,"Dtrk3EtaErr[Dsize]/D");
  dnt->Branch("Dtrk4EtaErr",Dtrk4EtaErr,"Dtrk4EtaErr[Dsize]/D");
  dnt->Branch("Dtrk1PhiErr",Dtrk1PhiErr,"Dtrk1PhiErr[Dsize]/D");
  dnt->Branch("Dtrk2PhiErr",Dtrk2PhiErr,"Dtrk2PhiErr[Dsize]/D");
  dnt->Branch("Dtrk3PhiErr",Dtrk3PhiErr,"Dtrk3PhiErr[Dsize]/D");
  dnt->Branch("Dtrk4PhiErr",Dtrk4PhiErr,"Dtrk4PhiErr[Dsize]/D");
  dnt->Branch("Dtrk1Y",Dtrk1Y,"Dtrk1Y[Dsize]/D");
  dnt->Branch("Dtrk2Y",Dtrk2Y,"Dtrk2Y[Dsize]/D");
  dnt->Branch("Dtrk3Y",Dtrk3Y,"Dtrk3Y[Dsize]/D");
  dnt->Branch("Dtrk4Y",Dtrk4Y,"Dtrk4Y[Dsize]/D");
  dnt->Branch("Dtrk1Dxy",Dtrk1Dxy,"Dtrk1Dxy[Dsize]/D");
  dnt->Branch("Dtrk2Dxy",Dtrk2Dxy,"Dtrk2Dxy[Dsize]/D");
  dnt->Branch("Dtrk3Dxy",Dtrk3Dxy,"Dtrk3Dxy[Dsize]/D");
  dnt->Branch("Dtrk4Dxy",Dtrk4Dxy,"Dtrk4Dxy[Dsize]/D");
  dnt->Branch("Dtrk1D0Err",Dtrk1D0Err,"Dtrk1D0Err[Dsize]/D");
  dnt->Branch("Dtrk2D0Err",Dtrk2D0Err,"Dtrk2D0Err[Dsize]/D");
  dnt->Branch("Dtrk3D0Err",Dtrk3D0Err,"Dtrk3D0Err[Dsize]/D");
  dnt->Branch("Dtrk4D0Err",Dtrk4D0Err,"Dtrk4D0Err[Dsize]/D");
  dnt->Branch("Dtrk1PixelHit",Dtrk1PixelHit,"Dtrk1PixelHit[Dsize]/D");
  dnt->Branch("Dtrk2PixelHit",Dtrk2PixelHit,"Dtrk2PixelHit[Dsize]/D");
  dnt->Branch("Dtrk3PixelHit",Dtrk3PixelHit,"Dtrk3PixelHit[Dsize]/D");
  dnt->Branch("Dtrk4PixelHit",Dtrk4PixelHit,"Dtrk4PixelHit[Dsize]/D");
  dnt->Branch("Dtrk1StripHit",Dtrk1StripHit,"Dtrk1StripHit[Dsize]/D");
  dnt->Branch("Dtrk2StripHit",Dtrk2StripHit,"Dtrk2StripHit[Dsize]/D");
  dnt->Branch("Dtrk3StripHit",Dtrk3StripHit,"Dtrk3StripHit[Dsize]/D");
  dnt->Branch("Dtrk4StripHit",Dtrk4StripHit,"Dtrk4StripHit[Dsize]/D");
  dnt->Branch("Dtrk1nStripLayer",Dtrk1nStripLayer,"Dtrk1nStripLayer[Dsize]/D");
  dnt->Branch("Dtrk2nStripLayer",Dtrk2nStripLayer,"Dtrk2nStripLayer[Dsize]/D");
  dnt->Branch("Dtrk3nStripLayer",Dtrk3nStripLayer,"Dtrk3nStripLayer[Dsize]/D");
  dnt->Branch("Dtrk4nStripLayer",Dtrk4nStripLayer,"Dtrk4nStripLayer[Dsize]/D");
  dnt->Branch("Dtrk1nPixelLayer",Dtrk1nPixelLayer,"Dtrk1nPixelLayer[Dsize]/D");
  dnt->Branch("Dtrk2nPixelLayer",Dtrk2nPixelLayer,"Dtrk2nPixelLayer[Dsize]/D");
  dnt->Branch("Dtrk3nPixelLayer",Dtrk3nPixelLayer,"Dtrk3nPixelLayer[Dsize]/D");
  dnt->Branch("Dtrk4nPixelLayer",Dtrk4nPixelLayer,"Dtrk4nPixelLayer[Dsize]/D");
  dnt->Branch("Dtrk1Chi2ndf",Dtrk1Chi2ndf,"Dtrk1Chi2ndf[Dsize]/D");
  dnt->Branch("Dtrk2Chi2ndf",Dtrk2Chi2ndf,"Dtrk2Chi2ndf[Dsize]/D");
  dnt->Branch("Dtrk3Chi2ndf",Dtrk3Chi2ndf,"Dtrk3Chi2ndf[Dsize]/D");
  dnt->Branch("Dtrk4Chi2ndf",Dtrk4Chi2ndf,"Dtrk4Chi2ndf[Dsize]/D");
  dnt->Branch("Dtrk1MassHypo",Dtrk1MassHypo,"Dtrk1MassHypo[Dsize]/D");
  dnt->Branch("Dtrk2MassHypo",Dtrk2MassHypo,"Dtrk2MassHypo[Dsize]/D");
  dnt->Branch("Dtrk3MassHypo",Dtrk3MassHypo,"Dtrk3MassHypo[Dsize]/D");
  dnt->Branch("Dtrk4MassHypo",Dtrk4MassHypo,"Dtrk4MassHypo[Dsize]/D");
  dnt->Branch("Dtrkminpt",Dtrkminpt,"Dtrkminpt[Dsize]/D");
  dnt->Branch("Dtrkmaxpt",Dtrkmaxpt,"Dtrkmaxpt[Dsize]/D");
  dnt->Branch("Dtrkminptindex",Dtrkminptindex,"Dtrkminptindex[Dsize]/I");
  dnt->Branch("Dtrkmaxptindex",Dtrkmaxptindex,"Dtrkmaxptindex[Dsize]/I");
  dnt->Branch("Dtrk1MVAVal",Dtrk1MVAVal,"Dtrk1MVAVal[Dsize]/D");
  dnt->Branch("Dtrk2MVAVal",Dtrk2MVAVal,"Dtrk2MVAVal[Dsize]/D");
  dnt->Branch("Dtrk3MVAVal",Dtrk3MVAVal,"Dtrk3MVAVal[Dsize]/D");
  dnt->Branch("Dtrk4MVAVal",Dtrk4MVAVal,"Dtrk4MVAVal[Dsize]/D");
  dnt->Branch("Dtrk1Algo",Dtrk1Algo,"Dtrk1Algo[Dsize]/I");
  dnt->Branch("Dtrk2Algo",Dtrk2Algo,"Dtrk2Algo[Dsize]/I");
  dnt->Branch("Dtrk3Algo",Dtrk3Algo,"Dtrk3Algo[Dsize]/I");
  dnt->Branch("Dtrk4Algo",Dtrk4Algo,"Dtrk4Algo[Dsize]/I");
  dnt->Branch("Dtrk1highPurity",Dtrk1highPurity,"Dtrk1highPurity[Dsize]/O");
  dnt->Branch("Dtrk2highPurity",Dtrk2highPurity,"Dtrk2highPurity[Dsize]/O");
  dnt->Branch("Dtrk3highPurity",Dtrk3highPurity,"Dtrk3highPurity[Dsize]/O");
  dnt->Branch("Dtrk4highPurity",Dtrk4highPurity,"Dtrk4highPurity[Dsize]/O");
  dnt->Branch("Dtrk1Quality",Dtrk1Quality,"Dtrk1Quality[Dsize]/I");
  dnt->Branch("Dtrk2Quality",Dtrk2Quality,"Dtrk2Quality[Dsize]/I");
  dnt->Branch("Dtrk3Quality",Dtrk3Quality,"Dtrk3Quality[Dsize]/I");
  dnt->Branch("Dtrk4Quality",Dtrk4Quality,"Dtrk4Quality[Dsize]/I");
  //DInfo.tktkResInfo
  dnt->Branch("DtktkResmass",DtktkResmass,"DtktkResmass[Dsize]/D");
  dnt->Branch("DtktkRespt",DtktkRespt,"DtktkRespt[Dsize]/D");
  dnt->Branch("DtktkReseta",DtktkReseta,"DtktkReseta[Dsize]/D");
  dnt->Branch("DtktkResphi",DtktkResphi,"DtktkResphi[Dsize]/D");
  //DInfo.genInfo
  dnt->Branch("Dgen",Dgen,"Dgen[Dsize]/D");
  dnt->Branch("DgenIndex",DgenIndex,"DgenIndex[Dsize]/I");
  dnt->Branch("DgennDa",DgennDa,"DgennDa[Dsize]/I");
  dnt->Branch("Dgenpt",Dgenpt,"Dgenpt[Dsize]/D");
  dnt->Branch("Dgeneta",Dgeneta,"Dgeneta[Dsize]/D");
  dnt->Branch("Dgenphi",Dgenphi,"Dgenphi[Dsize]/D");
  dnt->Branch("Dgeny",Dgeny,"Dgeny[Dsize]/D");
}

//GenInfo
Int_t      Gsize;
Double_t   Gy[MAX_GEN];
Double_t   Geta[MAX_GEN];
Double_t   Gphi[MAX_GEN];
Double_t   Gpt[MAX_GEN];
Double_t   GpdgId[MAX_GEN];
Int_t      GisSignal[MAX_GEN];

void buildGenBranch(TTree* nt)
{
  nt->Branch("Gsize",&Gsize);
  nt->Branch("Gy",Gy,"Gy[Gsize]/D");
  nt->Branch("Geta",Geta,"Geta[Gsize]/D");
  nt->Branch("Gphi",Gphi,"Gphi[Gsize]/D");
  nt->Branch("Gpt",Gpt,"Gpt[Gsize]/D");
  nt->Branch("GpdgId",GpdgId,"GpdgId[Gsize]/D");
  nt->Branch("GisSignal",GisSignal,"GisSignal[Gsize]/I");
}


//EvtInfo
Int_t           EvtInfo_RunNo;
Int_t           EvtInfo_EvtNo;
Int_t           EvtInfo_BxNo;
Int_t           EvtInfo_LumiNo;
Int_t           EvtInfo_Orbit;
Bool_t          EvtInfo_McFlag;
Int_t           EvtInfo_nBX;
Int_t           EvtInfo_BXPU[MAX_BX];
Int_t           EvtInfo_nPU[MAX_BX];
Float_t         EvtInfo_trueIT[MAX_BX];
Double_t        EvtInfo_PVx;
Double_t        EvtInfo_PVy;
Double_t        EvtInfo_PVz;
Double_t        EvtInfo_PVxE;
Double_t        EvtInfo_PVyE;
Double_t        EvtInfo_PVzE;
Double_t        EvtInfo_PVnchi2;
Double_t        EvtInfo_PVchi2;
Double_t        EvtInfo_BSx;
Double_t        EvtInfo_BSy;
Double_t        EvtInfo_BSz;
Double_t        EvtInfo_BSxErr;
Double_t        EvtInfo_BSyErr;
Double_t        EvtInfo_BSzErr;
Double_t        EvtInfo_BSdxdz;
Double_t        EvtInfo_BSdydz;
Double_t        EvtInfo_BSdxdzErr;
Double_t        EvtInfo_BSdydzErr;
Double_t        EvtInfo_BSWidthX;
Double_t        EvtInfo_BSWidthXErr;
Double_t        EvtInfo_BSWidthY;
Double_t        EvtInfo_BSWidthYErr;
//TrackInfo
Int_t           TrackInfo_size;
Int_t           TrackInfo_index[MAX_TRACK];
Int_t           TrackInfo_handle_index[MAX_TRACK];
Int_t           TrackInfo_charge[MAX_TRACK];
Double_t        TrackInfo_pt[MAX_TRACK];
Double_t        TrackInfo_eta[MAX_TRACK];
Double_t        TrackInfo_phi[MAX_TRACK];
Double_t        TrackInfo_ptErr[MAX_TRACK];
Double_t        TrackInfo_etaErr[MAX_TRACK];
Double_t        TrackInfo_phiErr[MAX_TRACK];
Int_t           TrackInfo_striphit[MAX_TRACK];
Int_t           TrackInfo_pixelhit[MAX_TRACK];
Int_t           TrackInfo_nStripLayer[MAX_TRACK];
Int_t           TrackInfo_nPixelLayer[MAX_TRACK];
Int_t           TrackInfo_fpbarrelhit[MAX_TRACK];
Int_t           TrackInfo_fpendcaphit[MAX_TRACK];
Double_t        TrackInfo_chi2[MAX_TRACK];
Double_t        TrackInfo_ndf[MAX_TRACK];
Double_t        TrackInfo_d0[MAX_TRACK];
Double_t        TrackInfo_d0error[MAX_TRACK];
Double_t        TrackInfo_dzPV[MAX_TRACK];
Double_t        TrackInfo_dxyPV[MAX_TRACK];
Int_t           TrackInfo_geninfo_index[MAX_TRACK];
Int_t           TrackInfo_trackQuality[MAX_TRACK];
Bool_t          TrackInfo_highPurity[MAX_TRACK];
Double_t        TrackInfo_trkMVAVal[MAX_TRACK];
Int_t           TrackInfo_trkAlgo[MAX_TRACK];

//DInfo
Int_t           DInfo_size;
Int_t           DInfo_index[MAX_XB];
Int_t           DInfo_type[MAX_XB];    
Double_t        DInfo_b4fit_mass[MAX_XB];
Double_t        DInfo_b4fit_pt[MAX_XB];
Double_t        DInfo_b4fit_eta[MAX_XB];
Double_t        DInfo_b4fit_phi[MAX_XB];
Double_t        DInfo_tktkRes_mass[MAX_XB];
Double_t        DInfo_tktkRes_pt[MAX_XB];
Double_t        DInfo_tktkRes_eta[MAX_XB];
Double_t        DInfo_tktkRes_phi[MAX_XB];
Double_t        DInfo_tktkRes_px[MAX_XB];
Double_t        DInfo_tktkRes_py[MAX_XB];
Double_t        DInfo_tktkRes_pz[MAX_XB];
Double_t        DInfo_tktkRes_vtxX[MAX_XB];
Double_t        DInfo_tktkRes_vtxY[MAX_XB];
Double_t        DInfo_tktkRes_vtxZ[MAX_XB];
Double_t        DInfo_tktkRes_vtxXErr[MAX_XB];
Double_t        DInfo_tktkRes_vtxYErr[MAX_XB];
Double_t        DInfo_tktkRes_vtxZErr[MAX_XB];
Double_t        DInfo_tktkRes_vtxdof[MAX_XB];
Double_t        DInfo_tktkRes_vtxchi2[MAX_XB];
Double_t        DInfo_tktkRes_rftk1_pt[MAX_XB];
Double_t        DInfo_tktkRes_rftk1_eta[MAX_XB];
Double_t        DInfo_tktkRes_rftk1_phi[MAX_XB];
Double_t        DInfo_tktkRes_rftk2_pt[MAX_XB];
Double_t        DInfo_tktkRes_rftk2_eta[MAX_XB];
Double_t        DInfo_tktkRes_rftk2_phi[MAX_XB];
Double_t        DInfo_mass[MAX_XB];
Double_t        DInfo_pt[MAX_XB];
Double_t        DInfo_eta[MAX_XB];
Double_t        DInfo_phi[MAX_XB];
Double_t        DInfo_px[MAX_XB];
Double_t        DInfo_py[MAX_XB];
Double_t        DInfo_pz[MAX_XB];
Double_t        DInfo_alpha[MAX_XB];
Double_t        DInfo_svpvDistance[MAX_XB];
Double_t        DInfo_svpvDisErr[MAX_XB];
Double_t        DInfo_svpvDistance_2D[MAX_XB];
Double_t        DInfo_svpvDisErr_2D[MAX_XB];
Double_t        DInfo_MaxDoca[MAX_XB];
Double_t        DInfo_vtxX[MAX_XB];
Double_t        DInfo_vtxY[MAX_XB];
Double_t        DInfo_vtxZ[MAX_XB];
Double_t        DInfo_vtxXErr[MAX_XB];
Double_t        DInfo_vtxYErr[MAX_XB];
Double_t        DInfo_vtxZErr[MAX_XB];
Double_t        DInfo_vtxYXErr[MAX_XB];
Double_t        DInfo_vtxZXErr[MAX_XB];
Double_t        DInfo_vtxZYErr[MAX_XB];
Double_t        DInfo_vtxdof[MAX_XB];
Double_t        DInfo_vtxchi2[MAX_XB];
Double_t        DInfo_rftk1_px[MAX_XB];
Double_t        DInfo_rftk1_py[MAX_XB];
Double_t        DInfo_rftk1_pz[MAX_XB];
Double_t        DInfo_rftk2_px[MAX_XB];
Double_t        DInfo_rftk2_py[MAX_XB];
Double_t        DInfo_rftk2_pz[MAX_XB];
Double_t        DInfo_rftk3_px[MAX_XB];
Double_t        DInfo_rftk3_py[MAX_XB];
Double_t        DInfo_rftk3_pz[MAX_XB];
Double_t        DInfo_rftk4_px[MAX_XB];
Double_t        DInfo_rftk4_py[MAX_XB];
Double_t        DInfo_rftk4_pz[MAX_XB];
Double_t        DInfo_rftk1_pt[MAX_XB];
Double_t        DInfo_rftk1_eta[MAX_XB];
Double_t        DInfo_rftk1_phi[MAX_XB];
Double_t        DInfo_rftk2_pt[MAX_XB];
Double_t        DInfo_rftk2_eta[MAX_XB];
Double_t        DInfo_rftk2_phi[MAX_XB];
Double_t        DInfo_rftk3_pt[MAX_XB];
Double_t        DInfo_rftk3_eta[MAX_XB];
Double_t        DInfo_rftk3_phi[MAX_XB];
Double_t        DInfo_rftk4_pt[MAX_XB];
Double_t        DInfo_rftk4_eta[MAX_XB];
Double_t        DInfo_rftk4_phi[MAX_XB];
Int_t           DInfo_rftk1_index[MAX_XB];
Int_t           DInfo_rftk2_index[MAX_XB];
Int_t           DInfo_rftk3_index[MAX_XB];
Int_t           DInfo_rftk4_index[MAX_XB];
Int_t           DInfo_rftk1_MassHypo[MAX_XB];
Int_t           DInfo_rftk2_MassHypo[MAX_XB];
Int_t           DInfo_rftk3_MassHypo[MAX_XB];
Int_t           DInfo_rftk4_MassHypo[MAX_XB];
//GenInfo
Int_t           GenInfo_size;
Int_t           GenInfo_index[MAX_GEN];
Int_t           GenInfo_handle_index[MAX_GEN];
Double_t        GenInfo_pt[MAX_GEN];
Double_t        GenInfo_eta[MAX_GEN];
Double_t        GenInfo_phi[MAX_GEN];
Double_t        GenInfo_mass[MAX_GEN];
Int_t           GenInfo_pdgId[MAX_GEN];
Int_t           GenInfo_status[MAX_GEN];
Int_t           GenInfo_nMo[MAX_GEN];
Int_t           GenInfo_nDa[MAX_GEN];
Int_t           GenInfo_mo1[MAX_GEN];
Int_t           GenInfo_mo2[MAX_GEN];
Int_t           GenInfo_da1[MAX_GEN];
Int_t           GenInfo_da2[MAX_GEN];

void setDBranch(TTree *root)
{
  //EvtInfo
  root->SetBranchAddress("EvtInfo.RunNo",&EvtInfo_RunNo);
  root->SetBranchAddress("EvtInfo.EvtNo",&EvtInfo_EvtNo);
  root->SetBranchAddress("EvtInfo.BxNo",&EvtInfo_BxNo);
  root->SetBranchAddress("EvtInfo.LumiNo",&EvtInfo_LumiNo);
  root->SetBranchAddress("EvtInfo.Orbit",&EvtInfo_Orbit);
  root->SetBranchAddress("EvtInfo.McFlag",&EvtInfo_McFlag);
  root->SetBranchAddress("EvtInfo.nBX",&EvtInfo_nBX);
  root->SetBranchAddress("EvtInfo.BXPU",EvtInfo_BXPU);
  root->SetBranchAddress("EvtInfo.nPU",EvtInfo_nPU);
  root->SetBranchAddress("EvtInfo.trueIT",EvtInfo_trueIT);
  root->SetBranchAddress("EvtInfo.PVx",&EvtInfo_PVx);
  root->SetBranchAddress("EvtInfo.PVy",&EvtInfo_PVy);
  root->SetBranchAddress("EvtInfo.PVz",&EvtInfo_PVz);
  root->SetBranchAddress("EvtInfo.PVxE",&EvtInfo_PVxE);
  root->SetBranchAddress("EvtInfo.PVyE",&EvtInfo_PVyE);
  root->SetBranchAddress("EvtInfo.PVzE",&EvtInfo_PVzE);
  root->SetBranchAddress("EvtInfo.PVnchi2",&EvtInfo_PVnchi2);
  root->SetBranchAddress("EvtInfo.PVchi2",&EvtInfo_PVchi2);
  root->SetBranchAddress("EvtInfo.BSx",&EvtInfo_BSx);
  root->SetBranchAddress("EvtInfo.BSy",&EvtInfo_BSy);
  root->SetBranchAddress("EvtInfo.BSz",&EvtInfo_BSz);
  root->SetBranchAddress("EvtInfo.BSxErr",&EvtInfo_BSxErr);
  root->SetBranchAddress("EvtInfo.BSyErr",&EvtInfo_BSyErr);
  root->SetBranchAddress("EvtInfo.BSzErr",&EvtInfo_BSzErr);
  root->SetBranchAddress("EvtInfo.BSdxdz",&EvtInfo_BSdxdz);
  root->SetBranchAddress("EvtInfo.BSdydz",&EvtInfo_BSdydz);
  root->SetBranchAddress("EvtInfo.BSdxdzErr",&EvtInfo_BSdxdzErr);
  root->SetBranchAddress("EvtInfo.BSdydzErr",&EvtInfo_BSdydzErr);
  root->SetBranchAddress("EvtInfo.BSWidthX",&EvtInfo_BSWidthX);
  root->SetBranchAddress("EvtInfo.BSWidthXErr",&EvtInfo_BSWidthXErr);
  root->SetBranchAddress("EvtInfo.BSWidthY",&EvtInfo_BSWidthY);
  root->SetBranchAddress("EvtInfo.BSWidthYErr",&EvtInfo_BSWidthYErr);
  //TrackInfo
  root->SetBranchAddress("TrackInfo.size",&TrackInfo_size);
  root->SetBranchAddress("TrackInfo.index",TrackInfo_index);
  root->SetBranchAddress("TrackInfo.handle_index",TrackInfo_handle_index);
  root->SetBranchAddress("TrackInfo.charge",TrackInfo_charge);
  root->SetBranchAddress("TrackInfo.pt",TrackInfo_pt);
  root->SetBranchAddress("TrackInfo.eta",TrackInfo_eta);
  root->SetBranchAddress("TrackInfo.phi",TrackInfo_phi);
  root->SetBranchAddress("TrackInfo.ptErr",TrackInfo_ptErr);
  root->SetBranchAddress("TrackInfo.etaErr",TrackInfo_etaErr);
  root->SetBranchAddress("TrackInfo.phiErr",TrackInfo_phiErr);
  root->SetBranchAddress("TrackInfo.striphit",TrackInfo_striphit);
  root->SetBranchAddress("TrackInfo.pixelhit",TrackInfo_pixelhit);
  root->SetBranchAddress("TrackInfo.nStripLayer",TrackInfo_nStripLayer);
  root->SetBranchAddress("TrackInfo.nPixelLayer",TrackInfo_nPixelLayer);
  root->SetBranchAddress("TrackInfo.fpbarrelhit",TrackInfo_fpbarrelhit);
  root->SetBranchAddress("TrackInfo.fpendcaphit",TrackInfo_fpendcaphit);
  root->SetBranchAddress("TrackInfo.chi2",TrackInfo_chi2);
  root->SetBranchAddress("TrackInfo.ndf",TrackInfo_ndf);
  root->SetBranchAddress("TrackInfo.d0",TrackInfo_d0);
  root->SetBranchAddress("TrackInfo.d0error",TrackInfo_d0error);
  root->SetBranchAddress("TrackInfo.dzPV",TrackInfo_dzPV);
  root->SetBranchAddress("TrackInfo.dxyPV",TrackInfo_dxyPV);
  root->SetBranchAddress("TrackInfo.geninfo_index",TrackInfo_geninfo_index);
  root->SetBranchAddress("TrackInfo.trackQuality",TrackInfo_trackQuality);
  root->SetBranchAddress("TrackInfo.highPurity",TrackInfo_highPurity);
  root->SetBranchAddress("TrackInfo.trkMVAVal",TrackInfo_trkMVAVal);
  root->SetBranchAddress("TrackInfo.trkAlgo",TrackInfo_trkAlgo);
  //DInfo
  root->SetBranchAddress("DInfo.size",&DInfo_size);
  root->SetBranchAddress("DInfo.index",DInfo_index);
  root->SetBranchAddress("DInfo.type",DInfo_type);
  root->SetBranchAddress("DInfo.tktkRes_mass",DInfo_tktkRes_mass);
  root->SetBranchAddress("DInfo.tktkRes_pt",DInfo_tktkRes_pt);
  root->SetBranchAddress("DInfo.tktkRes_eta",DInfo_tktkRes_eta);
  root->SetBranchAddress("DInfo.tktkRes_phi",DInfo_tktkRes_phi);
  root->SetBranchAddress("DInfo.tktkRes_px",DInfo_tktkRes_px);
  root->SetBranchAddress("DInfo.tktkRes_py",DInfo_tktkRes_py);
  root->SetBranchAddress("DInfo.tktkRes_pz",DInfo_tktkRes_pz);
  root->SetBranchAddress("DInfo.tktkRes_vtxX",DInfo_tktkRes_vtxX);
  root->SetBranchAddress("DInfo.tktkRes_vtxY",DInfo_tktkRes_vtxY);
  root->SetBranchAddress("DInfo.tktkRes_vtxZ",DInfo_tktkRes_vtxZ);
  root->SetBranchAddress("DInfo.tktkRes_vtxXErr",DInfo_tktkRes_vtxXErr);
  root->SetBranchAddress("DInfo.tktkRes_vtxYErr",DInfo_tktkRes_vtxYErr);
  root->SetBranchAddress("DInfo.tktkRes_vtxZErr",DInfo_tktkRes_vtxZErr);
  root->SetBranchAddress("DInfo.tktkRes_vtxdof",DInfo_tktkRes_vtxdof);
  root->SetBranchAddress("DInfo.tktkRes_vtxchi2",DInfo_tktkRes_vtxchi2);
  root->SetBranchAddress("DInfo.tktkRes_rftk1_pt",DInfo_tktkRes_rftk1_pt);
  root->SetBranchAddress("DInfo.tktkRes_rftk1_eta",DInfo_tktkRes_rftk1_eta);
  root->SetBranchAddress("DInfo.tktkRes_rftk1_phi",DInfo_tktkRes_rftk1_phi);
  root->SetBranchAddress("DInfo.tktkRes_rftk2_pt",DInfo_tktkRes_rftk2_pt);
  root->SetBranchAddress("DInfo.tktkRes_rftk2_eta",DInfo_tktkRes_rftk2_eta);
  root->SetBranchAddress("DInfo.tktkRes_rftk2_phi",DInfo_tktkRes_rftk2_phi);
  root->SetBranchAddress("DInfo.mass",DInfo_mass);
  root->SetBranchAddress("DInfo.pt",DInfo_pt);
  root->SetBranchAddress("DInfo.eta",DInfo_eta);
  root->SetBranchAddress("DInfo.phi",DInfo_phi);
  root->SetBranchAddress("DInfo.px",DInfo_px);
  root->SetBranchAddress("DInfo.py",DInfo_py);
  root->SetBranchAddress("DInfo.pz",DInfo_pz);
  root->SetBranchAddress("DInfo.alpha",DInfo_alpha);
  root->SetBranchAddress("DInfo.svpvDistance",DInfo_svpvDistance);
  root->SetBranchAddress("DInfo.svpvDisErr",DInfo_svpvDisErr);
  root->SetBranchAddress("DInfo.svpvDistance_2D",DInfo_svpvDistance_2D);
  root->SetBranchAddress("DInfo.svpvDisErr_2D",DInfo_svpvDisErr_2D);
  root->SetBranchAddress("DInfo.MaxDoca",DInfo_MaxDoca);
  root->SetBranchAddress("DInfo.vtxX",DInfo_vtxX);
  root->SetBranchAddress("DInfo.vtxY",DInfo_vtxY);
  root->SetBranchAddress("DInfo.vtxZ",DInfo_vtxZ);
  root->SetBranchAddress("DInfo.vtxXErr",DInfo_vtxXErr);
  root->SetBranchAddress("DInfo.vtxYErr",DInfo_vtxYErr);
  root->SetBranchAddress("DInfo.vtxZErr",DInfo_vtxZErr);
  root->SetBranchAddress("DInfo.vtxYXErr",DInfo_vtxYXErr);
  root->SetBranchAddress("DInfo.vtxZXErr",DInfo_vtxZXErr);
  root->SetBranchAddress("DInfo.vtxZYErr",DInfo_vtxZYErr);
  root->SetBranchAddress("DInfo.vtxdof",DInfo_vtxdof);
  root->SetBranchAddress("DInfo.vtxchi2",DInfo_vtxchi2);
  root->SetBranchAddress("DInfo.rftk1_px",DInfo_rftk1_px);
  root->SetBranchAddress("DInfo.rftk1_py",DInfo_rftk1_py);
  root->SetBranchAddress("DInfo.rftk1_pz",DInfo_rftk1_pz);
  root->SetBranchAddress("DInfo.rftk2_px",DInfo_rftk2_px);
  root->SetBranchAddress("DInfo.rftk2_py",DInfo_rftk2_py);
  root->SetBranchAddress("DInfo.rftk2_pz",DInfo_rftk2_pz);
  root->SetBranchAddress("DInfo.rftk3_px",DInfo_rftk3_px);
  root->SetBranchAddress("DInfo.rftk3_py",DInfo_rftk3_py);
  root->SetBranchAddress("DInfo.rftk3_pz",DInfo_rftk3_pz);
  root->SetBranchAddress("DInfo.rftk4_px",DInfo_rftk4_px);
  root->SetBranchAddress("DInfo.rftk4_py",DInfo_rftk4_py);
  root->SetBranchAddress("DInfo.rftk4_pz",DInfo_rftk4_pz);
  root->SetBranchAddress("DInfo.rftk1_pt",DInfo_rftk1_pt);
  root->SetBranchAddress("DInfo.rftk1_eta",DInfo_rftk1_eta);
  root->SetBranchAddress("DInfo.rftk1_phi",DInfo_rftk1_phi);
  root->SetBranchAddress("DInfo.rftk2_pt",DInfo_rftk2_pt);
  root->SetBranchAddress("DInfo.rftk2_eta",DInfo_rftk2_eta);
  root->SetBranchAddress("DInfo.rftk2_phi",DInfo_rftk2_phi);
  root->SetBranchAddress("DInfo.rftk3_pt",DInfo_rftk3_pt);
  root->SetBranchAddress("DInfo.rftk3_eta",DInfo_rftk3_eta);
  root->SetBranchAddress("DInfo.rftk3_phi",DInfo_rftk3_phi);
  root->SetBranchAddress("DInfo.rftk4_pt",DInfo_rftk4_pt);
  root->SetBranchAddress("DInfo.rftk4_eta",DInfo_rftk4_eta);
  root->SetBranchAddress("DInfo.rftk4_phi",DInfo_rftk4_phi);
  root->SetBranchAddress("DInfo.rftk1_index",DInfo_rftk1_index);
  root->SetBranchAddress("DInfo.rftk2_index",DInfo_rftk2_index);
  root->SetBranchAddress("DInfo.rftk3_index",DInfo_rftk3_index);
  root->SetBranchAddress("DInfo.rftk4_index",DInfo_rftk4_index);
  root->SetBranchAddress("DInfo.rftk1_MassHypo",DInfo_rftk1_MassHypo);
  root->SetBranchAddress("DInfo.rftk2_MassHypo",DInfo_rftk2_MassHypo);
  root->SetBranchAddress("DInfo.rftk3_MassHypo",DInfo_rftk3_MassHypo);
  root->SetBranchAddress("DInfo.rftk4_MassHypo",DInfo_rftk4_MassHypo);
  //GenInfo
  root->SetBranchAddress("GenInfo.size",&GenInfo_size);
  root->SetBranchAddress("GenInfo.index",GenInfo_index);
  root->SetBranchAddress("GenInfo.handle_index",GenInfo_handle_index);
  root->SetBranchAddress("GenInfo.pt",GenInfo_pt);
  root->SetBranchAddress("GenInfo.eta",GenInfo_eta);
  root->SetBranchAddress("GenInfo.phi",GenInfo_phi);
  root->SetBranchAddress("GenInfo.mass",GenInfo_mass);
  root->SetBranchAddress("GenInfo.pdgId",GenInfo_pdgId);
  root->SetBranchAddress("GenInfo.status",GenInfo_status);
  root->SetBranchAddress("GenInfo.nMo",GenInfo_nMo);
  root->SetBranchAddress("GenInfo.nDa",GenInfo_nDa);
  root->SetBranchAddress("GenInfo.mo1",GenInfo_mo1);
  root->SetBranchAddress("GenInfo.mo2",GenInfo_mo2);
  root->SetBranchAddress("GenInfo.da1",GenInfo_da1);
  root->SetBranchAddress("GenInfo.da2",GenInfo_da2);
}

//HltInfo
Int_t           Df_HLT_Run;
ULong64_t       Df_HLT_Event;
Int_t           Df_HLT_LumiBlock;
void setHltTreeBranch(TTree* hltroot)
{
  hltroot->SetBranchAddress("Run",&Df_HLT_Run);
  hltroot->SetBranchAddress("Event",&Df_HLT_Event);
  hltroot->SetBranchAddress("LumiBlock",&Df_HLT_LumiBlock);
}

//hiEvtInfo
Int_t           Df_HiTree_Run;
Int_t           Df_HiTree_Evt;
Int_t           Df_HiTree_Lumi;
void setHiTreeBranch(TTree* hitreeroot)
{
  hitreeroot->SetBranchAddress("run",&Df_HiTree_Run);
  hitreeroot->SetBranchAddress("evt",&Df_HiTree_Evt);
  hitreeroot->SetBranchAddress("lumi",&Df_HiTree_Lumi);
}
