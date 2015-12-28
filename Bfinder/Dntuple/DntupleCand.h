// vim:set ts=4 sw=4 fdm=marker et:
using namespace std;

#ifndef _DNTUPLECAND_H_
#define _DNTUPLECAND_H_
#include "format.h"

class DntupleBranchesCand
{//{{{
 public:
  //EvtInfo
  int     RunNo;
  int     EvtNo;
  int     LumiNo;
  float   PVx;
  float   PVy;
  float   PVz;
  float   PVxE;
  float   PVyE;
  float   PVzE;
  float   PVnchi2;
  float   PVchi2;
  float   BSx;
  float   BSy;
  float   BSz;
  float   BSxErr;
  float   BSyErr;
  float   BSzErr;
  float   BSdxdz;
  float   BSdydz;
  float   BSdxdzErr;
  float   BSdydzErr;
  float   BSWidthX;
  float   BSWidthXErr;
  float   BSWidthY;
  float   BSWidthYErr;
  //DInfo
  int     Dindex;
  int     Dtype;
  float   Dmass;
  float   Dpt;
  float   Deta;
  float   Dphi;
  float   Dy;
  float   DvtxX;
  float   DvtxY;
  float   Dd0;
  float   Dd0Err;
  float   Ddxyz;
  float   DdxyzErr;
  float   Dchi2ndf;
  float   Dchi2cl;
  float   Ddtheta;
  float   Dlxy;
  float   Dalpha;
  float   DsvpvDistance;
  float   DsvpvDisErr;
  float   DsvpvDistance_2D;
  float   DsvpvDisErr_2D;
  float   DtktkRes_chi2ndf;
  float   DtktkRes_chi2cl;
  float   DtktkRes_alpha;
  float   DtktkRes_svpvDistance;
  float   DtktkRes_svpvDisErr;
  float   DlxyBS;
  float   DlxyBSErr;
  float   DMaxDoca;
  bool    Dmaxpt;
  bool    Dmaxprob;
  bool    DmaxptMatched;
  bool    DmaxprobMatched;
  //DInfo.trkInfo
  int     Dtrk1Idx;
  int     Dtrk2Idx;
  int     Dtrk3Idx;
  int     Dtrk4Idx;
  float   Dtrk1Pt;
  float   Dtrk2Pt;
  float   Dtrk3Pt;
  float   Dtrk4Pt;
  float   Dtrk1Eta;
  float   Dtrk2Eta;
  float   Dtrk3Eta;
  float   Dtrk4Eta;
  float   Dtrk1Phi;
  float   Dtrk2Phi;
  float   Dtrk3Phi;
  float   Dtrk4Phi;
  float   Dtrk1PtErr;
  float   Dtrk2PtErr;
  float   Dtrk3PtErr;
  float   Dtrk4PtErr;
  float   Dtrk1EtaErr;
  float   Dtrk2EtaErr;
  float   Dtrk3EtaErr;
  float   Dtrk4EtaErr;
  float   Dtrk1PhiErr;
  float   Dtrk2PhiErr;
  float   Dtrk3PhiErr;
  float   Dtrk4PhiErr;
  float   Dtrk1Y;
  float   Dtrk2Y;
  float   Dtrk3Y;
  float   Dtrk4Y;
  float   Dtrk1Dxy;
  float   Dtrk2Dxy;
  float   Dtrk3Dxy;
  float   Dtrk4Dxy;
  float   Dtrk1D0Err;
  float   Dtrk2D0Err;
  float   Dtrk3D0Err;
  float   Dtrk4D0Err;
  float   Dtrk1PixelHit;
  float   Dtrk2PixelHit;
  float   Dtrk3PixelHit;
  float   Dtrk4PixelHit;
  float   Dtrk1StripHit;
  float   Dtrk2StripHit;
  float   Dtrk3StripHit;
  float   Dtrk4StripHit;
  float   Dtrk1nStripLayer;
  float   Dtrk2nStripLayer;
  float   Dtrk3nStripLayer;
  float   Dtrk4nStripLayer;
  float   Dtrk1nPixelLayer;
  float   Dtrk2nPixelLayer;
  float   Dtrk3nPixelLayer;
  float   Dtrk4nPixelLayer;
  float   Dtrk1Chi2ndf;
  float   Dtrk2Chi2ndf;
  float   Dtrk3Chi2ndf;
  float   Dtrk4Chi2ndf;
  float   Dtrk1MassHypo;
  float   Dtrk2MassHypo;
  float   Dtrk3MassHypo;
  float   Dtrk4MassHypo;
  float   Dtrk1MVAVal;
  float   Dtrk2MVAVal;
  float   Dtrk3MVAVal;
  float   Dtrk4MVAVal;
  int     Dtrk1Algo;
  int     Dtrk2Algo;
  int     Dtrk3Algo;
  int     Dtrk4Algo;
  bool    Dtrk1highPurity;
  bool    Dtrk2highPurity;
  bool    Dtrk3highPurity;
  bool    Dtrk4highPurity;
  int     Dtrk1Quality;
  int     Dtrk2Quality;
  int     Dtrk3Quality;
  int     Dtrk4Quality;
  //DInfo.tktkResInfo
  float   DtktkResmass;
  float   DtktkRespt;
  float   DtktkReseta;
  float   DtktkResphi;
  float   DRestrk1Pt;
  float   DRestrk1Eta;
  float   DRestrk1Phi;
  float   DRestrk1Y;
  float   DRestrk1Dxy;
  float   DRestrk1D0Err;
  float   DRestrk2Pt;
  float   DRestrk2Eta;
  float   DRestrk2Phi;
  float   DRestrk2Y;
  float   DRestrk2Dxy;
  float   DRestrk2D0Err;
  float   DRestrk3Pt;
  float   DRestrk3Eta;
  float   DRestrk3Phi;
  float   DRestrk3Y;
  float   DRestrk3Dxy;
  float   DRestrk3D0Err;
  float   DRestrk4Pt;
  float   DRestrk4Eta;
  float   DRestrk4Phi;
  float   DRestrk4Y;
  float   DRestrk4Dxy;
  float   DRestrk4D0Err;
  //DInfo.genInfo
  float   Dgen;
  int     DgennDa;
  int     DgenIndex;
  float   Dgenpt;
  float   Dgeneta;
  float   Dgenphi;
  float   Dgeny;
  
  void buildDBranch(TTree* dnt)
  {
    //EvtInfo
    dnt->Branch("RunNo",&RunNo);
    dnt->Branch("EvtNo",&EvtNo);
    dnt->Branch("LumiNo",&LumiNo);
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
    dnt->Branch("Dindex",&Dindex);
    dnt->Branch("Dtype",&Dtype);
    dnt->Branch("Dmass",&Dmass);
    dnt->Branch("Dpt",&Dpt);
    dnt->Branch("Deta",&Deta);
    dnt->Branch("Dphi",&Dphi);
    dnt->Branch("Dy",&Dy);
    dnt->Branch("DvtxX",&DvtxX);
    dnt->Branch("DvtxY",&DvtxY);
    dnt->Branch("Dd0",&Dd0);
    dnt->Branch("Dd0Err",&Dd0Err);
    dnt->Branch("Ddxyz",&Ddxyz);
    dnt->Branch("DdxyzErr",&DdxyzErr);
    dnt->Branch("Dchi2ndf",&Dchi2ndf);
    dnt->Branch("Dchi2cl",&Dchi2cl);
    dnt->Branch("Ddtheta",&Ddtheta);
    dnt->Branch("Dlxy",&Dlxy);
    dnt->Branch("Dalpha",&Dalpha);
    dnt->Branch("DsvpvDistance",&DsvpvDistance);
    dnt->Branch("DsvpvDisErr",&DsvpvDisErr);
    dnt->Branch("DsvpvDistance_2D",&DsvpvDistance_2D);
    dnt->Branch("DsvpvDisErr_2D",&DsvpvDisErr_2D);
    dnt->Branch("DtktkRes_chi2ndf",&DtktkRes_chi2ndf);
    dnt->Branch("DtktkRes_chi2cl",&DtktkRes_chi2cl);
    dnt->Branch("DtktkRes_alpha",&DtktkRes_alpha);
    dnt->Branch("DtktkRes_svpvDistance",&DtktkRes_svpvDistance);
    dnt->Branch("DtktkRes_svpvDisErr",&DtktkRes_svpvDisErr);
    dnt->Branch("DlxyBS",&DlxyBS);
    dnt->Branch("DlxyBSErr",&DlxyBSErr);
    dnt->Branch("DMaxDoca",&DMaxDoca);
    dnt->Branch("Dmaxpt",&Dmaxpt);
    dnt->Branch("Dmaxprob",&Dmaxprob);
    dnt->Branch("DmaxptMatched",&DmaxptMatched);
    dnt->Branch("DmaxprobMatched",&DmaxprobMatched);
    //DInfo.trkInfo
    dnt->Branch("Dtrk1Idx",&Dtrk1Idx);
    dnt->Branch("Dtrk2Idx",&Dtrk2Idx);
    dnt->Branch("Dtrk3Idx",&Dtrk3Idx);
    dnt->Branch("Dtrk4Idx",&Dtrk4Idx);
    dnt->Branch("Dtrk1Pt",&Dtrk1Pt);
    dnt->Branch("Dtrk2Pt",&Dtrk2Pt);
    dnt->Branch("Dtrk3Pt",&Dtrk3Pt);
    dnt->Branch("Dtrk4Pt",&Dtrk4Pt);
    dnt->Branch("Dtrk1Eta",&Dtrk1Eta);
    dnt->Branch("Dtrk2Eta",&Dtrk2Eta);
    dnt->Branch("Dtrk3Eta",&Dtrk3Eta);
    dnt->Branch("Dtrk4Eta",&Dtrk4Eta);
    dnt->Branch("Dtrk1Phi",&Dtrk1Phi);
    dnt->Branch("Dtrk2Phi",&Dtrk2Phi);
    dnt->Branch("Dtrk3Phi",&Dtrk3Phi);
    dnt->Branch("Dtrk4Phi",&Dtrk4Phi);
    dnt->Branch("Dtrk1PtErr",&Dtrk1PtErr);
    dnt->Branch("Dtrk2PtErr",&Dtrk2PtErr);
    dnt->Branch("Dtrk3PtErr",&Dtrk3PtErr);
    dnt->Branch("Dtrk4PtErr",&Dtrk4PtErr);
    dnt->Branch("Dtrk1EtaErr",&Dtrk1EtaErr);
    dnt->Branch("Dtrk2EtaErr",&Dtrk2EtaErr);
    dnt->Branch("Dtrk3EtaErr",&Dtrk3EtaErr);
    dnt->Branch("Dtrk4EtaErr",&Dtrk4EtaErr);
    dnt->Branch("Dtrk1PhiErr",&Dtrk1PhiErr);
    dnt->Branch("Dtrk2PhiErr",&Dtrk2PhiErr);
    dnt->Branch("Dtrk3PhiErr",&Dtrk3PhiErr);
    dnt->Branch("Dtrk4PhiErr",&Dtrk4PhiErr);
    dnt->Branch("Dtrk1Y",&Dtrk1Y);
    dnt->Branch("Dtrk2Y",&Dtrk2Y);
    dnt->Branch("Dtrk3Y",&Dtrk3Y);
    dnt->Branch("Dtrk4Y",&Dtrk4Y);
    dnt->Branch("Dtrk1Dxy",&Dtrk1Dxy);
    dnt->Branch("Dtrk2Dxy",&Dtrk2Dxy);
    dnt->Branch("Dtrk3Dxy",&Dtrk3Dxy);
    dnt->Branch("Dtrk4Dxy",&Dtrk4Dxy);
    dnt->Branch("Dtrk1D0Err",&Dtrk1D0Err);
    dnt->Branch("Dtrk2D0Err",&Dtrk2D0Err);
    dnt->Branch("Dtrk3D0Err",&Dtrk3D0Err);
    dnt->Branch("Dtrk4D0Err",&Dtrk4D0Err);
    dnt->Branch("Dtrk1PixelHit",&Dtrk1PixelHit);
    dnt->Branch("Dtrk2PixelHit",&Dtrk2PixelHit);
    dnt->Branch("Dtrk3PixelHit",&Dtrk3PixelHit);
    dnt->Branch("Dtrk4PixelHit",&Dtrk4PixelHit);
    dnt->Branch("Dtrk1StripHit",&Dtrk1StripHit);
    dnt->Branch("Dtrk2StripHit",&Dtrk2StripHit);
    dnt->Branch("Dtrk3StripHit",&Dtrk3StripHit);
    dnt->Branch("Dtrk4StripHit",&Dtrk4StripHit);
    dnt->Branch("Dtrk1nStripLayer",&Dtrk1nStripLayer);
    dnt->Branch("Dtrk2nStripLayer",&Dtrk2nStripLayer);
    dnt->Branch("Dtrk3nStripLayer",&Dtrk3nStripLayer);
    dnt->Branch("Dtrk4nStripLayer",&Dtrk4nStripLayer);
    dnt->Branch("Dtrk1nPixelLayer",&Dtrk1nPixelLayer);
    dnt->Branch("Dtrk2nPixelLayer",&Dtrk2nPixelLayer);
    dnt->Branch("Dtrk3nPixelLayer",&Dtrk3nPixelLayer);
    dnt->Branch("Dtrk4nPixelLayer",&Dtrk4nPixelLayer);
    dnt->Branch("Dtrk1Chi2ndf",&Dtrk1Chi2ndf);
    dnt->Branch("Dtrk2Chi2ndf",&Dtrk2Chi2ndf);
    dnt->Branch("Dtrk3Chi2ndf",&Dtrk3Chi2ndf);
    dnt->Branch("Dtrk4Chi2ndf",&Dtrk4Chi2ndf);
    dnt->Branch("Dtrk1MassHypo",&Dtrk1MassHypo);
    dnt->Branch("Dtrk2MassHypo",&Dtrk2MassHypo);
    dnt->Branch("Dtrk3MassHypo",&Dtrk3MassHypo);
    dnt->Branch("Dtrk4MassHypo",&Dtrk4MassHypo);
    dnt->Branch("Dtrk1MVAVal",&Dtrk1MVAVal);
    dnt->Branch("Dtrk2MVAVal",&Dtrk2MVAVal);
    dnt->Branch("Dtrk3MVAVal",&Dtrk3MVAVal);
    dnt->Branch("Dtrk4MVAVal",&Dtrk4MVAVal);
    dnt->Branch("Dtrk1Algo",&Dtrk1Algo);
    dnt->Branch("Dtrk2Algo",&Dtrk2Algo);
    dnt->Branch("Dtrk3Algo",&Dtrk3Algo);
    dnt->Branch("Dtrk4Algo",&Dtrk4Algo);
    dnt->Branch("Dtrk1highPurity",&Dtrk1highPurity);
    dnt->Branch("Dtrk2highPurity",&Dtrk2highPurity);
    dnt->Branch("Dtrk3highPurity",&Dtrk3highPurity);
    dnt->Branch("Dtrk4highPurity",&Dtrk4highPurity);
    dnt->Branch("Dtrk1Quality",&Dtrk1Quality);
    dnt->Branch("Dtrk2Quality",&Dtrk2Quality);
    dnt->Branch("Dtrk3Quality",&Dtrk3Quality);
    dnt->Branch("Dtrk4Quality",&Dtrk4Quality);
    //DInfo.tktkResInfo
    dnt->Branch("DtktkResmass",&DtktkResmass);
    dnt->Branch("DtktkRespt",&DtktkRespt);
    dnt->Branch("DtktkReseta",&DtktkReseta);
    dnt->Branch("DtktkResphi",&DtktkResphi);
    dnt->Branch("DRestrk1Pt",&DRestrk1Pt);
    dnt->Branch("DRestrk1Eta",&DRestrk1Eta);
    dnt->Branch("DRestrk1Phi",&DRestrk1Phi);
    dnt->Branch("DRestrk1Y",&DRestrk1Y);
    dnt->Branch("DRestrk1Dxy",&DRestrk1Dxy);
    dnt->Branch("DRestrk1D0Err",&DRestrk1D0Err);
    dnt->Branch("DRestrk2Pt",&DRestrk2Pt);
    dnt->Branch("DRestrk2Eta",&DRestrk2Eta);
    dnt->Branch("DRestrk2Phi",&DRestrk2Phi);
    dnt->Branch("DRestrk2Y",&DRestrk2Y);
    dnt->Branch("DRestrk2Dxy",&DRestrk2Dxy);
    dnt->Branch("DRestrk2D0Err",&DRestrk2D0Err);
    dnt->Branch("DRestrk3Pt",&DRestrk3Pt);
    dnt->Branch("DRestrk3Eta",&DRestrk3Eta);
    dnt->Branch("DRestrk3Phi",&DRestrk3Phi);
    dnt->Branch("DRestrk3Y",&DRestrk3Y);
    dnt->Branch("DRestrk3Dxy",&DRestrk3Dxy);
    dnt->Branch("DRestrk3D0Err",&DRestrk3D0Err);
    dnt->Branch("DRestrk4Pt",&DRestrk4Pt);
    dnt->Branch("DRestrk4Eta",&DRestrk4Eta);
    dnt->Branch("DRestrk4Phi",&DRestrk4Phi);
    dnt->Branch("DRestrk4Y",&DRestrk4Y);
    dnt->Branch("DRestrk4Dxy",&DRestrk4Dxy);
    dnt->Branch("DRestrk4D0Err",&DRestrk4D0Err);
    //DInfo.genInfo
    dnt->Branch("Dgen",&Dgen);
    dnt->Branch("DgenIndex",&DgenIndex);
    dnt->Branch("DgennDa",&DgennDa);
    dnt->Branch("Dgenpt",&Dgenpt);
    dnt->Branch("Dgeneta",&Dgeneta);
    dnt->Branch("Dgenphi",&Dgenphi);
    dnt->Branch("Dgeny",&Dgeny);
  }
  
  //GenInfo
  float   Gy;
  float   Geta;
  float   Gphi;
  float   Gpt;
  float   GpdgId;
  int     GisSignal;
  
  void buildGenBranch(TTree* nt)
  {
    nt->Branch("Gy",&Gy);
    nt->Branch("Geta",&Geta);
    nt->Branch("Gphi",&Gphi);
    nt->Branch("Gpt",&Gpt);
    nt->Branch("GpdgId",&GpdgId);
    nt->Branch("GisSignal",&GisSignal);
  }
  
  void makeDNtuple(int isDchannel[], bool REAL, bool skim, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, TrackInfoBranches *TrackInfo, DInfoBranches *DInfo, GenInfoBranches *GenInfo, TTree* ntD1, TTree* ntD2, TTree* ntD3, TTree* ntD4, TTree* ntD5, TTree* ntD6)
  {//{{{
    TVector3* bP = new TVector3;
    TVector3* bVtx = new TVector3;
    TLorentzVector* b4P = new TLorentzVector;
    for(int t=0;t<12;t++)
      {
        if(isDchannel[t]==1)
          {
            for(int j=0;j<DInfo->size;j++)
              {
                if(skim)
                  {
                    if(DInfo->alpha[j]>0.13) continue;
                    if((DInfo->pt[j]>=13.&&(DInfo->svpvDistance[j]/DInfo->svpvDisErr[j])<2.5)||
                       (DInfo->pt[j]>=5.5&&DInfo->pt[j]<13.&&(DInfo->svpvDistance[j]/DInfo->svpvDisErr[j])<4.)||
                       (DInfo->pt[j]<5.5&&(DInfo->svpvDistance[j]/DInfo->svpvDisErr[j])<5.)) continue;
                  }
                if(DInfo->type[j]==(t+1))
                  {
                    fillTreeEvt(EvtInfo);
                    fillDTree(bP,bVtx,b4P,j,REAL,EvtInfo,VtxInfo,TrackInfo,DInfo,GenInfo);
                    if(t==1||t==0)       ntD1->Fill();
                    else if(t==3||t==2)  ntD2->Fill();
                    else if(t==5||t==4)  ntD3->Fill();
                    else if(t==7||t==6)  ntD4->Fill();
                    else if(t==9||t==8)  ntD5->Fill();
                    else if(t==11||t==10) ntD6->Fill();
                  }
              }
          }
      }
  }//}}}
  
  void fillDGenTree(TTree* ntGen, GenInfoBranches *GenInfo, bool gskim=true)
  {
    TLorentzVector* bGen = new TLorentzVector;
    int gt=0,sigtype=0;
    for(int j=0;j<GenInfo->size;j++)
      {
        if(TMath::Abs(GenInfo->pdgId[j])!=DZERO_PDGID&&
           TMath::Abs(GenInfo->pdgId[j])!=DPLUS_PDGID&&
           TMath::Abs(GenInfo->pdgId[j])!=DSUBS_PDGID&&
           TMath::Abs(GenInfo->pdgId[j])!=DSTAR_PDGID&&gskim) continue;
        Gpt = GenInfo->pt[j];
        Geta = GenInfo->eta[j];
        Gphi = GenInfo->phi[j];
        GpdgId = GenInfo->pdgId[j];
        bGen->SetPtEtaPhiM(GenInfo->pt[j],GenInfo->eta[j],GenInfo->phi[j],GenInfo->mass[j]);
        Gy = bGen->Rapidity();
        sigtype=0;
        for(gt=1;gt<13;gt++)
          {
            if(isDsignalGen(gt,j, GenInfo))
              {
                sigtype=gt;
                break;
              }
          }
        GisSignal = sigtype;
        ntGen->Fill();
      }
  }

  double findMass(int particlePdgId)
  {
    if(TMath::Abs(particlePdgId)==211) return PION_MASS;
    if(TMath::Abs(particlePdgId)==321) return KAON_MASS;
    else
      {
        cout<<"ERROR: find particle mass falied >> Particle pdgId: "<<particlePdgId<<endl;
        return 0;
      }
  }

  int findPdgid(Double_t tkmass)
  {
    if(TMath::Abs(tkmass-KAON_MASS)<0.1) return KAON_PDGID;
    else if(TMath::Abs(tkmass-PION_MASS)<0.1) return PION_PDGID;
    else
      {
        cout<<"ERROR: find particle pdgid falied >> Particle mass: "<<tkmass<<endl;
        return 0;
      }
  }
  
  void fillTreeEvt(EvtInfoBranches *EvtInfo)
  {
    //Event Info
    RunNo = EvtInfo->RunNo;
    EvtNo = EvtInfo->EvtNo;
    LumiNo = EvtInfo->LumiNo;
    PVx = EvtInfo->PVx;
    PVy = EvtInfo->PVy;
    PVz = EvtInfo->PVz;
    PVxE = EvtInfo->PVxE;
    PVyE = EvtInfo->PVyE;
    PVzE = EvtInfo->PVzE;
    PVnchi2 = EvtInfo->PVnchi2;
    PVchi2 = EvtInfo->PVchi2;
    BSx = EvtInfo->BSx;
    BSy = EvtInfo->BSy;
    BSz = EvtInfo->BSz;
    BSxErr = EvtInfo->BSxErr;
    BSyErr = EvtInfo->BSyErr;
    BSzErr = EvtInfo->BSzErr;
    BSdxdz = EvtInfo->BSdxdz;
    BSdydz = EvtInfo->BSdydz;
    BSdxdzErr = EvtInfo->BSdxdzErr;
    BSdydzErr = EvtInfo->BSdydzErr;
    BSWidthX = EvtInfo->BSWidthX;
    BSWidthXErr = EvtInfo->BSWidthXErr;
    BSWidthY = EvtInfo->BSWidthY;
    BSWidthYErr = EvtInfo->BSWidthYErr;
  }

  void fillDTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, int j, bool REAL, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, TrackInfoBranches *TrackInfo, DInfoBranches *DInfo, GenInfoBranches *GenInfo)
  {
    //DInfo
    bP->SetPtEtaPhi(DInfo->pt[j],DInfo->eta[j]*0,DInfo->phi[j]);
    bVtx->SetXYZ(DInfo->vtxX[j]-EvtInfo->PVx,
                 DInfo->vtxY[j]-EvtInfo->PVy,
                 DInfo->vtxZ[j]*0-EvtInfo->PVz*0);
    b4P->SetPtEtaPhiM(DInfo->pt[j],DInfo->eta[j],DInfo->phi[j],DInfo->mass[j]);
    Dtype = DInfo->type[j];
    Dmass = DInfo->mass[j];
    Dpt = DInfo->pt[j];
    Deta = DInfo->eta[j];
    Dphi = DInfo->phi[j];
    Dy = b4P->Rapidity();
    DvtxX = DInfo->vtxX[j] - EvtInfo->PVx;
    DvtxY = DInfo->vtxY[j] - EvtInfo->PVy;
    Dd0 = TMath::Sqrt((DInfo->vtxX[j]-EvtInfo->PVx)*(DInfo->vtxX[j]-EvtInfo->PVx)+(DInfo->vtxY[j]-EvtInfo->PVy)*(DInfo->vtxY[j]-EvtInfo->PVy));
    Dd0Err = TMath::Sqrt(DInfo->vtxXErr[j]*DInfo->vtxXErr[j]+DInfo->vtxYErr[j]*DInfo->vtxYErr[j]);
    Ddxyz = TMath::Sqrt((DInfo->vtxX[j]-EvtInfo->PVx)*(DInfo->vtxX[j]-EvtInfo->PVx)+(DInfo->vtxY[j]-EvtInfo->PVy)*(DInfo->vtxY[j]-EvtInfo->PVy)+(DInfo->vtxZ[j]-EvtInfo->PVz)*(DInfo->vtxZ[j]-EvtInfo->PVz));
    DdxyzErr = TMath::Sqrt(DInfo->vtxXErr[j]*DInfo->vtxXErr[j]+DInfo->vtxYErr[j]*DInfo->vtxYErr[j]+DInfo->vtxZErr[j]*DInfo->vtxZErr[j]);
    Dchi2ndf = DInfo->vtxchi2[j]/DInfo->vtxdof[j];
    Dchi2cl = TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]);
    Ddtheta = bP->Angle(*bVtx);
    Dlxy = ((DInfo->vtxX[j]-EvtInfo->PVx)*b4P->Px() + (DInfo->vtxY[j]-EvtInfo->PVy)*b4P->Py())/DInfo->pt[j];
    Dalpha = DInfo->alpha[j];
    DsvpvDistance = DInfo->svpvDistance[j];
    DsvpvDisErr = DInfo->svpvDisErr[j];
    DsvpvDistance_2D = DInfo->svpvDistance_2D[j];
    DsvpvDisErr_2D = DInfo->svpvDisErr_2D[j];
    DtktkRes_chi2ndf = DInfo->tktkRes_vtxchi2[j]/DInfo->tktkRes_vtxdof[j];
    DtktkRes_chi2cl = TMath::Prob(DInfo->tktkRes_vtxchi2[j], DInfo->tktkRes_vtxdof[j]);
    DtktkRes_alpha = DInfo->tktkRes_alpha[j];
    DtktkRes_svpvDistance = DInfo->tktkRes_svpvDistance[j];
    DtktkRes_svpvDisErr = DInfo->tktkRes_svpvDisErr[j];
    float r2lxyBS = (DInfo->vtxX[j]-EvtInfo->BSx+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz) * (DInfo->vtxX[j]-EvtInfo->BSx+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz)
      + (DInfo->vtxY[j]-EvtInfo->BSy+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz) * (DInfo->vtxY[j]-EvtInfo->BSy+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz);
    float xlxyBS = DInfo->vtxX[j]-EvtInfo->BSx + (DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz;
    float ylxyBS = DInfo->vtxY[j]-EvtInfo->BSy + (DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz;
    DlxyBS = TMath::Sqrt(r2lxyBS);
    DlxyBSErr = (1./r2lxyBS) * ((xlxyBS*xlxyBS)*DInfo->vtxXErr[j] + (2*xlxyBS*ylxyBS)*DInfo->vtxYXErr[j] + (ylxyBS*ylxyBS)*DInfo->vtxYErr[j]);
    DMaxDoca = DInfo->MaxDoca[j];
    Dmaxpt = false;
    Dmaxprob = false;
    DmaxptMatched = false;
    DmaxprobMatched = false;
    
    //DInfo.trkInfo
    float trk1mass,trk2mass,trk3mass,trk4mass;
    if(DInfo->type[j]==1||DInfo->type[j]==2||DInfo->type[j]==3||DInfo->type[j]==4||DInfo->type[j]==5||DInfo->type[j]==6)
      {
        Dtrk1Idx = DInfo->rftk1_index[j];
        Dtrk1Pt = TrackInfo->pt[DInfo->rftk1_index[j]];
        Dtrk1Eta = TrackInfo->eta[DInfo->rftk1_index[j]];
        Dtrk1Phi = TrackInfo->phi[DInfo->rftk1_index[j]];
        Dtrk1PtErr = TrackInfo->ptErr[DInfo->rftk1_index[j]];
        Dtrk1EtaErr = TrackInfo->etaErr[DInfo->rftk1_index[j]];
        Dtrk1PhiErr = TrackInfo->phiErr[DInfo->rftk1_index[j]];
        trk1mass = findMass(DInfo->rftk1_MassHypo[j]);
        b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk1_index[j]],TrackInfo->eta[DInfo->rftk1_index[j]],TrackInfo->phi[DInfo->rftk1_index[j]],trk1mass);
        Dtrk1Y = b4P->Rapidity();
        Dtrk1Dxy = TrackInfo->dxyPV[DInfo->rftk1_index[j]];
        Dtrk1D0Err = TrackInfo->d0error[DInfo->rftk1_index[j]];
        Dtrk1PixelHit = TrackInfo->pixelhit[DInfo->rftk1_index[j]];
        Dtrk1StripHit = TrackInfo->striphit[DInfo->rftk1_index[j]];
        Dtrk1nPixelLayer = TrackInfo->nPixelLayer[DInfo->rftk1_index[j]];
        Dtrk1nStripLayer = TrackInfo->nStripLayer[DInfo->rftk1_index[j]];
        Dtrk1MassHypo = DInfo->rftk1_MassHypo[j]*TrackInfo->charge[DInfo->rftk1_index[j]];
        Dtrk1Chi2ndf = TrackInfo->chi2[DInfo->rftk1_index[j]]/TrackInfo->ndf[DInfo->rftk1_index[j]];
        Dtrk1MVAVal = TrackInfo->trkMVAVal[DInfo->rftk1_index[j]];
        Dtrk1Algo = TrackInfo->trkAlgo[DInfo->rftk1_index[j]];
        Dtrk1highPurity = TrackInfo->highPurity[DInfo->rftk1_index[j]];
        Dtrk1Quality = TrackInfo->trackQuality[DInfo->rftk1_index[j]];

        Dtrk2Idx = DInfo->rftk2_index[j];
        Dtrk2Pt = TrackInfo->pt[DInfo->rftk2_index[j]];
        Dtrk2Eta = TrackInfo->eta[DInfo->rftk2_index[j]];
        Dtrk2Phi = TrackInfo->phi[DInfo->rftk2_index[j]];
        Dtrk2PtErr = TrackInfo->ptErr[DInfo->rftk2_index[j]];
        Dtrk2EtaErr = TrackInfo->etaErr[DInfo->rftk2_index[j]];
        Dtrk2PhiErr = TrackInfo->phiErr[DInfo->rftk2_index[j]];
        trk2mass = findMass(DInfo->rftk2_MassHypo[j]);
        b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk2_index[j]],TrackInfo->eta[DInfo->rftk2_index[j]],TrackInfo->phi[DInfo->rftk2_index[j]],trk2mass);
        Dtrk2Y = b4P->Rapidity();
        Dtrk2Dxy = TrackInfo->dxyPV[DInfo->rftk2_index[j]];
        Dtrk2D0Err = TrackInfo->d0error[DInfo->rftk2_index[j]];
        Dtrk2PixelHit = TrackInfo->pixelhit[DInfo->rftk2_index[j]];
        Dtrk2StripHit = TrackInfo->striphit[DInfo->rftk2_index[j]];
        Dtrk2nPixelLayer = TrackInfo->nPixelLayer[DInfo->rftk2_index[j]];
        Dtrk2nStripLayer = TrackInfo->nStripLayer[DInfo->rftk2_index[j]];
        Dtrk2Chi2ndf = TrackInfo->chi2[DInfo->rftk2_index[j]]/TrackInfo->ndf[DInfo->rftk2_index[j]];
        Dtrk2MassHypo = DInfo->rftk2_MassHypo[j]*TrackInfo->charge[DInfo->rftk2_index[j]];
        Dtrk2MVAVal = TrackInfo->trkMVAVal[DInfo->rftk2_index[j]];
        Dtrk2Algo = TrackInfo->trkAlgo[DInfo->rftk2_index[j]];
        Dtrk2highPurity = TrackInfo->highPurity[DInfo->rftk2_index[j]];
        Dtrk2Quality = TrackInfo->trackQuality[DInfo->rftk2_index[j]];

        if(DInfo->type[j]==1||DInfo->type[j]==2)
          {
            Dtrk3Idx = -1;
            Dtrk3Pt = -1;
            Dtrk3Eta = -20;
            Dtrk3Phi = -20;
            Dtrk3PtErr = 0;
            Dtrk3EtaErr = 0;
            Dtrk3PhiErr = 0;
            Dtrk3Y = -1;
            Dtrk3Dxy = -1;
            Dtrk3D0Err = -1;
            Dtrk3PixelHit = -1;
            Dtrk3StripHit = -1;
            Dtrk3nPixelLayer = -1;
            Dtrk3nStripLayer = -1;
            Dtrk3Chi2ndf = -1;
            Dtrk3MassHypo = 0;
            Dtrk3MVAVal = -100;
            Dtrk3Algo = 0;
            Dtrk3Quality = 0;
            Dtrk3highPurity = false;
            Dtrk4Idx = -1;
            Dtrk4Pt = -1;
            Dtrk4Eta = -20;
            Dtrk4Phi = -20;
            Dtrk4PtErr = 0;
            Dtrk4EtaErr = 0;
            Dtrk4PhiErr = 0;
            Dtrk4Y = -1;
            Dtrk4Dxy = -1;
            Dtrk4D0Err = -1;
            Dtrk4PixelHit = -1;
            Dtrk4StripHit = -1;
            Dtrk4nPixelLayer = -1;
            Dtrk4nStripLayer = -1;
            Dtrk4Chi2ndf = -1;
            Dtrk4MassHypo = 0;
            Dtrk4MVAVal = -100;
            Dtrk4Algo = 0;
            Dtrk4Quality = 0;
            Dtrk4highPurity = false;
          
            DtktkResmass = -1;
            DtktkRespt = -1;
            DtktkReseta = -20;
            DtktkResphi = -20;
          
            DRestrk1Pt = -1;
            DRestrk1Eta = -20;
            DRestrk1Phi = -20;
            DRestrk1Y = -1;
            DRestrk1Dxy = -1;
            DRestrk1D0Err = -1;
            DRestrk2Pt = -1;
            DRestrk2Eta = -20;
            DRestrk2Phi = -20;
            DRestrk2Y = -1;
            DRestrk2Dxy = -1;
            DRestrk2D0Err = -1;
            DRestrk3Pt = -1;
            DRestrk3Eta = -20;
            DRestrk3Phi = -20;
            DRestrk3Y = -1;
            DRestrk3Dxy = -1;
            DRestrk3D0Err = -1;
            DRestrk4Pt = -1;
            DRestrk4Eta = -20;
            DRestrk4Phi = -20;
            DRestrk4Y = -1;
            DRestrk4Dxy = -1;
            DRestrk4D0Err = -1;
          }
        else if(DInfo->type[j]==3||DInfo->type[j]==4)
          {
            Dtrk3Idx = DInfo->rftk3_index[j];
            Dtrk3Pt = TrackInfo->pt[DInfo->rftk3_index[j]];
            Dtrk3Eta = TrackInfo->eta[DInfo->rftk3_index[j]];
            Dtrk3Phi = TrackInfo->phi[DInfo->rftk3_index[j]];
            Dtrk3PtErr = TrackInfo->ptErr[DInfo->rftk3_index[j]];
            Dtrk3EtaErr = TrackInfo->etaErr[DInfo->rftk3_index[j]];
            Dtrk3PhiErr = TrackInfo->phiErr[DInfo->rftk3_index[j]];
            trk3mass = findMass(DInfo->rftk3_MassHypo[j]);
            b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk3_index[j]],TrackInfo->eta[DInfo->rftk3_index[j]],TrackInfo->phi[DInfo->rftk3_index[j]],trk3mass);
            Dtrk3Y = b4P->Rapidity();
            Dtrk3Dxy = TrackInfo->dxyPV[DInfo->rftk3_index[j]];
            Dtrk3D0Err = TrackInfo->d0error[DInfo->rftk3_index[j]];
            Dtrk3PixelHit = TrackInfo->pixelhit[DInfo->rftk3_index[j]];
            Dtrk3StripHit = TrackInfo->striphit[DInfo->rftk3_index[j]];
            Dtrk3nPixelLayer = TrackInfo->nPixelLayer[DInfo->rftk3_index[j]];
            Dtrk3nStripLayer = TrackInfo->nStripLayer[DInfo->rftk3_index[j]];
            Dtrk3Chi2ndf = TrackInfo->chi2[DInfo->rftk3_index[j]]/TrackInfo->ndf[DInfo->rftk3_index[j]];
            Dtrk3MassHypo = DInfo->rftk3_MassHypo[j]*TrackInfo->charge[DInfo->rftk3_index[j]];
            Dtrk3MVAVal = TrackInfo->trkMVAVal[DInfo->rftk3_index[j]];
            Dtrk3Algo = TrackInfo->trkAlgo[DInfo->rftk3_index[j]];
            Dtrk3highPurity = TrackInfo->highPurity[DInfo->rftk3_index[j]];
            Dtrk3Quality = TrackInfo->trackQuality[DInfo->rftk3_index[j]];
            Dtrk4Idx = -1;
            Dtrk4Pt = -1;
            Dtrk4Eta = -20;
            Dtrk4Phi = -20;
            Dtrk4PtErr = 0;
            Dtrk4EtaErr = 0;
            Dtrk4PhiErr = 0;
            Dtrk4Y = -1;
            Dtrk4Dxy = -1;
            Dtrk4D0Err = -1;
            Dtrk4PixelHit = -1;
            Dtrk4StripHit = -1;
            Dtrk4nPixelLayer = -1;
            Dtrk4nStripLayer = -1;
            Dtrk4Chi2ndf = -1;
            Dtrk4MassHypo = 0;
            Dtrk4MVAVal = -100;
            Dtrk4Algo = 0;
            Dtrk4Quality = 0;
            Dtrk4highPurity = false;

            DtktkResmass = -1;
            DtktkRespt = -1;
            DtktkReseta = -20;
            DtktkResphi = -20;

            DRestrk1Pt = -1;
            DRestrk1Eta = -20;
            DRestrk1Phi = -20;
            DRestrk1Y = -1;
            DRestrk1Dxy = -1;
            DRestrk1D0Err = -1;
            DRestrk2Pt = -1;
            DRestrk2Eta = -20;
            DRestrk2Phi = -20;
            DRestrk2Y = -1;
            DRestrk2Dxy = -1;
            DRestrk2D0Err = -1;
            DRestrk3Pt = -1;
            DRestrk3Eta = -20;
            DRestrk3Phi = -20;
            DRestrk3Y = -1;
            DRestrk3Dxy = -1;
            DRestrk3D0Err = -1;
            DRestrk4Pt = -1;
            DRestrk4Eta = -20;
            DRestrk4Phi = -20;
            DRestrk4Y = -1;
            DRestrk4Dxy = -1;
            DRestrk4D0Err = -1;
          }
        else if(DInfo->type[j]==5||DInfo->type[j]==6)
          {
            Dtrk3Idx = DInfo->rftk3_index[j];
            Dtrk4Idx = DInfo->rftk4_index[j];
            Dtrk3Pt = TrackInfo->pt[DInfo->rftk3_index[j]];
            Dtrk4Pt = TrackInfo->pt[DInfo->rftk4_index[j]];
            Dtrk3Eta = TrackInfo->eta[DInfo->rftk3_index[j]];
            Dtrk4Eta = TrackInfo->eta[DInfo->rftk4_index[j]];
            Dtrk3Phi = TrackInfo->phi[DInfo->rftk3_index[j]];
            Dtrk4Phi = TrackInfo->phi[DInfo->rftk4_index[j]];
            Dtrk3PtErr = TrackInfo->ptErr[DInfo->rftk3_index[j]];
            Dtrk4PtErr = TrackInfo->ptErr[DInfo->rftk4_index[j]];
            Dtrk3EtaErr = TrackInfo->etaErr[DInfo->rftk3_index[j]];
            Dtrk4EtaErr = TrackInfo->etaErr[DInfo->rftk4_index[j]];
            Dtrk3PhiErr = TrackInfo->phiErr[DInfo->rftk3_index[j]];
            Dtrk4PhiErr = TrackInfo->phiErr[DInfo->rftk4_index[j]];
            trk3mass = findMass(DInfo->rftk3_MassHypo[j]);
            trk4mass = findMass(DInfo->rftk4_MassHypo[j]);
            b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk3_index[j]],TrackInfo->eta[DInfo->rftk3_index[j]],TrackInfo->phi[DInfo->rftk3_index[j]],trk3mass);
            Dtrk3Y = b4P->Rapidity();
            b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk4_index[j]],TrackInfo->eta[DInfo->rftk4_index[j]],TrackInfo->phi[DInfo->rftk4_index[j]],trk4mass);
            Dtrk4Y = b4P->Rapidity();
            Dtrk3Dxy = TrackInfo->dxyPV[DInfo->rftk3_index[j]];
            Dtrk4Dxy = TrackInfo->dxyPV[DInfo->rftk4_index[j]];
            Dtrk3D0Err = TrackInfo->d0error[DInfo->rftk3_index[j]];
            Dtrk4D0Err = TrackInfo->d0error[DInfo->rftk4_index[j]];
            Dtrk3PixelHit = TrackInfo->pixelhit[DInfo->rftk3_index[j]];
            Dtrk4PixelHit = TrackInfo->pixelhit[DInfo->rftk4_index[j]];
            Dtrk3StripHit = TrackInfo->striphit[DInfo->rftk3_index[j]];
            Dtrk4StripHit = TrackInfo->striphit[DInfo->rftk4_index[j]];
            Dtrk3nPixelLayer = TrackInfo->nPixelLayer[DInfo->rftk3_index[j]];
            Dtrk4nPixelLayer = TrackInfo->nPixelLayer[DInfo->rftk4_index[j]];
            Dtrk3nStripLayer = TrackInfo->nStripLayer[DInfo->rftk3_index[j]];
            Dtrk4nStripLayer = TrackInfo->nStripLayer[DInfo->rftk4_index[j]];
            Dtrk3Chi2ndf = TrackInfo->chi2[DInfo->rftk3_index[j]]/TrackInfo->ndf[DInfo->rftk3_index[j]];
            Dtrk4Chi2ndf = TrackInfo->chi2[DInfo->rftk4_index[j]]/TrackInfo->ndf[DInfo->rftk4_index[j]];
            Dtrk3MassHypo = DInfo->rftk3_MassHypo[j]*TrackInfo->charge[DInfo->rftk3_index[j]];
            Dtrk4MassHypo = DInfo->rftk4_MassHypo[j]*TrackInfo->charge[DInfo->rftk4_index[j]];
            Dtrk3MVAVal = TrackInfo->trkMVAVal[DInfo->rftk3_index[j]];
            Dtrk4MVAVal = TrackInfo->trkMVAVal[DInfo->rftk4_index[j]];
            Dtrk3Algo = TrackInfo->trkAlgo[DInfo->rftk3_index[j]];
            Dtrk4Algo = TrackInfo->trkAlgo[DInfo->rftk4_index[j]];
            Dtrk3highPurity = TrackInfo->highPurity[DInfo->rftk3_index[j]];
            Dtrk4highPurity = TrackInfo->highPurity[DInfo->rftk4_index[j]];
            Dtrk3Quality = TrackInfo->trackQuality[DInfo->rftk3_index[j]];
            Dtrk4Quality = TrackInfo->trackQuality[DInfo->rftk4_index[j]];

            DtktkResmass = -1;
            DtktkRespt = -1;
            DtktkReseta = -20;
            DtktkResphi = -20;

            DRestrk1Pt = -1;
            DRestrk1Eta = -20;
            DRestrk1Phi = -20;
            DRestrk1Y = -1;
            DRestrk1Dxy = -1;
            DRestrk1D0Err = -1;
            DRestrk2Pt = -1;
            DRestrk2Eta = -20;
            DRestrk2Phi = -20;
            DRestrk2Y = -1;
            DRestrk2Dxy = -1;
            DRestrk2D0Err = -1;
            DRestrk3Pt = -1;
            DRestrk3Eta = -20;
            DRestrk3Phi = -20;
            DRestrk3Y = -1;
            DRestrk3Dxy = -1;
            DRestrk3D0Err = -1;
            DRestrk4Pt = -1;
            DRestrk4Eta = -20;
            DRestrk4Phi = -20;
            DRestrk4Y = -1;
            DRestrk4Dxy = -1;
            DRestrk4D0Err = -1;
          }
      }
    else if(DInfo->type[j]==7||DInfo->type[j]==8||DInfo->type[j]==9||DInfo->type[j]==10||DInfo->type[j]==11||DInfo->type[j]==12)
      {
        Dtrk1Idx = DInfo->rftk2_index[j];
        Dtrk1Pt = TrackInfo->pt[DInfo->rftk2_index[j]];
        Dtrk1Eta = TrackInfo->eta[DInfo->rftk2_index[j]];
        Dtrk1Phi = TrackInfo->phi[DInfo->rftk2_index[j]];
        Dtrk1PtErr = TrackInfo->ptErr[DInfo->rftk2_index[j]];
        Dtrk1EtaErr = TrackInfo->etaErr[DInfo->rftk2_index[j]];
        Dtrk1PhiErr = TrackInfo->phiErr[DInfo->rftk2_index[j]];
        b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk2_index[j]],TrackInfo->eta[DInfo->rftk2_index[j]],TrackInfo->phi[DInfo->rftk2_index[j]],PION_MASS);
        Dtrk1Y = b4P->Rapidity();
        Dtrk1Dxy = TrackInfo->dxyPV[DInfo->rftk2_index[j]];
        Dtrk1D0Err = TrackInfo->d0error[DInfo->rftk2_index[j]];
        Dtrk1PixelHit = TrackInfo->pixelhit[DInfo->rftk2_index[j]];
        Dtrk1StripHit = TrackInfo->striphit[DInfo->rftk2_index[j]];
        Dtrk1nPixelLayer = TrackInfo->nPixelLayer[DInfo->rftk2_index[j]];
        Dtrk1nStripLayer = TrackInfo->nStripLayer[DInfo->rftk2_index[j]];
        Dtrk1Chi2ndf = TrackInfo->chi2[DInfo->rftk2_index[j]]/TrackInfo->ndf[DInfo->rftk2_index[j]];
        Dtrk1MassHypo = DInfo->rftk2_MassHypo[j]*TrackInfo->charge[DInfo->rftk2_index[j]];
        Dtrk1MVAVal = TrackInfo->trkMVAVal[DInfo->rftk2_index[j]];
        Dtrk1Algo = TrackInfo->trkAlgo[DInfo->rftk2_index[j]];
        Dtrk1highPurity = TrackInfo->highPurity[DInfo->rftk2_index[j]];
        Dtrk1Quality = TrackInfo->trackQuality[DInfo->rftk2_index[j]];

        Dtrk2Idx = -1;
        Dtrk2Pt = -1;
        Dtrk2Eta = -20;
        Dtrk2Phi = -20;
        Dtrk2PtErr = 0;
        Dtrk2EtaErr = 0;
        Dtrk2PhiErr = 0;
        Dtrk2Y = -1;
        Dtrk2Dxy = -1;
        Dtrk2D0Err = -1;
        Dtrk2PixelHit = -1;
        Dtrk2StripHit = -1;
        Dtrk2nPixelLayer = -1;
        Dtrk2nStripLayer = -1;
        Dtrk2Chi2ndf = -1;
        Dtrk2MassHypo = 0;
        Dtrk2MVAVal = -100;
        Dtrk2Algo = 0;
        Dtrk2Quality = 0;
        Dtrk2highPurity = false;
        Dtrk3Idx = -1;
        Dtrk3Pt = -1;
        Dtrk3Eta = -20;
        Dtrk3Phi = -20;
        Dtrk3PtErr = 0;
        Dtrk3EtaErr = 0;
        Dtrk3PhiErr = 0;
        Dtrk3Y = -1;
        Dtrk3Dxy = -1;
        Dtrk3D0Err = -1;
        Dtrk3PixelHit = -1;
        Dtrk3StripHit = -1;
        Dtrk3nPixelLayer = -1;
        Dtrk3nStripLayer = -1;
        Dtrk3Chi2ndf = -1;
        Dtrk3MassHypo = 0;
        Dtrk3MVAVal = -100;
        Dtrk3Algo = 0;
        Dtrk3Quality = 0;
        Dtrk3highPurity = false;
        Dtrk4Idx = -1;
        Dtrk4Pt = -1;
        Dtrk4Eta = -20;
        Dtrk4Phi = -20;
        Dtrk4PtErr = 0;
        Dtrk4EtaErr = 0;
        Dtrk4PhiErr = 0;
        Dtrk4Y = -1;
        Dtrk4Dxy = -1;
        Dtrk4D0Err = -1;
        Dtrk4PixelHit = -1;
        Dtrk4StripHit = -1;
        Dtrk4nPixelLayer = -1;
        Dtrk4nStripLayer = -1;
        Dtrk4Chi2ndf = -1;
        Dtrk4MassHypo = 0;
        Dtrk4MVAVal = -100;
        Dtrk4Algo = 0;
        Dtrk4Quality = 0;
        Dtrk4highPurity = false;

        DtktkResmass = DInfo->tktkRes_mass[j];
        DtktkRespt = DInfo->tktkRes_pt[j];
        DtktkReseta = DInfo->tktkRes_eta[j];
        DtktkResphi = DInfo->tktkRes_phi[j];

        DRestrk1Pt = TrackInfo->pt[DInfo->tktkRes_rftk1_index[j]];
        DRestrk1Eta = TrackInfo->eta[DInfo->tktkRes_rftk1_index[j]];
        DRestrk1Phi = TrackInfo->phi[DInfo->tktkRes_rftk1_index[j]];
        b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->tktkRes_rftk1_index[j]],TrackInfo->eta[DInfo->tktkRes_rftk1_index[j]],TrackInfo->phi[DInfo->tktkRes_rftk1_index[j]],DInfo->tktkRes_rftk1_mass[j]);
        DRestrk1Y = b4P->Rapidity();
        DRestrk1Dxy = TrackInfo->dxyPV[DInfo->tktkRes_rftk1_index[j]];
        DRestrk1D0Err = TrackInfo->d0error[DInfo->tktkRes_rftk1_index[j]];
        DRestrk2Pt = TrackInfo->pt[DInfo->tktkRes_rftk2_index[j]];
        DRestrk2Eta = TrackInfo->eta[DInfo->tktkRes_rftk2_index[j]];
        DRestrk2Phi = TrackInfo->phi[DInfo->tktkRes_rftk2_index[j]];
        b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->tktkRes_rftk2_index[j]],TrackInfo->eta[DInfo->tktkRes_rftk2_index[j]],TrackInfo->phi[DInfo->tktkRes_rftk2_index[j]],DInfo->tktkRes_rftk2_mass[j]);
        DRestrk2Y = b4P->Rapidity();
        DRestrk2Dxy = TrackInfo->dxyPV[DInfo->tktkRes_rftk2_index[j]];
        DRestrk2D0Err = TrackInfo->d0error[DInfo->tktkRes_rftk2_index[j]];
        DRestrk3Pt = -1;
        DRestrk3Eta = -20;
        DRestrk3Phi = -20;
        DRestrk3Y = -1;
        DRestrk3Dxy = -1;
        DRestrk3D0Err = -1;
        DRestrk4Pt = -1;
        DRestrk4Eta = -20;
        DRestrk4Phi = -20;
        DRestrk4Y = -1;
        DRestrk4Dxy = -1;
        DRestrk4D0Err = -1;
        if(DInfo->type[j]==11||DInfo->type[j]==12)
          {
            DRestrk3Pt = TrackInfo->pt[DInfo->tktkRes_rftk3_index[j]];
            DRestrk3Eta = TrackInfo->eta[DInfo->tktkRes_rftk3_index[j]];
            DRestrk3Phi = TrackInfo->phi[DInfo->tktkRes_rftk3_index[j]];
            b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->tktkRes_rftk3_index[j]],TrackInfo->eta[DInfo->tktkRes_rftk3_index[j]],TrackInfo->phi[DInfo->tktkRes_rftk3_index[j]],DInfo->tktkRes_rftk3_mass[j]);
            DRestrk3Y = b4P->Rapidity();
            DRestrk3Dxy = TrackInfo->dxyPV[DInfo->tktkRes_rftk3_index[j]];
            DRestrk3D0Err = TrackInfo->d0error[DInfo->tktkRes_rftk3_index[j]];
            DRestrk4Pt = TrackInfo->pt[DInfo->tktkRes_rftk4_index[j]];
            DRestrk4Eta = TrackInfo->eta[DInfo->tktkRes_rftk4_index[j]];
            DRestrk4Phi = TrackInfo->phi[DInfo->tktkRes_rftk4_index[j]];
            b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->tktkRes_rftk4_index[j]],TrackInfo->eta[DInfo->tktkRes_rftk4_index[j]],TrackInfo->phi[DInfo->tktkRes_rftk4_index[j]],DInfo->tktkRes_rftk4_mass[j]);
            DRestrk4Y = b4P->Rapidity();
            DRestrk4Dxy = TrackInfo->dxyPV[DInfo->tktkRes_rftk4_index[j]];
            DRestrk4D0Err = TrackInfo->d0error[DInfo->tktkRes_rftk4_index[j]];
          }
      }
    
    int DpdgId=0,RpdgId=0;
    if(DInfo->type[j]==1||DInfo->type[j]==2||DInfo->type[j]==5||DInfo->type[j]==6) DpdgId=DZERO_PDGID;
    else if(DInfo->type[j]==3||DInfo->type[j]==4) DpdgId=DPLUS_PDGID;
    else if(DInfo->type[j]==7||DInfo->type[j]==8) DpdgId=DSUBS_PDGID;
    else if(DInfo->type[j]==9||DInfo->type[j]==10||DInfo->type[j]==11||DInfo->type[j]==12) DpdgId=DSTAR_PDGID;
    if(DInfo->type[j]==7||DInfo->type[j]==8) RpdgId=PHI_PDGID;
    else if(DInfo->type[j]==9||DInfo->type[j]==10||DInfo->type[j]==11||DInfo->type[j]==12) RpdgId=DZERO_PDGID;
    Dgen = 0;//gen init
    DgenIndex = -1;
    DgennDa = -1;
    Dgenpt = -1;
    Dgeneta = -20;
    Dgenphi = -20;
    Dgeny = -1;
    if(!REAL)
      {
        if(DInfo->type[j]==1||DInfo->type[j]==2)
          {
            if(DInfo->rftk1_index[j]>-1 && DInfo->rftk2_index[j]>-1)
              {
                if(TrackInfo->geninfo_index[DInfo->rftk1_index[j]]>-1 && 
                   TrackInfo->geninfo_index[DInfo->rftk2_index[j]]>-1)
                  {
                    if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]>-1 && 
                       GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1)
                      {
                        if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]])==DpdgId && 
                           GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])
                          {
                            if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j] && 
                               TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j])
                              {
                                Dgen = 23333;
                              }
                            else if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk1_MassHypo[j] && 
                                    TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk2_MassHypo[j])
                              {
                                Dgen = 23344;
                              }
                          }
                      }
                  }
              }
          }
        else if(DInfo->type[j]==3||DInfo->type[j]==4)
          {
            if(DInfo->rftk1_index[j]>-1 && DInfo->rftk2_index[j]>-1 && DInfo->rftk3_index[j]>-1)
              {
                if(TrackInfo->geninfo_index[DInfo->rftk1_index[j]]>-1 && 
                   TrackInfo->geninfo_index[DInfo->rftk2_index[j]]>-1 && 
                   TrackInfo->geninfo_index[DInfo->rftk3_index[j]]>-1)
                  {
                    if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]>-1 && 
                       GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1 && 
                       GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]]>-1)
                      {
                        if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]])==DpdgId && 
                           GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]] &&
                           GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])
                          {
                            if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j] && 
                               TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j] &&
                               TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j])
                              {
                                Dgen = 23333;
                              }
                            else if((TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk2_MassHypo[j]&&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk1_MassHypo[j]&&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j]&&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==PION_PDGID) ||
                                    (TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk3_MassHypo[j]&&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk1_MassHypo[j]&&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j]&&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID) ||
                                    (TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk3_MassHypo[j]&&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk2_MassHypo[j]&&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j]&&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==PION_PDGID))
                              {
                                Dgen = 23344;
                              }
                          }
                      }
                  }
              }
          }
        else if(DInfo->type[j]==5||DInfo->type[j]==6)
          {
            if(DInfo->rftk1_index[j]>-1 && DInfo->rftk2_index[j]>-1 && DInfo->rftk3_index[j]>-1 && DInfo->rftk4_index[j]>-1)
              {
                if(TrackInfo->geninfo_index[DInfo->rftk1_index[j]]>-1 && 
                   TrackInfo->geninfo_index[DInfo->rftk2_index[j]]>-1 && 
                   TrackInfo->geninfo_index[DInfo->rftk3_index[j]]>-1 && 
                   TrackInfo->geninfo_index[DInfo->rftk4_index[j]]>-1)
                  {
                    if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]>-1 && 
                       GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1 && 
                       GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]]>-1 && 
                       GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]]>-1)
                      {
                        if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]])==DpdgId && 
                           GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]] &&
                           GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]] &&
                           GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])
                          {
                            if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j] && 
                               TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j] &&
                               TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j] &&
                               TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk4_MassHypo[j])
                              {
                                Dgen = 23333;
                              }
                            else if((TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk2_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk1_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk4_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==PION_PDGID &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==PION_PDGID) ||
                                    (TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk3_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk1_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk4_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==PION_PDGID) ||
                                    (TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk4_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk1_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==PION_PDGID) ||
                                    (TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk3_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk2_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk4_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==PION_PDGID &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==PION_PDGID) ||
                                    (TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk4_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk2_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==PION_PDGID &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==PION_PDGID) ||
                                    (TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk4_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk3_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j] &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==PION_PDGID &&
                                     TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID))
                              {
                                Dgen = 23344;
                              }
                          }
                      }
                  }
              }
          }
        int dGenIdxRes = -1;
        if(DInfo->type[j]==7||DInfo->type[j]==8||DInfo->type[j]==9||DInfo->type[j]==10)
          {
            if(DInfo->tktkRes_rftk1_index[j]>-1 && DInfo->tktkRes_rftk2_index[j]>-1)
              {
                if(TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]>-1 && 
                   TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]>-1)
                  {
                    if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]>-1 && 
                       GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]]>-1)
                      {
                        if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]])==RpdgId && 
                           GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])
                          {
                            if(GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]]>-1 &&
                               GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]]]>-1)
                              {
                                if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]]])==DpdgId)
                                  {
                                    if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) && 
                                       TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]))
                                      {
                                        Dgen = 3333;
                                        dGenIdxRes = GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]];
                                      }
                                    if((DInfo->type[j]==9||DInfo->type[j]==10) &&
                                       TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) && 
                                       TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]))
                                      {
                                        Dgen = 3344;
                                        dGenIdxRes = GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]];
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
            if(DInfo->rftk2_index[j]>-1)
              {
                if(TrackInfo->geninfo_index[DInfo->rftk2_index[j]]>-1)
                  {
                    if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1)
                      {
                        if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]])==DpdgId &&
                           GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]==dGenIdxRes)
                          {
                            if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID)
                              {
                                Dgen+=20000;
                              }
                          }
                      }
                  }
              }
          }
        dGenIdxRes = -1;
        if(DInfo->type[j]==11||DInfo->type[j]==12)
          {
            if(DInfo->tktkRes_rftk1_index[j]>-1 && DInfo->tktkRes_rftk2_index[j]>-1 && DInfo->tktkRes_rftk3_index[j]>-1 && DInfo->tktkRes_rftk4_index[j]>-1)
              {
                if(TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]>-1 && 
                   TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]>-1 && 
                   TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]>-1 && 
                   TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]>-1)
                  {
                    if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]>-1 && 
                       GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]]>-1 &&
                       GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]]>-1 &&
                       GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]]>-1)
                      {
                        if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]])==RpdgId && 
                           GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]] &&
                           GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]] &&
                           GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])
                          {
                            if(GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]]>-1)
                              {
                                if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]]])==DpdgId)
                                  {
                                    if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) && 
                                       TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) &&
                                       TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) &&
                                       TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]))
                                      {
                                        Dgen = 3333;
                                        dGenIdxRes = GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]];
                                      }
                                    else if((TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) && 
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==PION_PDGID &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==PION_PDGID) ||
                                            (TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) && 
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==PION_PDGID &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==PION_PDGID) ||
                                            (TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]) && 
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==PION_PDGID &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==PION_PDGID) ||
                                            (TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) && 
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==PION_PDGID &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==PION_PDGID) ||
                                            (TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]) && 
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==PION_PDGID &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==PION_PDGID) ||
                                            (TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]) && 
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==PION_PDGID &&
                                             TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==PION_PDGID))
                                      {
                                        Dgen = 3344;
                                        dGenIdxRes = GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]];
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
            if(DInfo->rftk2_index[j]>-1)
              {
                if(TrackInfo->geninfo_index[DInfo->rftk2_index[j]]>-1)
                  {
                    if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1)
                      {
                        if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]])==DpdgId &&
                           GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]==dGenIdxRes)
                          {
                            if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID)
                              {
                                Dgen+=20000;
                              }
                          }
                      }
                  }
              }
          }

        if(Dgen==23333||Dgen==23344)
          {
            DgenIndex = dGenIdxRes;
            if((DInfo->type[j]==1||DInfo->type[j]==2)&&GenInfo->nDa[DgenIndex]>2) Dgen=41000;
            DgennDa = GenInfo->nDa[DgenIndex];
            Dgenpt = GenInfo->pt[DgenIndex];
            Dgeneta = GenInfo->eta[DgenIndex];
            Dgenphi = GenInfo->phi[DgenIndex];
            b4P->SetXYZM(GenInfo->pt[DgenIndex]*cos(GenInfo->phi[DgenIndex]),
                         GenInfo->pt[DgenIndex]*sin(GenInfo->phi[DgenIndex]),
                         GenInfo->pt[DgenIndex]*sinh(GenInfo->eta[DgenIndex]),
                         GenInfo->mass[DgenIndex]);
            Dgeny = b4P->Rapidity();
          }
      }//if(!real)
  }//fillDtree
  
  bool isDsignalGen(int dmesontype, int j, GenInfoBranches *GenInfo)
  {
    bool flag=false;
    if(dmesontype==1||dmesontype==2)
      {
        if(TMath::Abs(GenInfo->pdgId[j])==DZERO_PDGID&&GenInfo->nDa[j]==2&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1)
          {
            if((((GenInfo->pdgId[GenInfo->da1[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID)||
                 (GenInfo->pdgId[GenInfo->da2[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID))&&dmesontype==1) ||
               (((GenInfo->pdgId[GenInfo->da1[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID)||
                 (GenInfo->pdgId[GenInfo->da2[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID))&&dmesontype==2))
              {
                flag=true;
              }
          }
      }
    if(dmesontype==3||dmesontype==4)
      {
        if(TMath::Abs(GenInfo->pdgId[j])==DPLUS_PDGID&&GenInfo->nDa[j]==3&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1&&GenInfo->da3[j]!=-1)
          {
            if((((GenInfo->pdgId[GenInfo->da1[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID)||
                 (GenInfo->pdgId[GenInfo->da2[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID)||
                 (GenInfo->pdgId[GenInfo->da3[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID))&&dmesontype==4) ||
               (((GenInfo->pdgId[GenInfo->da1[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID)||
                 (GenInfo->pdgId[GenInfo->da2[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID)||
                 (GenInfo->pdgId[GenInfo->da3[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID))&&dmesontype==3))
              {
                flag=true;
              }
          }
      }
    if(dmesontype==5||dmesontype==6)
      {
        if(TMath::Abs(GenInfo->pdgId[j])==DZERO_PDGID&&GenInfo->nDa[j]==4&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1&&GenInfo->da3[j]!=-1&&GenInfo->da4[j]!=-1)
          {
            if((((GenInfo->pdgId[GenInfo->da1[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[j]])==PION_PDGID)||
                 (GenInfo->pdgId[GenInfo->da2[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[j]])==PION_PDGID)||
                 (GenInfo->pdgId[GenInfo->da3[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[j]])==PION_PDGID)||
                 (GenInfo->pdgId[GenInfo->da4[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID))&&dmesontype==6) ||
               (((GenInfo->pdgId[GenInfo->da1[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[j]])==PION_PDGID)||
                 (GenInfo->pdgId[GenInfo->da2[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[j]])==PION_PDGID)||
                 (GenInfo->pdgId[GenInfo->da3[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[j]])==PION_PDGID)||
                 (GenInfo->pdgId[GenInfo->da4[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID))&&dmesontype==5))
              {
                flag=true;
              }
          }
      }
    if(dmesontype==7||dmesontype==8)
      {
        if(TMath::Abs(GenInfo->pdgId[j])==DSUBS_PDGID&&GenInfo->nDa[j]==2&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1)
          {
            if(TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PHI_PDGID)
              {
                if(GenInfo->nDa[GenInfo->da1[j]]==2&&GenInfo->da1[GenInfo->da1[j]]!=-1&&GenInfo->da2[GenInfo->da1[j]]!=-1)
                  {
                    if(TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==KAON_PDGID)
                      {
                        if((GenInfo->pdgId[GenInfo->da2[j]]==PION_PDGID&&dmesontype==7) ||
                           (GenInfo->pdgId[GenInfo->da2[j]]==(0-PION_PDGID)&&dmesontype==8))
                          flag=true;                      
                      }
                  }
              }
          }
      }
    if(dmesontype==9||dmesontype==10)
      {
        if(TMath::Abs(GenInfo->pdgId[j])==DSTAR_PDGID&&GenInfo->nDa[j]==2&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1)
          {
            if(TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==DZERO_PDGID)
              {
                if(GenInfo->nDa[GenInfo->da1[j]]==2&&GenInfo->da1[GenInfo->da1[j]]!=-1&&GenInfo->da2[GenInfo->da1[j]]!=-1)
                  {
                    if((((GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID)||
                         (GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID))&&GenInfo->pdgId[GenInfo->da2[j]]==(0-PION_PDGID)&&dmesontype==10) ||
                       (((GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID)||
                         (GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID))&&GenInfo->pdgId[GenInfo->da2[j]]==PION_PDGID&&dmesontype==9))
                      {
                        flag=true;                      
                      }
                  }
              }
          }
      }
    if(dmesontype==11||dmesontype==12)
      {
        if(TMath::Abs(GenInfo->pdgId[j])==DSTAR_PDGID&&GenInfo->nDa[j]==2&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1)
          {
            if(TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==DZERO_PDGID)
              {
                if(GenInfo->nDa[GenInfo->da1[j]]==4&&GenInfo->da1[GenInfo->da1[j]]!=-1&&GenInfo->da2[GenInfo->da1[j]]!=-1)
                  {
                    if((((GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]])==PION_PDGID)||
                         (GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]])==PION_PDGID)||
                         (GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]])==PION_PDGID)||
                         (GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]])==PION_PDGID))&&GenInfo->pdgId[GenInfo->da2[j]]==(0-PION_PDGID)&&dmesontype==12) ||
                       (((GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]])==PION_PDGID)||
                         (GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]])==PION_PDGID)||
                         (GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]])==PION_PDGID)||
                         (GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]])==PION_PDGID))&&GenInfo->pdgId[GenInfo->da2[j]]==PION_PDGID&&dmesontype==11))
                      {
                        flag=true;                      
                      }
                  }
              }
          }
      }
    return flag;
  }

};

#endif
