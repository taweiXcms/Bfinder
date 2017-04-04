// vim:set ts=4 sw=4 fdm=marker et:
#ifndef _BNTUPLE_H_
#define _BNTUPLE_H_
#include "format.h"

class BntupleBranches
{//{{{
 public:
  float tk1mass[7] = {KAON_MASS, PION_MASS, PION_MASS,   KAON_MASS,  KAON_MASS,  KAON_MASS, PION_MASS};
  float tk2mass[7] = {0,         0,         PION_MASS,   PION_MASS,  PION_MASS,  KAON_MASS, PION_MASS};
  float midmass[7] = {0,         0,         KSHORT_MASS, KSTAR_MASS, KSTAR_MASS, PHI_MASS,  0};

  //EvtInfo
  int      RunNo;
  int      EvtNo;
  int      LumiNo;
  int      Bsize;
  int      Jsize;
  float    PVx;
  float    PVy;
  float    PVz;
  float    PVxE;
  float    PVyE;
  float    PVzE;
  float    PVnchi2;
  float    PVchi2;
  float    BSx;
  float    BSy;
  float    BSz;
  float    BSxErr;
  float    BSyErr;
  float    BSzErr;
  float    BSdxdz;
  float    BSdydz;
  float    BSdxdzErr;
  float    BSdydzErr;
  float    BSWidthX;
  float    BSWidthXErr;
  float    BSWidthY;
  float    BSWidthYErr;
  
  //BInfo
  int       Bindex[MAX_XB];
  int       Btype[MAX_XB];
  float     Bmass[MAX_XB];
  float     Bmass_unfitted[MAX_XB];
  float     Bpt[MAX_XB];
  float     Beta[MAX_XB];
  float     Bphi[MAX_XB];
  float     By[MAX_XB];
  float     BvtxX[MAX_XB];
  float     BvtxY[MAX_XB];
  float     Bd0[MAX_XB];
  float     Bd0Err[MAX_XB];
  float     Bdxyz[MAX_XB];
  float     BdxyzErr[MAX_XB];
  float     Bchi2ndf[MAX_XB];
  float     Bchi2cl[MAX_XB];
  float     Bdtheta[MAX_XB];
  float     Blxy[MAX_XB];
  float     BlxyBS[MAX_XB];
  float     BlxyBSErr[MAX_XB];
  float     BMaxDoca[MAX_XB];
  float     Balpha[MAX_XB];
  float     BsvpvDistance[MAX_XB];
  float     BsvpvDisErr[MAX_XB];
  float     BsvpvDistance_2D[MAX_XB];
  float     BsvpvDisErr_2D[MAX_XB];
  int       Bisbestchi2[MAX_XB];
  
  //BInfo.muonInfo
  float     Bmu1pt[MAX_XB];
  float     Bmu2pt[MAX_XB];
  float     Bmu1p[MAX_XB];
  float     Bmu2p[MAX_XB];
  float     Bmu1eta[MAX_XB];
  float     Bmu2eta[MAX_XB];
  float     Bmu1phi[MAX_XB];
  float     Bmu2phi[MAX_XB];
  float     Bmu1y[MAX_XB];
  float     Bmu2y[MAX_XB];
  float     Bmu1dzPV[MAX_XB];
  float     Bmu2dzPV[MAX_XB];
  float     Bmu1dxyPV[MAX_XB];
  float     Bmu2dxyPV[MAX_XB];
  float     Bmu1normchi2[MAX_XB];
  float     Bmu2normchi2[MAX_XB];
  float     Bmu1Chi2ndf[MAX_XB];
  float     Bmu2Chi2ndf[MAX_XB];
  int       Bmu1muqual[MAX_XB];
  int       Bmu2muqual[MAX_XB];
  bool      Bmu1TrackerMuArbitrated[MAX_XB];
  bool      Bmu2TrackerMuArbitrated[MAX_XB];
  bool      Bmu1isTrackerMuon[MAX_XB];
  bool      Bmu2isTrackerMuon[MAX_XB];
  bool      Bmu1isGlobalMuon[MAX_XB];
  bool      Bmu2isGlobalMuon[MAX_XB];
  bool      Bmu1TMOneStationTight[MAX_XB];
  bool      Bmu2TMOneStationTight[MAX_XB];
  int       Bmu1striphit[MAX_XB];
  int       Bmu2striphit[MAX_XB];
  int       Bmu1pixelhit[MAX_XB];
  int       Bmu2pixelhit[MAX_XB];
  int       Bmu1trackerhit[MAX_XB];
  int       Bmu2trackerhit[MAX_XB];
  int       Bmu1InPixelLayer[MAX_XB];
  int       Bmu2InPixelLayer[MAX_XB];
  int       Bmu1InStripLayer[MAX_XB];
  int       Bmu2InStripLayer[MAX_XB];
  int       Bmu1InTrackerLayer[MAX_XB];
  int       Bmu2InTrackerLayer[MAX_XB];
  int       Bmu1TrkQuality[MAX_XB];
  int       Bmu2TrkQuality[MAX_XB];
  float     Bmu1TrgMatchFilterE[MAX_XB];
  float     Bmu1TrgMatchFilterPt[MAX_XB];
  float     Bmu1TrgMatchFilterEta[MAX_XB];
  float     Bmu1TrgMatchFilterPhi[MAX_XB];
  float     Bmu2TrgMatchFilterE[MAX_XB];
  float     Bmu2TrgMatchFilterPt[MAX_XB];
  float     Bmu2TrgMatchFilterEta[MAX_XB];
  float     Bmu2TrgMatchFilterPhi[MAX_XB];
  
  //BInfo.mumuInfo
  float      Bmumumass[MAX_XB];
  float      Bmumueta[MAX_XB];
  float      Bmumuphi[MAX_XB];
  float      Bmumuy[MAX_XB];
  float      Bmumupt[MAX_XB];
  
  //BInfo.ujInfo
  float     Bujmass[MAX_XB];
  float     BujvProb[MAX_XB];
  float     Bujpt[MAX_XB];
  float     Bujeta[MAX_XB];
  float     Bujphi[MAX_XB];
  float     Bujy[MAX_XB];
  float     Bujlxy[MAX_XB];
  
  //BInfo.trkInfo
  int       Btrk1Idx[MAX_XB];
  int       Btrk2Idx[MAX_XB];
  float     Btrk1Pt[MAX_XB];
  float     Btrk2Pt[MAX_XB];
  float     Btrk1Eta[MAX_XB];
  float     Btrk2Eta[MAX_XB];
  float     Btrk1Phi[MAX_XB];
  float     Btrk2Phi[MAX_XB];
  float     Btrk1PtErr[MAX_XB];
  float     Btrk2PtErr[MAX_XB];
  float     Btrk1EtaErr[MAX_XB];
  float     Btrk2EtaErr[MAX_XB];
  float     Btrk1PhiErr[MAX_XB];
  float     Btrk2PhiErr[MAX_XB];
  float     Btrk1Y[MAX_XB];
  float     Btrk2Y[MAX_XB];
  float     Btrk1Dz[MAX_XB];
  float     Btrk2Dz[MAX_XB];
  float     Btrk1Dxy[MAX_XB];
  float     Btrk2Dxy[MAX_XB];
  float     Btrk1D0[MAX_XB];
  float     Btrk2D0[MAX_XB];
  float     Btrk1D0Err[MAX_XB];
  float     Btrk2D0Err[MAX_XB];
  float     Btrk1PixelHit[MAX_XB];
  float     Btrk2PixelHit[MAX_XB];
  float     Btrk1StripHit[MAX_XB];
  float     Btrk2StripHit[MAX_XB];
  float     Btrk1nPixelLayer[MAX_XB];
  float     Btrk2nPixelLayer[MAX_XB];
  float     Btrk1nStripLayer[MAX_XB];
  float     Btrk2nStripLayer[MAX_XB];
  float     Btrk1Chi2ndf[MAX_XB];
  float     Btrk2Chi2ndf[MAX_XB];
  float     Btrk1MVAVal[MAX_XB];
  float     Btrk2MVAVal[MAX_XB];
  int       Btrk1Algo[MAX_XB];
  int       Btrk2Algo[MAX_XB];
  bool      Btrk1highPurity[MAX_XB];
  bool      Btrk2highPurity[MAX_XB];
  int       Btrk1Quality[MAX_XB];
  int       Btrk2Quality[MAX_XB];
  //BInfo.tktkInfo
  float     Btktkmass[MAX_XB];
  float     BtktkmassKK[MAX_XB];
  float     BtktkvProb[MAX_XB];
  float     Btktkpt[MAX_XB];
  float     Btktketa[MAX_XB];
  float     Btktkphi[MAX_XB];
  float     Btktky[MAX_XB];
  float     Bdoubletmass[MAX_XB];
  float     Bdoubletpt[MAX_XB];
  float     Bdoubleteta[MAX_XB];
  float     Bdoubletphi[MAX_XB];
  float     Bdoublety[MAX_XB];
  
  //BInfo.genInfo
  float     Bgen[MAX_XB];
  int       BgenIndex[MAX_XB];
  float     Bgenpt[MAX_XB];
  float     Bgeneta[MAX_XB];
  float     Bgenphi[MAX_XB];
  float     Bgeny[MAX_XB];
  
  /*
    int      kstar[MAX_XB]; 
    float     Btrk1MassHypo[MAX_XB];
    float     Btrk2MassHypo[MAX_XB];
    float     Btrkminpt[MAX_XB];
    float     Btrkmaxpt[MAX_XB];
    int        Btrkminptindex[MAX_XB];
    int        Btrkmaxptindex[MAX_XB];
  */

  int      Jindex[MAX_XB];
  float    Jmass[MAX_XB];
  float    Jmass_unfitted[MAX_XB];
  float    Jpt[MAX_XB];
  float    Jeta[MAX_XB];
  float    Jphi[MAX_XB];
  float    Jy[MAX_XB];
  float    JvtxX[MAX_XB];
  float    JvtxY[MAX_XB];
  float    Jd0[MAX_XB];
  float    Jdxyz[MAX_XB];
  float    Jchi2ndf[MAX_XB];
  float    Jchi2cl[MAX_XB];
  float    Jdtheta[MAX_XB];
  float    Jlxy[MAX_XB];
  float    JlxyBS[MAX_XB];
  
  //JInfo.muonInfo
  float    Jmu1pt[MAX_XB];
  float    Jmu2pt[MAX_XB];
  float    Jmu1p[MAX_XB];
  float    Jmu2p[MAX_XB];
  float    Jmu1eta[MAX_XB];
  float    Jmu2eta[MAX_XB];
  float    Jmu1phi[MAX_XB];
  float    Jmu2phi[MAX_XB];
  float    Jmu1y[MAX_XB];
  float    Jmu2y[MAX_XB];
  float    Jmu1dzPV[MAX_XB];
  float    Jmu2dzPV[MAX_XB];
  float    Jmu1dxyPV[MAX_XB];
  float    Jmu2dxyPV[MAX_XB];
  float    Jmu1normchi2[MAX_XB];
  float    Jmu2normchi2[MAX_XB];
  float    Jmu1Chi2ndf[MAX_XB];
  float    Jmu2Chi2ndf[MAX_XB];
  int      Jmu1muqual[MAX_XB];
  int      Jmu2muqual[MAX_XB];
  bool     Jmu1TrackerMuArbitrated[MAX_XB];
  bool     Jmu2TrackerMuArbitrated[MAX_XB];
  bool     Jmu1isTrackerMuon[MAX_XB];
  bool     Jmu2isTrackerMuon[MAX_XB];
  bool     Jmu1isGlobalMuon[MAX_XB];
  bool     Jmu2isGlobalMuon[MAX_XB];
  bool     Jmu1TMOneStationTight[MAX_XB];
  bool     Jmu2TMOneStationTight[MAX_XB];
  int      Jmu1striphit[MAX_XB];
  int      Jmu2striphit[MAX_XB];
  int      Jmu1pixelhit[MAX_XB];
  int      Jmu2pixelhit[MAX_XB];
  int      Jmu1trackerhit[MAX_XB];
  int      Jmu2trackerhit[MAX_XB];
  int      Jmu1InPixelLayer[MAX_XB];
  int      Jmu2InPixelLayer[MAX_XB];
  int      Jmu1InStripLayer[MAX_XB];
  int      Jmu2InStripLayer[MAX_XB];
  int      Jmu1InTrackerLayer[MAX_XB];
  int      Jmu2InTrackerLayer[MAX_XB];
  int      Jmu1TrkQuality[MAX_XB];
  int      Jmu2TrkQuality[MAX_XB];
  float    Jmu1TrgMatchFilterE[MAX_XB];
  float    Jmu1TrgMatchFilterPt[MAX_XB];
  float    Jmu1TrgMatchFilterEta[MAX_XB];
  float    Jmu1TrgMatchFilterPhi[MAX_XB];
  float    Jmu2TrgMatchFilterE[MAX_XB];
  float    Jmu2TrgMatchFilterPt[MAX_XB];
  float    Jmu2TrgMatchFilterEta[MAX_XB];
  float    Jmu2TrgMatchFilterPhi[MAX_XB];
  
  float    Jgen[MAX_XB];
  int      JgenIndex[MAX_XB];
  float    Jgenpt[MAX_XB];
  float    Jgeneta[MAX_XB];
  float    Jgenphi[MAX_XB];
  float    Jgeny[MAX_XB];

  void buildBranch(TTree* nt, Bool_t isJpsi=false)
  {
    //EvtInfo
    nt->Branch("RunNo",&RunNo);
    nt->Branch("EvtNo",&EvtNo);
    nt->Branch("LumiNo",&LumiNo);
    nt->Branch("PVx",&PVx);
    nt->Branch("PVy",&PVy);
    nt->Branch("PVz",&PVz);
    nt->Branch("PVxE",&PVxE);
    nt->Branch("PVyE",&PVyE);
    nt->Branch("PVzE",&PVzE);
    nt->Branch("PVnchi2",&PVnchi2);
    nt->Branch("PVchi2",&PVchi2);
    nt->Branch("BSx",&BSx);
    nt->Branch("BSy",&BSy);
    nt->Branch("BSz",&BSz);
    nt->Branch("BSxErr",&BSxErr);
    nt->Branch("BSyErr",&BSyErr);
    nt->Branch("BSzErr",&BSzErr);
    nt->Branch("BSdxdz",&BSdxdz);
    nt->Branch("BSdydz",&BSdydz);
    nt->Branch("BSdxdzErr",&BSdxdzErr);
    nt->Branch("BSdydzErr",&BSdydzErr);
    nt->Branch("BSWidthX",&BSWidthX);
    nt->Branch("BSWidthXErr",&BSWidthXErr);
    nt->Branch("BSWidthY",&BSWidthY);
    nt->Branch("BSWidthYErr",&BSWidthYErr);
    
    if(isJpsi)
      {
        nt->Branch("Jsize",&Jsize);
        nt->Branch("Jindex",Jindex,"Jindex[Jsize]/I");
        nt->Branch("Jmass",Jmass,"Jmass[Jsize]/F");
        nt->Branch("Jmass_unfitted",Jmass_unfitted,"Jmass_unfitted[Jsize]/F");
        nt->Branch("Jpt",Jpt,"Jpt[Jsize]/F");
        nt->Branch("Jeta",Jeta,"Jeta[Jsize]/F");
        nt->Branch("Jphi",Jphi,"Jphi[Jsize]/F");
        nt->Branch("Jy",Jy,"Jy[Jsize]/F");
        nt->Branch("JvtxX",JvtxX,"JvtxX[Jsize]/F");
        nt->Branch("JvtxY",JvtxY,"JvtxY[Jsize]/F");
        nt->Branch("Jd0",Jd0,"Jd0[Jsize]/F");
        nt->Branch("Jdxyz",Jdxyz,"Jdxyz[Jsize]/F");
        nt->Branch("Jchi2ndf",Jchi2ndf,"Jchi2ndf[Jsize]/F");
        nt->Branch("Jchi2cl",Jchi2cl,"Jchi2cl[Jsize]/F");
        nt->Branch("Jdtheta",Jdtheta,"Jdtheta[Jsize]/F");
        nt->Branch("Jlxy",Jlxy,"Jlxy[Jsize]/F");
        nt->Branch("JlxyBS",JlxyBS,"JlxyBS[Jsize]/F");
        
        nt->Branch("Jmu1pt",Jmu1pt,"Jmu1pt[Jsize]/F");
        nt->Branch("Jmu2pt",Jmu2pt,"Jmu2pt[Jsize]/F");
        nt->Branch("Jmu1p",Jmu1p,"Jmu1p[Jsize]/F");
        nt->Branch("Jmu2p",Jmu2p,"Jmu2p[Jsize]/F");
        nt->Branch("Jmu1eta",Jmu1eta,"Jmu1eta[Jsize]/F");
        nt->Branch("Jmu2eta",Jmu2eta,"Jmu2eta[Jsize]/F");
        nt->Branch("Jmu1phi",Jmu1phi,"Jmu1phi[Jsize]/F");
        nt->Branch("Jmu2phi",Jmu2phi,"Jmu2phi[Jsize]/F");
        nt->Branch("Jmu1y",Jmu1y,"Jmu1y[Jsize]/F");
        nt->Branch("Jmu2y",Jmu2y,"Jmu2y[Jsize]/F");
        nt->Branch("Jmu1dzPV",Jmu1dzPV,"Jmu1dzPV[Jsize]/F");
        nt->Branch("Jmu2dzPV",Jmu2dzPV,"Jmu2dzPV[Jsize]/F");
        nt->Branch("Jmu1dxyPV",Jmu1dxyPV,"Jmu1dxyPV[Jsize]/F");
        nt->Branch("Jmu2dxyPV",Jmu2dxyPV,"Jmu2dxyPV[Jsize]/F");
        nt->Branch("Jmu1normchi2",Jmu1normchi2,"Jmu1normchi2[Jsize]/F");
        nt->Branch("Jmu2normchi2",Jmu2normchi2,"Jmu2normchi2[Jsize]/F");
        nt->Branch("Jmu1Chi2ndf",Jmu1Chi2ndf,"Jmu1Chi2ndf[Jsize]/F");
        nt->Branch("Jmu2Chi2ndf",Jmu2Chi2ndf,"Jmu2Chi2ndf[Jsize]/F");
        nt->Branch("Jmu1muqual",Jmu1muqual,"Jmu1muqual[Jsize]/I");
        nt->Branch("Jmu2muqual",Jmu2muqual,"Jmu2muqual[Jsize]/I");
        nt->Branch("Jmu1TrackerMuArbitrated",Jmu1TrackerMuArbitrated,"Jmu1TrackerMuArbitrated[Jsize]/O");
        nt->Branch("Jmu2TrackerMuArbitrated",Jmu2TrackerMuArbitrated,"Jmu2TrackerMuArbitrated[Jsize]/O");
        nt->Branch("Jmu1isTrackerMuon",Jmu1isTrackerMuon,"Jmu1isTrackerMuon[Jsize]/O");
        nt->Branch("Jmu2isTrackerMuon",Jmu2isTrackerMuon,"Jmu2isTrackerMuon[Jsize]/O");
        nt->Branch("Jmu1isGlobalMuon",Jmu1isGlobalMuon,"Jmu1isGlobalMuon[Jsize]/O");
        nt->Branch("Jmu2isGlobalMuon",Jmu2isGlobalMuon,"Jmu2isGlobalMuon[Jsize]/O");
        nt->Branch("Jmu1TMOneStationTight",Jmu1TMOneStationTight,"Jmu1TMOneStationTight[Jsize]/O");
        nt->Branch("Jmu2TMOneStationTight",Jmu2TMOneStationTight,"Jmu2TMOneStationTight[Jsize]/O");
        nt->Branch("Jmu1striphit",Jmu1striphit,"Jmu1striphit[Jsize]/I");
        nt->Branch("Jmu2striphit",Jmu2striphit,"Jmu2striphit[Jsize]/I");
        nt->Branch("Jmu1pixelhit",Jmu1pixelhit,"Jmu1pixelhit[Jsize]/I");
        nt->Branch("Jmu2pixelhit",Jmu2pixelhit,"Jmu2pixelhit[Jsize]/I");
        nt->Branch("Jmu1trackerhit",Jmu1trackerhit,"Jmu1trackerhit[Jsize]/I");
        nt->Branch("Jmu2trackerhit",Jmu2trackerhit,"Jmu2trackerhit[Jsize]/I");
        nt->Branch("Jmu1InPixelLayer",Jmu1InPixelLayer,"Jmu1InPixelLayer[Jsize]/I");
        nt->Branch("Jmu2InPixelLayer",Jmu2InPixelLayer,"Jmu2InPixelLayer[Jsize]/I");
        nt->Branch("Jmu1InStripLayer",Jmu1InStripLayer,"Jmu1InStripLayer[Jsize]/I");
        nt->Branch("Jmu2InStripLayer",Jmu2InStripLayer,"Jmu2InStripLayer[Jsize]/I");
        nt->Branch("Jmu1InTrackerLayer",Jmu1InTrackerLayer,"Jmu1InTrackerLayer[Jsize]/I");
        nt->Branch("Jmu2InTrackerLayer",Jmu2InTrackerLayer,"Jmu2InTrackerLayer[Jsize]/I");
        nt->Branch("Jmu1TrkQuality",Jmu1TrkQuality,"Jmu1TrkQuality[Jsize]/I");
        nt->Branch("Jmu2TrkQuality",Jmu2TrkQuality,"Jmu2TrkQuality[Jsize]/I");
        nt->Branch("Jmu1TrgMatchFilterE"  ,Bmu1TrgMatchFilterE,  "Bmu1TrgMatchFilterE[Jsize]/F");
        nt->Branch("Jmu1TrgMatchFilterPt" ,Bmu1TrgMatchFilterPt, "Bmu1TrgMatchFilterPt[Jsize]/F");
        nt->Branch("Jmu1TrgMatchFilterEta",Jmu1TrgMatchFilterEta,"Jmu1TrgMatchFilterEta[Jsize]/F");
        nt->Branch("Jmu1TrgMatchFilterPhi",Jmu1TrgMatchFilterPhi,"Jmu1TrgMatchFilterPhi[Jsize]/F");
        nt->Branch("Jmu2TrgMatchFilterE"  ,Bmu2TrgMatchFilterE,  "Bmu2TrgMatchFilterE[Jsize]/F");
        nt->Branch("Jmu2TrgMatchFilterPt" ,Bmu2TrgMatchFilterPt, "Bmu2TrgMatchFilterPt[Jsize]/F");
        nt->Branch("Jmu2TrgMatchFilterEta",Jmu2TrgMatchFilterEta,"Jmu2TrgMatchFilterEta[Jsize]/F");
        nt->Branch("Jmu2TrgMatchFilterPhi",Jmu2TrgMatchFilterPhi,"Jmu2TrgMatchFilterPhi[Jsize]/F");

        nt->Branch("Jgen",Jgen,"Jgen[Jsize]/F");
        nt->Branch("JgenIndex",JgenIndex,"JgenIndex[Jsize]/I");
        nt->Branch("Jgenpt",Jgenpt,"Jgenpt[Jsize]/F");
        nt->Branch("Jgeny",Jgeny,"Jgeny[Jsize]/F");
        nt->Branch("Jgeneta",Jgeneta,"Jgeneta[Jsize]/F");
        nt->Branch("Jgenphi",Jgenphi,"Jgenphi[Jsize]/F");
      }
    else
      {
        nt->Branch("Bsize",&Bsize);
        //BInfo
        nt->Branch("Bindex",Bindex,"Bindex[Bsize]/I");
        nt->Branch("Btype",Btype,"Btype[Bsize]/I");
        nt->Branch("Bmass",Bmass,"Bmass[Bsize]/F");
        nt->Branch("Bmass_unfitted",Bmass_unfitted,"Bmass_unfitted[Bsize]/F");
        nt->Branch("Bpt",Bpt,"Bpt[Bsize]/F");
        nt->Branch("Beta",Beta,"Beta[Bsize]/F");
        nt->Branch("Bphi",Bphi,"Bphi[Bsize]/F");
        nt->Branch("By",By,"By[Bsize]/F");
        nt->Branch("BvtxX",BvtxX,"BvtxX[Bsize]/F");
        nt->Branch("BvtxY",BvtxY,"BvtxY[Bsize]/F");
        nt->Branch("Bd0",Bd0,"Bd0[Bsize]/F");
        nt->Branch("Bd0Err",Bd0Err,"Bd0Err[Bsize]/F");
        nt->Branch("Bdxyz",Bdxyz,"Bdxyz[Bsize]/F");
        nt->Branch("BdxyzErr",BdxyzErr,"BdxyzErr[Bsize]/F");
        nt->Branch("Bchi2ndf",Bchi2ndf,"Bchi2ndf[Bsize]/F");
        nt->Branch("Bchi2cl",Bchi2cl,"Bchi2cl[Bsize]/F");
        nt->Branch("Bdtheta",Bdtheta,"Bdtheta[Bsize]/F");
        nt->Branch("Blxy",Blxy,"Blxy[Bsize]/F");
        nt->Branch("BlxyBS",BlxyBS,"BlxyBS[Bsize]/F");
        nt->Branch("BlxyBSErr",BlxyBSErr,"BlxyBSErr[Bsize]/F");
        nt->Branch("Balpha",Balpha,"Balpha[Bsize]/F");
        nt->Branch("BsvpvDistance",BsvpvDistance,"BsvpvDistance[Bsize]/F");
        nt->Branch("BsvpvDisErr",BsvpvDisErr,"BsvpvDisErr[Bsize]/F");
        nt->Branch("BsvpvDistance_2D",BsvpvDistance_2D,"BsvpvDistance_2D[Bsize]/F");
        nt->Branch("BsvpvDisErr_2D",BsvpvDisErr_2D,"BsvpvDisErr_2D[Bsize]/F");
        nt->Branch("BMaxDoca",BMaxDoca,"BMaxDoca[Bsize]/F");
        nt->Branch("Bisbestchi2",Bisbestchi2,"Bisbestchi2[Bsize]/I");
        
        //BInfo.trkInfo
        nt->Branch("Btrk1Idx",Btrk1Idx,"Btrk1Idx[Bsize]/I");
        nt->Branch("Btrk2Idx",Btrk2Idx,"Btrk2Idx[Bsize]/I");
        nt->Branch("Btrk1Pt",Btrk1Pt,"Btrk1Pt[Bsize]/F");
        nt->Branch("Btrk2Pt",Btrk2Pt,"Btrk2Pt[Bsize]/F");
        nt->Branch("Btrk1Eta",Btrk1Eta,"Btrk1Eta[Bsize]/F");  
        nt->Branch("Btrk2Eta",Btrk2Eta,"Btrk2Eta[Bsize]/F");  
        nt->Branch("Btrk1Phi",Btrk1Phi,"Btrk1Phi[Bsize]/F");  
        nt->Branch("Btrk2Phi",Btrk2Phi,"Btrk2Phi[Bsize]/F");  
        nt->Branch("Btrk1PtErr",Btrk1PtErr,"Btrk1PtErr[Bsize]/F");  
        nt->Branch("Btrk2PtErr",Btrk2PtErr,"Btrk2PtErr[Bsize]/F");
        nt->Branch("Btrk1EtaErr",Btrk1EtaErr,"Btrk1EtaErr[Bsize]/F");
        nt->Branch("Btrk2EtaErr",Btrk2EtaErr,"Btrk2EtaErr[Bsize]/F");
        nt->Branch("Btrk1PhiErr",Btrk1PhiErr,"Btrk1PhiErr[Bsize]/F");
        nt->Branch("Btrk2PhiErr",Btrk2PhiErr,"Btrk2PhiErr[Bsize]/F");
        nt->Branch("Btrk1Y",Btrk1Y,"Btrk1Y[Bsize]/F");  
        nt->Branch("Btrk2Y",Btrk2Y,"Btrk2Y[Bsize]/F");  
        nt->Branch("Btrk1Dz",Btrk1Dz,"Btrk1Dz[Bsize]/F");
        nt->Branch("Btrk2Dz",Btrk2Dz,"Btrk2Dz[Bsize]/F");
        nt->Branch("Btrk1Dxy",Btrk1Dxy,"Btrk1Dxy[Bsize]/F");
        nt->Branch("Btrk2Dxy",Btrk2Dxy,"Btrk2Dxy[Bsize]/F");
        nt->Branch("Btrk1D0",Btrk1D0,"Btrk1D0[Bsize]/F");
        nt->Branch("Btrk2D0",Btrk2D0,"Btrk2D0[Bsize]/F");
        nt->Branch("Btrk1D0Err",Btrk1D0Err,"Btrk1D0Err[Bsize]/F");
        nt->Branch("Btrk2D0Err",Btrk2D0Err,"Btrk2D0Err[Bsize]/F");
        nt->Branch("Btrk1PixelHit",Btrk1PixelHit,"Btrk1PixelHit[Bsize]/F");
        nt->Branch("Btrk2PixelHit",Btrk2PixelHit,"Btrk2PixelHit[Bsize]/F");
        nt->Branch("Btrk1StripHit",Btrk1StripHit,"Btrk1StripHit[Bsize]/F");
        nt->Branch("Btrk2StripHit",Btrk2StripHit,"Btrk2StripHit[Bsize]/F");
        nt->Branch("Btrk1nPixelLayer",Btrk1nPixelLayer,"Btrk1nPixelLayer[Bsize]/F");
        nt->Branch("Btrk2nPixelLayer",Btrk2nPixelLayer,"Btrk2nPixelLayer[Bsize]/F");
        nt->Branch("Btrk1nStripLayer",Btrk1nStripLayer,"Btrk1nStripLayer[Bsize]/F");
        nt->Branch("Btrk2nStripLayer",Btrk2nStripLayer,"Btrk2nStripLayer[Bsize]/F");
        nt->Branch("Btrk1Chi2ndf",Btrk1Chi2ndf,"Btrk1Chi2ndf[Bsize]/F");
        nt->Branch("Btrk2Chi2ndf",Btrk2Chi2ndf,"Btrk2Chi2ndf[Bsize]/F");
        nt->Branch("Btrk1MVAVal",Btrk1MVAVal,"Btrk1MVAVal[Bsize]/F");
        nt->Branch("Btrk2MVAVal",Btrk2MVAVal,"Btrk2MVAVal[Bsize]/F");
        nt->Branch("Btrk1Algo",Btrk1Algo,"Btrk1Algo[Bsize]/I");
        nt->Branch("Btrk2Algo",Btrk2Algo,"Btrk2Algo[Bsize]/I");
        nt->Branch("Btrk1highPurity",Btrk1highPurity,"Btrk1highPurity[Bsize]/O");
        nt->Branch("Btrk2highPurity",Btrk2highPurity,"Btrk2highPurity[Bsize]/O");
        nt->Branch("Btrk1Quality",Btrk1Quality,"Btrk1Quality[Bsize]/I");
        nt->Branch("Btrk2Quality",Btrk2Quality,"Btrk2Quality[Bsize]/I");
        
        //BInfo.tktkInfo
        nt->Branch("Btktkmass",Btktkmass,"Btktkmass[Bsize]/F");
        nt->Branch("BtktkmassKK",BtktkmassKK,"BtktkmassKK[Bsize]/F");
        nt->Branch("BtktkvProb",BtktkvProb,"BtktkvProb[Bsize]/F");
        nt->Branch("Btktkpt",Btktkpt,"Btktkpt[Bsize]/F");
        nt->Branch("Btktketa",Btktketa,"Btktketa[Bsize]/F");
        nt->Branch("Btktkphi",Btktkphi,"Btktkphi[Bsize]/F");
        nt->Branch("Btktky",Btktky,"Btktky[Bsize]/F");
        nt->Branch("Bdoubletmass",Bdoubletmass,"Bdoubletmass[Bsize]/F");
        nt->Branch("Bdoubletpt",Bdoubletpt,"Bdoubletpt[Bsize]/F");
        nt->Branch("Bdoubleteta",Bdoubleteta,"Bdoubleteta[Bsize]/F");  
        nt->Branch("Bdoubletphi",Bdoubletphi,"Bdoubletphi[Bsize]/F");  
        nt->Branch("Bdoublety",Bdoublety,"Bdoublety[Bsize]/F");
        
        //BInfo.muonInfo
        nt->Branch("Bmu1pt",Bmu1pt,"Bmu1pt[Bsize]/F");
        nt->Branch("Bmu2pt",Bmu2pt,"Bmu2pt[Bsize]/F");
        nt->Branch("Bmu1p",Bmu1p,"Bmu1p[Bsize]/F");
        nt->Branch("Bmu2p",Bmu2p,"Bmu2p[Bsize]/F");
        nt->Branch("Bmu1eta",Bmu1eta,"Bmu1eta[Bsize]/F");
        nt->Branch("Bmu2eta",Bmu2eta,"Bmu2eta[Bsize]/F");
        nt->Branch("Bmu1phi",Bmu1phi,"Bmu1phi[Bsize]/F");
        nt->Branch("Bmu2phi",Bmu2phi,"Bmu2phi[Bsize]/F");
        nt->Branch("Bmu1y",Bmu1y,"Bmu1y[Bsize]/F");
        nt->Branch("Bmu2y",Bmu2y,"Bmu2y[Bsize]/F");
        nt->Branch("Bmu1dzPV",Bmu1dzPV,"Bmu1dzPV[Bsize]/F");
        nt->Branch("Bmu2dzPV",Bmu2dzPV,"Bmu2dzPV[Bsize]/F");
        nt->Branch("Bmu1dxyPV",Bmu1dxyPV,"Bmu1dxyPV[Bsize]/F");
        nt->Branch("Bmu2dxyPV",Bmu2dxyPV,"Bmu2dxyPV[Bsize]/F");
        nt->Branch("Bmu1normchi2",Bmu1normchi2,"Bmu1normchi2[Bsize]/F");
        nt->Branch("Bmu2normchi2",Bmu2normchi2,"Bmu2normchi2[Bsize]/F");
        nt->Branch("Bmu1Chi2ndf",Bmu1Chi2ndf,"Bmu1Chi2ndf[Bsize]/F");
        nt->Branch("Bmu2Chi2ndf",Bmu2Chi2ndf,"Bmu2Chi2ndf[Bsize]/F");
        nt->Branch("Bmu1muqual",Bmu1muqual,"Bmu1muqual[Bsize]/I");
        nt->Branch("Bmu2muqual",Bmu2muqual,"Bmu2muqual[Bsize]/I");
        nt->Branch("Bmu1TrackerMuArbitrated",Bmu1TrackerMuArbitrated,"Bmu1TrackerMuArbitrated[Bsize]/O");
        nt->Branch("Bmu2TrackerMuArbitrated",Bmu2TrackerMuArbitrated,"Bmu2TrackerMuArbitrated[Bsize]/O");
        nt->Branch("Bmu1isTrackerMuon",Bmu1isTrackerMuon,"Bmu1isTrackerMuon[Bsize]/O");
        nt->Branch("Bmu2isTrackerMuon",Bmu2isTrackerMuon,"Bmu2isTrackerMuon[Bsize]/O");
        nt->Branch("Bmu1isGlobalMuon",Bmu1isGlobalMuon,"Bmu1isGlobalMuon[Bsize]/O");
        nt->Branch("Bmu2isGlobalMuon",Bmu2isGlobalMuon,"Bmu2isGlobalMuon[Bsize]/O");
        nt->Branch("Bmu1TMOneStationTight",Bmu1TMOneStationTight,"Bmu1TMOneStationTight[Bsize]/O");
        nt->Branch("Bmu2TMOneStationTight",Bmu2TMOneStationTight,"Bmu2TMOneStationTight[Bsize]/O");
        nt->Branch("Bmu1striphit",Bmu1striphit,"Bmu1striphit[Bsize]/I");
        nt->Branch("Bmu2striphit",Bmu2striphit,"Bmu2striphit[Bsize]/I");
        nt->Branch("Bmu1pixelhit",Bmu1pixelhit,"Bmu1pixelhit[Bsize]/I");
        nt->Branch("Bmu2pixelhit",Bmu2pixelhit,"Bmu2pixelhit[Bsize]/I");
        nt->Branch("Bmu1trackerhit",Bmu1trackerhit,"Bmu1trackerhit[Bsize]/I");
        nt->Branch("Bmu2trackerhit",Bmu2trackerhit,"Bmu2trackerhit[Bsize]/I");
        nt->Branch("Bmu1InPixelLayer",Bmu1InPixelLayer,"Bmu1InPixelLayer[Bsize]/I");
        nt->Branch("Bmu2InPixelLayer",Bmu2InPixelLayer,"Bmu2InPixelLayer[Bsize]/I");
        nt->Branch("Bmu1InStripLayer",Bmu1InStripLayer,"Bmu1InStripLayer[Bsize]/I");
        nt->Branch("Bmu2InStripLayer",Bmu2InStripLayer,"Bmu2InStripLayer[Bsize]/I");
        nt->Branch("Bmu1InTrackerLayer",Bmu1InTrackerLayer,"Bmu1InTrackerLayer[Bsize]/I");
        nt->Branch("Bmu2InTrackerLayer",Bmu2InTrackerLayer,"Bmu2InTrackerLayer[Bsize]/I");
        nt->Branch("Bmu1TrkQuality",Bmu1TrkQuality,"Bmu1TrkQuality[Bsize]/I");
        nt->Branch("Bmu2TrkQuality",Bmu2TrkQuality,"Bmu2TrkQuality[Bsize]/I");
        nt->Branch("Bmu1TrgMatchFilterE"  ,Bmu1TrgMatchFilterE,  "Bmu1TrgMatchFilterE[Bsize]/F");
        nt->Branch("Bmu1TrgMatchFilterPt" ,Bmu1TrgMatchFilterPt, "Bmu1TrgMatchFilterPt[Bsize]/F");
        nt->Branch("Bmu1TrgMatchFilterEta",Bmu1TrgMatchFilterEta,"Bmu1TrgMatchFilterEta[Bsize]/F");
        nt->Branch("Bmu1TrgMatchFilterPhi",Bmu1TrgMatchFilterPhi,"Bmu1TrgMatchFilterPhi[Bsize]/F");
        nt->Branch("Bmu2TrgMatchFilterE"  ,Bmu2TrgMatchFilterE,  "Bmu2TrgMatchFilterE[Bsize]/F");
        nt->Branch("Bmu2TrgMatchFilterPt" ,Bmu2TrgMatchFilterPt, "Bmu2TrgMatchFilterPt[Bsize]/F");
        nt->Branch("Bmu2TrgMatchFilterEta",Bmu2TrgMatchFilterEta,"Bmu2TrgMatchFilterEta[Bsize]/F");
        nt->Branch("Bmu2TrgMatchFilterPhi",Bmu2TrgMatchFilterPhi,"Bmu2TrgMatchFilterPhi[Bsize]/F");
        nt->Branch("Bmumumass",Bmumumass,"Bmumumass[Bsize]/F");
        nt->Branch("Bmumueta",Bmumueta,"Bmumueta[Bsize]/F");
        nt->Branch("Bmumuphi",Bmumuphi,"Bmumuphi[Bsize]/F");
        nt->Branch("Bmumuy",Bmumuy,"Bmumuy[Bsize]/F");
        nt->Branch("Bmumupt",Bmumupt,"Bmumupt[Bsize]/F");
        nt->Branch("Bujmass",Bujmass,"Bujmass[Bsize]/F");
        nt->Branch("BujvProb",BujvProb,"BujvProb[Bsize]/F");
        nt->Branch("Bujpt",Bujpt,"Bujpt[Bsize]/F");
        nt->Branch("Bujeta",Bujeta,"Bujeta[Bsize]/F");
        nt->Branch("Bujphi",Bujphi,"Bujphi[Bsize]/F");
        nt->Branch("Bujy",Bujy,"Bujy[Bsize]/F");
        nt->Branch("Bujlxy",Bujlxy,"Bujlxy[Bsize]/F");
        
        //BInfo.genInfo
        nt->Branch("Bgen",Bgen,"Bgen[Bsize]/F");
        nt->Branch("BgenIndex",BgenIndex,"BgenIndex[Bsize]/I");
        nt->Branch("Bgenpt",Bgenpt,"Bgenpt[Bsize]/F");
        nt->Branch("Bgeny",Bgeny,"Bgeny[Bsize]/F");
        nt->Branch("Bgeneta",Bgeneta,"Bgeneta[Bsize]/F");
        nt->Branch("Bgenphi",Bgenphi,"Bgenphi[Bsize]/F");
        
        /*
          nt->Branch("Btrk1MassHypo",Btrk1MassHypo,"Btrk1MassHypo[Bsize]/F");
          nt->Branch("Btrk2MassHypo",Btrk2MassHypo,"Btrk2MassHypo[Bsize]/F");
          nt->Branch("Btrkminpt",Btrkminpt,"Btrkminpt[Bsize]/F");
          nt->Branch("Btrkmaxpt",Btrkmaxpt,"Btrkmaxpt[Bsize]/F");
          nt->Branch("Btrkminptindex",Btrkminptindex,"Btrkminptindex[Bsize]/I");
          nt->Branch("Btrkmaxptindex",Btrkmaxptindex,"Btrkmaxptindex[Bsize]/I");
        */
      }
  }
  
  int      Gsize;
  float   Gy[MAX_GEN];
  float   Geta[MAX_GEN];
  float   Gphi[MAX_GEN];
  float   Gpt[MAX_GEN];
  int     GpdgId[MAX_GEN];
  int     GcollisionId[MAX_GEN];
  int     GisSignal[MAX_GEN];
  float   Gmu1pt[MAX_GEN];
  float   Gmu2pt[MAX_GEN];
  float   Gmu1p[MAX_GEN];
  float   Gmu2p[MAX_GEN];
  float   Gmu1eta[MAX_GEN];
  float   Gmu2eta[MAX_GEN];
  float   Gmu1phi[MAX_GEN];
  float   Gmu2phi[MAX_GEN];
  float   Gtk1pt[MAX_GEN];
  float   Gtk2pt[MAX_GEN];
  float   Gtk1eta[MAX_GEN];
  float   Gtk2eta[MAX_GEN];
  float   Gtk1phi[MAX_GEN];
  float   Gtk2phi[MAX_GEN];
  
  void buildGenBranch(TTree* nt)
  {
    nt->Branch("Gsize",&Gsize);
    nt->Branch("Gy",Gy,"Gy[Gsize]/F");
    nt->Branch("Geta",Geta,"Geta[Gsize]/F");
    nt->Branch("Gphi",Gphi,"Gphi[Gsize]/F");
    nt->Branch("Gpt",Gpt,"Gpt[Gsize]/F");
    nt->Branch("GpdgId",GpdgId,"GpdgId[Gsize]/I");
    nt->Branch("GcollisionId",GcollisionId,"GcollisionId[Gsize]/I");
    nt->Branch("GisSignal",GisSignal,"GisSignal[Gsize]/I");
    nt->Branch("Gmu1eta",Gmu1eta,"Gmu1eta[Gsize]/F");
    nt->Branch("Gmu1phi",Gmu1phi,"Gmu1phi[Gsize]/F");
    nt->Branch("Gmu1pt",Gmu1pt,"Gmu1pt[Gsize]/F");
    nt->Branch("Gmu1p",Gmu1p,"Gmu1p[Gsize]/F");
    nt->Branch("Gmu2eta",Gmu2eta,"Gmu2eta[Gsize]/F");
    nt->Branch("Gmu2phi",Gmu2phi,"Gmu2phi[Gsize]/F");
    nt->Branch("Gmu2pt",Gmu2pt,"Gmu2pt[Gsize]/F");
    nt->Branch("Gmu2p",Gmu2p,"Gmu2p[Gsize]/F");
    nt->Branch("Gtk1pt",Gtk1pt,"Gtk1pt[Gsize]/F");
    nt->Branch("Gtk1eta",Gtk1eta,"Gtk1eta[Gsize]/F");
    nt->Branch("Gtk1phi",Gtk1phi,"Gtk1phi[Gsize]/F");
    nt->Branch("Gtk2pt",Gtk2pt,"Gtk2pt[Gsize]/F");
    nt->Branch("Gtk2eta",Gtk2eta,"Gtk2eta[Gsize]/F");
    nt->Branch("Gtk2phi",Gtk2phi,"Gtk2phi[Gsize]/F");
  }
  
  void makeNtuple(int ifchannel[], bool REAL, bool skim, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, MuonInfoBranches *MuonInfo, TrackInfoBranches *TrackInfo, BInfoBranches *BInfo, GenInfoBranches *GenInfo, TTree* nt0, TTree* nt1, TTree* nt2, TTree* nt3, TTree* nt5, TTree* nt6, TTree* nt7)
  {//{{{
    TVector3* bP = new TVector3;
    TVector3* bVtx = new TVector3;
    TLorentzVector* b4P = new TLorentzVector;
    fillTreeEvt(EvtInfo);
    int Btypesize[8]={0,0,0,0,0,0,0,0};
    for(int t=0;t<7;t++)
      {
        int tidx = t-1;
        if(t!=4)
          {
            tidx = t;
            Bsize = 0;
          }
        if(ifchannel[t]==1)
          {
            int bestindex=-1;
            for(int j=0;j<BInfo->size;j++)
              {
                if(skim)
                  {
                    //if(BInfo->pt[j]<3.) continue;
                    //if(BInfo->pt[j]<10.) continue;
                    if(!( (BInfo->pt[j] > 7. && BInfo->pt[j] < 10. && BInfo->svpvDistance[j]/BInfo->svpvDisErr[j] > 5.5) || (BInfo->pt[j] > 10. && BInfo->svpvDistance[j]/BInfo->svpvDisErr[j] > 3.5) )) continue;
                    if(BInfo->mass[j]<5. || BInfo->mass[j]>6.) continue;
                  }
                if(BInfo->type[j]==(t+1))
                  {
                    fillTree(bP,bVtx,b4P,j,Btypesize[tidx],tk1mass[t],tk2mass[t],REAL, EvtInfo, VtxInfo, MuonInfo, TrackInfo, BInfo, GenInfo, bestindex);
                    Btypesize[tidx]++;
                  }
              }
            if(t==0)      nt0->Fill();
            else if(t==1) nt1->Fill();
            else if(t==2) nt2->Fill();
            else if(t==4) nt3->Fill();
            else if(t==5) nt5->Fill();
            else if(t==6) nt6->Fill();
          }
      }
    Jsize = 0;
    if(ifchannel[7]==1)
      {
        for(int j=0;j<BInfo->uj_size;j++)
          {
            if(skim)
              {
                ;
              }
            fillJpsiTree(bP, bVtx, b4P, j, Btypesize[7], REAL, EvtInfo, VtxInfo, MuonInfo, TrackInfo, BInfo, GenInfo);
            Btypesize[7]++;
          }
        nt7->Fill();
      }
  }//}}}
  
  void fillGenTree(TTree* ntGen, GenInfoBranches *GenInfo, bool gskim=true)
  {//{{{
    TLorentzVector* bGen = new TLorentzVector;
    int gt=0,sigtype=0;
    int gsize=0;
    Gsize = 0;
    for(int j=0;j<GenInfo->size;j++)
      {
        if((TMath::Abs(GenInfo->pdgId[j])!=BPLUS_PDGID&&TMath::Abs(GenInfo->pdgId[j])!=BZERO_PDGID&&TMath::Abs(GenInfo->pdgId[j])!=BSUBS_PDGID&&TMath::Abs(GenInfo->pdgId[j])!=JPSI_PDGID) && gskim) continue;
        Gsize = gsize+1;
        Gpt[gsize] = GenInfo->pt[j];
        Geta[gsize] = GenInfo->eta[j];
        Gphi[gsize] = GenInfo->phi[j];
        GpdgId[gsize] = GenInfo->pdgId[j];
        GcollisionId[gsize] = GenInfo->collisionId[j];
        bGen->SetPtEtaPhiM(GenInfo->pt[j],GenInfo->eta[j],GenInfo->phi[j],GenInfo->mass[j]);
        Gy[gsize] = bGen->Rapidity();
        sigtype=0;
        for(gt=1;gt<10;gt++)
          {
            if(signalGen(gt,j,GenInfo))
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
            if(sigtype==8 || sigtype==9)
              {
                Gmu1pt[gsize] = GenInfo->pt[GenInfo->da1[j]];
                Gmu1eta[gsize] = GenInfo->eta[GenInfo->da1[j]];
                Gmu1phi[gsize] = GenInfo->phi[GenInfo->da1[j]];
                Gmu1p[gsize] = Gmu1pt[gsize]*cosh(Gmu1eta[gsize]);
                Gmu2pt[gsize] = GenInfo->pt[GenInfo->da2[j]];
                Gmu2eta[gsize] = GenInfo->eta[GenInfo->da2[j]];
                Gmu2phi[gsize] = GenInfo->phi[GenInfo->da2[j]];
                Gmu2p[gsize] = Gmu2pt[gsize]*cosh(Gmu2eta[gsize]);
              }
            else
              {
                Gmu1pt[gsize] = GenInfo->pt[GenInfo->da1[GenInfo->da1[j]]];
                Gmu1eta[gsize] = GenInfo->eta[GenInfo->da1[GenInfo->da1[j]]];
                Gmu1phi[gsize] = GenInfo->phi[GenInfo->da1[GenInfo->da1[j]]];
                Gmu1p[gsize] = Gmu1pt[gsize]*cosh(Gmu1eta[gsize]);
                Gmu2pt[gsize] = GenInfo->pt[GenInfo->da2[GenInfo->da1[j]]];
                Gmu2eta[gsize] = GenInfo->eta[GenInfo->da2[GenInfo->da1[j]]];
                Gmu2phi[gsize] = GenInfo->phi[GenInfo->da2[GenInfo->da1[j]]];
                Gmu2p[gsize] = Gmu2pt[gsize]*cosh(Gmu2eta[gsize]);
                if(sigtype==1||sigtype==2)
                  {
                    Gtk1pt[gsize] = GenInfo->pt[GenInfo->da2[j]];
                    Gtk1eta[gsize] = GenInfo->eta[GenInfo->da2[j]];
                    Gtk1phi[gsize] = GenInfo->phi[GenInfo->da2[j]];
                  }
                else
                  {
                    Gtk1pt[gsize] = GenInfo->pt[GenInfo->da1[GenInfo->da2[j]]];
                    Gtk1eta[gsize] = GenInfo->eta[GenInfo->da1[GenInfo->da2[j]]];
                    Gtk1phi[gsize] = GenInfo->phi[GenInfo->da1[GenInfo->da2[j]]];
                    Gtk2pt[gsize] = GenInfo->pt[GenInfo->da2[GenInfo->da2[j]]];
                    Gtk2eta[gsize] = GenInfo->eta[GenInfo->da2[GenInfo->da2[j]]];
                    Gtk2phi[gsize] = GenInfo->phi[GenInfo->da2[GenInfo->da2[j]]];
                  }
              }
          }
        gsize++;
      }
    ntGen->Fill();
  }//}}}

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
  
  void fillTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, int j, int typesize, float track_mass1, float track_mass2, bool REAL, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, MuonInfoBranches *MuonInfo, TrackInfoBranches *TrackInfo, BInfoBranches *BInfo, GenInfoBranches *GenInfo, int &bestindex)
  {//{{{
    //Event Info
    Bsize = typesize+1;
    bP->SetPtEtaPhi(BInfo->pt[j],BInfo->eta[j]*0,BInfo->phi[j]);
    bVtx->SetXYZ(BInfo->vtxX[j]-EvtInfo->PVx,
                 BInfo->vtxY[j]-EvtInfo->PVy,
                 BInfo->vtxZ[j]*0-EvtInfo->PVz*0);
    b4P->SetPtEtaPhiM(BInfo->pt[j],BInfo->eta[j],BInfo->phi[j],BInfo->mass[j]);
    
    Bindex[typesize] = typesize;
    Btype[typesize] = BInfo->type[j];
    Bmass[typesize] = BInfo->mass[j];
    Bmass_unfitted[typesize] = BInfo->unfitted_mass[j];
    Bpt[typesize] = BInfo->pt[j];
    Beta[typesize] = BInfo->eta[j];
    Bphi[typesize] = BInfo->phi[j];
    By[typesize] = b4P->Rapidity();
    BvtxX[typesize] = BInfo->vtxX[j] - EvtInfo->PVx;
    BvtxY[typesize] = BInfo->vtxY[j] - EvtInfo->PVy;
    Bd0[typesize] = TMath::Sqrt((BInfo->vtxX[j]-EvtInfo->PVx)*(BInfo->vtxX[j]-EvtInfo->PVx)+(BInfo->vtxY[j]-EvtInfo->PVy)*(BInfo->vtxY[j]-EvtInfo->PVy));
    Bd0Err[typesize] = TMath::Sqrt(BInfo->vtxXErr[j]*BInfo->vtxXErr[j]+BInfo->vtxYErr[j]*BInfo->vtxYErr[j]);
    Bdxyz[typesize] = TMath::Sqrt((BInfo->vtxX[j]-EvtInfo->PVx)*(BInfo->vtxX[j]-EvtInfo->PVx)+(BInfo->vtxY[j]-EvtInfo->PVy)*(BInfo->vtxY[j]-EvtInfo->PVy)+(BInfo->vtxZ[j]-EvtInfo->PVz)*(BInfo->vtxZ[j]-EvtInfo->PVz));
    BdxyzErr[typesize] = TMath::Sqrt(BInfo->vtxXErr[j]*BInfo->vtxXErr[j]+BInfo->vtxYErr[j]*BInfo->vtxYErr[j]+BInfo->vtxZErr[j]*BInfo->vtxZErr[j]);
    Bchi2ndf[typesize] = BInfo->vtxchi2[j]/BInfo->vtxdof[j];
    Bchi2cl[typesize] = TMath::Prob(BInfo->vtxchi2[j],BInfo->vtxdof[j]);
    Bdtheta[typesize] = bP->Angle(*bVtx);
    Blxy[typesize] = ((BInfo->vtxX[j]-EvtInfo->PVx)*b4P->Px() + (BInfo->vtxY[j]-EvtInfo->PVy)*b4P->Py())/BInfo->pt[j];
    float r2lxyBS = (BInfo->vtxX[j]-EvtInfo->BSx+(BInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz) * (BInfo->vtxX[j]-EvtInfo->BSx+(BInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz)
      + (BInfo->vtxY[j]-EvtInfo->BSy+(BInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz) * (BInfo->vtxY[j]-EvtInfo->BSy+(BInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz);
    float xlxyBS = BInfo->vtxX[j]-EvtInfo->BSx + (BInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz;
    float ylxyBS = BInfo->vtxY[j]-EvtInfo->BSy + (BInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz;
    BlxyBS[typesize] = TMath::Sqrt(r2lxyBS);
    BlxyBSErr[typesize] = (1./r2lxyBS) * ((xlxyBS*xlxyBS)*BInfo->vtxXErr[j] + (2*xlxyBS*ylxyBS)*BInfo->vtxYXErr[j] + (ylxyBS*ylxyBS)*BInfo->vtxYErr[j]);
    Balpha[typesize] = BInfo->alpha[j];
    BsvpvDistance[typesize] = BInfo->svpvDistance[j];
    BsvpvDisErr[typesize] = BInfo->svpvDisErr[j];
    BsvpvDistance_2D[typesize] = BInfo->svpvDistance_2D[j];
    BsvpvDisErr_2D[typesize] = BInfo->svpvDisErr_2D[j];
    BMaxDoca[typesize] = BInfo->MaxDoca[j];
    
    b4P->SetPtEtaPhiM(MuonInfo->pt[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]],MuonInfo->eta[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]],MuonInfo->phi[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]],MUON_MASS);
    Bmu1pt[typesize] = MuonInfo->pt[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu1eta[typesize] = MuonInfo->eta[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu1phi[typesize] = MuonInfo->phi[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu1p[typesize] = b4P->P();
    Bmu1y[typesize] = b4P->Rapidity();
    Float_t mu1px,mu1py,mu1pz,mu1E;
    mu1px = b4P->Px();
    mu1py = b4P->Py();
    mu1pz = b4P->Pz();
    mu1E = b4P->E();

    b4P->SetPtEtaPhiM(MuonInfo->pt[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]],MuonInfo->eta[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]],MuonInfo->phi[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]],MUON_MASS);
    Bmu2pt[typesize] = MuonInfo->pt[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu2eta[typesize] = MuonInfo->eta[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu2phi[typesize] = MuonInfo->phi[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu2p[typesize] = b4P->P();
    Bmu2y[typesize] = b4P->Rapidity();
    Float_t mu2px,mu2py,mu2pz,mu2E;
    mu2px = b4P->Px();
    mu2py = b4P->Py();
    mu2pz = b4P->Pz();
    mu2E = b4P->E();

    Bmu1dzPV[typesize] = MuonInfo->dzPV[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2dzPV[typesize] = MuonInfo->dzPV[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1dxyPV[typesize] = MuonInfo->dxyPV[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2dxyPV[typesize] = MuonInfo->dxyPV[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1normchi2[typesize] = MuonInfo->normchi2[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2normchi2[typesize] = MuonInfo->normchi2[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1Chi2ndf[typesize] = MuonInfo->i_chi2[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]/MuonInfo->i_ndf[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2Chi2ndf[typesize] = MuonInfo->i_chi2[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]/MuonInfo->i_ndf[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1muqual[typesize] = MuonInfo->muqual[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2muqual[typesize] = MuonInfo->muqual[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1TrackerMuArbitrated[typesize] = MuonInfo->TrackerMuonArbitrated[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2TrackerMuArbitrated[typesize] = MuonInfo->TrackerMuonArbitrated[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1isTrackerMuon[typesize] = MuonInfo->isTrackerMuon[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2isTrackerMuon[typesize] = MuonInfo->isTrackerMuon[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1isGlobalMuon[typesize] = MuonInfo->isGlobalMuon[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2isGlobalMuon[typesize] = MuonInfo->isGlobalMuon[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1TMOneStationTight[typesize] = MuonInfo->TMOneStationTight[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2TMOneStationTight[typesize] = MuonInfo->TMOneStationTight[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1striphit[typesize] = MuonInfo->i_striphit[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu1pixelhit[typesize] = MuonInfo->i_pixelhit[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu1trackerhit[typesize] = Bmu1striphit[typesize]+Bmu1pixelhit[typesize];
    Bmu2striphit[typesize] = MuonInfo->i_striphit[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu2pixelhit[typesize] = MuonInfo->i_pixelhit[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu2trackerhit[typesize] = Bmu2striphit[typesize]+Bmu2pixelhit[typesize];
    Bmu1InPixelLayer[typesize] = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2InPixelLayer[typesize] = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1InStripLayer[typesize] = MuonInfo->i_nStripLayer[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2InStripLayer[typesize] = MuonInfo->i_nStripLayer[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1InTrackerLayer[typesize] = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]] + MuonInfo->i_nStripLayer[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2InTrackerLayer[typesize] = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]] + MuonInfo->i_nStripLayer[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1TrkQuality[typesize] = MuonInfo->innerTrackQuality[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2TrkQuality[typesize] = MuonInfo->innerTrackQuality[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];

    Bmu1TrgMatchFilterE[typesize] = MuonInfo->MuTrgMatchFilterTrgObjE[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]*MuonInfo->MuTrgMatchFilterSize+0];
    Bmu1TrgMatchFilterPt[typesize] = MuonInfo->MuTrgMatchFilterTrgObjPt[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]*MuonInfo->MuTrgMatchFilterSize+0];
    Bmu1TrgMatchFilterEta[typesize] = MuonInfo->MuTrgMatchFilterTrgObjEta[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]*MuonInfo->MuTrgMatchFilterSize+0];
    Bmu1TrgMatchFilterPhi[typesize] = MuonInfo->MuTrgMatchFilterTrgObjPhi[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]*MuonInfo->MuTrgMatchFilterSize+0];
    Bmu2TrgMatchFilterE[typesize] = MuonInfo->MuTrgMatchFilterTrgObjE[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]*MuonInfo->MuTrgMatchFilterSize+0];
    Bmu2TrgMatchFilterPt[typesize] = MuonInfo->MuTrgMatchFilterTrgObjPt[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]*MuonInfo->MuTrgMatchFilterSize+0];
    Bmu2TrgMatchFilterEta[typesize] = MuonInfo->MuTrgMatchFilterTrgObjEta[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]*MuonInfo->MuTrgMatchFilterSize+0];
    Bmu2TrgMatchFilterPhi[typesize] = MuonInfo->MuTrgMatchFilterTrgObjPhi[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]*MuonInfo->MuTrgMatchFilterSize+0];

    b4P->SetPxPyPzE(mu1px+mu2px,
                    mu1py+mu2py,
                    mu1pz+mu2pz,
                    mu1E+mu2E);
    Bmumumass[typesize] = b4P->Mag();
    Bmumueta[typesize] = b4P->Eta();
    Bmumuphi[typesize] = b4P->Phi();
    Bmumuy[typesize] = b4P->Rapidity();
    Bmumupt[typesize] = b4P->Pt();
    Bujmass[typesize] = BInfo->uj_mass[BInfo->rfuj_index[j]];
    BujvProb[typesize] = TMath::Prob(BInfo->uj_vtxchi2[BInfo->rfuj_index[j]],BInfo->uj_vtxdof[BInfo->rfuj_index[j]]);
    Bujpt[typesize] = BInfo->uj_pt[BInfo->rfuj_index[j]];
    Bujeta[typesize] = BInfo->uj_eta[BInfo->rfuj_index[j]];
    Bujphi[typesize] = BInfo->uj_phi[BInfo->rfuj_index[j]];
    b4P->SetPtEtaPhiM(BInfo->uj_pt[BInfo->rfuj_index[j]],
                      BInfo->uj_eta[BInfo->rfuj_index[j]],
                      BInfo->uj_phi[BInfo->rfuj_index[j]],
                      BInfo->uj_mass[BInfo->rfuj_index[j]]);
    Bujy[typesize] = b4P->Rapidity();
    Bujlxy[typesize] = ((BInfo->uj_vtxX[BInfo->rfuj_index[j]]-EvtInfo->PVx)*b4P->Px() + (BInfo->uj_vtxY[BInfo->rfuj_index[j]]-EvtInfo->PVy)*b4P->Py())/BInfo->uj_pt[BInfo->rfuj_index[j]];
    /*
      if(MuonInfo->muqual[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]&16) mu1TrackerMuArbitrated[typesize] = 1;
      else mu1TrackerMuArbitrated[typesize] = 0;
      if(MuonInfo->muqual[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]&4096) mu1TMOneStationTight[typesize] = 1;
      else mu1TMOneStationTight[typesize] = 0;
      if(MuonInfo->muqual[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]&16) mu2TrackerMuArbitrated[typesize] = 1;
      else mu2TrackerMuArbitrated[typesize] = 0;
      if(MuonInfo->muqual[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]&4096) mu2TMOneStationTight[typesize] = 1;
      else mu2TMOneStationTight[typesize] = 0;
    */
    
    float tk1px,tk1py,tk1pz,tk1E;
    float tk2px,tk2py,tk2pz,tk2E;
    Btrk1Idx[typesize] = BInfo->rftk1_index[j];
    Btrk2Idx[typesize] = BInfo->rftk2_index[j];
    if(BInfo->type[j]==1 || BInfo->type[j]==2)
      {
        b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk1_index[j]],TrackInfo->eta[BInfo->rftk1_index[j]],TrackInfo->phi[BInfo->rftk1_index[j]],track_mass1);
        Btrk1Pt[typesize] = TrackInfo->pt[BInfo->rftk1_index[j]];
        Btrk1Eta[typesize] = TrackInfo->eta[BInfo->rftk1_index[j]];
        Btrk1Phi[typesize] = TrackInfo->phi[BInfo->rftk1_index[j]];
        Btrk1PtErr[typesize] = TrackInfo->ptErr[BInfo->rftk1_index[j]];
        Btrk1EtaErr[typesize] = TrackInfo->etaErr[BInfo->rftk1_index[j]];
        Btrk1PhiErr[typesize] = TrackInfo->phiErr[BInfo->rftk1_index[j]];
        Btrk1Y[typesize] = b4P->Rapidity();
        Btrk1Dz[typesize] = TrackInfo->dzPV[BInfo->rftk1_index[j]];
        Btrk1Dxy[typesize] = TrackInfo->dxyPV[BInfo->rftk1_index[j]];
        Btrk1D0[typesize] = TrackInfo->d0[BInfo->rftk1_index[j]];
        Btrk1D0Err[typesize] = TrackInfo->d0error[BInfo->rftk1_index[j]];
        Btrk1PixelHit[typesize] = TrackInfo->pixelhit[BInfo->rftk1_index[j]];
        Btrk1StripHit[typesize] = TrackInfo->striphit[BInfo->rftk1_index[j]];
        Btrk1nPixelLayer[typesize] = TrackInfo->nPixelLayer[BInfo->rftk1_index[j]];
        Btrk1nStripLayer[typesize] = TrackInfo->nStripLayer[BInfo->rftk1_index[j]];
        Btrk1Chi2ndf[typesize] = TrackInfo->chi2[BInfo->rftk1_index[j]]/TrackInfo->ndf[BInfo->rftk1_index[j]];
        Btrk1MVAVal[typesize] = TrackInfo->trkMVAVal[BInfo->rftk1_index[j]];
        Btrk1Algo[typesize] = TrackInfo->trkAlgo[BInfo->rftk1_index[j]];
        Btrk1highPurity[typesize] = TrackInfo->highPurity[BInfo->rftk1_index[j]];
        Btrk1Quality[typesize] = TrackInfo->trackQuality[BInfo->rftk1_index[j]];
        Btrk2Pt[typesize] = -1;
        Btrk2Eta[typesize] = -20;
        Btrk2Phi[typesize] = -20;
        Btrk2PtErr[typesize] = 0;
        Btrk2EtaErr[typesize] = 0;
        Btrk2PhiErr[typesize] = 0;
        Btrk2Y[typesize] = -1;
        Btrk2Dz[typesize] = -1;
        Btrk2Dxy[typesize] = -1;
        Btrk2D0[typesize] = -1;
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
        //get best chi2 index
        Bisbestchi2[typesize] = -1;
        if(fabs(By[typesize])<2.4){// someselection
          if(bestindex == -1) {
            bestindex = typesize;
            Bisbestchi2[typesize] = 1;
          }
          else if(Bchi2cl[typesize] > Bchi2cl[bestindex]){
            Bisbestchi2[bestindex] = -1;
            bestindex = typesize;
            Bisbestchi2[typesize] = 1;
          }
        }
      }  
    else if(BInfo->type[j]==5)
      {
        b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk2_index[j]],TrackInfo->eta[BInfo->rftk2_index[j]],TrackInfo->phi[BInfo->rftk2_index[j]],track_mass1);
        Btrk1Pt[typesize] = TrackInfo->pt[BInfo->rftk2_index[j]];
        Btrk1Eta[typesize] = TrackInfo->eta[BInfo->rftk2_index[j]];
        Btrk1Phi[typesize] = TrackInfo->phi[BInfo->rftk2_index[j]];
        Btrk1PtErr[typesize] = TrackInfo->ptErr[BInfo->rftk2_index[j]];
        Btrk1EtaErr[typesize] = TrackInfo->etaErr[BInfo->rftk2_index[j]];
        Btrk1PhiErr[typesize] = TrackInfo->phiErr[BInfo->rftk2_index[j]];
        Btrk1Y[typesize] = b4P->Rapidity();
        Btrk1Dz[typesize] = TrackInfo->dzPV[BInfo->rftk2_index[j]];
        Btrk1Dxy[typesize] = TrackInfo->dxyPV[BInfo->rftk2_index[j]];
        Btrk1D0[typesize] = TrackInfo->d0[BInfo->rftk2_index[j]];
        Btrk1D0Err[typesize] = TrackInfo->d0error[BInfo->rftk2_index[j]];
        Btrk1PixelHit[typesize] = TrackInfo->pixelhit[BInfo->rftk2_index[j]];
        Btrk1StripHit[typesize] = TrackInfo->striphit[BInfo->rftk2_index[j]];
        Btrk1nPixelLayer[typesize] = TrackInfo->nPixelLayer[BInfo->rftk2_index[j]];
        Btrk1nStripLayer[typesize] = TrackInfo->nStripLayer[BInfo->rftk2_index[j]];
        Btrk1Chi2ndf[typesize] = TrackInfo->chi2[BInfo->rftk2_index[j]]/TrackInfo->ndf[BInfo->rftk2_index[j]];
        Btrk1MVAVal[typesize] = TrackInfo->trkMVAVal[BInfo->rftk2_index[j]];
        Btrk1Algo[typesize] = TrackInfo->trkAlgo[BInfo->rftk2_index[j]];
        Btrk1highPurity[typesize] = TrackInfo->highPurity[BInfo->rftk2_index[j]];
        Btrk1Quality[typesize] = TrackInfo->trackQuality[BInfo->rftk2_index[j]];
        tk1px = b4P->Px();
        tk1py = b4P->Py();
        tk1pz = b4P->Pz();
        tk1E = b4P->E();
        b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk1_index[j]],TrackInfo->eta[BInfo->rftk1_index[j]],TrackInfo->phi[BInfo->rftk1_index[j]],track_mass2);
        Btrk2Pt[typesize] = TrackInfo->pt[BInfo->rftk1_index[j]];
        Btrk2Eta[typesize] = TrackInfo->eta[BInfo->rftk1_index[j]];
        Btrk2Phi[typesize] = TrackInfo->phi[BInfo->rftk1_index[j]];
        Btrk2PtErr[typesize] = TrackInfo->ptErr[BInfo->rftk1_index[j]];
        Btrk2EtaErr[typesize] = TrackInfo->etaErr[BInfo->rftk1_index[j]];
        Btrk2PhiErr[typesize] = TrackInfo->phiErr[BInfo->rftk1_index[j]];
        Btrk2Y[typesize] = b4P->Rapidity();
        Btrk2Dz[typesize] = TrackInfo->dzPV[BInfo->rftk1_index[j]];
        Btrk2Dxy[typesize] = TrackInfo->dxyPV[BInfo->rftk1_index[j]];
        Btrk2D0[typesize] = TrackInfo->d0error[BInfo->rftk1_index[j]];
        Btrk2D0Err[typesize] = TrackInfo->d0[BInfo->rftk1_index[j]];
        Btrk2PixelHit[typesize] = TrackInfo->pixelhit[BInfo->rftk1_index[j]];
        Btrk2StripHit[typesize] = TrackInfo->striphit[BInfo->rftk1_index[j]];
        Btrk2nPixelLayer[typesize] = TrackInfo->nPixelLayer[BInfo->rftk1_index[j]];
        Btrk2nStripLayer[typesize] = TrackInfo->nStripLayer[BInfo->rftk1_index[j]];
        Btrk2Chi2ndf[typesize] = TrackInfo->chi2[BInfo->rftk1_index[j]]/TrackInfo->ndf[BInfo->rftk1_index[j]];
        Btrk2MVAVal[typesize] = TrackInfo->trkMVAVal[BInfo->rftk1_index[j]];
        Btrk2Algo[typesize] = TrackInfo->trkAlgo[BInfo->rftk1_index[j]];
        Btrk2highPurity[typesize] = TrackInfo->highPurity[BInfo->rftk1_index[j]];
        Btrk2Quality[typesize] = TrackInfo->trackQuality[BInfo->rftk1_index[j]];
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
        BtktkvProb[typesize] = TMath::Prob(BInfo->tktk_vtxchi2[j],BInfo->tktk_vtxdof[j]);

        Bdoubletmass[typesize] = BInfo->tktk_mass[j];
        Bdoubletpt[typesize] = BInfo->tktk_pt[j];
        Bdoubleteta[typesize] = BInfo->tktk_eta[j];
        Bdoubletphi[typesize] = BInfo->tktk_phi[j];
        b4P->SetPtEtaPhiM(BInfo->tktk_pt[j],BInfo->tktk_eta[j],BInfo->tktk_phi[j],BInfo->tktk_mass[j]);
        Bdoublety[typesize] = b4P->Rapidity();
        
        b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk1_index[j]],TrackInfo->eta[BInfo->rftk1_index[j]],TrackInfo->phi[BInfo->rftk1_index[j]],KAON_MASS);
        float tk1EK = b4P->E();
        b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk2_index[j]],TrackInfo->eta[BInfo->rftk2_index[j]],TrackInfo->phi[BInfo->rftk2_index[j]],KAON_MASS);
        float tk2EK = b4P->E();
        b4P->SetPxPyPzE(tk1px+tk2px,
                        tk1py+tk2py,
                        tk1pz+tk2pz,
                        tk1EK+tk2EK);
        BtktkmassKK[typesize] = b4P->Mag();
      }
    else
      {
        b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk1_index[j]],TrackInfo->eta[BInfo->rftk1_index[j]],TrackInfo->phi[BInfo->rftk1_index[j]],track_mass1);
        Btrk1Pt[typesize] = TrackInfo->pt[BInfo->rftk1_index[j]];
        Btrk1Eta[typesize] = TrackInfo->eta[BInfo->rftk1_index[j]];
        Btrk1Phi[typesize] = TrackInfo->phi[BInfo->rftk1_index[j]];
        Btrk1PtErr[typesize] = TrackInfo->ptErr[BInfo->rftk1_index[j]];
        Btrk1EtaErr[typesize] = TrackInfo->etaErr[BInfo->rftk1_index[j]];
        Btrk1PhiErr[typesize] = TrackInfo->phiErr[BInfo->rftk1_index[j]];
        Btrk1Y[typesize] = b4P->Rapidity();
        Btrk1Dz[typesize] = TrackInfo->dzPV[BInfo->rftk1_index[j]];
        Btrk1Dxy[typesize] = TrackInfo->dxyPV[BInfo->rftk1_index[j]];
        Btrk1D0[typesize] = TrackInfo->d0[BInfo->rftk1_index[j]];
        Btrk1D0Err[typesize] = TrackInfo->d0error[BInfo->rftk1_index[j]];
        Btrk1PixelHit[typesize] = TrackInfo->pixelhit[BInfo->rftk1_index[j]];
        Btrk1StripHit[typesize] = TrackInfo->striphit[BInfo->rftk1_index[j]];
        Btrk1nPixelLayer[typesize] = TrackInfo->nPixelLayer[BInfo->rftk1_index[j]];
        Btrk1nStripLayer[typesize] = TrackInfo->nStripLayer[BInfo->rftk1_index[j]];
        Btrk1Chi2ndf[typesize] = TrackInfo->chi2[BInfo->rftk1_index[j]]/TrackInfo->ndf[BInfo->rftk1_index[j]];
        Btrk1MVAVal[typesize] = TrackInfo->trkMVAVal[BInfo->rftk1_index[j]];
        Btrk1Algo[typesize] = TrackInfo->trkAlgo[BInfo->rftk1_index[j]];
        Btrk1highPurity[typesize] = TrackInfo->highPurity[BInfo->rftk1_index[j]];
        Btrk1Quality[typesize] = TrackInfo->trackQuality[BInfo->rftk1_index[j]];
        tk1px = b4P->Px();
        tk1py = b4P->Py();
        tk1pz = b4P->Pz();
        tk1E = b4P->E();
        b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk2_index[j]],TrackInfo->eta[BInfo->rftk2_index[j]],TrackInfo->phi[BInfo->rftk2_index[j]],track_mass2);
        Btrk2Pt[typesize] = TrackInfo->pt[BInfo->rftk2_index[j]];
        Btrk2Eta[typesize] = TrackInfo->eta[BInfo->rftk2_index[j]];
        Btrk2Phi[typesize] = TrackInfo->phi[BInfo->rftk2_index[j]];
        Btrk2PtErr[typesize] = TrackInfo->ptErr[BInfo->rftk2_index[j]];
        Btrk2EtaErr[typesize] = TrackInfo->etaErr[BInfo->rftk2_index[j]];
        Btrk2PhiErr[typesize] = TrackInfo->phiErr[BInfo->rftk2_index[j]];
        Btrk2Y[typesize] = b4P->Rapidity();
        Btrk2Dz[typesize] = TrackInfo->dzPV[BInfo->rftk2_index[j]];
        Btrk2Dxy[typesize] = TrackInfo->dxyPV[BInfo->rftk2_index[j]];
        Btrk2D0[typesize] = TrackInfo->d0[BInfo->rftk2_index[j]];
        Btrk2D0Err[typesize] = TrackInfo->d0error[BInfo->rftk2_index[j]];
        Btrk2PixelHit[typesize] = TrackInfo->pixelhit[BInfo->rftk2_index[j]];
        Btrk2StripHit[typesize] = TrackInfo->striphit[BInfo->rftk2_index[j]];
        Btrk2nPixelLayer[typesize] = TrackInfo->nPixelLayer[BInfo->rftk2_index[j]];
        Btrk2nStripLayer[typesize] = TrackInfo->nStripLayer[BInfo->rftk2_index[j]];
        Btrk2Chi2ndf[typesize] = TrackInfo->chi2[BInfo->rftk2_index[j]]/TrackInfo->ndf[BInfo->rftk2_index[j]];
        Btrk2MVAVal[typesize] = TrackInfo->trkMVAVal[BInfo->rftk2_index[j]];
        Btrk2Algo[typesize] = TrackInfo->trkAlgo[BInfo->rftk2_index[j]];
        Btrk2highPurity[typesize] = TrackInfo->highPurity[BInfo->rftk2_index[j]];
        Btrk2Quality[typesize] = TrackInfo->trackQuality[BInfo->rftk2_index[j]];
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
        BtktkvProb[typesize] = TMath::Prob(BInfo->tktk_vtxchi2[j],BInfo->tktk_vtxdof[j]);
        Bdoubletmass[typesize] = BInfo->tktk_mass[j];
        Bdoubletpt[typesize] = BInfo->tktk_pt[j];
        Bdoubleteta[typesize] = BInfo->tktk_eta[j];
        Bdoubletphi[typesize] = BInfo->tktk_phi[j];
        b4P->SetPtEtaPhiM(BInfo->tktk_pt[j],BInfo->tktk_eta[j],BInfo->tktk_phi[j],BInfo->tktk_mass[j]);
        Bdoublety[typesize] = b4P->Rapidity();
        
        b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk1_index[j]],TrackInfo->eta[BInfo->rftk1_index[j]],TrackInfo->phi[BInfo->rftk1_index[j]],KAON_MASS);
        float tk1EK = b4P->E();
        b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk2_index[j]],TrackInfo->eta[BInfo->rftk2_index[j]],TrackInfo->phi[BInfo->rftk2_index[j]],KAON_MASS);
        float tk2EK = b4P->E();
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
        //int ujGenIdxMu1=-1;
        //int ujGenIdxMu2=-1;
        
        float BId,MId,tk1Id,tk2Id;
        //tk1:positive, tk2:negtive
        if(BInfo->type[j]==1)
          {
            BId = 521;//B+-
            MId = -1;
            tk1Id = 321;//K+-
            tk2Id = -1;
          }
        if(BInfo->type[j]==2)
          {
            BId = 521;//B+-
            MId = -1;
            tk1Id = 211;//pi+-
            tk2Id = -1;
          }
        if(BInfo->type[j]==3)
          {
            BId = 511;//B0
            MId = 310;//Ks
            tk1Id = 211;//pi+
            tk2Id = 211;//pi-
          }
        if(BInfo->type[j]==4)
          {
            BId = 511;//B0
            MId = 313;//K*0
            tk1Id = 321;//K+
            tk2Id = 211;//pi-
          }
        if(BInfo->type[j]==5)
          {
            BId = 511;//B0
            MId = 313;//K*0
            tk1Id = 211;//pi+
            tk2Id = 321;//K-
          }
        if(BInfo->type[j]==6)
          {
            BId = 531;//Bs
            MId = 333;//phi
            tk1Id = 321;//K+
            tk2Id = 321;//K-
          }

        int twoTks,kStar,flagkstar=0;
        if(BInfo->type[j]==1 || BInfo->type[j]==2) twoTks=0;
        else twoTks=1;
        if(BInfo->type[j]==4 || BInfo->type[j]==5) kStar=1;
        else kStar=0;
        //int nonprompt=0,prompt=0;
        
        // tk1
        if(TrackInfo->geninfo_index[BInfo->rftk1_index[j]]>-1)
          {
            int level =0;
            if(abs(GenInfo->pdgId[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]])==tk1Id)
              {
                level = 1;
                if(GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]]>-1)
                  {
                    if(!twoTks)//one trk channel
                      {
                        mGenIdxTk1=0;
                        if(abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]]])==BId)
                          {
                            level = 3;
                            bGenIdxTk1=GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]];
                          }		  
                      }
                    else//two trk channel
                      {
                        if(abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]]])==MId)
                          {
                            level = 2;
                            if(GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]]]>-1)
                              {
                                if(abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]]]])==BId)
                                  {
                                    level = 3;
                                    bGenIdxTk1=GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]]];
                                  }
                              }
                            mGenIdxTk1=GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]];
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
            if(TrackInfo->geninfo_index[BInfo->rftk2_index[j]]>-1)
              {
                int level =0;
                if(abs(GenInfo->pdgId[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]])==tk2Id)
                  {
                    level = 1;
                    if(GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]]>-1)
                      {
                        if(abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]]])==MId)
                          {
                            level = 2;
                            if(GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]]]>-1)
                              {
                                if(abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]]]])==BId)
                                  {
                                    level = 3;
                                    bGenIdxTk2 = GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]]];
                                  }
                              }
                            mGenIdxTk2 = GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]];
                          }
                      }
                  }
                Bgen[typesize]+=(level*10);
              }
          }
        
        //mu1
        //cout<<MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]<<endl;
        if(MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]>-1)
          {
            int level =0;
            if(abs(GenInfo->pdgId[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]])==13)
              {
                level=1;
                if(GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]]>-1)
                  {
                    if(GenInfo->pdgId[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]]]==443)
                      {
                        //ujGenIdxMu1 = GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]];
                        level=2;
                        if(GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]]]>-1)
                          {
                            if(abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]]]])==BId)
                              {
                                //nonprompt=1;
                                level = 3;
                                bGenIdxMu1=GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]]];
                                flagkstar++;///////////////////////////////////////////////=1
                              }
                          }
                        else 
                          {
                            //prompt=1;
                          }
                      } 
                  }
              }
            Bgen[typesize]+=(level*100);
          }
        
        //mu2
        if(MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]>-1)
          {  
            int level =0;
            if(abs(GenInfo->pdgId[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]])==13)
              {
                level = 1;
                if(GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]]>-1)
                  {
                    if(GenInfo->pdgId[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]]]==443)
                      {
                        //ujGenIdxMu2 = GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]];
                        level = 2;
                        if(GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]]]>-1)
                          {
                            if(abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]]]])==BId)
                              {
                                level = 3;
                                bGenIdxMu2=GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]]];
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
            if(TrackInfo->geninfo_index[BInfo->rftk1_index[j]]>-1)
              {
                if(abs(GenInfo->pdgId[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]])==tk2Id)
                  {
                    if(GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]]>-1)
                      {
                        if(abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]]])==MId)
                          {
                            if(GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]]]>-1)
                              {
                                if(abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]]]])==BId)
                                  {
                                    flagkstar++;//////////////////////////////////////////////=3
                                    bGenIdxTk1=GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]]];
                                  }
                              }
                            mGenIdxTk1=GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk1_index[j]]];
                          }
                      }
                  }
              }
            
            //tk2
            if(TrackInfo->geninfo_index[BInfo->rftk2_index[j]]>-1)
              {
                if(abs(GenInfo->pdgId[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]])==tk1Id)
                  {
                    if(GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]]>-1)
                      {
                        if(abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]]])==MId)
                          {
                            if(GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]]]>-1)
                              {
                                if(abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]]]])==BId)
                                  {
                                    flagkstar++;////////////////////////////////////////////////////=4
                                    bGenIdxTk2 = GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]]];
                                  }
                              }
                            mGenIdxTk2 = GenInfo->mo1[TrackInfo->geninfo_index[BInfo->rftk2_index[j]]];
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
            Bgenpt[typesize] = GenInfo->pt[tgenIndex];
            Bgeneta[typesize] = GenInfo->eta[tgenIndex];
            Bgenphi[typesize] = GenInfo->phi[tgenIndex];
            b4P->SetXYZM(GenInfo->pt[tgenIndex]*cos(GenInfo->phi[tgenIndex]),
                         GenInfo->pt[tgenIndex]*sin(GenInfo->phi[tgenIndex]),
                         GenInfo->pt[tgenIndex]*sinh(GenInfo->eta[tgenIndex]),
                         GenInfo->mass[tgenIndex]);
            Bgeny[typesize] = b4P->Rapidity();
          }
      }
  }//}}}

  void fillJpsiTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, int j, int typesize, bool REAL, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, MuonInfoBranches *MuonInfo, TrackInfoBranches *TrackInfo, BInfoBranches *BInfo, GenInfoBranches *GenInfo)
  {
    Jsize = typesize+1;
    bP->SetPtEtaPhi(BInfo->uj_pt[j],BInfo->uj_eta[j]*0,BInfo->uj_phi[j]);
    bVtx->SetXYZ(BInfo->uj_vtxX[j]-EvtInfo->PVx,
                 BInfo->uj_vtxY[j]-EvtInfo->PVy,
                 BInfo->uj_vtxZ[j]*0-EvtInfo->PVz*0);
    b4P->SetPtEtaPhiM(BInfo->uj_pt[j],BInfo->uj_eta[j],BInfo->uj_phi[j],BInfo->uj_mass[j]);
    
    Jindex[typesize] = typesize;
    Jmass[typesize] = BInfo->uj_mass[j];
    Jpt[typesize] = BInfo->uj_pt[j];
    Jeta[typesize] = BInfo->uj_eta[j];
    Jphi[typesize] = BInfo->uj_phi[j];
    Jy[typesize] = b4P->Rapidity();
    JvtxX[typesize] = BInfo->uj_vtxX[j] - EvtInfo->PVx;
    JvtxY[typesize] = BInfo->uj_vtxY[j] - EvtInfo->PVy;
    Jd0[typesize] = TMath::Sqrt((BInfo->uj_vtxX[j]-EvtInfo->PVx)*(BInfo->uj_vtxX[j]-EvtInfo->PVx)+(BInfo->uj_vtxY[j]-EvtInfo->PVy)*(BInfo->uj_vtxY[j]-EvtInfo->PVy));
    Jdxyz[typesize] = TMath::Sqrt((BInfo->uj_vtxX[j]-EvtInfo->PVx)*(BInfo->uj_vtxX[j]-EvtInfo->PVx)+(BInfo->uj_vtxY[j]-EvtInfo->PVy)*(BInfo->uj_vtxY[j]-EvtInfo->PVy)+(BInfo->uj_vtxZ[j]-EvtInfo->PVz)*(BInfo->uj_vtxZ[j]-EvtInfo->PVz));
    Jchi2ndf[typesize] = BInfo->uj_vtxchi2[j]/BInfo->uj_vtxdof[j];
    Jchi2cl[typesize] = TMath::Prob(BInfo->uj_vtxchi2[j],BInfo->uj_vtxdof[j]);
    Jdtheta[typesize] = bP->Angle(*bVtx);
    Jlxy[typesize] = ((BInfo->uj_vtxX[j]-EvtInfo->PVx)*b4P->Px() + (BInfo->uj_vtxY[j]-EvtInfo->PVy)*b4P->Py())/BInfo->uj_pt[j];
    float r2lxyBS = (BInfo->uj_vtxX[j]-EvtInfo->BSx+(BInfo->uj_vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz) * (BInfo->uj_vtxX[j]-EvtInfo->BSx+(BInfo->uj_vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz)
      + (BInfo->uj_vtxY[j]-EvtInfo->BSy+(BInfo->uj_vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz) * (BInfo->uj_vtxY[j]-EvtInfo->BSy+(BInfo->uj_vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz);
    JlxyBS[typesize] = TMath::Sqrt(r2lxyBS);
    
    b4P->SetPtEtaPhiM(MuonInfo->pt[BInfo->uj_rfmu1_index[j]],MuonInfo->eta[BInfo->uj_rfmu1_index[j]],MuonInfo->phi[BInfo->uj_rfmu1_index[j]],MUON_MASS);
    Jmu1pt[typesize] = MuonInfo->pt[BInfo->uj_rfmu1_index[j]];
    Jmu1eta[typesize] = MuonInfo->eta[BInfo->uj_rfmu1_index[j]];
    Jmu1phi[typesize] = MuonInfo->phi[BInfo->uj_rfmu1_index[j]];
    Jmu1p[typesize] = b4P->P();
    Jmu1y[typesize] = b4P->Rapidity();
    Float_t mu1px,mu1py,mu1pz,mu1E;
    mu1px = b4P->Px();
    mu1py = b4P->Py();
    mu1pz = b4P->Pz();
    mu1E = b4P->E();
    
    b4P->SetPtEtaPhiM(MuonInfo->pt[BInfo->uj_rfmu2_index[j]],MuonInfo->eta[BInfo->uj_rfmu2_index[j]],MuonInfo->phi[BInfo->uj_rfmu2_index[j]],MUON_MASS);
    Jmu2pt[typesize] = MuonInfo->pt[BInfo->uj_rfmu2_index[j]];
    Jmu2eta[typesize] = MuonInfo->eta[BInfo->uj_rfmu2_index[j]];
    Jmu2phi[typesize] = MuonInfo->phi[BInfo->uj_rfmu2_index[j]];
    Jmu2p[typesize] = b4P->P();
    Jmu2y[typesize] = b4P->Rapidity();
    Float_t mu2px,mu2py,mu2pz,mu2E;
    mu2px = b4P->Px();
    mu2py = b4P->Py();
    mu2pz = b4P->Pz();
    mu2E = b4P->E();

    Jmu1dzPV[typesize] = MuonInfo->dzPV[BInfo->uj_rfmu1_index[j]];
    Jmu2dzPV[typesize] = MuonInfo->dzPV[BInfo->uj_rfmu2_index[j]];
    Jmu1dxyPV[typesize] = MuonInfo->dxyPV[BInfo->uj_rfmu1_index[j]];
    Jmu2dxyPV[typesize] = MuonInfo->dxyPV[BInfo->uj_rfmu2_index[j]];
    Jmu1normchi2[typesize] = MuonInfo->normchi2[BInfo->uj_rfmu1_index[j]];
    Jmu2normchi2[typesize] = MuonInfo->normchi2[BInfo->uj_rfmu2_index[j]];
    Jmu1Chi2ndf[typesize] = MuonInfo->i_chi2[BInfo->uj_rfmu1_index[j]]/MuonInfo->i_ndf[BInfo->uj_rfmu1_index[j]];
    Jmu2Chi2ndf[typesize] = MuonInfo->i_chi2[BInfo->uj_rfmu2_index[j]]/MuonInfo->i_ndf[BInfo->uj_rfmu2_index[j]];
    Jmu1muqual[typesize] = MuonInfo->muqual[BInfo->uj_rfmu1_index[j]];
    Jmu2muqual[typesize] = MuonInfo->muqual[BInfo->uj_rfmu2_index[j]];
    Jmu1TrackerMuArbitrated[typesize] = MuonInfo->TrackerMuonArbitrated[BInfo->uj_rfmu1_index[j]];
    Jmu2TrackerMuArbitrated[typesize] = MuonInfo->TrackerMuonArbitrated[BInfo->uj_rfmu2_index[j]];
    Jmu1isTrackerMuon[typesize] = MuonInfo->isTrackerMuon[BInfo->uj_rfmu1_index[j]];
    Jmu2isTrackerMuon[typesize] = MuonInfo->isTrackerMuon[BInfo->uj_rfmu2_index[j]];
    Jmu1isGlobalMuon[typesize] = MuonInfo->isGlobalMuon[BInfo->uj_rfmu1_index[j]];
    Jmu2isGlobalMuon[typesize] = MuonInfo->isGlobalMuon[BInfo->uj_rfmu2_index[j]];
    Jmu1TMOneStationTight[typesize] = MuonInfo->TMOneStationTight[BInfo->uj_rfmu1_index[j]];
    Jmu2TMOneStationTight[typesize] = MuonInfo->TMOneStationTight[BInfo->uj_rfmu2_index[j]];
    Jmu1striphit[typesize] = MuonInfo->i_striphit[BInfo->uj_rfmu1_index[j]];
    Jmu1pixelhit[typesize] = MuonInfo->i_pixelhit[BInfo->uj_rfmu1_index[j]];
    Jmu1trackerhit[typesize] = Jmu1striphit[typesize]+Bmu1pixelhit[typesize];
    Jmu2striphit[typesize] = MuonInfo->i_striphit[BInfo->uj_rfmu2_index[j]];
    Jmu2pixelhit[typesize] = MuonInfo->i_pixelhit[BInfo->uj_rfmu2_index[j]];
    Jmu2trackerhit[typesize] = Bmu2striphit[typesize]+Bmu2pixelhit[typesize];
    Jmu1InPixelLayer[typesize] = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu1_index[j]];
    Jmu2InPixelLayer[typesize] = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu2_index[j]];
    Jmu1InStripLayer[typesize] = MuonInfo->i_nStripLayer[BInfo->uj_rfmu1_index[j]];
    Jmu2InStripLayer[typesize] = MuonInfo->i_nStripLayer[BInfo->uj_rfmu2_index[j]];
    Jmu1InTrackerLayer[typesize] = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu1_index[j]] + MuonInfo->i_nStripLayer[BInfo->uj_rfmu1_index[j]];
    Jmu2InTrackerLayer[typesize] = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu2_index[j]] + MuonInfo->i_nStripLayer[BInfo->uj_rfmu2_index[j]];
    Jmu1TrkQuality[typesize] = MuonInfo->innerTrackQuality[BInfo->uj_rfmu1_index[j]];
    Jmu2TrkQuality[typesize] = MuonInfo->innerTrackQuality[BInfo->uj_rfmu2_index[j]];

    Jmu1TrgMatchFilterE[typesize] = MuonInfo->MuTrgMatchFilterTrgObjE[BInfo->uj_rfmu1_index[j]*MuonInfo->MuTrgMatchFilterSize+0];
    Jmu1TrgMatchFilterPt[typesize] = MuonInfo->MuTrgMatchFilterTrgObjPt[BInfo->uj_rfmu1_index[j]*MuonInfo->MuTrgMatchFilterSize+0];
    Jmu1TrgMatchFilterEta[typesize] = MuonInfo->MuTrgMatchFilterTrgObjEta[BInfo->uj_rfmu1_index[j]*MuonInfo->MuTrgMatchFilterSize+0];
    Jmu1TrgMatchFilterPhi[typesize] = MuonInfo->MuTrgMatchFilterTrgObjPhi[BInfo->uj_rfmu1_index[j]*MuonInfo->MuTrgMatchFilterSize+0];
    Jmu2TrgMatchFilterE[typesize] = MuonInfo->MuTrgMatchFilterTrgObjE[BInfo->uj_rfmu2_index[j]*MuonInfo->MuTrgMatchFilterSize+0];
    Jmu2TrgMatchFilterPt[typesize] = MuonInfo->MuTrgMatchFilterTrgObjPt[BInfo->uj_rfmu2_index[j]*MuonInfo->MuTrgMatchFilterSize+0];
    Jmu2TrgMatchFilterEta[typesize] = MuonInfo->MuTrgMatchFilterTrgObjEta[BInfo->uj_rfmu2_index[j]*MuonInfo->MuTrgMatchFilterSize+0];
    Jmu2TrgMatchFilterPhi[typesize] = MuonInfo->MuTrgMatchFilterTrgObjPhi[BInfo->uj_rfmu2_index[j]*MuonInfo->MuTrgMatchFilterSize+0];

    b4P->SetPxPyPzE(mu1px+mu2px,
                    mu1py+mu2py,
                    mu1pz+mu2pz,
                    mu1E+mu2E);
    Jmass_unfitted[typesize] = b4P->Mag();

    if(!REAL)
      {
        Jgen[typesize] = 0;
        Jgen[typesize]+=3300;
        JgenIndex[typesize] = -1;
        Jgenpt[typesize] = -1;
        Jgeneta[typesize] = -20;
        Jgenphi[typesize] = -20;
        Jgeny[typesize] = -1;
        int bGenIdxMu1=-1;
        int bGenIdxMu2=-1;
        int ujGenIdxMu1=-1;
        int ujGenIdxMu2=-1;

        if(MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[j]]>-1)
          {
            int level =0;
            if(abs(GenInfo->pdgId[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[j]]])==13)
              {
                level=1;
                if(GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[j]]]>-1)
                  {
                    if(GenInfo->pdgId[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[j]]]]==443)
                      {
                        ujGenIdxMu1 = GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[j]]];
                        level = 3;
                        if(GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[j]]]]>-1)
                          {
                            int bPDG = abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[j]]]]]);
                            bPDG /= 100;
                            bPDG = bPDG%10;
                            if(bPDG==5)
                              {
                                level = 4;
                                bGenIdxMu1 = GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu1_index[j]]]];
                              }
                          }
                      } 
                  }
              }
            Jgen[typesize]+=(level);
          }
        if(MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[j]]>-1)
          {  
            int level =0;
            if(abs(GenInfo->pdgId[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[j]]])==13)
              {
                level = 1;
                if(GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[j]]]>-1)
                  {
                    if(GenInfo->pdgId[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[j]]]]==443)
                      {
                        ujGenIdxMu2 = GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[j]]];
                        level = 3;
                        if(GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[j]]]]>-1)
                          {
                            int bPDG = abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[j]]]]]);
                            bPDG /= 100;
                            bPDG = bPDG%10;
                            if(bPDG==5)
                              {
                                level = 4;
                                bGenIdxMu2 = GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[BInfo->uj_rfmu2_index[j]]]];
                              }
                          }
                      }
                  }
              }
            Jgen[typesize]+=(level*10);
          }

        int level=0;
        if(ujGenIdxMu1!=-1 && ujGenIdxMu1==ujGenIdxMu2)
          {
            level = 2;
            JgenIndex[typesize] = ujGenIdxMu1;
            if(Jgen[typesize]==3344 && bGenIdxMu1!=bGenIdxMu2) level = 1;
          }
        Jgen[typesize]+=(level*10000);
        if(Jgen[typesize]==23333 || Jgen[typesize]==23344)
          {
            int tgenIndex = JgenIndex[typesize];
            Jgenpt[typesize] = GenInfo->pt[tgenIndex];
            Jgeneta[typesize] = GenInfo->eta[tgenIndex];
            Jgenphi[typesize] = GenInfo->phi[tgenIndex];
            b4P->SetXYZM(GenInfo->pt[tgenIndex]*cos(GenInfo->phi[tgenIndex]),
                         GenInfo->pt[tgenIndex]*sin(GenInfo->phi[tgenIndex]),
                         GenInfo->pt[tgenIndex]*sinh(GenInfo->eta[tgenIndex]),
                         GenInfo->mass[tgenIndex]);
            Jgeny[typesize] = b4P->Rapidity();
          }
      }//!REAL

  }

  
  bool signalGen(int Btype, int j, GenInfoBranches *GenInfo)
  {//{{{

    if(Btype==8 || Btype==9)
      {
        int flag=0;
        if(abs(GenInfo->pdgId[j])==JPSI_PDGID&&GenInfo->nDa[j]==2&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1)
          {
            if(abs(GenInfo->pdgId[GenInfo->da1[j]])==13&&abs(GenInfo->pdgId[GenInfo->da2[j]])==13)
              {
                int bPDG = 0;
                if(GenInfo->mo1[j]!=-1)
                  {
                    bPDG = abs(GenInfo->pdgId[GenInfo->mo1[j]]);
                    bPDG /= 100;
                    bPDG = bPDG%10;
                  }
                if(bPDG==5)
                  {
                    if(Btype==9) flag++;
                  }
                else
                  {
                    if(Btype==8) flag++;
                  }
              }
          }
        return flag;
      }
    else
      {
        float BId,MId,tk1Id,tk2Id;
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
        if (abs(GenInfo->pdgId[j])==BId&&GenInfo->nDa[j]==2&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1)
          {
            if (abs(GenInfo->pdgId[GenInfo->da1[j]])==443)//jpsi
              {
                if(GenInfo->da1[GenInfo->da1[j]]!=-1&&GenInfo->da2[GenInfo->da1[j]]!=-1)
                  {
                    if(abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==13&&abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==13)
                      {
                        if(!twoTks)
                          {
                            if(abs(GenInfo->pdgId[GenInfo->da2[j]])==tk1Id) flag++;
                          }
                        else
                          {
                            if (abs(GenInfo->pdgId[GenInfo->da2[j]])==MId) 
                              {
                                if(GenInfo->da1[GenInfo->da2[j]]!=-1 && GenInfo->da2[GenInfo->da2[j]]!=-1)
                                  {
                                    if(GenInfo->pdgId[GenInfo->da1[GenInfo->da2[j]]]==tk1Id && GenInfo->pdgId[GenInfo->da2[GenInfo->da2[j]]]==tk2Id) flag++;
                                  }
                              }
                          }
                      }
                  }
              }
          }
        return flag;
      }

  }//}}}
  
};//}}}

#endif
