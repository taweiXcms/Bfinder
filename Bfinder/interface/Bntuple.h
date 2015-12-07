// vim:set ts=4 sw=4 fdm=marker et:
#ifndef _BNTUPLE_H_
#define _BNTUPLE_H_
#include "format.h"

class BntupleBranches
{//{{{
 public:
  double tk1mass[7] = {KAON_MASS, PION_MASS, PION_MASS,   KAON_MASS,  KAON_MASS,  KAON_MASS, PION_MASS};
  double tk2mass[7] = {0,         0,         PION_MASS,   PION_MASS,  PION_MASS,  KAON_MASS, PION_MASS};
  double midmass[7] = {0,         0,         KSHORT_MASS, KSTAR_MASS, KSTAR_MASS, PHI_MASS,  0};

  //EvtInfo
  int      RunNo;
  int      EvtNo;
  int      LumiNo;
  int      Bsize;
  double   PVx;
  double   PVy;
  double   PVz;
  double   PVxE;
  double   PVyE;
  double   PVzE;
  double   PVnchi2;
  double   PVchi2;
  double   BSx;
  double   BSy;
  double   BSz;
  double   BSxErr;
  double   BSyErr;
  double   BSzErr;
  double   BSdxdz;
  double   BSdydz;
  double   BSdxdzErr;
  double   BSdydzErr;
  double   BSWidthX;
  double   BSWidthXErr;
  double   BSWidthY;
  double   BSWidthYErr;
  
  //DInfo
  int      Bindex[MAX_XB];
  int      Btype[MAX_XB];
  double   Bmass[MAX_XB];
  double   Bpt[MAX_XB];
  double   Beta[MAX_XB];
  double   Bphi[MAX_XB];
  double   By[MAX_XB];
  double   BvtxX[MAX_XB];
  double   BvtxY[MAX_XB];
  double   Bd0[MAX_XB];
  double   Bd0Err[MAX_XB];
  double   Bdxyz[MAX_XB];
  double   BdxyzErr[MAX_XB];
  double   Bchi2ndf[MAX_XB];
  double   Bchi2cl[MAX_XB];
  double   Bdtheta[MAX_XB];
  double   Blxy[MAX_XB];
  double   BlxyBS[MAX_XB];
  double   BlxyBSErr[MAX_XB];
  bool     Bmaxpt[MAX_XB];
  bool     Bmaxprob[MAX_XB];
  bool     Bbesttktkmass[MAX_XB];
  bool     BmaxptMatched[MAX_XB];
  bool     BmaxprobMatched[MAX_XB];
  bool     BbesttktkmassMatched[MAX_XB];
  
  //BInfo.muonInfo
  double   Bmu1pt[MAX_XB];
  double   Bmu2pt[MAX_XB];
  double   Bmu1p[MAX_XB];
  double   Bmu2p[MAX_XB];
  double   Bmu1eta[MAX_XB];
  double   Bmu2eta[MAX_XB];
  double   Bmu1phi[MAX_XB];
  double   Bmu2phi[MAX_XB];
  double   Bmu1y[MAX_XB];
  double   Bmu2y[MAX_XB];
  double   Bmu1dzPV[MAX_XB];
  double   Bmu2dzPV[MAX_XB];
  double   Bmu1dxyPV[MAX_XB];
  double   Bmu2dxyPV[MAX_XB];
  double   Bmu1normchi2[MAX_XB];
  double   Bmu2normchi2[MAX_XB];
  double   Bmu1Chi2ndf[MAX_XB];
  double   Bmu2Chi2ndf[MAX_XB];
  int      Bmu1muqual[MAX_XB];
  int      Bmu2muqual[MAX_XB];
  bool     Bmu1TrackerMuArbitrated[MAX_XB];
  bool     Bmu2TrackerMuArbitrated[MAX_XB];
  bool     Bmu1isTrackerMuon[MAX_XB];
  bool     Bmu2isTrackerMuon[MAX_XB];
  bool     Bmu1isGlobalMuon[MAX_XB];
  bool     Bmu2isGlobalMuon[MAX_XB];
  bool     Bmu1TMOneStationTight[MAX_XB];
  bool     Bmu2TMOneStationTight[MAX_XB];
  int      Bmu1InPixelLayer[MAX_XB];
  int      Bmu2InPixelLayer[MAX_XB];
  int      Bmu1InStripLayer[MAX_XB];
  int      Bmu2InStripLayer[MAX_XB];
  int      Bmu1InTrackerLayer[MAX_XB];
  int      Bmu2InTrackerLayer[MAX_XB];
  
  //BInfo.mumuInfo
  double   Bmumumass[MAX_XB];
  double   Bmumueta[MAX_XB];
  double   Bmumuphi[MAX_XB];
  double   Bmumuy[MAX_XB];
  double   Bmumupt[MAX_XB];
  
  //BInfo.ujInfo
  double   Bujmass[MAX_XB];
  double   BujvProb[MAX_XB];
  double   Bujpt[MAX_XB];
  double   Bujeta[MAX_XB];
  double   Bujphi[MAX_XB];
  double   Bujy[MAX_XB];
  double   Bujlxy[MAX_XB];
  
  //BInfo.trkInfo
  int      Btrk1Idx[MAX_XB];
  int      Btrk2Idx[MAX_XB];
  double   Btrk1Pt[MAX_XB];
  double   Btrk2Pt[MAX_XB];
  double   Btrk1Eta[MAX_XB];
  double   Btrk2Eta[MAX_XB];
  double   Btrk1Phi[MAX_XB];
  double   Btrk2Phi[MAX_XB];
  double   Btrk1PtErr[MAX_XB];
  double   Btrk2PtErr[MAX_XB];
  double   Btrk1EtaErr[MAX_XB];
  double   Btrk2EtaErr[MAX_XB];
  double   Btrk1PhiErr[MAX_XB];
  double   Btrk2PhiErr[MAX_XB];
  double   Btrk1Y[MAX_XB];
  double   Btrk2Y[MAX_XB];
  double   Btrk1Dxy[MAX_XB];
  double   Btrk2Dxy[MAX_XB];
  double   Btrk1D0Err[MAX_XB];
  double   Btrk2D0Err[MAX_XB];
  double   Btrk1PixelHit[MAX_XB];
  double   Btrk2PixelHit[MAX_XB];
  double   Btrk1StripHit[MAX_XB];
  double   Btrk2StripHit[MAX_XB];
  double   Btrk1nPixelLayer[MAX_XB];
  double   Btrk2nPixelLayer[MAX_XB];
  double   Btrk1nStripLayer[MAX_XB];
  double   Btrk2nStripLayer[MAX_XB];
  double   Btrk1Chi2ndf[MAX_XB];
  double   Btrk2Chi2ndf[MAX_XB];
  double   Btrk1MVAVal[MAX_XB];
  double   Btrk2MVAVal[MAX_XB];
  int      Btrk1Algo[MAX_XB];
  int      Btrk2Algo[MAX_XB];
  bool     Btrk1highPurity[MAX_XB];
  bool     Btrk2highPurity[MAX_XB];
  int      Btrk1Quality[MAX_XB];
  int      Btrk2Quality[MAX_XB];
  //BInfo.tktkInfo
  double   Btktkmass[MAX_XB];
  double   BtktkmassKK[MAX_XB];
  double   BtktkvProb[MAX_XB];
  double   Btktkpt[MAX_XB];
  double   Btktketa[MAX_XB];
  double   Btktkphi[MAX_XB];
  double   Btktky[MAX_XB];
  double   Bdoubletmass[MAX_XB];
  double   Bdoubletpt[MAX_XB];
  double   Bdoubleteta[MAX_XB];
  double   Bdoubletphi[MAX_XB];
  double   Bdoublety[MAX_XB];
  
  //BInfo.genInfo
  double   Bgen[MAX_XB];
  int      BgenIndex[MAX_XB];
  double   Bgenpt[MAX_XB];
  double   Bgeneta[MAX_XB];
  double   Bgenphi[MAX_XB];
  double   Bgeny[MAX_XB];
  
  /*
    double   BMaxDoca[MAX_XB];
    double   Balpha[MAX_XB];
    double   BsvpvDistance[MAX_XB];
    double   BsvpvDisErr[MAX_XB];
    double   BsvpvDistance_2D[MAX_XB];
    double   BsvpvDisErr_2D[MAX_XB];
    int      kstar[MAX_XB]; 
    double   Btrk1MassHypo[MAX_XB];
    double   Btrk2MassHypo[MAX_XB];
    double   Btrkminpt[MAX_XB];
    double   Btrkmaxpt[MAX_XB];
    int      Btrkminptindex[MAX_XB];
    int      Btrkmaxptindex[MAX_XB];
  */
  
  void buildBranch(TTree* nt)
  {
    //EvtInfo
    nt->Branch("RunNo",&RunNo);
    nt->Branch("EvtNo",&EvtNo);
    nt->Branch("LumiNo",&LumiNo);
    nt->Branch("Bsize",&Bsize);
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
    
    //BInfo
    nt->Branch("Bindex",Bindex,"Bindex[Bsize]/I");
    nt->Branch("Btype",Btype,"Btype[Bsize]/I");
    nt->Branch("Bmass",Bmass,"Bmass[Bsize]/D");
    nt->Branch("Bpt",Bpt,"Bpt[Bsize]/D");
    nt->Branch("Beta",Beta,"Beta[Bsize]/D");
    nt->Branch("Bphi",Bphi,"Bphi[Bsize]/D");
    nt->Branch("By",By,"By[Bsize]/D");
    nt->Branch("BvtxX",BvtxX,"BvtxX[Bsize]/D");
    nt->Branch("BvtxY",BvtxY,"BvtxY[Bsize]/D");
    nt->Branch("Bd0",Bd0,"Bd0[Bsize]/D");
    nt->Branch("Bd0Err",Bd0Err,"Bd0Err[Bsize]/D");
    nt->Branch("Bdxyz",Bdxyz,"Bdxyz[Bsize]/D");
    nt->Branch("BdxyzErr",BdxyzErr,"BdxyzErr[Bsize]/D");
    nt->Branch("Bchi2ndf",Bchi2ndf,"Bchi2ndf[Bsize]/D");
    nt->Branch("Bchi2cl",Bchi2cl,"Bchi2cl[Bsize]/D");
    nt->Branch("Bdtheta",Bdtheta,"Bdtheta[Bsize]/D");
    nt->Branch("Blxy",Blxy,"Blxy[Bsize]/D");
    nt->Branch("BlxyBS",BlxyBS,"BlxyBS[Bsize]/D");
    nt->Branch("BlxyBSErr",BlxyBSErr,"BlxyBSErr[Bsize]/D");
    nt->Branch("Bmaxpt",Bmaxpt,"Bmaxpt[Bsize]/O");
    nt->Branch("Bmaxprob",Bmaxprob,"Bmaxprob[Bsize]/O");
    nt->Branch("Bbesttktkmass",Bbesttktkmass,"Bbesttktkmass[Bsize]/O");
    nt->Branch("BmaxptMatched",BmaxptMatched,"BmaxptMatched[Bsize]/O");
    nt->Branch("BmaxprobMatched",BmaxprobMatched,"BmaxprobMatched[Bsize]/O");
    nt->Branch("BbesttktkmassMatched",BbesttktkmassMatched,"BbesttktkmassMatched[Bsize]/O");
    
    //BInfo.trkInfo
    nt->Branch("Btrk1Idx",Btrk1Idx,"Btrk1Idx[Bsize]/I");
    nt->Branch("Btrk2Idx",Btrk2Idx,"Btrk2Idx[Bsize]/I");
    nt->Branch("Btrk1Pt",Btrk1Pt,"Btrk1Pt[Bsize]/D");
    nt->Branch("Btrk2Pt",Btrk2Pt,"Btrk2Pt[Bsize]/D");
    nt->Branch("Btrk1Eta",Btrk1Eta,"Btrk1Eta[Bsize]/D");  
    nt->Branch("Btrk2Eta",Btrk2Eta,"Btrk2Eta[Bsize]/D");  
    nt->Branch("Btrk1Phi",Btrk1Phi,"Btrk1Phi[Bsize]/D");  
    nt->Branch("Btrk2Phi",Btrk2Phi,"Btrk2Phi[Bsize]/D");  
    nt->Branch("Btrk1PtErr",Btrk1PtErr,"Btrk1PtErr[Bsize]/D");  
    nt->Branch("Btrk2PtErr",Btrk2PtErr,"Btrk2PtErr[Bsize]/D");
    nt->Branch("Btrk1EtaErr",Btrk1EtaErr,"Btrk1EtaErr[Bsize]/D");
    nt->Branch("Btrk2EtaErr",Btrk2EtaErr,"Btrk2EtaErr[Bsize]/D");
    nt->Branch("Btrk1PhiErr",Btrk1PhiErr,"Btrk1PhiErr[Bsize]/D");
    nt->Branch("Btrk2PhiErr",Btrk2PhiErr,"Btrk2PhiErr[Bsize]/D");
    nt->Branch("Btrk1Y",Btrk1Y,"Btrk1Y[Bsize]/D");  
    nt->Branch("Btrk2Y",Btrk2Y,"Btrk2Y[Bsize]/D");  
    nt->Branch("Btrk1Dxy",Btrk1Dxy,"Btrk1Dxy[Bsize]/D");
    nt->Branch("Btrk2Dxy",Btrk2Dxy,"Btrk2Dxy[Bsize]/D");
    nt->Branch("Btrk1D0Err",Btrk1D0Err,"Btrk1D0Err[Bsize]/D");
    nt->Branch("Btrk2D0Err",Btrk2D0Err,"Btrk2D0Err[Bsize]/D");
    nt->Branch("Btrk1PixelHit",Btrk1PixelHit,"Btrk1PixelHit[Bsize]/D");
    nt->Branch("Btrk2PixelHit",Btrk2PixelHit,"Btrk2PixelHit[Bsize]/D");
    nt->Branch("Btrk1StripHit",Btrk1StripHit,"Btrk1StripHit[Bsize]/D");
    nt->Branch("Btrk2StripHit",Btrk2StripHit,"Btrk2StripHit[Bsize]/D");
    nt->Branch("Btrk1nPixelLayer",Btrk1nPixelLayer,"Btrk1nPixelLayer[Bsize]/D");
    nt->Branch("Btrk2nPixelLayer",Btrk2nPixelLayer,"Btrk2nPixelLayer[Bsize]/D");
    nt->Branch("Btrk1nStripLayer",Btrk1nStripLayer,"Btrk1nStripLayer[Bsize]/D");
    nt->Branch("Btrk2nStripLayer",Btrk2nStripLayer,"Btrk2nStripLayer[Bsize]/D");
    nt->Branch("Btrk1Chi2ndf",Btrk1Chi2ndf,"Btrk1Chi2ndf[Bsize]/D");
    nt->Branch("Btrk2Chi2ndf",Btrk2Chi2ndf,"Btrk2Chi2ndf[Bsize]/D");
    nt->Branch("Btrk1MVAVal",Btrk1MVAVal,"Btrk1MVAVal[Bsize]/D");
    nt->Branch("Btrk2MVAVal",Btrk2MVAVal,"Btrk2MVAVal[Bsize]/D");
    nt->Branch("Btrk1Algo",Btrk1Algo,"Btrk1Algo[Bsize]/I");
    nt->Branch("Btrk2Algo",Btrk2Algo,"Btrk2Algo[Bsize]/I");
    nt->Branch("Btrk1highPurity",Btrk1highPurity,"Btrk1highPurity[Bsize]/O");
    nt->Branch("Btrk2highPurity",Btrk2highPurity,"Btrk2highPurity[Bsize]/O");
    nt->Branch("Btrk1Quality",Btrk1Quality,"Btrk1Quality[Bsize]/I");
    nt->Branch("Btrk2Quality",Btrk2Quality,"Btrk2Quality[Bsize]/I");
    
    //BInfo.tktkInfo
    nt->Branch("Btktkmass",Btktkmass,"Btktkmass[Bsize]/D");
    nt->Branch("BtktkmassKK",BtktkmassKK,"BtktkmassKK[Bsize]/D");
    nt->Branch("BtktkvProb",BtktkvProb,"BtktkvProb[Bsize]/D");
    nt->Branch("Btktkpt",Btktkpt,"Btktkpt[Bsize]/D");
    nt->Branch("Btktketa",Btktketa,"Btktketa[Bsize]/D");
    nt->Branch("Btktkphi",Btktkphi,"Btktkphi[Bsize]/D");
    nt->Branch("Btktky",Btktky,"Btktky[Bsize]/D");
    nt->Branch("Bdoubletmass",Bdoubletmass,"Bdoubletmass[Bsize]/D");
    nt->Branch("Bdoubletpt",Bdoubletpt,"Bdoubletpt[Bsize]/D");
    nt->Branch("Bdoubleteta",Bdoubleteta,"Bdoubleteta[Bsize]/D");  
    nt->Branch("Bdoubletphi",Bdoubletphi,"Bdoubletphi[Bsize]/D");  
    nt->Branch("Bdoublety",Bdoublety,"Bdoublety[Bsize]/D");
    
    //BInfo.muonInfo
    nt->Branch("Bmu1pt",Bmu1pt,"Bmu1pt[Bsize]/D");
    nt->Branch("Bmu2pt",Bmu2pt,"Bmu2pt[Bsize]/D");
    nt->Branch("Bmu1p",Bmu1p,"Bmu1p[Bsize]/D");
    nt->Branch("Bmu2p",Bmu2p,"Bmu2p[Bsize]/D");
    nt->Branch("Bmu1eta",Bmu1eta,"Bmu1eta[Bsize]/D");
    nt->Branch("Bmu2eta",Bmu2eta,"Bmu2eta[Bsize]/D");
    nt->Branch("Bmu1phi",Bmu1phi,"Bmu1phi[Bsize]/D");
    nt->Branch("Bmu2phi",Bmu2phi,"Bmu2phi[Bsize]/D");
    nt->Branch("Bmu1y",Bmu1y,"Bmu1y[Bsize]/D");
    nt->Branch("Bmu2y",Bmu2y,"Bmu2y[Bsize]/D");
    nt->Branch("Bmu1dzPV",Bmu1dzPV,"Bmu1dzPV[Bsize]/D");
    nt->Branch("Bmu2dzPV",Bmu2dzPV,"Bmu2dzPV[Bsize]/D");
    nt->Branch("Bmu1dxyPV",Bmu1dxyPV,"Bmu1dxyPV[Bsize]/D");
    nt->Branch("Bmu2dxyPV",Bmu2dxyPV,"Bmu2dxyPV[Bsize]/D");
    nt->Branch("Bmu1normchi2",Bmu1normchi2,"Bmu1normchi2[Bsize]/D");
    nt->Branch("Bmu2normchi2",Bmu2normchi2,"Bmu2normchi2[Bsize]/D");
    nt->Branch("Bmu1Chi2ndf",Bmu1Chi2ndf,"Bmu1Chi2ndf[Bsize]/D");
    nt->Branch("Bmu2Chi2ndf",Bmu2Chi2ndf,"Bmu2Chi2ndf[Bsize]/D");
    nt->Branch("Bmu1muqual",Bmu1muqual,"Bmu1muqual[Bsize]/I");
    nt->Branch("Bmu1muqual",Bmu1muqual,"Bmu1muqual[Bsize]/I");
    nt->Branch("Bmu1TrackerMuArbitrated",Bmu1TrackerMuArbitrated,"Bmu1TrackerMuArbitrated[Bsize]/O");
    nt->Branch("Bmu2TrackerMuArbitrated",Bmu2TrackerMuArbitrated,"Bmu2TrackerMuArbitrated[Bsize]/O");
    nt->Branch("Bmu1isTrackerMuon",Bmu1isTrackerMuon,"Bmu1isTrackerMuon[Bsize]/O");
    nt->Branch("Bmu2isTrackerMuon",Bmu2isTrackerMuon,"Bmu2isTrackerMuon[Bsize]/O");
    nt->Branch("Bmu1isGlobalMuon",Bmu1isGlobalMuon,"Bmu1isGlobalMuon[Bsize]/O");
    nt->Branch("Bmu2isGlobalMuon",Bmu2isGlobalMuon,"Bmu2isGlobalMuon[Bsize]/O");
    nt->Branch("Bmu1TMOneStationTight",Bmu1TMOneStationTight,"Bmu1TMOneStationTight[Bsize]/O");
    nt->Branch("Bmu2TMOneStationTight",Bmu2TMOneStationTight,"Bmu2TMOneStationTight[Bsize]/O");
    nt->Branch("Bmu1InPixelLayer",Bmu1InPixelLayer,"Bmu1InPixelLayer[Bsize]/I");
    nt->Branch("Bmu2InPixelLayer",Bmu2InPixelLayer,"Bmu2InPixelLayer[Bsize]/I");
    nt->Branch("Bmu1InStripLayer",Bmu1InStripLayer,"Bmu1InStripLayer[Bsize]/I");
    nt->Branch("Bmu2InStripLayer",Bmu2InStripLayer,"Bmu2InStripLayer[Bsize]/I");
    nt->Branch("Bmu1InTrackerLayer",Bmu1InTrackerLayer,"Bmu1InTrackerLayer[Bsize]/I");
    nt->Branch("Bmu2InTrackerLayer",Bmu2InTrackerLayer,"Bmu2InTrackerLayer[Bsize]/I");
    nt->Branch("Bmumumass",Bmumumass,"Bmumumass[Bsize]/D");
    nt->Branch("Bmumueta",Bmumueta,"Bmumueta[Bsize]/D");
    nt->Branch("Bmumuphi",Bmumuphi,"Bmumuphi[Bsize]/D");
    nt->Branch("Bmumuy",Bmumuy,"Bmumuy[Bsize]/D");
    nt->Branch("Bmumupt",Bmumupt,"Bmumupt[Bsize]/D");
    nt->Branch("Bujmass",Bujmass,"Bujmass[Bsize]/D");
    nt->Branch("BujvProb",BujvProb,"BujvProb[Bsize]/D");
    nt->Branch("Bujpt",Bujpt,"Bujpt[Bsize]/D");
    nt->Branch("Bujeta",Bujeta,"Bujeta[Bsize]/D");
    nt->Branch("Bujphi",Bujphi,"Bujphi[Bsize]/D");
    nt->Branch("Bujy",Bujy,"Bujy[Bsize]/D");
    nt->Branch("Bujlxy",Bujlxy,"Bujlxy[Bsize]/D");
    
    //BInfo.genInfo
    nt->Branch("Bgen",Bgen,"Bgen[Bsize]/D");
    nt->Branch("BgenIndex",BgenIndex,"BgenIndex[Bsize]/I");
    nt->Branch("Bgenpt",Bgenpt,"Bgenpt[Bsize]/D");
    nt->Branch("Bgeny",Bgeny,"Bgeny[Bsize]/D");
    nt->Branch("Bgeneta",Bgeneta,"Bgeneta[Bsize]/D");
    nt->Branch("Bgenphi",Bgenphi,"Bgenphi[Bsize]/D");
    
    /*
      nt->Branch("Balpha",Balpha,"Balpha[Bsize]/D");
      nt->Branch("BsvpvDistance",BsvpvDistance,"BsvpvDistance[Bsize]/D");
      nt->Branch("BsvpvDisErr",BsvpvDisErr,"BsvpvDisErr[Bsize]/D");
      nt->Branch("BsvpvDistance_2D",BsvpvDistance_2D,"BsvpvDistance_2D[Bsize]/D");
      nt->Branch("BsvpvDisErr_2D",BsvpvDisErr_2D,"BsvpvDisErr_2D[Bsize]/D");
      nt->Branch("BMaxDoca",BMaxDoca,"BMaxDoca[Bsize]/D");
      nt->Branch("Btrk1MassHypo",Btrk1MassHypo,"Btrk1MassHypo[Bsize]/D");
      nt->Branch("Btrk2MassHypo",Btrk2MassHypo,"Btrk2MassHypo[Bsize]/D");
      nt->Branch("Btrkminpt",Btrkminpt,"Btrkminpt[Bsize]/D");
      nt->Branch("Btrkmaxpt",Btrkmaxpt,"Btrkmaxpt[Bsize]/D");
      nt->Branch("Btrkminptindex",Btrkminptindex,"Btrkminptindex[Bsize]/I");
      nt->Branch("Btrkmaxptindex",Btrkmaxptindex,"Btrkmaxptindex[Bsize]/I");
    */
  }
  
  int      Gsize;
  double   Gy[MAX_GEN];
  double   Geta[MAX_GEN];
  double   Gphi[MAX_GEN];
  double   Gpt[MAX_GEN];
  double   GpdgId[MAX_GEN];
  double   GisSignal[MAX_GEN];
  double   Gmu1pt[MAX_GEN];
  double   Gmu2pt[MAX_GEN];
  double   Gmu1p[MAX_GEN];
  double   Gmu2p[MAX_GEN];
  double   Gmu1eta[MAX_GEN];
  double   Gmu2eta[MAX_GEN];
  double   Gmu1phi[MAX_GEN];
  double   Gmu2phi[MAX_GEN];
  double   Gtk1pt[MAX_GEN];
  double   Gtk2pt[MAX_GEN];
  double   Gtk1eta[MAX_GEN];
  double   Gtk2eta[MAX_GEN];
  double   Gtk1phi[MAX_GEN];
  double   Gtk2phi[MAX_GEN];
  
  void buildGenBranch(TTree* nt)
  {
    nt->Branch("Gsize",&Gsize);
    nt->Branch("Gy",Gy,"Gy[Gsize]/D");
    nt->Branch("Geta",Geta,"Geta[Gsize]/D");
    nt->Branch("Gphi",Gphi,"Gphi[Gsize]/D");
    nt->Branch("Gpt",Gpt,"Gpt[Gsize]/D");
    nt->Branch("GpdgId",GpdgId,"GpdgId[Gsize]/D");
    nt->Branch("GisSignal",GisSignal,"GisSignal[Gsize]/D");
    nt->Branch("Gmu1eta",Gmu1eta,"Gmu1eta[Gsize]/D");
    nt->Branch("Gmu1phi",Gmu1phi,"Gmu1phi[Gsize]/D");
    nt->Branch("Gmu1pt",Gmu1pt,"Gmu1pt[Gsize]/D");
    nt->Branch("Gmu1p",Gmu1p,"Gmu1p[Gsize]/D");
    nt->Branch("Gmu2eta",Gmu2eta,"Gmu2eta[Gsize]/D");
    nt->Branch("Gmu2phi",Gmu2phi,"Gmu2phi[Gsize]/D");
    nt->Branch("Gmu2pt",Gmu2pt,"Gmu2pt[Gsize]/D");
    nt->Branch("Gmu2p",Gmu2p,"Gmu2p[Gsize]/D");
    nt->Branch("Gtk1pt",Gtk1pt,"Gtk1pt[Gsize]/D");
    nt->Branch("Gtk1eta",Gtk1eta,"Gtk1eta[Gsize]/D");
    nt->Branch("Gtk1phi",Gtk1phi,"Gtk1phi[Gsize]/D");
    nt->Branch("Gtk2pt",Gtk2pt,"Gtk2pt[Gsize]/D");
    nt->Branch("Gtk2eta",Gtk2eta,"Gtk2eta[Gsize]/D");
    nt->Branch("Gtk2phi",Gtk2phi,"Gtk2phi[Gsize]/D");
  }
  
  void makeNtuple(int ifchannel[], bool REAL, bool skim, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, MuonInfoBranches *MuonInfo, TrackInfoBranches *TrackInfo, BInfoBranches *BInfo, GenInfoBranches *GenInfo, TTree* nt0, TTree* nt1, TTree* nt2, TTree* nt3, TTree* nt5, TTree* nt6)
  {//{{{
    TVector3* bP = new TVector3;
    TVector3* bVtx = new TVector3;
    TLorentzVector* b4P = new TLorentzVector;
    //TLorentzVector* b4Pout = new TLorentzVector;
    int Btypesize[7]={0,0,0,0,0,0,0};
    int ptflag=-1,ptMatchedflag=-1,probflag=-1,probMatchedflag=-1,tktkflag=-1,tktkMatchedflag=-1;
    double pttem=0,ptMatchedtem=0,probtem=0,probMatchedtem=0,tktktem=0,tktkMatchedtem=0;
    for(int t=0;t<7;t++)
      {
        int tidx = t-1;
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
            for(int j=0;j<BInfo->size;j++)
              {
                if(skim)
                  {
                    if(BInfo->pt[j]<3.) continue;
                  }
                if(BInfo->type[j]==(t+1))
                  {
                    fillTree(bP,bVtx,b4P,j,Btypesize[tidx],tk1mass[t],tk2mass[t],REAL, EvtInfo, VtxInfo, MuonInfo, TrackInfo, BInfo, GenInfo);
                    if(BInfo->pt[j]>pttem)
                      {
                        ptflag = Btypesize[tidx];
                        pttem = BInfo->pt[j];
                      }
                    if(TMath::Prob(BInfo->vtxchi2[j],BInfo->vtxdof[j])>probtem)
                      {
                        probflag = Btypesize[tidx];
                        probtem = TMath::Prob(BInfo->vtxchi2[j],BInfo->vtxdof[j]);
                      }
                    if(BInfo->type[j]>2&&BInfo->type[j]<7)
                      {
                        if(TMath::Abs(BInfo->tktk_mass[j]-midmass[t])<tktktem)
                          {
                            tktkflag = Btypesize[tidx];
                            tktktem = TMath::Abs(BInfo->tktk_mass[j]-midmass[t]);
                          }
                      }
                    if((!REAL&&(Bgen[Btypesize[tidx]]==23333||Bgen[Btypesize[tidx]]==41000))||REAL)//////////////////////////////
                      {
                        if(BInfo->pt[j]>ptMatchedtem)
                          {
                            ptMatchedflag = Btypesize[tidx];
                            ptMatchedtem = BInfo->pt[j];
                          }
                        if(TMath::Prob(BInfo->vtxchi2[j],BInfo->vtxdof[j])>probMatchedtem)
                          {
                            probMatchedflag = Btypesize[tidx];
                            probMatchedtem = TMath::Prob(BInfo->vtxchi2[j],BInfo->vtxdof[j]);
                          }
                        if(BInfo->type[j]>2&&BInfo->type[j]<7)
                          {
                            if(TMath::Abs(BInfo->tktk_mass[j]-midmass[t])<tktkMatchedtem)
                              {
                                tktkMatchedflag = Btypesize[tidx];
                                tktkMatchedtem = TMath::Abs(BInfo->tktk_mass[j]-midmass[t]);
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
  }//}}}
  
  void fillGenTree(TTree* ntGen, GenInfoBranches *GenInfo, bool gskim=true)
  {//{{{
    TLorentzVector* bGen = new TLorentzVector;
    int gt=0,sigtype=0;
    int gsize=0;
    Gsize = 0;
    for(int j=0;j<GenInfo->size;j++)
      {
        if((TMath::Abs(GenInfo->pdgId[j])!=BPLUS_PDGID&&TMath::Abs(GenInfo->pdgId[j])!=BZERO_PDGID&&TMath::Abs(GenInfo->pdgId[j])!=BSUBS_PDGID) && gskim) continue;
        Gsize = gsize+1;
        Gpt[gsize] = GenInfo->pt[j];
        Geta[gsize] = GenInfo->eta[j];
        Gphi[gsize] = GenInfo->phi[j];
        GpdgId[gsize] = GenInfo->pdgId[j];
        bGen->SetPtEtaPhiM(GenInfo->pt[j],GenInfo->eta[j],GenInfo->phi[j],GenInfo->mass[j]);
        Gy[gsize] = bGen->Rapidity();
        sigtype=0;
        for(gt=1;gt<8;gt++)
          {
            if(signalGen(gt,j, GenInfo))
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
        gsize++;
      }
    ntGen->Fill();
  }//}}}
  
  void fillTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, int j, int typesize, double track_mass1, double track_mass2, bool REAL, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, MuonInfoBranches *MuonInfo, TrackInfoBranches *TrackInfo, BInfoBranches *BInfo, GenInfoBranches *GenInfo)
  {//{{{
    //Event Info
    RunNo = EvtInfo->RunNo;
    EvtNo = EvtInfo->EvtNo;
    LumiNo = EvtInfo->LumiNo;
    Bsize = typesize+1;
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
    bP->SetPtEtaPhi(BInfo->pt[j],BInfo->eta[j]*0,BInfo->phi[j]);
    bVtx->SetXYZ(BInfo->vtxX[j]-EvtInfo->PVx,
                 BInfo->vtxY[j]-EvtInfo->PVy,
                 BInfo->vtxZ[j]*0-EvtInfo->PVz*0);
    b4P->SetPtEtaPhiM(BInfo->pt[j],BInfo->eta[j],BInfo->phi[j],BInfo->mass[j]);
    
    Bindex[typesize] = typesize;
    Btype[typesize] = BInfo->type[j];
    Bmass[typesize] = BInfo->mass[j];
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
    double r2lxyBS = (BInfo->vtxX[j]-EvtInfo->BSx+(BInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz) * (BInfo->vtxX[j]-EvtInfo->BSx+(BInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz)
      + (BInfo->vtxY[j]-EvtInfo->BSy+(BInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz) * (BInfo->vtxY[j]-EvtInfo->BSy+(BInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz);
    double xlxyBS = BInfo->vtxX[j]-EvtInfo->BSx + (BInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz;
    double ylxyBS = BInfo->vtxY[j]-EvtInfo->BSy + (BInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz;
    BlxyBS[typesize] = TMath::Sqrt(r2lxyBS);
    //BlxyBSErr[typesize] = 0;
    BlxyBSErr[typesize] = (1./r2lxyBS) * ((xlxyBS*xlxyBS)*BInfo->vtxXErr[j] + (2*xlxyBS*ylxyBS)*BInfo->vtxYXErr[j] + (ylxyBS*ylxyBS)*BInfo->vtxYErr[j]);
    Bmaxpt[typesize] = false;
    Bmaxprob[typesize] = false;
    Bbesttktkmass[typesize] = false;
    BmaxptMatched[typesize] = false;
    BmaxprobMatched[typesize] = false;
    BbesttktkmassMatched[typesize] = false;
    
    b4P->SetPtEtaPhiM(MuonInfo->pt[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]],MuonInfo->eta[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]],MuonInfo->phi[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]],MUON_MASS);
    Bmu1pt[typesize] = MuonInfo->pt[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu1eta[typesize] = MuonInfo->eta[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu1phi[typesize] = MuonInfo->phi[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu1p[typesize] = b4P->P();
    Bmu1y[typesize] = b4P->Rapidity();
    Double_t mu1px,mu1py,mu1pz,mu1E;
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
    Double_t mu2px,mu2py,mu2pz,mu2E;
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
    Bmu1InPixelLayer[typesize] = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2InPixelLayer[typesize] = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1InStripLayer[typesize] = MuonInfo->i_nStripLayer[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2InStripLayer[typesize] = MuonInfo->i_nStripLayer[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
    Bmu1InTrackerLayer[typesize] = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]] + MuonInfo->i_nStripLayer[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    Bmu2InTrackerLayer[typesize] = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]] + MuonInfo->i_nStripLayer[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
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
    
    double tk1px,tk1py,tk1pz,tk1E;
    double tk2px,tk2py,tk2pz,tk2E;
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
        Btrk1Dxy[typesize] = TrackInfo->dxyPV[BInfo->rftk1_index[j]];
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
        Btrk1Dxy[typesize] = TrackInfo->dxyPV[BInfo->rftk2_index[j]];
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
        Btrk2Dxy[typesize] = TrackInfo->dxyPV[BInfo->rftk1_index[j]];
        Btrk2D0Err[typesize] = TrackInfo->d0error[BInfo->rftk1_index[j]];
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
        double tk1EK = b4P->E();
        b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk2_index[j]],TrackInfo->eta[BInfo->rftk2_index[j]],TrackInfo->phi[BInfo->rftk2_index[j]],KAON_MASS);
        double tk2EK = b4P->E();
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
        Btrk1Dxy[typesize] = TrackInfo->dxyPV[BInfo->rftk1_index[j]];
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
        Btrk2Dxy[typesize] = TrackInfo->dxyPV[BInfo->rftk2_index[j]];
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
        double tk1EK = b4P->E();
        b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk2_index[j]],TrackInfo->eta[BInfo->rftk2_index[j]],TrackInfo->phi[BInfo->rftk2_index[j]],KAON_MASS);
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
        //int ujGenIdxMu1=-1;
        //int ujGenIdxMu2=-1;
        
        double BId,MId,tk1Id,tk2Id;
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
  
  bool signalGen(int Btype, int j, GenInfoBranches *GenInfo)
  {//{{{
    double BId,MId,tk1Id,tk2Id;
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
  }//}}}
  
};//}}}

#endif
