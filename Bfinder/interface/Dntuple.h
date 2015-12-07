// vim:set ts=4 sw=4 fdm=marker et:
#ifndef _DNTUPLE_H_
#define _DNTUPLE_H_
#include "format.h"

class DntupleBranches{//{{{
    public:
        //EvtInfo
        int      RunNo;
        int      EvtNo;
        int      Bsize;
        int      Dsize;
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
        int      Dindex[MAX_XB];
        int      Dtype[MAX_XB];
        double   Dmass[MAX_XB];
        double   Dpt[MAX_XB];
        double   Deta[MAX_XB];
        double   Dphi[MAX_XB];
        double   Dy[MAX_XB];
        double   DvtxX[MAX_XB];
        double   DvtxY[MAX_XB];
        double   Dd0[MAX_XB];
        double   Dd0Err[MAX_XB];
        double   Ddxyz[MAX_XB];
        double   DdxyzErr[MAX_XB];
        double   Dchi2ndf[MAX_XB];
        double   Dchi2cl[MAX_XB];
        double   Ddtheta[MAX_XB];
        double   Dlxy[MAX_XB];
        double   Dalpha[MAX_XB];
        double   DsvpvDistance[MAX_XB];
        double   DsvpvDisErr[MAX_XB];
        double   DsvpvDistance_2D[MAX_XB];
        double   DsvpvDisErr_2D[MAX_XB];
        double   DlxyBS[MAX_XB];
        double   DlxyBSErr[MAX_XB];
        double   DMaxDoca[MAX_XB];
        bool     Dmaxpt[MAX_XB];
        bool     Dmaxprob[MAX_XB];
        bool     DmaxptMatched[MAX_XB];
        bool     DmaxprobMatched[MAX_XB];

        //DInfo.trkInfo
        int      Dtrk1Idx[MAX_XB];
        int      Dtrk2Idx[MAX_XB];
        int      Dtrk3Idx[MAX_XB];
        int      Dtrk4Idx[MAX_XB];
        double   Dtrk1Pt[MAX_XB];
        double   Dtrk2Pt[MAX_XB];
        double   Dtrk3Pt[MAX_XB];
        double   Dtrk4Pt[MAX_XB];
        double   Dtrk1Eta[MAX_XB];
        double   Dtrk2Eta[MAX_XB];
        double   Dtrk3Eta[MAX_XB];
        double   Dtrk4Eta[MAX_XB];
        double   Dtrk1Phi[MAX_XB];
        double   Dtrk2Phi[MAX_XB];
        double   Dtrk3Phi[MAX_XB];
        double   Dtrk4Phi[MAX_XB];
        double   Dtrk1PtErr[MAX_XB];
        double   Dtrk2PtErr[MAX_XB];
        double   Dtrk3PtErr[MAX_XB];
        double   Dtrk4PtErr[MAX_XB];
        double   Dtrk1EtaErr[MAX_XB];
        double   Dtrk2EtaErr[MAX_XB];
        double   Dtrk3EtaErr[MAX_XB];
        double   Dtrk4EtaErr[MAX_XB];
        double   Dtrk1PhiErr[MAX_XB];
        double   Dtrk2PhiErr[MAX_XB];
        double   Dtrk3PhiErr[MAX_XB];
        double   Dtrk4PhiErr[MAX_XB];
        double   Dtrk1Y[MAX_XB];
        double   Dtrk2Y[MAX_XB];
        double   Dtrk3Y[MAX_XB];
        double   Dtrk4Y[MAX_XB];
        double   Dtrk1Dxy[MAX_XB];
        double   Dtrk2Dxy[MAX_XB];
        double   Dtrk3Dxy[MAX_XB];
        double   Dtrk4Dxy[MAX_XB];
        double   Dtrk1D0Err[MAX_XB];
        double   Dtrk2D0Err[MAX_XB];
        double   Dtrk3D0Err[MAX_XB];
        double   Dtrk4D0Err[MAX_XB];
        double   Dtrk1PixelHit[MAX_XB];
        double   Dtrk2PixelHit[MAX_XB];
        double   Dtrk3PixelHit[MAX_XB];
        double   Dtrk4PixelHit[MAX_XB];
        double   Dtrk1StripHit[MAX_XB];
        double   Dtrk2StripHit[MAX_XB];
        double   Dtrk3StripHit[MAX_XB];
        double   Dtrk4StripHit[MAX_XB];
        double   Dtrk1nStripLayer[MAX_XB];
        double   Dtrk2nStripLayer[MAX_XB];
        double   Dtrk3nStripLayer[MAX_XB];
        double   Dtrk4nStripLayer[MAX_XB];
        double   Dtrk1nPixelLayer[MAX_XB];
        double   Dtrk2nPixelLayer[MAX_XB];
        double   Dtrk3nPixelLayer[MAX_XB];
        double   Dtrk4nPixelLayer[MAX_XB];
        double   Dtrk1Chi2ndf[MAX_XB];
        double   Dtrk2Chi2ndf[MAX_XB];
        double   Dtrk3Chi2ndf[MAX_XB];
        double   Dtrk4Chi2ndf[MAX_XB];
        double   Dtrk1MassHypo[MAX_XB];
        double   Dtrk2MassHypo[MAX_XB];
        double   Dtrk3MassHypo[MAX_XB];
        double   Dtrk4MassHypo[MAX_XB];
        double   Dtrkminpt[MAX_XB];
        double   Dtrkmaxpt[MAX_XB];
        int      Dtrkminptindex[MAX_XB];
        int      Dtrkmaxptindex[MAX_XB];
        double   Dtrk1MVAVal[MAX_XB];
        double   Dtrk2MVAVal[MAX_XB];
        double   Dtrk3MVAVal[MAX_XB];
        double   Dtrk4MVAVal[MAX_XB];
        int      Dtrk1Algo[MAX_XB];
        int      Dtrk2Algo[MAX_XB];
        int      Dtrk3Algo[MAX_XB];
        int      Dtrk4Algo[MAX_XB];
        bool     Dtrk1highPurity[MAX_XB];
        bool     Dtrk2highPurity[MAX_XB];
        bool     Dtrk3highPurity[MAX_XB];
        bool     Dtrk4highPurity[MAX_XB];
        int      Dtrk1Quality[MAX_XB];
        int      Dtrk2Quality[MAX_XB];
        int      Dtrk3Quality[MAX_XB];
        int      Dtrk4Quality[MAX_XB];
        //DInfo.tktkResInfo
        double   DtktkResmass[MAX_XB];
        double   DtktkRespt[MAX_XB];
        double   DtktkReseta[MAX_XB];
        double   DtktkResphi[MAX_XB];
        //DInfo.genInfo
        double   Dgen[MAX_XB];
        int      DgennDa[MAX_XB];
        int      DgenIndex[MAX_XB];
        double   Dgenpt[MAX_XB];
        double   Dgeneta[MAX_XB];
        double   Dgenphi[MAX_XB];
        double   Dgeny[MAX_XB];

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
        int      Gsize;
        double   Gy[MAX_GEN];
        double   Geta[MAX_GEN];
        double   Gphi[MAX_GEN];
        double   Gpt[MAX_GEN];
        double   GpdgId[MAX_GEN];
        int      GisSignal[MAX_GEN];

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

        void makeDNtuple(int isDchannel[], bool REAL, bool skim, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, TrackInfoBranches *TrackInfo, DInfoBranches *DInfo, GenInfoBranches *GenInfo, TTree* ntD1, TTree* ntD2, TTree* ntD3)
        {//{{{
            TVector3* bP = new TVector3;
            TVector3* bVtx = new TVector3;
            TLorentzVector* b4P = new TLorentzVector;
            int Dtypesize[3]={0,0,0};
            //int Ndbc=0;
            int ptflag=-1,ptMatchedflag=-1,probflag=-1,probMatchedflag=-1;
            double pttem=0,ptMatchedtem=0,probtem=0,probMatchedtem=0;
            for(int t=0;t<6;t++)
            {
                if(t%2==0)
                {
                    Dsize = 0;
                    ptflag = -1;
                    pttem = 0;
                    ptMatchedflag = -1;
                    ptMatchedtem = 0;
                    probflag = -1;
                    probtem = 0;
                    probMatchedflag = -1;
                    probMatchedtem = 0;
                }
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
                            //fillDTree(bP,bVtx,b4P,j,Dtypesize[t/2],REAL);
                            fillDTree(bP,bVtx,b4P,j,Dtypesize[t/2],REAL, EvtInfo, VtxInfo, TrackInfo, DInfo, GenInfo);
                            if(DInfo->pt[j]>pttem)
                            {
                                ptflag = Dtypesize[t/2];
                                pttem = DInfo->pt[j];
                            }
                            if(TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j])>probtem)
                            {
                                probflag = Dtypesize[t/2];
                                probtem = TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]);
                            }
                            if(((!REAL&&(Dgen[Dtypesize[t/2]]==23333||Dgen[Dtypesize[t/2]]==23344))||REAL)&&Dtrk1Pt[Dtypesize[t/2]]>8.&&Dtrk2Pt[Dtypesize[t/2]]>8.&&Dchi2cl[Dtypesize[t/2]]>0.05&&(DsvpvDistance[Dtypesize[t/2]]/DsvpvDisErr[Dtypesize[t/2]])>2.5&&TMath::Cos(Dalpha[Dtypesize[t/2]])>0.95&&(Dtrk1highPurity[Dtypesize[t/2]]&&Dtrk2highPurity[Dtypesize[t/2]]))
                            {
                                if(DInfo->pt[j]>ptMatchedtem)
                                {
                                    ptMatchedflag = Dtypesize[t/2];
                                    ptMatchedtem = DInfo->pt[j];
                                }
                                if(TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j])>probMatchedtem)
                                {
                                    probMatchedflag = Dtypesize[t/2];
                                    probMatchedtem = TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]);
                                }
                            }
                            Dtypesize[t/2]++;
                        }
                    }
                    if(t%2==1)
                    {
                        if(ptflag>=0) Dmaxpt[ptflag] = true;
                        if(probflag>=0) Dmaxprob[probflag] = true;
                        if(ptMatchedflag>=0) DmaxptMatched[ptMatchedflag] = true;
                        if(probMatchedflag>=0) DmaxprobMatched[probMatchedflag] = true;
                    }	      
                    if(t==1)      ntD1->Fill();
                    else if(t==3) ntD2->Fill();
                    else if(t==5) ntD3->Fill();
                }
            }
        }//}}}

        void fillDGenTree(TTree* ntGen, GenInfoBranches *GenInfo, bool gskim=true)
        {//{{{
            TLorentzVector* bGen = new TLorentzVector;
            int gt=0,sigtype=0;
            int gsize=0;
            Gsize = 0;
            for(int j=0;j<GenInfo->size;j++)
            {
                if(TMath::Abs(GenInfo->pdgId[j])!=DZERO_PDGID&&gskim) continue;
                Gsize = gsize+1;
                Gpt[gsize] = GenInfo->pt[j];
                Geta[gsize] = GenInfo->eta[j];
                Gphi[gsize] = GenInfo->phi[j];
                GpdgId[gsize] = GenInfo->pdgId[j];
                bGen->SetPtEtaPhiM(GenInfo->pt[j],GenInfo->eta[j],GenInfo->phi[j],GenInfo->mass[j]);
                Gy[gsize] = bGen->Rapidity();
                sigtype=0;
                for(gt=1;gt<5;gt++)
                {
                    if(isDsignalGen(gt,j, GenInfo))
                    {
                        sigtype=gt;
                        break;
                    }
                }
                GisSignal[gsize] = sigtype;
                gsize++;
            }
            ntGen->Fill();
        }//}}}

        double findMass(int particlePdgId)
        {//{{{
            if(TMath::Abs(particlePdgId)==211) return PION_MASS;
            if(TMath::Abs(particlePdgId)==321) return KAON_MASS;
            else
            {
                std::cout<<"ERROR: find particle mass falied >> Particle pdgId: "<<particlePdgId<<std::endl;
                return 0;
            }
        }//}}}

        void fillDTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, int j, int typesize, bool REAL, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, TrackInfoBranches *TrackInfo, DInfoBranches *DInfo, GenInfoBranches *GenInfo)
        {//{{{
            //EvtInfo
            RunNo = EvtInfo->RunNo;
            EvtNo = EvtInfo->EvtNo;
            Dsize = typesize+1;
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

            //DInfo
            bP->SetPtEtaPhi(DInfo->pt[j],DInfo->eta[j]*0,DInfo->phi[j]);
            bVtx->SetXYZ(DInfo->vtxX[j]-EvtInfo->PVx,
                         DInfo->vtxY[j]-EvtInfo->PVy,
                         DInfo->vtxZ[j]*0-EvtInfo->PVz*0);
            b4P->SetPtEtaPhiM(DInfo->pt[j],DInfo->eta[j],DInfo->phi[j],DInfo->mass[j]);
            Dindex[typesize] = typesize;
            Dtype[typesize] = DInfo->type[j];
            Dmass[typesize] = DInfo->mass[j];
            Dpt[typesize] = DInfo->pt[j];
            Deta[typesize] = DInfo->eta[j];
            Dphi[typesize] = DInfo->phi[j];
            Dy[typesize] = b4P->Rapidity();
            DvtxX[typesize] = DInfo->vtxX[j] - EvtInfo->PVx;
            DvtxY[typesize] = DInfo->vtxY[j] - EvtInfo->PVy;
            Dd0[typesize] = TMath::Sqrt((DInfo->vtxX[j]-EvtInfo->PVx)*(DInfo->vtxX[j]-EvtInfo->PVx)+(DInfo->vtxY[j]-EvtInfo->PVy)*(DInfo->vtxY[j]-EvtInfo->PVy));
            Dd0Err[typesize] = TMath::Sqrt(DInfo->vtxXErr[j]*DInfo->vtxXErr[j]+DInfo->vtxYErr[j]*DInfo->vtxYErr[j]);
            Ddxyz[typesize] = TMath::Sqrt((DInfo->vtxX[j]-EvtInfo->PVx)*(DInfo->vtxX[j]-EvtInfo->PVx)+(DInfo->vtxY[j]-EvtInfo->PVy)*(DInfo->vtxY[j]-EvtInfo->PVy)+(DInfo->vtxZ[j]-EvtInfo->PVz)*(DInfo->vtxZ[j]-EvtInfo->PVz));
            DdxyzErr[typesize] = TMath::Sqrt(DInfo->vtxXErr[j]*DInfo->vtxXErr[j]+DInfo->vtxYErr[j]*DInfo->vtxYErr[j]+DInfo->vtxZErr[j]*DInfo->vtxZErr[j]);
            Dchi2ndf[typesize] = DInfo->vtxchi2[j]/DInfo->vtxdof[j];
            Dchi2cl[typesize] = TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]);
            Ddtheta[typesize] = bP->Angle(*bVtx);
            Dlxy[typesize] = ((DInfo->vtxX[j]-EvtInfo->PVx)*b4P->Px() + (DInfo->vtxY[j]-EvtInfo->PVy)*b4P->Py())/DInfo->pt[j];
            Dalpha[typesize] = DInfo->alpha[j];
            DsvpvDistance[typesize] = DInfo->svpvDistance[j];
            DsvpvDisErr[typesize] = DInfo->svpvDisErr[j];
            DsvpvDistance_2D[typesize] = DInfo->svpvDistance_2D[j];
            DsvpvDisErr_2D[typesize] = DInfo->svpvDisErr_2D[j];
            double r2lxyBS = (DInfo->vtxX[j]-EvtInfo->BSx+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz) * (DInfo->vtxX[j]-EvtInfo->BSx+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz)
                + (DInfo->vtxY[j]-EvtInfo->BSy+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz) * (DInfo->vtxY[j]-EvtInfo->BSy+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz);
            double xlxyBS = DInfo->vtxX[j]-EvtInfo->BSx + (DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz;
            double ylxyBS = DInfo->vtxY[j]-EvtInfo->BSy + (DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz;
            DlxyBS[typesize] = TMath::Sqrt(r2lxyBS);
            DlxyBSErr[typesize] = (1./r2lxyBS) * ((xlxyBS*xlxyBS)*DInfo->vtxXErr[j] + (2*xlxyBS*ylxyBS)*DInfo->vtxYXErr[j] + (ylxyBS*ylxyBS)*DInfo->vtxYErr[j]);
            DMaxDoca[typesize] = DInfo->MaxDoca[j];
            Dmaxpt[typesize] = false;
            Dmaxprob[typesize] = false;
            DmaxptMatched[typesize] = false;
            DmaxprobMatched[typesize] = false;

            //DInfo.trkInfo
            double trk1mass,trk2mass,trk3mass,trk4mass;
            Dtrk1Idx[typesize] = DInfo->rftk1_index[j];
            Dtrk2Idx[typesize] = DInfo->rftk2_index[j];
            Dtrk1Pt[typesize] = TrackInfo->pt[DInfo->rftk1_index[j]];
            Dtrk2Pt[typesize] = TrackInfo->pt[DInfo->rftk2_index[j]];
            Dtrk1Eta[typesize] = TrackInfo->eta[DInfo->rftk1_index[j]];
            Dtrk2Eta[typesize] = TrackInfo->eta[DInfo->rftk2_index[j]];
            Dtrk1Phi[typesize] = TrackInfo->phi[DInfo->rftk1_index[j]];
            Dtrk2Phi[typesize] = TrackInfo->phi[DInfo->rftk2_index[j]];
            Dtrk1PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk1_index[j]];
            Dtrk2PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk2_index[j]];
            Dtrk1EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk1_index[j]];
            Dtrk2EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk2_index[j]];
            Dtrk1PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk1_index[j]];
            Dtrk2PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk2_index[j]];
            trk1mass = findMass(DInfo->rftk1_MassHypo[j]);
            trk2mass = findMass(DInfo->rftk2_MassHypo[j]);
            b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk1_index[j]],TrackInfo->eta[DInfo->rftk1_index[j]],TrackInfo->phi[DInfo->rftk1_index[j]],trk1mass);
            Dtrk1Y[typesize] = b4P->Rapidity();
            b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk2_index[j]],TrackInfo->eta[DInfo->rftk2_index[j]],TrackInfo->phi[DInfo->rftk2_index[j]],trk2mass);
            Dtrk2Y[typesize] = b4P->Rapidity();
            Dtrk1Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk1_index[j]];
            Dtrk2Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk2_index[j]];
            Dtrk1D0Err[typesize] = TrackInfo->d0error[DInfo->rftk1_index[j]];
            Dtrk2D0Err[typesize] = TrackInfo->d0error[DInfo->rftk2_index[j]];
            Dtrk1PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk1_index[j]];
            Dtrk2PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk2_index[j]];
            Dtrk1StripHit[typesize] = TrackInfo->striphit[DInfo->rftk1_index[j]];
            Dtrk2StripHit[typesize] = TrackInfo->striphit[DInfo->rftk2_index[j]];
            Dtrk1nPixelLayer[typesize] = TrackInfo->nPixelLayer[DInfo->rftk1_index[j]];
            Dtrk2nPixelLayer[typesize] = TrackInfo->nPixelLayer[DInfo->rftk2_index[j]];
            Dtrk1nStripLayer[typesize] = TrackInfo->nStripLayer[DInfo->rftk1_index[j]];
            Dtrk2nStripLayer[typesize] = TrackInfo->nStripLayer[DInfo->rftk2_index[j]];
            Dtrk1Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk1_index[j]]/TrackInfo->ndf[DInfo->rftk1_index[j]];
            Dtrk2Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk2_index[j]]/TrackInfo->ndf[DInfo->rftk2_index[j]];
            Dtrk1MassHypo[typesize] = DInfo->rftk1_MassHypo[j]*TrackInfo->charge[DInfo->rftk1_index[j]];
            Dtrk2MassHypo[typesize] = DInfo->rftk2_MassHypo[j]*TrackInfo->charge[DInfo->rftk2_index[j]];
            Dtrk1MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk1_index[j]];
            Dtrk2MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk2_index[j]];
            Dtrk1Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk1_index[j]];
            Dtrk2Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk2_index[j]];
            Dtrk1highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk1_index[j]];
            Dtrk2highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk2_index[j]];
            Dtrk1Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk1_index[j]];
            Dtrk2Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk2_index[j]];
            Dtrkminpt[typesize] = (TrackInfo->pt[DInfo->rftk1_index[j]]<TrackInfo->pt[DInfo->rftk2_index[j]])?TrackInfo->pt[DInfo->rftk1_index[j]]:TrackInfo->pt[DInfo->rftk2_index[j]];
            Dtrkmaxpt[typesize] = (TrackInfo->pt[DInfo->rftk1_index[j]]>TrackInfo->pt[DInfo->rftk2_index[j]])?TrackInfo->pt[DInfo->rftk1_index[j]]:TrackInfo->pt[DInfo->rftk2_index[j]];
            Dtrkminptindex[typesize] = (TrackInfo->pt[DInfo->rftk1_index[j]]<TrackInfo->pt[DInfo->rftk2_index[j]])?1:2;
            Dtrkmaxptindex[typesize] = (TrackInfo->pt[DInfo->rftk1_index[j]]>TrackInfo->pt[DInfo->rftk2_index[j]])?1:2;
            if(DInfo->type[j]==1||DInfo->type[j]==2)
            {
                Dtrk3Idx[typesize] = -1;
                Dtrk4Idx[typesize] = -1;
                Dtrk3Pt[typesize] = -1;
                Dtrk4Pt[typesize] = -1;
                Dtrk3Eta[typesize] = -20;
                Dtrk4Eta[typesize] = -20;
                Dtrk3Phi[typesize] = -20;
                Dtrk4Phi[typesize] = -20;
                Dtrk3PtErr[typesize] = 0;
                Dtrk4PtErr[typesize] = 0;
                Dtrk3EtaErr[typesize] = 0;
                Dtrk4EtaErr[typesize] = 0;
                Dtrk3PhiErr[typesize] = 0;
                Dtrk4PhiErr[typesize] = 0;
                Dtrk3Y[typesize] = -1;
                Dtrk4Y[typesize] = -1;
                Dtrk3Dxy[typesize] = -1;
                Dtrk4Dxy[typesize] = -1;
                Dtrk3D0Err[typesize] = -1;
                Dtrk4D0Err[typesize] = -1;
                Dtrk3PixelHit[typesize] = -1;
                Dtrk4PixelHit[typesize] = -1;
                Dtrk3StripHit[typesize] = -1;
                Dtrk4StripHit[typesize] = -1;
                Dtrk1nPixelLayer[typesize] = -1;
                Dtrk2nPixelLayer[typesize] = -1;
                Dtrk1nStripLayer[typesize] = -1;
                Dtrk2nStripLayer[typesize] = -1;
                Dtrk3Chi2ndf[typesize] = -1;
                Dtrk4Chi2ndf[typesize] = -1;
                Dtrk3MassHypo[typesize] = 0;
                Dtrk4MassHypo[typesize] = 0;
                Dtrk3MVAVal[typesize] = -100;
                Dtrk4MVAVal[typesize] = -100;
                Dtrk3highPurity[typesize] = false;
                Dtrk4highPurity[typesize] = false;
                DtktkResmass[typesize] = -1;
                DtktkRespt[typesize] = -1;
                DtktkReseta[typesize] = -20;
                DtktkResphi[typesize] = -20;
            }
            else if(DInfo->type[j]==3||DInfo->type[j]==4)
            {
                Dtrk3Idx[typesize] = DInfo->rftk3_index[j];
                Dtrk3Pt[typesize] = TrackInfo->pt[DInfo->rftk3_index[j]];
                Dtrk3Eta[typesize] = TrackInfo->eta[DInfo->rftk3_index[j]];
                Dtrk3Phi[typesize] = TrackInfo->phi[DInfo->rftk3_index[j]];
                Dtrk3PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk3_index[j]];
                Dtrk3EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk3_index[j]];
                Dtrk3PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk3_index[j]];
                trk3mass = findMass(DInfo->rftk3_MassHypo[j]);
                b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk3_index[j]],TrackInfo->eta[DInfo->rftk3_index[j]],TrackInfo->phi[DInfo->rftk3_index[j]],trk3mass);
                Dtrk3Y[typesize] = b4P->Rapidity();
                Dtrk3Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk3_index[j]];
                Dtrk3D0Err[typesize] = TrackInfo->d0error[DInfo->rftk3_index[j]];
                Dtrk3PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk3_index[j]];
                Dtrk3StripHit[typesize] = TrackInfo->striphit[DInfo->rftk3_index[j]];
                Dtrk3nPixelLayer[typesize] = TrackInfo->nPixelLayer[DInfo->rftk3_index[j]];
                Dtrk3nStripLayer[typesize] = TrackInfo->nStripLayer[DInfo->rftk3_index[j]];
                Dtrk3Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk3_index[j]]/TrackInfo->ndf[DInfo->rftk3_index[j]];
                Dtrk3MassHypo[typesize] = DInfo->rftk3_MassHypo[j]*TrackInfo->charge[DInfo->rftk3_index[j]];
                Dtrk3MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk3_index[j]];
                Dtrk3Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk3_index[j]];
                Dtrk3highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk3_index[j]];
                Dtrk3Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk3_index[j]];
                Dtrk4Idx[typesize] = -1;
                Dtrk4Pt[typesize] = -1;
                Dtrk4Eta[typesize] = -20;
                Dtrk4Phi[typesize] = -20;
                Dtrk4PtErr[typesize] = 0;
                Dtrk4EtaErr[typesize] = 0;
                Dtrk4PhiErr[typesize] = 0;
                Dtrk4Y[typesize] = -1;
                Dtrk4Dxy[typesize] = -1;
                Dtrk4D0Err[typesize] = -1;
                Dtrk4PixelHit[typesize] = -1;
                Dtrk4StripHit[typesize] = -1;
                Dtrk4nPixelLayer[typesize] = -1;
                Dtrk4nStripLayer[typesize] = -1;
                Dtrk4Chi2ndf[typesize] = -1;
                Dtrk4MassHypo[typesize] = 0;
                Dtrk4MVAVal[typesize] = -100;
                Dtrk4Algo[typesize] = 0;
                Dtrk4Quality[typesize] = 0;
                Dtrk4highPurity[typesize] = false;
                DtktkResmass[typesize] = -1;
                DtktkRespt[typesize] = -1;
                DtktkReseta[typesize] = -20;
                DtktkResphi[typesize] = -20;
            }
            else if(DInfo->type[j]==5||DInfo->type[j]==6)
            {
                Dtrk3Idx[typesize] = DInfo->rftk3_index[j];
                Dtrk4Idx[typesize] = DInfo->rftk4_index[j];
                Dtrk3Pt[typesize] = TrackInfo->pt[DInfo->rftk3_index[j]];
                Dtrk4Pt[typesize] = TrackInfo->pt[DInfo->rftk4_index[j]];
                Dtrk3Eta[typesize] = TrackInfo->eta[DInfo->rftk3_index[j]];
                Dtrk4Eta[typesize] = TrackInfo->eta[DInfo->rftk4_index[j]];
                Dtrk3Phi[typesize] = TrackInfo->phi[DInfo->rftk3_index[j]];
                Dtrk4Phi[typesize] = TrackInfo->phi[DInfo->rftk4_index[j]];
                Dtrk3PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk3_index[j]];
                Dtrk4PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk4_index[j]];
                Dtrk3EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk3_index[j]];
                Dtrk4EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk4_index[j]];
                Dtrk3PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk3_index[j]];
                Dtrk4PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk4_index[j]];
                trk3mass = findMass(DInfo->rftk3_MassHypo[j]);
                trk4mass = findMass(DInfo->rftk4_MassHypo[j]);
                b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk3_index[j]],TrackInfo->eta[DInfo->rftk3_index[j]],TrackInfo->phi[DInfo->rftk3_index[j]],trk3mass);
                Dtrk3Y[typesize] = b4P->Rapidity();
                b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk4_index[j]],TrackInfo->eta[DInfo->rftk4_index[j]],TrackInfo->phi[DInfo->rftk4_index[j]],trk4mass);
                Dtrk4Y[typesize] = b4P->Rapidity();
                Dtrk3Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk3_index[j]];
                Dtrk4Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk4_index[j]];
                Dtrk3D0Err[typesize] = TrackInfo->d0error[DInfo->rftk3_index[j]];
                Dtrk4D0Err[typesize] = TrackInfo->d0error[DInfo->rftk4_index[j]];
                Dtrk3PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk3_index[j]];
                Dtrk4PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk4_index[j]];
                Dtrk3StripHit[typesize] = TrackInfo->striphit[DInfo->rftk3_index[j]];
                Dtrk4StripHit[typesize] = TrackInfo->striphit[DInfo->rftk4_index[j]];
                Dtrk3nPixelLayer[typesize] = TrackInfo->nPixelLayer[DInfo->rftk3_index[j]];
                Dtrk4nPixelLayer[typesize] = TrackInfo->nPixelLayer[DInfo->rftk4_index[j]];
                Dtrk3nStripLayer[typesize] = TrackInfo->nStripLayer[DInfo->rftk3_index[j]];
                Dtrk4nStripLayer[typesize] = TrackInfo->nStripLayer[DInfo->rftk4_index[j]];
                Dtrk3Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk3_index[j]]/TrackInfo->ndf[DInfo->rftk3_index[j]];
                Dtrk4Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk4_index[j]]/TrackInfo->ndf[DInfo->rftk4_index[j]];
                Dtrk3MassHypo[typesize] = DInfo->rftk3_MassHypo[j]*TrackInfo->charge[DInfo->rftk3_index[j]];
                Dtrk4MassHypo[typesize] = DInfo->rftk4_MassHypo[j]*TrackInfo->charge[DInfo->rftk4_index[j]];
                Dtrk3MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk3_index[j]];
                Dtrk4MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk4_index[j]];
                Dtrk3Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk3_index[j]];
                Dtrk4Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk4_index[j]];
                Dtrk3highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk3_index[j]];
                Dtrk4highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk4_index[j]];
                Dtrk3Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk3_index[j]];
                Dtrk4Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk4_index[j]];
                DtktkResmass[typesize] = -1;
                DtktkRespt[typesize] = -1;
                DtktkReseta[typesize] = -20;
                DtktkResphi[typesize] = -20;
            }
            else if(DInfo->type[j]==7||DInfo->type[j]==8)
            {
                Dtrk3Idx[typesize] = DInfo->rftk3_index[j];
                Dtrk4Idx[typesize] = DInfo->rftk4_index[j];
                Dtrk3Pt[typesize] = TrackInfo->pt[DInfo->rftk3_index[j]];
                Dtrk4Pt[typesize] = TrackInfo->pt[DInfo->rftk4_index[j]];
                Dtrk3Eta[typesize] = TrackInfo->eta[DInfo->rftk3_index[j]];
                Dtrk4Eta[typesize] = TrackInfo->eta[DInfo->rftk4_index[j]];
                Dtrk3Phi[typesize] = TrackInfo->phi[DInfo->rftk3_index[j]];
                Dtrk4Phi[typesize] = TrackInfo->phi[DInfo->rftk4_index[j]];
                Dtrk3PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk3_index[j]];
                Dtrk4PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk4_index[j]];
                Dtrk3EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk3_index[j]];
                Dtrk4EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk4_index[j]];
                Dtrk3PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk3_index[j]];
                Dtrk4PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk4_index[j]];
                trk3mass = findMass(DInfo->rftk3_MassHypo[j]);
                trk4mass = findMass(DInfo->rftk4_MassHypo[j]);
                b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk3_index[j]],TrackInfo->eta[DInfo->rftk3_index[j]],TrackInfo->phi[DInfo->rftk3_index[j]],trk3mass);
                Dtrk3Y[typesize] = b4P->Rapidity();
                b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk4_index[j]],TrackInfo->eta[DInfo->rftk4_index[j]],TrackInfo->phi[DInfo->rftk4_index[j]],trk4mass);
                Dtrk4Y[typesize] = b4P->Rapidity();
                Dtrk3Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk3_index[j]];
                Dtrk4Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk4_index[j]];
                Dtrk3D0Err[typesize] = TrackInfo->d0error[DInfo->rftk3_index[j]];
                Dtrk4D0Err[typesize] = TrackInfo->d0error[DInfo->rftk4_index[j]];
                Dtrk3PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk3_index[j]];
                Dtrk4PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk4_index[j]];
                Dtrk3StripHit[typesize] = TrackInfo->striphit[DInfo->rftk3_index[j]];
                Dtrk4StripHit[typesize] = TrackInfo->striphit[DInfo->rftk4_index[j]];
                Dtrk3Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk3_index[j]]/TrackInfo->ndf[DInfo->rftk3_index[j]];
                Dtrk4Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk4_index[j]]/TrackInfo->ndf[DInfo->rftk4_index[j]];
                Dtrk3MassHypo[typesize] = DInfo->rftk3_MassHypo[j]*TrackInfo->charge[DInfo->rftk3_index[j]];
                Dtrk4MassHypo[typesize] = DInfo->rftk4_MassHypo[j]*TrackInfo->charge[DInfo->rftk4_index[j]];
                Dtrk3MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk3_index[j]];
                Dtrk4MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk4_index[j]];
                Dtrk3Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk3_index[j]];
                Dtrk4Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk4_index[j]];
                Dtrk3highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk3_index[j]];
                Dtrk4highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk4_index[j]];
                Dtrk3Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk3_index[j]];
                Dtrk4Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk4_index[j]];
                DtktkResmass[typesize] = DInfo->tktkRes_mass[j];
                DtktkRespt[typesize] = DInfo->tktkRes_pt[j];
                DtktkReseta[typesize] = DInfo->tktkRes_eta[j];
                DtktkResphi[typesize] = DInfo->tktkRes_phi[j];
            }

            int hypo=-1,DpdgId=0,level=0;
            if(DInfo->type[j]==1||DInfo->type[j]==2||DInfo->type[j]==5) DpdgId=DZERO_PDGID;
            else if(DInfo->type[j]==3||DInfo->type[j]==4) DpdgId=DPLUS_PDGID;
            else if(DInfo->type[j]==6||DInfo->type[j]==7) DpdgId=DSUBS_PDGID;
            Dgen[typesize] = 0;//gen init
            DgenIndex[typesize] = -1;
            DgennDa[typesize] = -1;
            Dgenpt[typesize] = -1;
            Dgeneta[typesize] = -20;
            Dgenphi[typesize] = -20;
            Dgeny[typesize] = -1;
            //int rGenIdxTk1 = -1;
            //int rGenIdxTk2 = -1;
            int dGenIdxTk1 = -1;
            int dGenIdxTk2 = -1;
            int dGenIdxTk3 = -1;
            int dGenIdxTk4 = -1;
            //int dGenIdxRes = -1;
            if(!REAL)
            {
                if(DInfo->type[j]==1||DInfo->type[j]==2||DInfo->type[j]==3||DInfo->type[j]==4||DInfo->type[j]==5||DInfo->type[j]==6)
                {
                    if(TrackInfo->geninfo_index[DInfo->rftk1_index[j]]>-1)
                    {
                        level=0;
                        if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j])
                        {
                            hypo=-1;
                            if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==KAON_PDGID) hypo=0;
                            else if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==PION_PDGID) hypo=1;
                            if(hypo==0||hypo==1)
                            {
                                level=1;		      
                                if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]>-1)
                                {
                                    if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]])==DpdgId)
                                    {
                                        level=3;
                                        dGenIdxTk1=GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]];
                                    }
                                }
                            }
                        }
                        else if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk2_MassHypo[j] && (DInfo->type[j]==1||DInfo->type[j]==2))
                        {
                            hypo=-1;
                            if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==KAON_PDGID) hypo=0;
                            else if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==PION_PDGID) hypo=1;
                            if(hypo==0||hypo==1)
                            {
                                if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]>-1)
                                {
                                    if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]])==DpdgId)
                                    {
                                        level=4;
                                        dGenIdxTk1=GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]];
                                    }
                                }
                            }
                        }
                        Dgen[typesize]=level;
                    }
                    if(TrackInfo->geninfo_index[DInfo->rftk2_index[j]]>-1)
                    {
                        level=0;
                        if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j])
                        {
                            if((hypo==0&&TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID)||(hypo==1&&TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==KAON_PDGID))
                            {
                                if(hypo==1&&TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==KAON_PDGID) hypo=2;
                                level=1;
                                if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1)
                                {
                                    if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]])==DpdgId)
                                    {
                                        level=3;
                                        dGenIdxTk2=GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]];
                                    }
                                }
                            }
                            else if(hypo==1&&TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID&&(DInfo->type[j]==3||DInfo->type[j]==4||DInfo->type[j]==5||DInfo->type[j]==6))
                            {
                                hypo=3;
                                level=1;
                                if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1)
                                {
                                    if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]])==DpdgId)
                                    {
                                        level=3;
                                        dGenIdxTk2=GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]];
                                    }
                                }
                            }
                        }
                        else if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk1_MassHypo[j] && (DInfo->type[j]==1||DInfo->type[j]==2))
                        {
                            if((hypo==0&&TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID)||(hypo==1&&TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==KAON_PDGID))
                            {
                                if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1)
                                {
                                    if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]])==DpdgId)
                                    {
                                        level=4;
                                        dGenIdxTk2=GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]];
                                    }
                                }
                            }
                        }
                        Dgen[typesize]+=(level*10);
                    }
                    if(DInfo->type[j]==1||DInfo->type[j]==2)
                    {
                        Dgen[typesize]+=3300;
                    }
                    else
                    {
                        if(TrackInfo->geninfo_index[DInfo->rftk3_index[j]]>-1)
                        {
                            level=0;
                            if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j])
                            {
                                if(((hypo==0||hypo==2)&&TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==PION_PDGID)||(hypo==3&&TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==KAON_PDGID))
                                {
                                    if(hypo==3&&TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==KAON_PDGID) hypo=4;
                                    level=1;
                                    if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]]>-1)
                                    {
                                        if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]]])==DpdgId)
                                        {
                                            level=3;
                                            dGenIdxTk3=GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]];
                                        }
                                    }
                                    else if(hypo==3&&TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==PION_PDGID&&(DInfo->type[j]==5||DInfo->type[j]==6))
                                    {
                                        hypo=5;
                                        level=1;
                                        if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]]>-1)
                                        {
                                            if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]]])==DpdgId)
                                            {
                                                level=3;
                                                dGenIdxTk3=GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]];
                                            }
                                        }
                                    }
                                }
                            }
                            Dgen[typesize]+=(level*100);
                        }
                        if(DInfo->type[j]==3||DInfo->type[j]==4)
                        {
                            Dgen[typesize]+=3000;
                        }
                        else
                        {
                            if(TrackInfo->geninfo_index[DInfo->rftk4_index[j]]>-1)
                            {
                                level=0;
                                if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk4_MassHypo[j])
                                {
                                    if(((hypo==0||hypo==2||hypo==4)&&TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==PION_PDGID)||(hypo==5&&TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==KAON_PDGID))
                                    {
                                        level=1;
                                        if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]]>-1)
                                        {
                                            if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]]])==DpdgId)
                                            {
                                                level=3;
                                                dGenIdxTk4=GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]];
                                            }
                                        }
                                    }
                                }
                                Dgen[typesize]+=(level*1000);
                            }
                        }
                    }
                    if((Dgen[typesize]==3333||Dgen[typesize]==3344)&&dGenIdxTk1==dGenIdxTk2)
                    {
                        if(DInfo->type[j]==1||DInfo->type[j]==2)
                        {
                            Dgen[typesize]+=20000;
                        }
                        else if(dGenIdxTk1==dGenIdxTk3)
                        {
                            if(DInfo->type[j]==3||DInfo->type[j]==4)
                            {
                                Dgen[typesize]+=20000;
                            }
                            else if(dGenIdxTk1==dGenIdxTk4)
                            {
                                Dgen[typesize]+=20000;
                            }
                        }	    
                    }
                }//if(DInfo->type[j]==1||DInfo->type[j]==2||DInfo->type[j]==3||DInfo->type[j]==4||DInfo->type[j]==5)
                if(Dgen[typesize]==23333||Dgen[typesize]==23344)
                {
                    DgenIndex[typesize] = dGenIdxTk1;
                    if((DInfo->type[j]==1||DInfo->type[j]==2)&&GenInfo->nDa[DgenIndex[typesize]]>2) Dgen[typesize]=41000;
                    DgennDa[typesize] = GenInfo->nDa[DgenIndex[typesize]];
                    Dgenpt[typesize] = GenInfo->pt[DgenIndex[typesize]];
                    Dgeneta[typesize] = GenInfo->eta[DgenIndex[typesize]];
                    Dgenphi[typesize] = GenInfo->phi[DgenIndex[typesize]];
                    b4P->SetXYZM(GenInfo->pt[DgenIndex[typesize]]*cos(GenInfo->phi[DgenIndex[typesize]]),
                                 GenInfo->pt[DgenIndex[typesize]]*sin(GenInfo->phi[DgenIndex[typesize]]),
                                 GenInfo->pt[DgenIndex[typesize]]*sinh(GenInfo->eta[DgenIndex[typesize]]),
                                 GenInfo->mass[DgenIndex[typesize]]);
                    Dgeny[typesize] = b4P->Rapidity();
                }
            }//if(!real)
        }//fillDtree//}}}

        bool isDsignalGen(int dmesontype, int j, GenInfoBranches *GenInfo)
        {//{{{
            bool flag=false;
            if(dmesontype==1||dmesontype==2)
            {
                if(TMath::Abs(GenInfo->pdgId[j])==DZERO_PDGID&&GenInfo->nDa[j]==2&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1)
                {
                    if((TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID) || 
                            (TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==KAON_PDGID))
                    {
                        flag=true;
                    }
                }
            }
            return flag;
        }//}}}
};//}}}

#endif
