// vim:set ts=4 sw=4 fdm=marker et:
#ifndef _XBFRAMEFORMAT_H_
#define _XBFRAMEFORMAT_H_

//Note, when the array size gett too large, SetBranchAddress will fail, root will abort w/o error msg
#define MAX_XB 16384
#define MAX_MUON 4096
#define MAX_TRACK 8192
//#define MAX_XB 8192
//#define MAX_MUON 2048
//#define MAX_TRACK 4096
#define MAX_GEN 4096
#define MAX_BX 128
#define MAX_Vertices 4096
//#define N_TRIGGER_BOOKINGS 5842

#include <memory>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <list>
#include <iomanip>
#include <cmath>

//#include "Bfinder/Bfinder/interface/TriggerBooking.h"
#include <TLorentzVector.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <TString.h>
#include <TNtuple.h>
#include <TVector3.h>

#define ELECTRON_MASS 0.0005
#define MUON_MASS   0.10565837
#define PION_MASS   0.13957018
#define KAON_MASS   0.493677
#define KSHORT_MASS 0.497614
#define KSTAR_MASS  0.89594
#define PHI_MASS    1.019455
#define JPSI_MASS   3.096916
#define PSI2S_MASS  3.686109
#define PROTON_MASS 0.9383
#define D0_MASS 1.8648
#define DSTAR_MASS 2.01028

#define MUON_PDGID 13
#define PION_PDGID 211
#define KAON_PDGID 321
#define KSTAR_PDGID 313
#define PHI_PDGID 333
#define JPSI_PDGID 443
#define BZERO_PDGID 511
#define BPLUS_PDGID 521
#define BSUBS_PDGID 531
#define DZERO_PDGID 421
#define DPLUS_PDGID 411
#define DSUBS_PDGID 431

class EvtInfoBranches{ //{{{
	public:
		int	    RunNo;
	    int	    EvtNo;
	    int	    BxNo;
	    int	    LumiNo;
	    int	    Orbit;
	    bool	McFlag;
        int     nBX;
        int     BXPU[MAX_BX];
        int     nPU[MAX_BX];
        float   trueIT[MAX_BX];
        //int     trgCount;                   //number of successfully triggered HLT path in the booking.
        //int     nTrgBook;                   //N_TRIGGER_BOOKING
        //char    trgBook[N_TRIGGER_BOOKINGS];//status of booked triggers
		//int 	nHLT;                       //# of HLT of the event 
        //bool    hltBits[N_TRIGGER_BOOKINGS];//is HLT of the event acceptted?
		//std::vector<std::string> *hltnames;
		//bool	hltflag[N_TRIGGER_NAMES]; //status of HLT
		//int	nHLTm; //# of HLT hope to be matched
		//char	hltflagm[N_TRIGGER_NAMES]; //status of HLT hope to be matched
		double	PVx;
		double	PVy;
		double	PVz;
		double	PVxE;
		double	PVyE;
		double	PVzE;
		double	PVnchi2;
		double	PVchi2;
        double  BSx;
        double  BSy;
        double  BSz;
        double  BSxErr;
        double  BSyErr;
        double  BSzErr;
        double  BSdxdz;
        double  BSdydz;
        double  BSdxdzErr;
        double  BSdydzErr;
        double  BSWidthX;
        double  BSWidthXErr;
        double  BSWidthY;
        double  BSWidthYErr;
		//double	PVc2p;
		
		void regTree(TTree *root){//{{{
			root->Branch("EvtInfo.RunNo"        , &RunNo                     , "EvtInfo.RunNo/I"			);
			root->Branch("EvtInfo.EvtNo"        , &EvtNo                     , "EvtInfo.EvtNo/I"			);
			root->Branch("EvtInfo.BxNo"         , &BxNo                      , "EvtInfo.BxNo/I"			);
			root->Branch("EvtInfo.LumiNo"       , &LumiNo                    , "EvtInfo.LumiNo/I"			);
			root->Branch("EvtInfo.Orbit"        , &Orbit                     , "EvtInfo.Orbit/I"			);
			root->Branch("EvtInfo.McFlag"       , &McFlag                    , "EvtInfo.McFlag/O"			);
			root->Branch("EvtInfo.nBX"          , &nBX                       , "EvtInfo.nBX/I" 			);
			root->Branch("EvtInfo.BXPU"         , BXPU                       , "EvtInfo.BXPU[EvtInfo.nBX]/I");
			root->Branch("EvtInfo.nPU"          , nPU                        , "EvtInfo.nPU[EvtInfo.nBX]/I");
			root->Branch("EvtInfo.trueIT"       , trueIT                     , "EvtInfo.trueIT[EvtInfo.nBX]/F");
            //root->Branch("EvtInfo.trgCount"     , &trgCount                  , "EvtInfo.trgCount/I"       );
            //root->Branch("EvtInfo.nTrgBook"     , &nTrgBook                  , "EvtInfo.nTrgBook/I"       );
            //root->Branch("EvtInfo.trgBook"      , trgBook                    , "EvtInfo.trgBook[EvtInfo.nTrgBook]/B");//notice /B
			//root->Branch("EvtInfo.nHLT"         , &nHLT                      , "EvtInfo.nHLT/I"			);
            //root->Branch("EvtInfo.hltBits"      , hltBits                    , "EvtInfo.hltBits[EvtInfo.nHLT]/O");
			//root->Branch("EvtInfo.hltnames"   , "std::vector<std::string>" , &hltnames);
			//root->Branch("EvtInfo.hltflag"    , hltflag                    , "EvtInfo.hltflag[EvtInfo.nHLT]/O"	);
			//root->Branch("EvtInfo.nHLTm"      , &nHLTm                     , "EvtInfo.nHLTm/I"			);
			//root->Branch("EvtInfo.hltflagm"   , hltflagm                   , "EvtInfo.hltflagm[nHLTm]/O"		);
  		    root->Branch("EvtInfo.PVx"          , &PVx                       , "EvtInfo.PVx/D"			);
  		    root->Branch("EvtInfo.PVy"          , &PVy                       , "EvtInfo.PVy/D"			);
  		    root->Branch("EvtInfo.PVz"          , &PVz                       , "EvtInfo.PVz/D"			);
            root->Branch("EvtInfo.PVxE"         , &PVxE                      , "EvtInfo.PVxE/D"           );
            root->Branch("EvtInfo.PVyE"         , &PVyE                      , "EvtInfo.PVyE/D"           );
            root->Branch("EvtInfo.PVzE"         , &PVzE                      , "EvtInfo.PVzE/D"           );
  		    root->Branch("EvtInfo.PVnchi2"      , &PVnchi2                   , "EvtInfo.PVnchi2/D"		);
  		    root->Branch("EvtInfo.PVchi2"       , &PVchi2                    , "EvtInfo.PVchi2/D"			);
  		    root->Branch("EvtInfo.BSx"          , &BSx                       , "EvtInfo.BSx/D"			);
  		    root->Branch("EvtInfo.BSy"          , &BSy                       , "EvtInfo.BSy/D"			);
  		    root->Branch("EvtInfo.BSz"          , &BSz                       , "EvtInfo.BSz/D"			);
  		    root->Branch("EvtInfo.BSxErr"       , &BSxErr                    , "EvtInfo.BSxErr/D"			);
  		    root->Branch("EvtInfo.BSyErr"       , &BSyErr                    , "EvtInfo.BSyErr/D"			);
  		    root->Branch("EvtInfo.BSzErr"       , &BSzErr                    , "EvtInfo.BSzErr/D"			);
  		    root->Branch("EvtInfo.BSdxdz"       , &BSdxdz                    , "EvtInfo.BSdxdz/D"			);
  		    root->Branch("EvtInfo.BSdydz"       , &BSdydz                    , "EvtInfo.BSdydz/D"			);
  		    root->Branch("EvtInfo.BSdxdzErr"    , &BSdxdzErr                 , "EvtInfo.BSdxdzErr/D"	);
  		    root->Branch("EvtInfo.BSdydzErr"    , &BSdydzErr                 , "EvtInfo.BSdydzErr/D"	);
  		    root->Branch("EvtInfo.BSWidthX"     , &BSWidthX                  , "EvtInfo.BSWidthX/D"		);
  		    root->Branch("EvtInfo.BSWidthXErr"  , &BSWidthXErr               , "EvtInfo.BSWidthXErr/D"	);
  		    root->Branch("EvtInfo.BSWidthY"     , &BSWidthY                  , "EvtInfo.BSWidthY/D"		);
  		    root->Branch("EvtInfo.BSWidthYErr"  , &BSWidthYErr               , "EvtInfo.BSWidthYErr/D"	);
			//root->Branch("EvtInfo.PVc2p"      , &PVc2p                     , "EvtInfo.PVc2p/D"			);//
		}//}}}

    	void setbranchadd(TTree *root){ //{{{
            root->SetBranchAddress("EvtInfo.RunNo"          ,&RunNo	    );
            root->SetBranchAddress("EvtInfo.EvtNo"          ,&EvtNo       );
            root->SetBranchAddress("EvtInfo.BxNo"           ,&BxNo        );
            root->SetBranchAddress("EvtInfo.LumiNo"         ,&LumiNo      );
            root->SetBranchAddress("EvtInfo.Orbit"          ,&Orbit       );
            root->SetBranchAddress("EvtInfo.McFlag"         ,&McFlag      );
			root->SetBranchAddress("EvtInfo.nBX"            ,&nBX         );
			root->SetBranchAddress("EvtInfo.BXPU"           ,BXPU         );
			root->SetBranchAddress("EvtInfo.nPU"            ,nPU          );
			root->SetBranchAddress("EvtInfo.trueIT"         ,trueIT       );
            //root->SetBranchAddress("EvtInfo.trgCount"     ,&trgCount    );
            //root->SetBranchAddress("EvtInfo.nTrgBook"     ,&nTrgBook    );
            //root->SetBranchAddress("EvtInfo.trgBook"      ,trgBook      );
            //root->SetBranchAddress("EvtInfo.nHLT"         ,&nHLT	    );
            //root->SetBranchAddress("EvtInfo.hltBits"      ,hltBits      );
            //root->SetBranchAddress("EvtInfo.hltnames"     ,&hltnames	);
            //root->SetBranchAddress("EvtInfo.hltflag"      ,hltflag	);
            //root->SetBranchAddress("EvtInfo.nHLTm"	    ,&nHLTm	);
            //root->SetBranchAddress("EvtInfo.hltflagm"     ,hltflagm	);
            root->SetBranchAddress("EvtInfo.PVx"            ,&PVx		);
            root->SetBranchAddress("EvtInfo.PVy"            ,&PVy		);
            root->SetBranchAddress("EvtInfo.PVz"            ,&PVz		);
            root->SetBranchAddress("EvtInfo.PVxE"           ,&PVxE	);
            root->SetBranchAddress("EvtInfo.PVyE"           ,&PVyE	);
            root->SetBranchAddress("EvtInfo.PVzE"           ,&PVzE	);
            root->SetBranchAddress("EvtInfo.PVnchi2"        ,&PVnchi2	);
            root->SetBranchAddress("EvtInfo.PVchi2"         ,&PVchi2	);
            root->SetBranchAddress("EvtInfo.BSx"            ,&BSx	    );
            root->SetBranchAddress("EvtInfo.BSy"            ,&BSy	    );
            root->SetBranchAddress("EvtInfo.BSz"            ,&BSz	    );
            root->SetBranchAddress("EvtInfo.BSxErr"         ,&BSxErr  );
            root->SetBranchAddress("EvtInfo.BSyErr"         ,&BSyErr  );
            root->SetBranchAddress("EvtInfo.BSzErr"         ,&BSzErr  );
            root->SetBranchAddress("EvtInfo.BSdxdz"         ,&BSdxdz  );
            root->SetBranchAddress("EvtInfo.BSdydz"         ,&BSdydz  );
            root->SetBranchAddress("EvtInfo.BSdxdzErr"      ,&BSdxdzErr  );
            root->SetBranchAddress("EvtInfo.BSdydzErr"      ,&BSdydzErr  );
            root->SetBranchAddress("EvtInfo.BSWidthX"       ,&BSWidthX  );
            root->SetBranchAddress("EvtInfo.BSWidthXErr"    ,&BSWidthXErr  );
            root->SetBranchAddress("EvtInfo.BSWidthY"       ,&BSWidthY  );
            root->SetBranchAddress("EvtInfo.BSWidthYErr"    ,&BSWidthYErr  );
            //root->SetBranchAddress("EvtInfo.PVc2p"    ,&PVc2p	);
        } //}}}
}; //}}}

class VtxInfoBranches { //{{{
	public:
		int     Size;
		int     isValid[MAX_Vertices];
		bool    isFake[MAX_Vertices];
		float   Ndof[MAX_Vertices];
		float   NormalizedChi2[MAX_Vertices];
		float   Pt_Sum[MAX_Vertices];
		float   Pt_Sum2[MAX_Vertices];
		float   x[MAX_Vertices];
		float   y[MAX_Vertices];
		float   z[MAX_Vertices];

		void regTree(TTree *root) { //{{{
			root->Branch("VtxInfo.Size"	    , &Size	       , "VtxInfo.Size/I"	    );
			root->Branch("VtxInfo.isValid"  , &isValid[0]      , "VtxInfo.isValid[VtxInfo.Size]/I"	    );
			root->Branch("VtxInfo.isFake"   , &isFake[0]       , "VtxInfo.isFake[VtxInfo.Size]/O"	    ); 
			root->Branch("VtxInfo.Ndof"	    , &Ndof[0]	       , "VtxInfo.Ndof[VtxInfo.Size]/F"	    );
			root->Branch("VtxInfo.NormalizedChi2"	    , &NormalizedChi2[0]	       , "VtxInfo.NormalizedChi2[VtxInfo.Size]/F"	    );
			root->Branch("VtxInfo.Pt_Sum"	    , &Pt_Sum[0]	       , "VtxInfo.Pt_Sum[VtxInfo.Size]/F"	    );
			root->Branch("VtxInfo.Pt_Sum2"	    , &Pt_Sum2[0]	       , "VtxInfo.Pt_Sum2[VtxInfo.Size]/F"	    );
			root->Branch("VtxInfo.x"	    , &x[0]	       , "VtxInfo.x[VtxInfo.Size]/F"	    );
			root->Branch("VtxInfo.y"	    , &y[0]	       , "VtxInfo.y[VtxInfo.Size]/F"	    );
			root->Branch("VtxInfo.z"	    , &z[0]	       , "VtxInfo.z[VtxInfo.Size]/F"	    );
		} //}}}
	    
		void setbranchadd(TTree *root) { //{{{
			root->SetBranchAddress("VtxInfo.Size"        , &Size  	 );
			root->SetBranchAddress("VtxInfo.isValid"     , &isValid[0]  	 );
			root->SetBranchAddress("VtxInfo.isFake"      , &isFake[0]  	 );
			root->SetBranchAddress("VtxInfo.Ndof"        , &Ndof[0]  	 );
			root->SetBranchAddress("VtxInfo.NormalizedChi2"        , &NormalizedChi2[0]  	 );
			root->SetBranchAddress("VtxInfo.Pt_Sum"        , &Pt_Sum[0]  	 );
			root->SetBranchAddress("VtxInfo.Pt_Sum2"        , &Pt_Sum2[0]  	 );
			root->SetBranchAddress("VtxInfo.x"        , &x[0]  	 );
			root->SetBranchAddress("VtxInfo.y"        , &y[0]  	 );
			root->SetBranchAddress("VtxInfo.z"        , &z[0]  	 );
		} //}}}		    
};//}}}

 class MuonInfoBranches{//{{{
    public:
        int	    size;
        int     index        [ MAX_MUON];
        int	    handle_index [ MAX_MUON];
        int 	charge       [ MAX_MUON];
        double 	pt           [ MAX_MUON];
        double	eta          [ MAX_MUON];
        double 	phi          [ MAX_MUON];
        double 	ptErr        [ MAX_MUON];
        double	etaErr       [ MAX_MUON];
        double 	phiErr       [ MAX_MUON];
        bool    isTrackerMuon[ MAX_MUON];
        bool    isGlobalMuon [ MAX_MUON];
        int	    muqual       [ MAX_MUON];
        double  iso_trk      [ MAX_MUON];
        double  iso_ecal     [ MAX_MUON];
        double  iso_hcal     [ MAX_MUON];
        int     type         [ MAX_MUON];
        double  n_matches    [ MAX_MUON];
        bool    TMOneStationTight [MAX_MUON];
        bool    TrackerMuonArbitrated [MAX_MUON];
        int     geninfo_index[ MAX_MUON];
        bool    isNeededMuon[MAX_MUON];//for intermediate Bfinder usage, not stored in output
        bool    BfinderMuID [MAX_MUON];
        bool    SoftMuID [MAX_MUON];

        bool    isStandAloneMuon            [ MAX_MUON];
        int 	StandAloneMuon_charge       [ MAX_MUON];
        double 	StandAloneMuon_pt           [ MAX_MUON];
        double	StandAloneMuon_eta          [ MAX_MUON];
        double 	StandAloneMuon_phi          [ MAX_MUON];
        double	StandAloneMuon_d0           [ MAX_MUON];
        double	StandAloneMuon_dz           [ MAX_MUON];
        double	StandAloneMuon_dzPV         [ MAX_MUON];
        double	StandAloneMuon_dxyPV        [ MAX_MUON];

        bool    outerTrackisNonnull  [MAX_MUON];
        bool    innerTrackisNonnull  [MAX_MUON];
        bool    globalTrackisNonnull [MAX_MUON];
        int     innerTrackQuality    [MAX_MUON];
        double  normchi2     [ MAX_MUON];
        int	    i_striphit   [ MAX_MUON];
        int	    i_pixelhit   [ MAX_MUON];
        int     i_nStripLayer[ MAX_MUON];
        int     i_nPixelLayer[ MAX_MUON];
        double	i_chi2       [ MAX_MUON];
        double	i_ndf        [ MAX_MUON];
        int	    fpbarrelhit  [ MAX_MUON];
        int	    fpendcaphit  [ MAX_MUON];
        double	d0           [ MAX_MUON];
        double	dz           [ MAX_MUON];
        double	dzPV         [ MAX_MUON];
        double	dxyPV        [ MAX_MUON];
        double	g_chi2       [ MAX_MUON];
        double	g_ndf        [ MAX_MUON];
        int	    g_striphit   [ MAX_MUON];
        int	    g_pixelhit   [ MAX_MUON];
        int	    nmuhit       [ MAX_MUON];

        bool    isTriggered  [MAX_MUON];
        int	    MuTrgMatchPathSize;
        std::vector<std::vector<double > > *MuTrgMatchTrgObjE;
        std::vector<std::vector<double > > *MuTrgMatchTrgObjPt;
        std::vector<std::vector<double > > *MuTrgMatchTrgObjEta;
        std::vector<std::vector<double > > *MuTrgMatchTrgObjPhi;

        void regTree(TTree *root){//{{{
            root->Branch("MuonInfo.size"          , &size         , "MuonInfo.size/I"			);
            root->Branch("MuonInfo.index"         , index         , "MuonInfo.index[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.handle_index"  , handle_index  , "MuonInfo.handle_index[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.charge"        , charge        , "MuonInfo.charge[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.pt"            , pt            , "MuonInfo.pt[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.eta"           , eta           , "MuonInfo.eta[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.phi"           , phi           , "MuonInfo.phi[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.ptErr"        , ptErr         , "MuonInfo.ptErr[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.etaErr"        , etaErr        , "MuonInfo.etaErr[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.phiErr"        , phiErr        , "MuonInfo.phiErr[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.isTrackerMuon" , isTrackerMuon , "MuonInfo.isTrackerMuon[MuonInfo.size]/O");
            root->Branch("MuonInfo.isGlobalMuon"  , isGlobalMuon  , "MuonInfo.isGlobalMuon[MuonInfo.size]/O");
            root->Branch("MuonInfo.muqual"        , muqual        , "MuonInfo.muqual[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.iso_trk"       , iso_trk       , "MuonInfo.iso_trk[MuonInfo.size]/D");
            root->Branch("MuonInfo.iso_ecal"      , iso_ecal      , "MuonInfo.iso_ecal[MuonInfo.size]/D");
            root->Branch("MuonInfo.iso_hcal"      , iso_hcal      , "MuonInfo.iso_hcal[MuonInfo.size]/D");
            root->Branch("MuonInfo.type"          , type         ,  "MuonInfo.type[MuonInfo.size]/I"   );
            root->Branch("MuonInfo.n_matches"     , n_matches     , "MuonInfo.n_matches[MuonInfo.size]/D");
            root->Branch("MuonInfo.TMOneStationTight" ,TMOneStationTight, "MuonInfo.TMOneStationTight[MuonInfo.size]/O");
            root->Branch("MuonInfo.TrackerMuonArbitrated" ,TrackerMuonArbitrated, "MuonInfo.TrackerMuonArbitrated[MuonInfo.size]/O");
            root->Branch("MuonInfo.geninfo_index"    , geninfo_index    , "MuonInfo.geninfo_index[MuonInfo.size]/I");
            root->Branch("MuonInfo.BfinderMuID" ,BfinderMuID, "MuonInfo.BfinderMuID[MuonInfo.size]/O");
            root->Branch("MuonInfo.SoftMuID" ,SoftMuID, "MuonInfo.SoftMuID[MuonInfo.size]/O");

            root->Branch("MuonInfo.isStandAloneMuon"             , isStandAloneMuon             , "MuonInfo.isStandAloneMuon[MuonInfo.size]/O"		);
            root->Branch("MuonInfo.StandAloneMuon_charge"        , StandAloneMuon_charge        , "MuonInfo.StandAloneMuon_charge[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.StandAloneMuon_pt"            , StandAloneMuon_pt            , "MuonInfo.StandAloneMuon_pt[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.StandAloneMuon_eta"           , StandAloneMuon_eta           , "MuonInfo.StandAloneMuon_eta[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.StandAloneMuon_phi"           , StandAloneMuon_phi           , "MuonInfo.StandAloneMuon_phi[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.StandAloneMuon_d0"            , StandAloneMuon_d0            , "MuonInfo.StandAloneMuon_d0[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.StandAloneMuon_dz"            , StandAloneMuon_dz            , "MuonInfo.StandAloneMuon_dz[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.StandAloneMuon_dzPV"          , StandAloneMuon_dzPV          , "MuonInfo.StandAloneMuon_dzPV[MuonInfo.size]/D"   );
            root->Branch("MuonInfo.StandAloneMuon_dxyPV"         , StandAloneMuon_dxyPV         , "MuonInfo.StandAloneMuon_dxyPV[MuonInfo.size]/D"	);

            root->Branch("MuonInfo.outerTrackisNonnull" ,outerTrackisNonnull, "MuonInfo.outerTrackisNonnull[MuonInfo.size]/O");
            root->Branch("MuonInfo.innerTrackisNonnull" ,innerTrackisNonnull, "MuonInfo.innerTrackisNonnull[MuonInfo.size]/O");
            root->Branch("MuonInfo.globalTrackisNonnull" ,globalTrackisNonnull, "MuonInfo.globalTrackisNonnull[MuonInfo.size]/O");
            root->Branch("MuonInfo.innerTrackQuality"    , innerTrackQuality    , "MuonInfo.innerTrackQuality[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.normchi2"      , normchi2      , "MuonInfo.normchi2[MuonInfo.size]/D");
            root->Branch("MuonInfo.i_striphit"    , i_striphit    , "MuonInfo.i_striphit[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.i_pixelhit"    , i_pixelhit    , "MuonInfo.i_pixelhit[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.i_nStripLayer" , i_nStripLayer , "MuonInfo.i_nStripLayer[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.i_nPixelLayer" , i_nPixelLayer , "MuonInfo.i_nPixelLayer[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.i_chi2"        , i_chi2        , "MuonInfo.i_chi2[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.i_ndf"         , i_ndf         , "MuonInfo.i_ndf[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.fpbarrelhit"   , fpbarrelhit   , "MuonInfo.fpbarrelhit[MuonInfo.size]/I");
            root->Branch("MuonInfo.fpendcaphit"   , fpendcaphit   , "MuonInfo.fpendcaphit[MuonInfo.size]/I");
            root->Branch("MuonInfo.d0"            , d0            , "MuonInfo.d0[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.dz"            , dz            , "MuonInfo.dz[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.dzPV"          , dzPV          , "MuonInfo.dzPV[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.dxyPV"         , dxyPV         , "MuonInfo.dxyPV[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.g_chi2"        , g_chi2        , "MuonInfo.g_chi2[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.g_ndf"         , g_ndf         , "MuonInfo.g_ndf[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.g_striphit"    , g_striphit    , "MuonInfo.g_striphit[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.g_pixelhit"    , g_pixelhit    , "MuonInfo.g_pixelhit[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.nmuhit"        , nmuhit        , "MuonInfo.nmuhit[MuonInfo.size]/I"	);

            root->Branch("MuonInfo.isTriggered" ,isTriggered, "MuonInfo.isTriggered[MuonInfo.size]/O");
            root->Branch("MuonInfo.MuTrgMatchPathSize", &MuTrgMatchPathSize, "MuonInfo.MuTrgMatchPathSize/I");
            root->Branch("MuonInfo.MuTrgMatchTrgObjE", "std::vector<std::vector<double>>", &MuTrgMatchTrgObjE);
            root->Branch("MuonInfo.MuTrgMatchTrgObjPt","std::vector<std::vector<double>>", &MuTrgMatchTrgObjPt);
            root->Branch("MuonInfo.MuTrgMatchTrgObjEta", "std::vector<std::vector<double>>", &MuTrgMatchTrgObjEta);
            root->Branch("MuonInfo.MuTrgMatchTrgObjPhi", "std::vector<std::vector<double>>", &MuTrgMatchTrgObjPhi);
        }//}}}

        void setbranchadd(TTree *root){//{{{
            root->SetBranchAddress("MuonInfo.size"          , &size          );
            root->SetBranchAddress("MuonInfo.index"         , index          );
            root->SetBranchAddress("MuonInfo.handle_index"  , handle_index          );
            root->SetBranchAddress("MuonInfo.charge"        , charge         );
            root->SetBranchAddress("MuonInfo.pt"            , pt             );
            root->SetBranchAddress("MuonInfo.eta"           , eta            );
            root->SetBranchAddress("MuonInfo.phi"           , phi            );
            root->SetBranchAddress("MuonInfo.ptErr"         , ptErr          );
            root->SetBranchAddress("MuonInfo.etaErr"        , etaErr         );
            root->SetBranchAddress("MuonInfo.phiErr"        , phiErr         );
            root->SetBranchAddress("MuonInfo.isTrackerMuon" , isTrackerMuon);
            root->SetBranchAddress("MuonInfo.isGlobalMuon"  , isGlobalMuon);
            root->SetBranchAddress("MuonInfo.muqual"        , muqual         );
            root->SetBranchAddress("MuonInfo.iso_trk"       , iso_trk);
            root->SetBranchAddress("MuonInfo.iso_ecal"      , iso_ecal);
            root->SetBranchAddress("MuonInfo.iso_hcal"      , iso_hcal);
            root->SetBranchAddress("MuonInfo.type"          , type          );
            root->SetBranchAddress("MuonInfo.n_matches"     , n_matches);
            root->SetBranchAddress("MuonInfo.TMOneStationTight" , TMOneStationTight);
            root->SetBranchAddress("MuonInfo.TrackerMuonArbitrated" , TrackerMuonArbitrated);
            root->SetBranchAddress("MuonInfo.geninfo_index"    , geninfo_index);
            root->SetBranchAddress("MuonInfo.BfinderMuID" , BfinderMuID);
            root->SetBranchAddress("MuonInfo.SoftMuID" , SoftMuID);

            root->SetBranchAddress("MuonInfo.isStandAloneMuon"             , isStandAloneMuon              );
            root->SetBranchAddress("MuonInfo.StandAloneMuon_charge"        , StandAloneMuon_charge         );
            root->SetBranchAddress("MuonInfo.StandAloneMuon_pt"            , StandAloneMuon_pt             );
            root->SetBranchAddress("MuonInfo.StandAloneMuon_eta"           , StandAloneMuon_eta            );
            root->SetBranchAddress("MuonInfo.StandAloneMuon_phi"           , StandAloneMuon_phi            );
            root->SetBranchAddress("MuonInfo.StandAloneMuon_d0"            , StandAloneMuon_d0             );
            root->SetBranchAddress("MuonInfo.StandAloneMuon_dz"            , StandAloneMuon_dz             );
            root->SetBranchAddress("MuonInfo.StandAloneMuon_dzPV"          , StandAloneMuon_dzPV           );
            root->SetBranchAddress("MuonInfo.StandAloneMuon_dxyPV"         , StandAloneMuon_dxyPV          );

            root->SetBranchAddress("MuonInfo.outerTrackisNonnull"    , outerTrackisNonnull);
            root->SetBranchAddress("MuonInfo.innerTrackisNonnull"    , innerTrackisNonnull);
            root->SetBranchAddress("MuonInfo.globalTrackisNonnull"    , globalTrackisNonnull);
            root->SetBranchAddress("MuonInfo.innerTrackQuality"    , innerTrackQuality	);
            root->SetBranchAddress("MuonInfo.normchi2"      , normchi2);
            root->SetBranchAddress("MuonInfo.i_striphit"    , i_striphit	);
            root->SetBranchAddress("MuonInfo.i_pixelhit"    , i_pixelhit	);
            root->SetBranchAddress("MuonInfo.i_nStripLayer" , i_nStripLayer);
            root->SetBranchAddress("MuonInfo.i_nPixelLayer" , i_nPixelLayer);
            root->SetBranchAddress("MuonInfo.i_chi2"        , i_chi2		);
            root->SetBranchAddress("MuonInfo.i_ndf"         , i_ndf		);
            root->SetBranchAddress("MuonInfo.fpbarrelhit"   , fpbarrelhit    );
            root->SetBranchAddress("MuonInfo.fpendcaphit"   , fpendcaphit    );
            root->SetBranchAddress("MuonInfo.d0"            , d0             );
            root->SetBranchAddress("MuonInfo.dz"            , dz             );
            root->SetBranchAddress("MuonInfo.dzPV"          , dzPV             );
            root->SetBranchAddress("MuonInfo.dxyPV"         , dxyPV             );
            root->SetBranchAddress("MuonInfo.g_chi2"        , g_chi2		);
            root->SetBranchAddress("MuonInfo.g_ndf"         , g_ndf		);
            root->SetBranchAddress("MuonInfo.g_striphit"    , g_striphit	);
            root->SetBranchAddress("MuonInfo.g_pixelhit"    , g_pixelhit	);
            root->SetBranchAddress("MuonInfo.nmuhit"        , nmuhit		);

            root->SetBranchAddress("MuonInfo.isTriggered" , isTriggered);
            MuTrgMatchTrgObjE = new std::vector<std::vector<double > >();
            MuTrgMatchTrgObjPt= new std::vector<std::vector<double > >();                                                                                                                                   
            MuTrgMatchTrgObjEta= new std::vector<std::vector<double > >();
            MuTrgMatchTrgObjPhi= new std::vector<std::vector<double > >();
            root->SetBranchAddress("MuonInfo.MuTrgMatchPathSize", &MuTrgMatchPathSize);
            root->SetBranchAddress("MuonInfo.MuTrgMatchTrgObjE", &MuTrgMatchTrgObjE);
            root->SetBranchAddress("MuonInfo.MuTrgMatchTrgObjPt", &MuTrgMatchTrgObjPt);
            root->SetBranchAddress("MuonInfo.MuTrgMatchTrgObjEta", &MuTrgMatchTrgObjEta);
            root->SetBranchAddress("MuonInfo.MuTrgMatchTrgObjPhi", &MuTrgMatchTrgObjPhi);
        }//}}}

};//}}}

class TrackInfoBranches{//{{{
    public:
        int     size;
        int     index        [ MAX_TRACK];
        int     handle_index [ MAX_TRACK];
        int     charge       [ MAX_TRACK];
        double  pt           [ MAX_TRACK];
        double  eta          [ MAX_TRACK];
        double  phi          [ MAX_TRACK];
        double  ptErr        [ MAX_TRACK];
        double  etaErr       [ MAX_TRACK];
        double  phiErr       [ MAX_TRACK];
        //double  p            [ MAX_TRACK];
        int     striphit     [ MAX_TRACK];
        int     pixelhit     [ MAX_TRACK];
        int     nStripLayer  [ MAX_TRACK];
        int     nPixelLayer  [ MAX_TRACK];
        int	    fpbarrelhit  [ MAX_TRACK];
        int	    fpendcaphit  [ MAX_TRACK];
        double	chi2         [ MAX_TRACK];
        double	ndf          [ MAX_TRACK];
        double	d0           [ MAX_TRACK];
        double	d0error      [ MAX_TRACK];
        double	dzPV         [ MAX_TRACK];
        double	dxyPV        [ MAX_TRACK];
        int     geninfo_index[ MAX_TRACK];
        int     trackQuality [ MAX_TRACK];
        bool    highPurity   [ MAX_TRACK];
        double  trkMVAVal    [ MAX_TRACK];
        int     trkAlgo      [ MAX_TRACK];

        void regTree(TTree *root){//{{{
            root->Branch("TrackInfo.size"           ,&size		    ,"TrackInfo.size/I"			);
            root->Branch("TrackInfo.index"          ,index          ,"TrackInfo.index[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.handle_index"   ,handle_index   ,"TrackInfo.handle_index[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.charge"  	    ,charge         ,"TrackInfo.charge[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.pt"             ,pt             ,"TrackInfo.pt[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.eta"            ,eta            ,"TrackInfo.eta[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.phi"            ,phi            ,"TrackInfo.phi[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.ptErr"          ,ptErr          ,"TrackInfo.ptErr[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.etaErr"         ,etaErr         ,"TrackInfo.etaErr[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.phiErr"         ,phiErr         ,"TrackInfo.phiErr[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.striphit"	    ,striphit	    ,"TrackInfo.striphit[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.pixelhit"	    ,pixelhit	    ,"TrackInfo.pixelhit[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.nStripLayer"	,nStripLayer	,"TrackInfo.nStripLayer[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.nPixelLayer"	,nPixelLayer	,"TrackInfo.nPixelLayer[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.fpbarrelhit"	,fpbarrelhit	,"TrackInfo.fpbarrelhit[TrackInfo.size]/I");
            root->Branch("TrackInfo.fpendcaphit"	,fpendcaphit	,"TrackInfo.fpendcaphit[TrackInfo.size]/I");
            root->Branch("TrackInfo.chi2"		    ,chi2		    ,"TrackInfo.chi2[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.ndf"		    ,ndf		    ,"TrackInfo.ndf[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.d0"		        ,d0		        ,"TrackInfo.d0[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.d0error"	    ,d0error	    ,"TrackInfo.d0error[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.dzPV"           ,dzPV           ,"TrackInfo.dzPV[TrackInfo.size]/D"		);
            root->Branch("TrackInfo.dxyPV"          ,dxyPV          ,"TrackInfo.dxyPV[TrackInfo.size]/D"		);
            root->Branch("TrackInfo.geninfo_index"  ,geninfo_index  ,"TrackInfo.geninfo_index[TrackInfo.size]/I");
            root->Branch("TrackInfo.trackQuality"   ,trackQuality   ,"TrackInfo.trackQuality[TrackInfo.size]/I");
            root->Branch("TrackInfo.highPurity"     ,highPurity     ,"TrackInfo.highPurity[TrackInfo.size]/O");
            root->Branch("TrackInfo.trkMVAVal"      ,trkMVAVal      ,"TrackInfo.trkMVAVal[TrackInfo.size]/D");
            root->Branch("TrackInfo.trkAlgo"        ,trkAlgo        ,"TrackInfo.trkAlgo[TrackInfo.size]/I");
        }//}}}

        void setbranchadd(TTree *root){//{{{
            root->SetBranchAddress("TrackInfo.size"          , &size       );
            root->SetBranchAddress("TrackInfo.index"         , index       );
            root->SetBranchAddress("TrackInfo.handle_index"  , handle_index       );
            root->SetBranchAddress("TrackInfo.charge"        , charge      );
            root->SetBranchAddress("TrackInfo.pt"            , pt          );
            root->SetBranchAddress("TrackInfo.eta"           , eta         );
            root->SetBranchAddress("TrackInfo.phi"           , phi         );
            root->SetBranchAddress("TrackInfo.ptErr"         , ptErr       );
            root->SetBranchAddress("TrackInfo.etaErr"        , etaErr      );
            root->SetBranchAddress("TrackInfo.phiErr"        , phiErr      );
            root->SetBranchAddress("TrackInfo.striphit"      , striphit    );
            root->SetBranchAddress("TrackInfo.pixelhit"      , pixelhit    );
            root->SetBranchAddress("TrackInfo.nStripLayer"   , nStripLayer );
            root->SetBranchAddress("TrackInfo.nPixelLayer"   , nPixelLayer );
            root->SetBranchAddress("TrackInfo.fpbarrelhit"   , fpbarrelhit );
            root->SetBranchAddress("TrackInfo.fpendcaphit"   , fpendcaphit );
            root->SetBranchAddress("TrackInfo.chi2"          , chi2        );
            root->SetBranchAddress("TrackInfo.ndf"           , ndf         );
            root->SetBranchAddress("TrackInfo.d0"            , d0          );
            root->SetBranchAddress("TrackInfo.d0error"       , d0error     );
            root->SetBranchAddress("TrackInfo.dzPV"          , dzPV        );
            root->SetBranchAddress("TrackInfo.dxyPV"         , dxyPV       );
            root->SetBranchAddress("TrackInfo.geninfo_index" , geninfo_index  );
            root->SetBranchAddress("TrackInfo.trackQuality"  , trackQuality  );
            root->SetBranchAddress("TrackInfo.highPurity"    , highPurity  );
            root->SetBranchAddress("TrackInfo.trkMVAVal"     , trkMVAVal  );
            root->SetBranchAddress("TrackInfo.trkAlgo"       , trkAlgo  );
        }//}}}
};//}}}

class BInfoBranches{//{{{
public:
    int	    uj_size;
    int	    uj_index[MAX_XB];
    double  uj_mass[MAX_XB];
    double  uj_pt[MAX_XB];
    double  uj_eta[MAX_XB];
    double  uj_phi[MAX_XB];
    double  uj_px[MAX_XB];
    double  uj_py[MAX_XB];
    double  uj_pz[MAX_XB];
    double	uj_vtxX[MAX_XB];
    double  uj_vtxY[MAX_XB];
    double  uj_vtxZ[MAX_XB];
    double  uj_vtxXErr[MAX_XB];
    double  uj_vtxYErr[MAX_XB];
    double  uj_vtxZErr[MAX_XB];
    double  uj_vtxYXErr[MAX_XB];
    double  uj_vtxZXErr[MAX_XB];
    double  uj_vtxZYErr[MAX_XB];
    double	uj_vtxdof[MAX_XB];
    double	uj_vtxchi2[MAX_XB];
    int     uj_rfmu1_index[MAX_XB];
    int     uj_rfmu2_index[MAX_XB];
    
    double  uj_rfmu1_px[MAX_XB];//after vertexing
    double  uj_rfmu1_py[MAX_XB];
    double  uj_rfmu1_pz[MAX_XB];
    double  uj_rfmu2_px[MAX_XB];
    double  uj_rfmu2_py[MAX_XB];
    double  uj_rfmu2_pz[MAX_XB];
    
    int	    size;
    int	    index[MAX_XB];
    double	mass[MAX_XB];
    double	pt[MAX_XB];
    double	eta[MAX_XB];
    double	phi[MAX_XB];
    double	px[MAX_XB];
    double	py[MAX_XB];
    double	pz[MAX_XB];
    double  pxE[MAX_XB];
    double  pyE[MAX_XB];
    double  pzE[MAX_XB];
    double  alpha[MAX_XB];
    double  svpvDistance[MAX_XB];
    double  svpvDisErr[MAX_XB];
    double  svpvDistance_2D[MAX_XB];
    double  svpvDisErr_2D[MAX_XB];
    double  MaxDoca[MAX_XB];
    double  vtxX[MAX_XB];
    double  vtxY[MAX_XB];
    double  vtxZ[MAX_XB];
    double  vtxXErr[MAX_XB];
    double  vtxYErr[MAX_XB];
    double  vtxZErr[MAX_XB];
    double  vtxYXErr[MAX_XB];
    double  vtxZXErr[MAX_XB];
    double  vtxZYErr[MAX_XB];
    double	vtxdof[MAX_XB];
    double	vtxchi2[MAX_XB];
    int     rfuj_index[MAX_XB];
    int     rftk1_index[MAX_XB];
    int     rftk2_index[MAX_XB];
    int     type[MAX_XB];
    
    double  rfmu1_px[MAX_XB];
    double  rfmu1_py[MAX_XB];
    double  rfmu1_pz[MAX_XB];
    double  rfmu2_px[MAX_XB];
    double  rfmu2_py[MAX_XB];
    double  rfmu2_pz[MAX_XB];
    double  rftk1_px[MAX_XB];
    double  rftk1_py[MAX_XB];
    double  rftk1_pz[MAX_XB];
    double  rftk2_px[MAX_XB];
    double  rftk2_py[MAX_XB];
    double  rftk2_pz[MAX_XB];
    
    double	tktk_mass[MAX_XB];
    double  tktk_pt[MAX_XB];
    double  tktk_eta[MAX_XB];
    double  tktk_phi[MAX_XB];
    double  tktk_px[MAX_XB];
    double  tktk_py[MAX_XB];
    double  tktk_pz[MAX_XB];
    double	tktk_vtxX[MAX_XB];
    double  tktk_vtxY[MAX_XB];
    double  tktk_vtxZ[MAX_XB];
    double  tktk_vtxXErr[MAX_XB];
    double  tktk_vtxYErr[MAX_XB];
    double  tktk_vtxZErr[MAX_XB];
    double  tktk_vtxYXErr[MAX_XB];
    double  tktk_vtxZXErr[MAX_XB];
    double  tktk_vtxZYErr[MAX_XB];
    double	tktk_vtxdof[MAX_XB];
    double	tktk_vtxchi2[MAX_XB];
    double  tktk_rftk1_px[MAX_XB];
    double  tktk_rftk1_py[MAX_XB];
    double  tktk_rftk1_pz[MAX_XB];
    double  tktk_rftk2_px[MAX_XB];
    double  tktk_rftk2_py[MAX_XB];
    double  tktk_rftk2_pz[MAX_XB];
    
    void regTree(TTree *root){//{{{
        root->Branch("BInfo.uj_size"          , &uj_size       , "BInfo.uj_size/I"			);
        root->Branch("BInfo.uj_index"         , uj_index       , "BInfo.uj_index[BInfo.uj_size]/I"	);
        root->Branch("BInfo.uj_mass"          , uj_mass        , "BInfo.uj_mass[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_pt"            , uj_pt          , "BInfo.uj_pt[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_eta"            , uj_eta          , "BInfo.uj_eta[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_phi"            , uj_phi          , "BInfo.uj_phi[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_px"            , uj_px          , "BInfo.uj_px[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_py"            , uj_py          , "BInfo.uj_py[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_pz"            , uj_pz          , "BInfo.uj_pz[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_vtxX"          , uj_vtxX        , "BInfo.uj_vtxX[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_vtxY"          , uj_vtxY        , "BInfo.uj_vtxY[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_vtxZ"          , uj_vtxZ        , "BInfo.uj_vtxZ[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_vtxXErr"       , uj_vtxXErr     , "BInfo.uj_vtxXErr[BInfo.uj_size]/D"   );
        root->Branch("BInfo.uj_vtxYErr"       , uj_vtxYErr     , "BInfo.uj_vtxYErr[BInfo.uj_size]/D"   );
        root->Branch("BInfo.uj_vtxZErr"       , uj_vtxZErr     , "BInfo.uj_vtxZErr[BInfo.uj_size]/D"   );
        root->Branch("BInfo.uj_vtxYXErr"      , uj_vtxYXErr    , "BInfo.uj_vtxYXErr[BInfo.uj_size]/D"   );
        root->Branch("BInfo.uj_vtxZXErr"      , uj_vtxZXErr    , "BInfo.uj_vtxZXErr[BInfo.uj_size]/D"   );
        root->Branch("BInfo.uj_vtxZYErr"      , uj_vtxZYErr    , "BInfo.uj_vtxZYErr[BInfo.uj_size]/D"   );
        root->Branch("BInfo.uj_vtxdof"        , uj_vtxdof      , "BInfo.uj_vtxdof[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_vtxchi2"       , uj_vtxchi2     , "BInfo.uj_vtxchi2[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_rfmu1_index"   , uj_rfmu1_index , "BInfo.uj_rfmu1_index[BInfo.uj_size]/I"	);
        root->Branch("BInfo.uj_rfmu2_index"   , uj_rfmu2_index , "BInfo.uj_rfmu2_index[BInfo.uj_size]/I"	);
        
        root->Branch("BInfo.uj_rfmu1_px"      , uj_rfmu1_px    , "BInfo.uj_rfmu1_px[BInfo.uj_size]/D");
        root->Branch("BInfo.uj_rfmu1_py"      , uj_rfmu1_py    , "BInfo.uj_rfmu1_py[BInfo.uj_size]/D");
        root->Branch("BInfo.uj_rfmu1_pz"      , uj_rfmu1_pz    , "BInfo.uj_rfmu1_pz[BInfo.uj_size]/D");
        root->Branch("BInfo.uj_rfmu2_px"      , uj_rfmu2_px    , "BInfo.uj_rfmu2_px[BInfo.uj_size]/D");
        root->Branch("BInfo.uj_rfmu2_py"      , uj_rfmu2_py    , "BInfo.uj_rfmu2_py[BInfo.uj_size]/D");
        root->Branch("BInfo.uj_rfmu2_pz"      , uj_rfmu2_pz    , "BInfo.uj_rfmu2_pz[BInfo.uj_size]/D");
        
        root->Branch("BInfo.size"             , &size          , "BInfo.size/I"			);
        root->Branch("BInfo.index"            , index          , "BInfo.index[BInfo.size]/I"		);
        root->Branch("BInfo.mass"             , mass           , "BInfo.mass[BInfo.size]/D"		);
        root->Branch("BInfo.pt"               , pt             , "BInfo.pt[BInfo.size]/D"		);
        root->Branch("BInfo.eta"              , eta            , "BInfo.eta[BInfo.size]/D"		);
        root->Branch("BInfo.phi"              , phi            , "BInfo.phi[BInfo.size]/D"		);
        root->Branch("BInfo.px"               , px             , "BInfo.px[BInfo.size]/D"		);
        root->Branch("BInfo.py"               , py             , "BInfo.py[BInfo.size]/D"		);
        root->Branch("BInfo.pz"               , pz             , "BInfo.pz[BInfo.size]/D"		);
        root->Branch("BInfo.pxE"              , pxE            , "BInfo.pxE[BInfo.size]/D"            );
        root->Branch("BInfo.pyE"              , pyE            , "BInfo.pyE[BInfo.size]/D"            );
        root->Branch("BInfo.pzE"              , pzE            , "BInfo.pzE[BInfo.size]/D"            );
        root->Branch("BInfo.alpha"            , alpha          , "BInfo.alpha[BInfo.size]/D"	);
        root->Branch("BInfo.svpvDistance"     , svpvDistance   , "BInfo.svpvDistance[BInfo.size]/D"	);
        root->Branch("BInfo.svpvDisErr"       , svpvDisErr     , "BInfo.svpvDisErr[BInfo.size]/D"	);
        root->Branch("BInfo.svpvDistance_2D"  , svpvDistance_2D, "BInfo.svpvDistance_2D[BInfo.size]/D"	);
        root->Branch("BInfo.svpvDisErr_2D"    , svpvDisErr_2D  , "BInfo.svpvDisErr_2D[BInfo.size]/D"	);
        root->Branch("BInfo.MaxDoca"          , MaxDoca        , "BInfo.MaxDoca[BInfo.size]/D"	);
        root->Branch("BInfo.vtxX"             , vtxX           , "BInfo.vtxX[BInfo.size]/D"		);
        root->Branch("BInfo.vtxY"             , vtxY           , "BInfo.vtxY[BInfo.size]/D"		);
        root->Branch("BInfo.vtxZ"             , vtxZ           , "BInfo.vtxZ[BInfo.size]/D"		);
        root->Branch("BInfo.vtxXErr"          , vtxXErr        , "BInfo.vtxXErr[BInfo.size]/D"          );
        root->Branch("BInfo.vtxYErr"          , vtxYErr        , "BInfo.vtxYErr[BInfo.size]/D"          );
        root->Branch("BInfo.vtxZErr"          , vtxZErr        , "BInfo.vtxZErr[BInfo.size]/D"          );
        root->Branch("BInfo.vtxYXErr"         , vtxYXErr       , "BInfo.vtxYXErr[BInfo.size]/D"          );
        root->Branch("BInfo.vtxZXErr"         , vtxZXErr       , "BInfo.vtxZXErr[BInfo.size]/D"          );
        root->Branch("BInfo.vtxZYErr"         , vtxZYErr       , "BInfo.vtxZYErr[BInfo.size]/D"          );
        root->Branch("BInfo.vtxdof"           , vtxdof         , "BInfo.vtxdof[BInfo.size]/D"		);
        root->Branch("BInfo.vtxchi2"          , vtxchi2        , "BInfo.vtxchi2[BInfo.size]/D"	);
        root->Branch("BInfo.rfuj_index"       , rfuj_index     , "BInfo.rfuj_index[BInfo.size]/I");
        root->Branch("BInfo.rftk1_index"      , rftk1_index    , "BInfo.rftk1_index[BInfo.size]/I");
        root->Branch("BInfo.rftk2_index"      , rftk2_index    , "BInfo.rftk2_index[BInfo.size]/I");
        root->Branch("BInfo.type"             , type           , "BInfo.type[BInfo.size]/I"	);
        
        root->Branch("BInfo.rfmu1_px"         , rfmu1_px       , "BInfo.rfmu1_px[BInfo.size]/D"	);
        root->Branch("BInfo.rfmu1_py"         , rfmu1_py       , "BInfo.rfmu1_py[BInfo.size]/D"	);
        root->Branch("BInfo.rfmu1_pz"         , rfmu1_pz       , "BInfo.rfmu1_pz[BInfo.size]/D"	);
        root->Branch("BInfo.rfmu2_px"         , rfmu2_px       , "BInfo.rfmu2_px[BInfo.size]/D"	);
        root->Branch("BInfo.rfmu2_py"         , rfmu2_py       , "BInfo.rfmu2_py[BInfo.size]/D"	);
        root->Branch("BInfo.rfmu2_pz"         , rfmu2_pz       , "BInfo.rfmu2_pz[BInfo.size]/D"	);
        root->Branch("BInfo.rftk1_px"         , rftk1_px       , "BInfo.rftk1_px[BInfo.size]/D"     );
        root->Branch("BInfo.rftk1_py"         , rftk1_py       , "BInfo.rftk1_py[BInfo.size]/D"     );
        root->Branch("BInfo.rftk1_pz"         , rftk1_pz       , "BInfo.rftk1_pz[BInfo.size]/D"     );
        root->Branch("BInfo.rftk2_px"         , rftk2_px       , "BInfo.rftk2_px[BInfo.size]/D"     );
        root->Branch("BInfo.rftk2_py"         , rftk2_py       , "BInfo.rftk2_py[BInfo.size]/D"     );
        root->Branch("BInfo.rftk2_pz"         , rftk2_pz       , "BInfo.rftk2_pz[BInfo.size]/D"     );
        
        root->Branch("BInfo.tktk_mass"          , tktk_mass        , "BInfo.tktk_mass[BInfo.size]/D"     );
        root->Branch("BInfo.tktk_pt"            , tktk_pt          , "BInfo.tktk_pt[BInfo.size]/D"	);
        root->Branch("BInfo.tktk_eta"            , tktk_eta          , "BInfo.tktk_eta[BInfo.size]/D"	);
        root->Branch("BInfo.tktk_phi"            , tktk_phi          , "BInfo.tktk_phi[BInfo.size]/D"	);
        root->Branch("BInfo.tktk_px"            , tktk_px          , "BInfo.tktk_px[BInfo.size]/D"	);
        root->Branch("BInfo.tktk_py"            , tktk_py          , "BInfo.tktk_py[BInfo.size]/D"	);
        root->Branch("BInfo.tktk_pz"            , tktk_pz          , "BInfo.tktk_pz[BInfo.size]/D"	);
        root->Branch("BInfo.tktk_vtxX"          , tktk_vtxX        , "BInfo.tktk_vtxX[BInfo.size]/D"	);
        root->Branch("BInfo.tktk_vtxY"          , tktk_vtxY        , "BInfo.tktk_vtxY[BInfo.size]/D"	);
        root->Branch("BInfo.tktk_vtxZ"          , tktk_vtxZ        , "BInfo.tktk_vtxZ[BInfo.size]/D"	);
        root->Branch("BInfo.tktk_vtxXErr"       , tktk_vtxXErr     , "BInfo.tktk_vtxXErr[BInfo.size]/D"   );
        root->Branch("BInfo.tktk_vtxYErr"       , tktk_vtxYErr     , "BInfo.tktk_vtxYErr[BInfo.size]/D"   );
        root->Branch("BInfo.tktk_vtxZErr"       , tktk_vtxZErr     , "BInfo.tktk_vtxZErr[BInfo.size]/D"   );
        root->Branch("BInfo.tktk_vtxYXErr"      , tktk_vtxYXErr    , "BInfo.tktk_vtxYXErr[BInfo.size]/D"   );
        root->Branch("BInfo.tktk_vtxZXErr"      , tktk_vtxZXErr    , "BInfo.tktk_vtxZXErr[BInfo.size]/D"   );
        root->Branch("BInfo.tktk_vtxZYErr"      , tktk_vtxZYErr    , "BInfo.tktk_vtxZYErr[BInfo.size]/D"   );
        root->Branch("BInfo.tktk_vtxdof"        , tktk_vtxdof      , "BInfo.tktk_vtxdof[BInfo.size]/D"	);
        root->Branch("BInfo.tktk_vtxchi2"       , tktk_vtxchi2     , "BInfo.tktk_vtxchi2[BInfo.size]/D"	);
        root->Branch("BInfo.tktk_rftk1_px"         , tktk_rftk1_px       , "BInfo.tktk_rftk1_px[BInfo.size]/D"     );
        root->Branch("BInfo.tktk_rftk1_py"         , tktk_rftk1_py       , "BInfo.tktk_rftk1_py[BInfo.size]/D"     );
        root->Branch("BInfo.tktk_rftk1_pz"         , tktk_rftk1_pz       , "BInfo.tktk_rftk1_pz[BInfo.size]/D"     );
        root->Branch("BInfo.tktk_rftk2_px"         , tktk_rftk2_px       , "BInfo.tktk_rftk2_px[BInfo.size]/D"     );
        root->Branch("BInfo.tktk_rftk2_py"         , tktk_rftk2_py       , "BInfo.tktk_rftk2_py[BInfo.size]/D"     );
        root->Branch("BInfo.tktk_rftk2_pz"         , tktk_rftk2_pz       , "BInfo.tktk_rftk2_pz[BInfo.size]/D"     );
        
    }//}}}
    
    void setbranchadd(TTree *root){//{{{
        root->SetBranchAddress("BInfo.uj_size"		    ,&uj_size	);
        root->SetBranchAddress("BInfo.uj_size"		    ,&uj_size	    );
        root->SetBranchAddress("BInfo.uj_index"        ,uj_index   );
        root->SetBranchAddress("BInfo.uj_mass"         ,uj_mass   	);
        root->SetBranchAddress("BInfo.uj_pt"           ,uj_pt     	);
        root->SetBranchAddress("BInfo.uj_eta"           ,uj_eta     	);
        root->SetBranchAddress("BInfo.uj_phi"           ,uj_phi     	);
        root->SetBranchAddress("BInfo.uj_px"           ,uj_px     	);
        root->SetBranchAddress("BInfo.uj_py"           ,uj_py    	);
        root->SetBranchAddress("BInfo.uj_pz"           ,uj_pz   	);
        root->SetBranchAddress("BInfo.uj_vtxX"         ,uj_vtxX    );
        root->SetBranchAddress("BInfo.uj_vtxY"         ,uj_vtxY    );
        root->SetBranchAddress("BInfo.uj_vtxZ"         ,uj_vtxZ    );
        root->SetBranchAddress("BInfo.uj_vtxXErr"      ,uj_vtxXErr  );
        root->SetBranchAddress("BInfo.uj_vtxYErr"      ,uj_vtxYErr  );
        root->SetBranchAddress("BInfo.uj_vtxZErr"      ,uj_vtxZErr  );
        root->SetBranchAddress("BInfo.uj_vtxYXErr"     ,uj_vtxYXErr  );
        root->SetBranchAddress("BInfo.uj_vtxZXErr"     ,uj_vtxZXErr  );
        root->SetBranchAddress("BInfo.uj_vtxZYErr"     ,uj_vtxZYErr  );
        root->SetBranchAddress("BInfo.uj_vtxdof"	   ,uj_vtxdof	);
        root->SetBranchAddress("BInfo.uj_vtxchi2"      ,uj_vtxchi2 );
        root->SetBranchAddress("BInfo.uj_rfmu1_index"  ,uj_rfmu1_index );
        root->SetBranchAddress("BInfo.uj_rfmu2_index"  ,uj_rfmu2_index );
        
        root->SetBranchAddress("BInfo.size"            ,&size        );
        root->SetBranchAddress("BInfo.index"           ,index       	);
        root->SetBranchAddress("BInfo.mass"		    ,mass		);
        root->SetBranchAddress("BInfo.pt"		    ,pt		);
        root->SetBranchAddress("BInfo.eta"		    ,eta		);
        root->SetBranchAddress("BInfo.phi"		    ,phi		);
        root->SetBranchAddress("BInfo.px"              ,px         	);
        root->SetBranchAddress("BInfo.py"              ,py        	);
        root->SetBranchAddress("BInfo.pz"              ,pz           );
        root->SetBranchAddress("BInfo.pxE"             ,pxE          );
        root->SetBranchAddress("BInfo.pyE"             ,pyE          );
        root->SetBranchAddress("BInfo.pzE"             ,pzE         	);
        root->SetBranchAddress("BInfo.alpha"           ,alpha   	);
        root->SetBranchAddress("BInfo.svpvDistance"    ,svpvDistance   	);
        root->SetBranchAddress("BInfo.svpvDisErr"      ,svpvDisErr   	);
        root->SetBranchAddress("BInfo.svpvDistance_2D" ,svpvDistance_2D   	);
        root->SetBranchAddress("BInfo.svpvDisErr_2D"   ,svpvDisErr_2D   	);
        root->SetBranchAddress("BInfo.MaxDoca"         ,MaxDoca   	);
        root->SetBranchAddress("BInfo.vtxX"            ,vtxX       	);
        root->SetBranchAddress("BInfo.vtxY"            ,vtxY      	);
        root->SetBranchAddress("BInfo.vtxZ"            ,vtxZ     	);
        root->SetBranchAddress("BInfo.vtxXErr"         ,vtxXErr 	);
        root->SetBranchAddress("BInfo.vtxYErr"         ,vtxYErr     );
        root->SetBranchAddress("BInfo.vtxZErr"         ,vtxZErr    	);
        root->SetBranchAddress("BInfo.vtxYXErr"        ,vtxYXErr 	);
        root->SetBranchAddress("BInfo.vtxZXErr"        ,vtxZXErr    );
        root->SetBranchAddress("BInfo.vtxZYErr"        ,vtxZYErr  	);
        root->SetBranchAddress("BInfo.vtxdof"		    ,vtxdof		);
        root->SetBranchAddress("BInfo.vtxchi2"         ,vtxchi2   	);
        root->SetBranchAddress("BInfo.rfuj_index"      ,rfuj_index   	);
        root->SetBranchAddress("BInfo.rftk1_index"     ,rftk1_index   	);
        root->SetBranchAddress("BInfo.rftk2_index"     ,rftk2_index   	);
        root->SetBranchAddress("BInfo.type"            ,type   	);
        
        root->SetBranchAddress("BInfo.rfmu1_px"        ,rfmu1_px 	);
        root->SetBranchAddress("BInfo.rfmu1_py"        ,rfmu1_py	);
        root->SetBranchAddress("BInfo.rfmu1_pz"        ,rfmu1_pz     );
        root->SetBranchAddress("BInfo.rfmu2_px"        ,rfmu2_px     );
        root->SetBranchAddress("BInfo.rfmu2_py"        ,rfmu2_py    	);
        root->SetBranchAddress("BInfo.rfmu2_pz"        ,rfmu2_pz   	);
        root->SetBranchAddress("BInfo.rftk1_px"        ,rftk1_px 	);
        root->SetBranchAddress("BInfo.rftk1_py"        ,rftk1_py     );
        root->SetBranchAddress("BInfo.rftk1_pz"        ,rftk1_pz     );
        root->SetBranchAddress("BInfo.rftk2_px"        ,rftk2_px    	);
        root->SetBranchAddress("BInfo.rftk2_py"        ,rftk2_py   	);
        root->SetBranchAddress("BInfo.rftk2_pz"        ,rftk2_pz  	);
        
        root->SetBranchAddress("BInfo.tktk_mass"         ,tktk_mass    );
        root->SetBranchAddress("BInfo.tktk_pt"           ,tktk_pt     	);
        root->SetBranchAddress("BInfo.tktk_eta"           ,tktk_eta     	);
        root->SetBranchAddress("BInfo.tktk_phi"           ,tktk_phi     	);
        root->SetBranchAddress("BInfo.tktk_px"           ,tktk_px     	);
        root->SetBranchAddress("BInfo.tktk_py"           ,tktk_py    	);
        root->SetBranchAddress("BInfo.tktk_pz"           ,tktk_pz   	);
        root->SetBranchAddress("BInfo.tktk_vtxX"         ,tktk_vtxX    );
        root->SetBranchAddress("BInfo.tktk_vtxY"         ,tktk_vtxY    );
        root->SetBranchAddress("BInfo.tktk_vtxZ"         ,tktk_vtxZ    );
        root->SetBranchAddress("BInfo.tktk_vtxXErr"      ,tktk_vtxXErr  );
        root->SetBranchAddress("BInfo.tktk_vtxYErr"      ,tktk_vtxYErr  );
        root->SetBranchAddress("BInfo.tktk_vtxZErr"      ,tktk_vtxZErr  );
        root->SetBranchAddress("BInfo.tktk_vtxYXErr"     ,tktk_vtxYXErr );
        root->SetBranchAddress("BInfo.tktk_vtxZXErr"     ,tktk_vtxZXErr );
        root->SetBranchAddress("BInfo.tktk_vtxZYErr"     ,tktk_vtxZYErr );
        root->SetBranchAddress("BInfo.tktk_vtxdof"	      ,tktk_vtxdof	);
        root->SetBranchAddress("BInfo.tktk_vtxchi2"      ,tktk_vtxchi2 );
        root->SetBranchAddress("BInfo.tktk_rftk1_px"        ,tktk_rftk1_px 	);
        root->SetBranchAddress("BInfo.tktk_rftk1_py"        ,tktk_rftk1_py     );
        root->SetBranchAddress("BInfo.tktk_rftk1_pz"        ,tktk_rftk1_pz     );
        root->SetBranchAddress("BInfo.tktk_rftk2_px"        ,tktk_rftk2_px    	);
        root->SetBranchAddress("BInfo.tktk_rftk2_py"        ,tktk_rftk2_py   	);
        root->SetBranchAddress("BInfo.tktk_rftk2_pz"        ,tktk_rftk2_pz  	);
    }//}}}
};//}}}

class DInfoBranches{//{{{
public:
    
    int	    size;
    int	    index[MAX_XB];
    int     type[MAX_XB];

    double	tktkRes_mass[MAX_XB];
    double  tktkRes_pt[MAX_XB];
    double  tktkRes_eta[MAX_XB];
    double  tktkRes_phi[MAX_XB];
    double  tktkRes_px[MAX_XB];
    double  tktkRes_py[MAX_XB];
    double  tktkRes_pz[MAX_XB];
    double	tktkRes_vtxX[MAX_XB];
    double  tktkRes_vtxY[MAX_XB];
    double  tktkRes_vtxZ[MAX_XB];
    double  tktkRes_vtxXErr[MAX_XB];
    double  tktkRes_vtxYErr[MAX_XB];
    double  tktkRes_vtxZErr[MAX_XB];
    double  tktkRes_vtxYXErr[MAX_XB];
    double  tktkRes_vtxZXErr[MAX_XB];
    double  tktkRes_vtxZYErr[MAX_XB];
    double	tktkRes_vtxdof[MAX_XB];
    double	tktkRes_vtxchi2[MAX_XB];
    double  tktkRes_rftk1_mass[MAX_XB];
    double  tktkRes_rftk1_pt[MAX_XB];
    double  tktkRes_rftk1_eta[MAX_XB];
    double  tktkRes_rftk1_phi[MAX_XB];
    double  tktkRes_rftk2_mass[MAX_XB];
    double  tktkRes_rftk2_pt[MAX_XB];
    double  tktkRes_rftk2_eta[MAX_XB];
    double  tktkRes_rftk2_phi[MAX_XB];
    double  tktkRes_rftk3_mass[MAX_XB];
    double  tktkRes_rftk3_pt[MAX_XB];
    double  tktkRes_rftk3_eta[MAX_XB];
    double  tktkRes_rftk3_phi[MAX_XB];
    double  tktkRes_rftk4_mass[MAX_XB];
    double  tktkRes_rftk4_pt[MAX_XB];
    double  tktkRes_rftk4_eta[MAX_XB];
    double  tktkRes_rftk4_phi[MAX_XB];
    int     tktkRes_rftk1_index[MAX_XB];
    int     tktkRes_rftk2_index[MAX_XB];
    int     tktkRes_rftk3_index[MAX_XB];
    int     tktkRes_rftk4_index[MAX_XB];

    double	mass[MAX_XB];
    double  pt[MAX_XB];
    double  eta[MAX_XB];
    double  phi[MAX_XB];
    double  px[MAX_XB];
    double  py[MAX_XB];
    double  pz[MAX_XB];
    double  alpha[MAX_XB];
    double  svpvDistance[MAX_XB];
    double  svpvDisErr[MAX_XB];
    double  svpvDistance_2D[MAX_XB];
    double  svpvDisErr_2D[MAX_XB];
    double  MaxDoca[MAX_XB];
    double	vtxX[MAX_XB];
    double  vtxY[MAX_XB];
    double  vtxZ[MAX_XB];
    double  vtxXErr[MAX_XB];
    double  vtxYErr[MAX_XB];
    double  vtxZErr[MAX_XB];
    double  vtxYXErr[MAX_XB];
    double  vtxZXErr[MAX_XB];
    double  vtxZYErr[MAX_XB];
    double	vtxdof[MAX_XB];
    double	vtxchi2[MAX_XB];
    double  rftk1_px[MAX_XB];
    double  rftk1_py[MAX_XB];
    double  rftk1_pz[MAX_XB];
    double  rftk2_px[MAX_XB];
    double  rftk2_py[MAX_XB];
    double  rftk2_pz[MAX_XB];
    double  rftk3_px[MAX_XB];
    double  rftk3_py[MAX_XB];
    double  rftk3_pz[MAX_XB];
    double  rftk4_px[MAX_XB];
    double  rftk4_py[MAX_XB];
    double  rftk4_pz[MAX_XB];
    double  rftk5_px[MAX_XB];
    double  rftk5_py[MAX_XB];
    double  rftk5_pz[MAX_XB];
    double  rftk1_mass[MAX_XB];
    double  rftk1_pt[MAX_XB];
    double  rftk1_eta[MAX_XB];
    double  rftk1_phi[MAX_XB];
    double  rftk2_mass[MAX_XB];
    double  rftk2_pt[MAX_XB];
    double  rftk2_eta[MAX_XB];
    double  rftk2_phi[MAX_XB];
    double  rftk3_mass[MAX_XB];
    double  rftk3_pt[MAX_XB];
    double  rftk3_eta[MAX_XB];
    double  rftk3_phi[MAX_XB];
    double  rftk4_mass[MAX_XB];
    double  rftk4_pt[MAX_XB];
    double  rftk4_eta[MAX_XB];
    double  rftk4_phi[MAX_XB];
    double  rftk5_mass[MAX_XB];
    double  rftk5_pt[MAX_XB];
    double  rftk5_eta[MAX_XB];
    double  rftk5_phi[MAX_XB];
    int     rftk1_index[MAX_XB];
    int     rftk2_index[MAX_XB];
    int     rftk3_index[MAX_XB];
    int     rftk4_index[MAX_XB];
    int     rftk5_index[MAX_XB];
    int     rftk1_MassHypo[MAX_XB];
    int     rftk2_MassHypo[MAX_XB];
    int     rftk3_MassHypo[MAX_XB];
    int     rftk4_MassHypo[MAX_XB];
    int     rftk5_MassHypo[MAX_XB];
   
    void regTree(TTree *root){//{{{
        root->Branch("DInfo.size"             , &size          , "DInfo.size/I"			);
        root->Branch("DInfo.index"            , index          , "DInfo.index[DInfo.size]/I"		);
        root->Branch("DInfo.type"             , type           , "DInfo.type[DInfo.size]/I"	);
        
        root->Branch("DInfo.tktkRes_mass"          , tktkRes_mass              , "DInfo.tktkRes_mass[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_pt"            , tktkRes_pt                , "DInfo.tktkRes_pt[DInfo.size]/D"	);
        root->Branch("DInfo.tktkRes_eta"           , tktkRes_eta               , "DInfo.tktkRes_eta[DInfo.size]/D"	);
        root->Branch("DInfo.tktkRes_phi"           , tktkRes_phi               , "DInfo.tktkRes_phi[DInfo.size]/D"	);
        root->Branch("DInfo.tktkRes_px"            , tktkRes_px                , "DInfo.tktkRes_px[DInfo.size]/D"	);
        root->Branch("DInfo.tktkRes_py"            , tktkRes_py                , "DInfo.tktkRes_py[DInfo.size]/D"	);
        root->Branch("DInfo.tktkRes_pz"            , tktkRes_pz                , "DInfo.tktkRes_pz[DInfo.size]/D"	);
        root->Branch("DInfo.tktkRes_vtxX"          , tktkRes_vtxX              , "DInfo.tktkRes_vtxX[DInfo.size]/D"	);
        root->Branch("DInfo.tktkRes_vtxY"          , tktkRes_vtxY              , "DInfo.tktkRes_vtxY[DInfo.size]/D"	);
        root->Branch("DInfo.tktkRes_vtxZ"          , tktkRes_vtxZ              , "DInfo.tktkRes_vtxZ[DInfo.size]/D"	);
        root->Branch("DInfo.tktkRes_vtxXErr"       , tktkRes_vtxXErr           , "DInfo.tktkRes_vtxXErr[DInfo.size]/D"   );
        root->Branch("DInfo.tktkRes_vtxYErr"       , tktkRes_vtxYErr           , "DInfo.tktkRes_vtxYErr[DInfo.size]/D"   );
        root->Branch("DInfo.tktkRes_vtxZErr"       , tktkRes_vtxZErr           , "DInfo.tktkRes_vtxZErr[DInfo.size]/D"   );
        root->Branch("DInfo.tktkRes_vtxYXErr"      , tktkRes_vtxYXErr          , "DInfo.tktkRes_vtxYXErr[DInfo.size]/D"   );
        root->Branch("DInfo.tktkRes_vtxZXErr"      , tktkRes_vtxZXErr          , "DInfo.tktkRes_vtxZXErr[DInfo.size]/D"   );
        root->Branch("DInfo.tktkRes_vtxZYErr"      , tktkRes_vtxZYErr          , "DInfo.tktkRes_vtxZYErr[DInfo.size]/D"   );
        root->Branch("DInfo.tktkRes_vtxdof"        , tktkRes_vtxdof            , "DInfo.tktkRes_vtxdof[DInfo.size]/D"	);
        root->Branch("DInfo.tktkRes_vtxchi2"       , tktkRes_vtxchi2           , "DInfo.tktkRes_vtxchi2[DInfo.size]/D"	);
        root->Branch("DInfo.tktkRes_rftk1_mass"    , tktkRes_rftk1_mass        , "DInfo.tktkRes_rftk1_mass[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk1_pt"      , tktkRes_rftk1_pt          , "DInfo.tktkRes_rftk1_pt[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk1_eta"     , tktkRes_rftk1_eta         , "DInfo.tktkRes_rftk1_eta[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk1_phi"     , tktkRes_rftk1_phi         , "DInfo.tktkRes_rftk1_phi[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk2_mass"    , tktkRes_rftk2_mass        , "DInfo.tktkRes_rftk2_mass[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk2_pt"      , tktkRes_rftk2_pt          , "DInfo.tktkRes_rftk2_pt[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk2_eta"     , tktkRes_rftk2_eta         , "DInfo.tktkRes_rftk2_eta[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk2_phi"     , tktkRes_rftk2_phi         , "DInfo.tktkRes_rftk2_phi[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk3_mass"    , tktkRes_rftk3_mass        , "DInfo.tktkRes_rftk3_mass[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk3_pt"      , tktkRes_rftk3_pt          , "DInfo.tktkRes_rftk3_pt[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk3_eta"     , tktkRes_rftk3_eta         , "DInfo.tktkRes_rftk3_eta[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk3_phi"     , tktkRes_rftk3_phi         , "DInfo.tktkRes_rftk3_phi[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk4_mass"    , tktkRes_rftk4_mass        , "DInfo.tktkRes_rftk4_mass[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk4_pt"      , tktkRes_rftk4_pt          , "DInfo.tktkRes_rftk4_pt[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk4_eta"     , tktkRes_rftk4_eta         , "DInfo.tktkRes_rftk4_eta[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk4_phi"     , tktkRes_rftk4_phi         , "DInfo.tktkRes_rftk4_phi[DInfo.size]/D"     );
        root->Branch("DInfo.tktkRes_rftk1_index"   , tktkRes_rftk1_index       , "DInfo.tktkRes_rftk1_index[DInfo.size]/I");
        root->Branch("DInfo.tktkRes_rftk2_index"   , tktkRes_rftk2_index       , "DInfo.tktkRes_rftk2_index[DInfo.size]/I");
        root->Branch("DInfo.tktkRes_rftk3_index"   , tktkRes_rftk3_index       , "DInfo.tktkRes_rftk3_index[DInfo.size]/I");
        root->Branch("DInfo.tktkRes_rftk4_index"   , tktkRes_rftk4_index       , "DInfo.tktkRes_rftk4_index[DInfo.size]/I");

        root->Branch("DInfo.mass"             , mass              , "DInfo.mass[DInfo.size]/D"     );
        root->Branch("DInfo.pt"               , pt                , "DInfo.pt[DInfo.size]/D"	);
        root->Branch("DInfo.eta"              , eta               , "DInfo.eta[DInfo.size]/D"	);
        root->Branch("DInfo.phi"              , phi               , "DInfo.phi[DInfo.size]/D"	);
        root->Branch("DInfo.px"               , px                , "DInfo.px[DInfo.size]/D"	);
        root->Branch("DInfo.py"               , py                , "DInfo.py[DInfo.size]/D"	);
        root->Branch("DInfo.pz"               , pz                , "DInfo.pz[DInfo.size]/D"	);
        root->Branch("DInfo.alpha"            , alpha             , "DInfo.alpha[DInfo.size]/D"	);
        root->Branch("DInfo.svpvDistance"     , svpvDistance      , "DInfo.svpvDistance[DInfo.size]/D"	);
        root->Branch("DInfo.svpvDisErr"       , svpvDisErr        , "DInfo.svpvDisErr[DInfo.size]/D"	);
        root->Branch("DInfo.svpvDistance_2D"  , svpvDistance_2D   , "DInfo.svpvDistance_2D[DInfo.size]/D"	);
        root->Branch("DInfo.svpvDisErr_2D"    , svpvDisErr_2D     , "DInfo.svpvDisErr_2D[DInfo.size]/D"	);
        root->Branch("DInfo.MaxDoca"          , MaxDoca           , "DInfo.MaxDoca[DInfo.size]/D"	);
        root->Branch("DInfo.vtxX"             , vtxX              , "DInfo.vtxX[DInfo.size]/D"	);
        root->Branch("DInfo.vtxY"             , vtxY              , "DInfo.vtxY[DInfo.size]/D"	);
        root->Branch("DInfo.vtxZ"             , vtxZ              , "DInfo.vtxZ[DInfo.size]/D"	);
        root->Branch("DInfo.vtxXErr"          , vtxXErr           , "DInfo.vtxXErr[DInfo.size]/D"   );
        root->Branch("DInfo.vtxYErr"          , vtxYErr           , "DInfo.vtxYErr[DInfo.size]/D"   );
        root->Branch("DInfo.vtxZErr"          , vtxZErr           , "DInfo.vtxZErr[DInfo.size]/D"   );
        root->Branch("DInfo.vtxYXErr"         , vtxYXErr          , "DInfo.vtxYXErr[DInfo.size]/D"   );
        root->Branch("DInfo.vtxZXErr"         , vtxZXErr          , "DInfo.vtxZXErr[DInfo.size]/D"   );
        root->Branch("DInfo.vtxZYErr"         , vtxZYErr          , "DInfo.vtxZYErr[DInfo.size]/D"   );
        root->Branch("DInfo.vtxdof"           , vtxdof            , "DInfo.vtxdof[DInfo.size]/D"	);
        root->Branch("DInfo.vtxchi2"          , vtxchi2           , "DInfo.vtxchi2[DInfo.size]/D"	);
        root->Branch("DInfo.rftk1_px"         , rftk1_px          , "DInfo.rftk1_px[DInfo.size]/D"     );
        root->Branch("DInfo.rftk1_py"         , rftk1_py          , "DInfo.rftk1_py[DInfo.size]/D"     );
        root->Branch("DInfo.rftk1_pz"         , rftk1_pz          , "DInfo.rftk1_pz[DInfo.size]/D"     );
        root->Branch("DInfo.rftk2_px"         , rftk2_px          , "DInfo.rftk2_px[DInfo.size]/D"     );
        root->Branch("DInfo.rftk2_py"         , rftk2_py          , "DInfo.rftk2_py[DInfo.size]/D"     );
        root->Branch("DInfo.rftk2_pz"         , rftk2_pz          , "DInfo.rftk2_pz[DInfo.size]/D"     );
        root->Branch("DInfo.rftk3_px"         , rftk3_px          , "DInfo.rftk3_px[DInfo.size]/D"     );
        root->Branch("DInfo.rftk3_py"         , rftk3_py          , "DInfo.rftk3_py[DInfo.size]/D"     );
        root->Branch("DInfo.rftk3_pz"         , rftk3_pz          , "DInfo.rftk3_pz[DInfo.size]/D"     );
        root->Branch("DInfo.rftk4_px"         , rftk4_px          , "DInfo.rftk4_px[DInfo.size]/D"     );
        root->Branch("DInfo.rftk4_py"         , rftk4_py          , "DInfo.rftk4_py[DInfo.size]/D"     );
        root->Branch("DInfo.rftk4_pz"         , rftk4_pz          , "DInfo.rftk4_pz[DInfo.size]/D"     );
        root->Branch("DInfo.rftk5_px"         , rftk5_px          , "DInfo.rftk5_px[DInfo.size]/D"     );
        root->Branch("DInfo.rftk5_py"         , rftk5_py          , "DInfo.rftk5_py[DInfo.size]/D"     );
        root->Branch("DInfo.rftk5_pz"         , rftk5_pz          , "DInfo.rftk5_pz[DInfo.size]/D"     );

        root->Branch("DInfo.rftk1_mass"       , rftk1_mass        , "DInfo.rftk1_mass[DInfo.size]/D"     );
        root->Branch("DInfo.rftk1_pt"         , rftk1_pt          , "DInfo.rftk1_pt[DInfo.size]/D"     );
        root->Branch("DInfo.rftk1_eta"        , rftk1_eta         , "DInfo.rftk1_eta[DInfo.size]/D"     );
        root->Branch("DInfo.rftk1_phi"        , rftk1_phi         , "DInfo.rftk1_phi[DInfo.size]/D"     );
        root->Branch("DInfo.rftk2_mass"       , rftk2_mass        , "DInfo.rftk2_mass[DInfo.size]/D"     );
        root->Branch("DInfo.rftk2_pt"         , rftk2_pt          , "DInfo.rftk2_pt[DInfo.size]/D"     );
        root->Branch("DInfo.rftk2_eta"        , rftk2_eta         , "DInfo.rftk2_eta[DInfo.size]/D"     );
        root->Branch("DInfo.rftk2_phi"        , rftk2_phi         , "DInfo.rftk2_phi[DInfo.size]/D"     );
        root->Branch("DInfo.rftk3_mass"       , rftk3_mass        , "DInfo.rftk3_mass[DInfo.size]/D"     );
        root->Branch("DInfo.rftk3_pt"         , rftk3_pt          , "DInfo.rftk3_pt[DInfo.size]/D"     );
        root->Branch("DInfo.rftk3_eta"        , rftk3_eta         , "DInfo.rftk3_eta[DInfo.size]/D"     );
        root->Branch("DInfo.rftk3_phi"        , rftk3_phi         , "DInfo.rftk3_phi[DInfo.size]/D"     );
        root->Branch("DInfo.rftk4_mass"       , rftk4_mass        , "DInfo.rftk4_mass[DInfo.size]/D"     );
        root->Branch("DInfo.rftk4_pt"         , rftk4_pt          , "DInfo.rftk4_pt[DInfo.size]/D"     );
        root->Branch("DInfo.rftk4_eta"        , rftk4_eta         , "DInfo.rftk4_eta[DInfo.size]/D"     );
        root->Branch("DInfo.rftk4_phi"        , rftk4_phi         , "DInfo.rftk4_phi[DInfo.size]/D"     );
        root->Branch("DInfo.rftk5_mass"       , rftk5_mass        , "DInfo.rftk5_mass[DInfo.size]/D"     );
        root->Branch("DInfo.rftk5_pt"         , rftk5_pt          , "DInfo.rftk5_pt[DInfo.size]/D"     );
        root->Branch("DInfo.rftk5_eta"        , rftk5_eta         , "DInfo.rftk5_eta[DInfo.size]/D"     );
        root->Branch("DInfo.rftk5_phi"        , rftk5_phi         , "DInfo.rftk5_phi[DInfo.size]/D"     );
        root->Branch("DInfo.rftk1_index"      , rftk1_index       , "DInfo.rftk1_index[DInfo.size]/I");
        root->Branch("DInfo.rftk2_index"      , rftk2_index       , "DInfo.rftk2_index[DInfo.size]/I");
        root->Branch("DInfo.rftk3_index"      , rftk3_index       , "DInfo.rftk3_index[DInfo.size]/I");
        root->Branch("DInfo.rftk4_index"      , rftk4_index       , "DInfo.rftk4_index[DInfo.size]/I");
        root->Branch("DInfo.rftk5_index"      , rftk5_index       , "DInfo.rftk5_index[DInfo.size]/I");

        root->Branch("DInfo.rftk1_MassHypo"   , rftk1_MassHypo       , "DInfo.rftk1_MassHypo[DInfo.size]/I");
        root->Branch("DInfo.rftk2_MassHypo"   , rftk2_MassHypo       , "DInfo.rftk2_MassHypo[DInfo.size]/I");
        root->Branch("DInfo.rftk3_MassHypo"   , rftk3_MassHypo       , "DInfo.rftk3_MassHypo[DInfo.size]/I");
        root->Branch("DInfo.rftk4_MassHypo"   , rftk4_MassHypo       , "DInfo.rftk4_MassHypo[DInfo.size]/I");
        root->Branch("DInfo.rftk5_MassHypo"   , rftk5_MassHypo       , "DInfo.rftk5_MassHypo[DInfo.size]/I");
       
    }//}}}
    
    void setbranchadd(TTree *root){//{{{
        root->SetBranchAddress("DInfo.size"            ,&size        );
        root->SetBranchAddress("DInfo.index"           ,index       	);
        root->SetBranchAddress("DInfo.type"            ,type   	);
        
        root->SetBranchAddress("DInfo.tktkRes_mass"            ,tktkRes_mass    );
        root->SetBranchAddress("DInfo.tktkRes_pt"              ,tktkRes_pt     	);
        root->SetBranchAddress("DInfo.tktkRes_eta"             ,tktkRes_eta     	);
        root->SetBranchAddress("DInfo.tktkRes_phi"             ,tktkRes_phi     	);
        root->SetBranchAddress("DInfo.tktkRes_px"              ,tktkRes_px     	);
        root->SetBranchAddress("DInfo.tktkRes_py"              ,tktkRes_py    	);
        root->SetBranchAddress("DInfo.tktkRes_pz"              ,tktkRes_pz   	);
        root->SetBranchAddress("DInfo.tktkRes_vtxX"            ,tktkRes_vtxX    );
        root->SetBranchAddress("DInfo.tktkRes_vtxY"            ,tktkRes_vtxY    );
        root->SetBranchAddress("DInfo.tktkRes_vtxZ"            ,tktkRes_vtxZ    );
        root->SetBranchAddress("DInfo.tktkRes_vtxXErr"         ,tktkRes_vtxXErr   );
        root->SetBranchAddress("DInfo.tktkRes_vtxYErr"         ,tktkRes_vtxYErr   );
        root->SetBranchAddress("DInfo.tktkRes_vtxZErr"         ,tktkRes_vtxZErr   );
        root->SetBranchAddress("DInfo.tktkRes_vtxYXErr"        ,tktkRes_vtxYXErr   );
        root->SetBranchAddress("DInfo.tktkRes_vtxZXErr"        ,tktkRes_vtxZXErr   );
        root->SetBranchAddress("DInfo.tktkRes_vtxZYErr"        ,tktkRes_vtxZYErr   );
        root->SetBranchAddress("DInfo.tktkRes_vtxdof"          ,tktkRes_vtxdof	);
        root->SetBranchAddress("DInfo.tktkRes_vtxchi2"         ,tktkRes_vtxchi2 );
        root->SetBranchAddress("DInfo.tktkRes_rftk1_mass"      ,tktkRes_rftk1_mass  );
        root->SetBranchAddress("DInfo.tktkRes_rftk1_pt"        ,tktkRes_rftk1_pt  );
        root->SetBranchAddress("DInfo.tktkRes_rftk1_eta"       ,tktkRes_rftk1_eta  );
        root->SetBranchAddress("DInfo.tktkRes_rftk1_phi"       ,tktkRes_rftk1_phi  );
        root->SetBranchAddress("DInfo.tktkRes_rftk2_mass"      ,tktkRes_rftk2_mass  );
        root->SetBranchAddress("DInfo.tktkRes_rftk2_pt"        ,tktkRes_rftk2_pt  );
        root->SetBranchAddress("DInfo.tktkRes_rftk2_eta"       ,tktkRes_rftk2_eta  );
        root->SetBranchAddress("DInfo.tktkRes_rftk2_phi"       ,tktkRes_rftk2_phi  );
        root->SetBranchAddress("DInfo.tktkRes_rftk3_mass"      ,tktkRes_rftk3_mass  );
        root->SetBranchAddress("DInfo.tktkRes_rftk3_pt"        ,tktkRes_rftk3_pt  );
        root->SetBranchAddress("DInfo.tktkRes_rftk3_eta"       ,tktkRes_rftk3_eta  );
        root->SetBranchAddress("DInfo.tktkRes_rftk3_phi"       ,tktkRes_rftk3_phi  );
        root->SetBranchAddress("DInfo.tktkRes_rftk4_mass"      ,tktkRes_rftk4_mass  );
        root->SetBranchAddress("DInfo.tktkRes_rftk4_pt"        ,tktkRes_rftk4_pt  );
        root->SetBranchAddress("DInfo.tktkRes_rftk4_eta"       ,tktkRes_rftk4_eta  );
        root->SetBranchAddress("DInfo.tktkRes_rftk4_phi"       ,tktkRes_rftk4_phi  );
        root->SetBranchAddress("DInfo.tktkRes_rftk1_index"     ,tktkRes_rftk1_index   	);
        root->SetBranchAddress("DInfo.tktkRes_rftk2_index"     ,tktkRes_rftk2_index   	);
        root->SetBranchAddress("DInfo.tktkRes_rftk3_index"     ,tktkRes_rftk3_index   	);
        root->SetBranchAddress("DInfo.tktkRes_rftk4_index"     ,tktkRes_rftk4_index   	);

        root->SetBranchAddress("DInfo.mass"            ,mass    );
        root->SetBranchAddress("DInfo.pt"              ,pt     	);
        root->SetBranchAddress("DInfo.eta"             ,eta     	);
        root->SetBranchAddress("DInfo.phi"             ,phi     	);
        root->SetBranchAddress("DInfo.px"              ,px     	);
        root->SetBranchAddress("DInfo.py"              ,py    	);
        root->SetBranchAddress("DInfo.pz"              ,pz   	);
        root->SetBranchAddress("DInfo.alpha"           ,alpha   	);
        root->SetBranchAddress("DInfo.svpvDistance"    ,svpvDistance   	);
        root->SetBranchAddress("DInfo.svpvDisErr"      ,svpvDisErr   	);
        root->SetBranchAddress("DInfo.svpvDistance_2D" ,svpvDistance_2D   	);
        root->SetBranchAddress("DInfo.svpvDisErr_2D"   ,svpvDisErr_2D   	);
        root->SetBranchAddress("DInfo.MaxDoca"         ,MaxDoca   	);
        root->SetBranchAddress("DInfo.vtxX"            ,vtxX    );
        root->SetBranchAddress("DInfo.vtxY"            ,vtxY    );
        root->SetBranchAddress("DInfo.vtxZ"            ,vtxZ    );
        root->SetBranchAddress("DInfo.vtxXErr"         ,vtxXErr   );
        root->SetBranchAddress("DInfo.vtxYErr"         ,vtxYErr   );
        root->SetBranchAddress("DInfo.vtxZErr"         ,vtxZErr   );
        root->SetBranchAddress("DInfo.vtxYXErr"        ,vtxYXErr   );
        root->SetBranchAddress("DInfo.vtxZXErr"        ,vtxZXErr   );
        root->SetBranchAddress("DInfo.vtxZYErr"        ,vtxZYErr   );
        root->SetBranchAddress("DInfo.vtxdof"	       ,vtxdof	);
        root->SetBranchAddress("DInfo.vtxchi2"         ,vtxchi2 );
        root->SetBranchAddress("DInfo.rftk1_px"        ,rftk1_px 	);
        root->SetBranchAddress("DInfo.rftk1_py"        ,rftk1_py  );
        root->SetBranchAddress("DInfo.rftk1_pz"        ,rftk1_pz  );
        root->SetBranchAddress("DInfo.rftk2_px"        ,rftk2_px  );
        root->SetBranchAddress("DInfo.rftk2_py"        ,rftk2_py  );
        root->SetBranchAddress("DInfo.rftk2_pz"        ,rftk2_pz  );
        root->SetBranchAddress("DInfo.rftk3_px"        ,rftk3_px 	);
        root->SetBranchAddress("DInfo.rftk3_py"        ,rftk3_py  );
        root->SetBranchAddress("DInfo.rftk3_pz"        ,rftk3_pz  );
        root->SetBranchAddress("DInfo.rftk4_px"        ,rftk4_px  );
        root->SetBranchAddress("DInfo.rftk4_py"        ,rftk4_py  );
        root->SetBranchAddress("DInfo.rftk4_pz"        ,rftk4_pz  );
        root->SetBranchAddress("DInfo.rftk5_px"        ,rftk5_px  );
        root->SetBranchAddress("DInfo.rftk5_py"        ,rftk5_py  );
        root->SetBranchAddress("DInfo.rftk5_pz"        ,rftk5_pz  );

        root->SetBranchAddress("DInfo.rftk1_mass"      ,rftk1_mass);
        root->SetBranchAddress("DInfo.rftk1_pt"        ,rftk1_pt  );
        root->SetBranchAddress("DInfo.rftk1_eta"       ,rftk1_eta  );
        root->SetBranchAddress("DInfo.rftk1_phi"       ,rftk1_phi  );
        root->SetBranchAddress("DInfo.rftk2_mass"      ,rftk2_mass);
        root->SetBranchAddress("DInfo.rftk2_pt"        ,rftk2_pt  );
        root->SetBranchAddress("DInfo.rftk2_eta"       ,rftk2_eta  );
        root->SetBranchAddress("DInfo.rftk2_phi"       ,rftk2_phi  );
        root->SetBranchAddress("DInfo.rftk3_mass"      ,rftk3_mass);
        root->SetBranchAddress("DInfo.rftk3_pt"        ,rftk3_pt  );
        root->SetBranchAddress("DInfo.rftk3_eta"       ,rftk3_eta  );
        root->SetBranchAddress("DInfo.rftk3_phi"       ,rftk3_phi  );
        root->SetBranchAddress("DInfo.rftk4_mass"      ,rftk4_mass);
        root->SetBranchAddress("DInfo.rftk4_pt"        ,rftk4_pt  );
        root->SetBranchAddress("DInfo.rftk4_eta"       ,rftk4_eta  );
        root->SetBranchAddress("DInfo.rftk4_phi"       ,rftk4_phi  );
        root->SetBranchAddress("DInfo.rftk5_mass"      ,rftk5_mass);
        root->SetBranchAddress("DInfo.rftk5_pt"        ,rftk5_pt  );
        root->SetBranchAddress("DInfo.rftk5_eta"       ,rftk5_eta  );
        root->SetBranchAddress("DInfo.rftk5_phi"       ,rftk5_phi  );

        root->SetBranchAddress("DInfo.rftk1_index"     ,rftk1_index   	);
        root->SetBranchAddress("DInfo.rftk2_index"     ,rftk2_index   	);
        root->SetBranchAddress("DInfo.rftk3_index"     ,rftk3_index   	);
        root->SetBranchAddress("DInfo.rftk4_index"     ,rftk4_index   	);
        root->SetBranchAddress("DInfo.rftk5_index"     ,rftk5_index   	);

        root->SetBranchAddress("DInfo.rftk1_MassHypo"     ,rftk1_MassHypo   	);
        root->SetBranchAddress("DInfo.rftk2_MassHypo"     ,rftk2_MassHypo   	);
        root->SetBranchAddress("DInfo.rftk3_MassHypo"     ,rftk3_MassHypo   	);
        root->SetBranchAddress("DInfo.rftk4_MassHypo"     ,rftk4_MassHypo   	);
        root->SetBranchAddress("DInfo.rftk5_MassHypo"     ,rftk5_MassHypo   	);
    }//}}}
};//}}}

class GenInfoBranches{//{{{
    public:
        int     size;
        int     index       [MAX_GEN];
        int     handle_index[MAX_GEN];
        double  pt          [MAX_GEN];
        double  eta         [MAX_GEN];
        double  phi         [MAX_GEN];
        double  mass        [MAX_GEN];
        int     pdgId       [MAX_GEN];
        int     status      [MAX_GEN];
        int     nMo         [MAX_GEN];
        int     nDa         [MAX_GEN];
        int     mo1         [MAX_GEN];
        int     mo2         [MAX_GEN];
        int     da1         [MAX_GEN];
        int     da2         [MAX_GEN];
        int     da3         [MAX_GEN];
        int     da4         [MAX_GEN];

        void regTree(TTree *root){//{{{
            root->Branch("GenInfo.size"         ,&size          ,"GenInfo.size/I");
            root->Branch("GenInfo.index"        ,index          ,"GenInfo.index[GenInfo.size]/I");
            root->Branch("GenInfo.handle_index" ,handle_index   ,"GenInfo.handle_index[GenInfo.size]/I");
            root->Branch("GenInfo.pt"           ,pt             ,"GenInfo.pt[GenInfo.size]/D");
            root->Branch("GenInfo.eta"          ,eta            ,"GenInfo.eta[GenInfo.size]/D");
            root->Branch("GenInfo.phi"          ,phi            ,"GenInfo.phi[GenInfo.size]/D");
            root->Branch("GenInfo.mass"         ,mass           ,"GenInfo.mass[GenInfo.size]/D");
            root->Branch("GenInfo.pdgId"        ,pdgId          ,"GenInfo.pdgId[GenInfo.size]/I");
            root->Branch("GenInfo.status"       ,status         ,"GenInfo.status[GenInfo.size]/I");
            root->Branch("GenInfo.nMo"          ,nMo            ,"GenInfo.nMo[GenInfo.size]/I");
            root->Branch("GenInfo.nDa"          ,nDa            ,"GenInfo.nDa[GenInfo.size]/I");
            root->Branch("GenInfo.mo1"          ,mo1            ,"GenInfo.mo1[GenInfo.size]/I");
            root->Branch("GenInfo.mo2"          ,mo2            ,"GenInfo.mo2[GenInfo.size]/I");
            root->Branch("GenInfo.da1"          ,da1            ,"GenInfo.da1[GenInfo.size]/I");
            root->Branch("GenInfo.da2"          ,da2            ,"GenInfo.da2[GenInfo.size]/I");
            root->Branch("GenInfo.da3"          ,da3            ,"GenInfo.da3[GenInfo.size]/I");
            root->Branch("GenInfo.da4"          ,da4            ,"GenInfo.da4[GenInfo.size]/I");
        }//}}}

        void setbranchadd(TTree *root){//{{{
            root->SetBranchAddress("GenInfo.size"         ,&size          );
            root->SetBranchAddress("GenInfo.index"        ,index          );
            root->SetBranchAddress("GenInfo.handle_index" ,handle_index   );
            root->SetBranchAddress("GenInfo.pt"           ,pt             );
            root->SetBranchAddress("GenInfo.eta"          ,eta            );
            root->SetBranchAddress("GenInfo.phi"          ,phi            );
            root->SetBranchAddress("GenInfo.mass"         ,mass           );
            root->SetBranchAddress("GenInfo.pdgId"        ,pdgId          );
            root->SetBranchAddress("GenInfo.status"       ,status         );
            root->SetBranchAddress("GenInfo.nMo"          ,nMo            );
            root->SetBranchAddress("GenInfo.nDa"          ,nDa            );
            root->SetBranchAddress("GenInfo.mo1"          ,mo1            );
            root->SetBranchAddress("GenInfo.mo2"          ,mo2            );
            root->SetBranchAddress("GenInfo.da1"          ,da1            );
            root->SetBranchAddress("GenInfo.da2"          ,da2            );
            root->SetBranchAddress("GenInfo.da3"          ,da3            );
            root->SetBranchAddress("GenInfo.da4"          ,da4            );
        }//}}}
};//}}}

class BntupleBranches{//{{{
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
        
        //DInfo.trkInfo
        double   Btrk1Idx[MAX_XB];
        double   Btrk2Idx[MAX_XB];
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

        void makeNtuple(int ifchannel[], bool REAL, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, MuonInfoBranches *MuonInfo, TrackInfoBranches *TrackInfo, BInfoBranches *BInfo, GenInfoBranches *GenInfo, TTree* nt0, TTree* nt1, TTree* nt2, TTree* nt3, TTree* nt5, TTree* nt6){
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
        }

        void fillGenTree(TTree* ntGen, GenInfoBranches *GenInfo, bool gskim=true)
        {
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
        }


        void fillTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, int j, int typesize, double track_mass1, double track_mass2, bool REAL, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, MuonInfoBranches *MuonInfo, TrackInfoBranches *TrackInfo, BInfoBranches *BInfo, GenInfoBranches *GenInfo)
        {
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
            bP->SetXYZ(BInfo->px[j],BInfo->py[j],BInfo->pz[j]*0);
            bVtx->SetXYZ(BInfo->vtxX[j]-EvtInfo->PVx,
                    BInfo->vtxY[j]-EvtInfo->PVy,
                    BInfo->vtxZ[j]*0-EvtInfo->PVz*0);
            b4P->SetXYZM(BInfo->px[j],BInfo->py[j],BInfo->pz[j],BInfo->mass[j]);

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
            Blxy[typesize] = ((BInfo->vtxX[j]-EvtInfo->PVx)*BInfo->px[j] + (BInfo->vtxY[j]-EvtInfo->PVy)*BInfo->py[j])/BInfo->pt[j];
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

            double mu1px,mu1py,mu1pz,mu1E;
            mu1px = MuonInfo->pt[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]*cos(MuonInfo->phi[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]);
            mu1py = MuonInfo->pt[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]*sin(MuonInfo->phi[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]);
            mu1pz = MuonInfo->pt[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]*sinh(MuonInfo->eta[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]);
            b4P->SetXYZM(mu1px,mu1py,mu1pz,MUON_MASS);
            mu1E = b4P->E();
            Bmu1pt[typesize] = b4P->Pt();
            Bmu1p[typesize] = b4P->P();
            Bmu1eta[typesize] = b4P->Eta();
            Bmu1phi[typesize] = b4P->Phi();
            Bmu1y[typesize] = b4P->Rapidity();
            double mu2px,mu2py,mu2pz,mu2E;
            mu2px = MuonInfo->pt[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]*cos(MuonInfo->phi[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]);
            mu2py = MuonInfo->pt[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]*sin(MuonInfo->phi[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]);
            mu2pz = MuonInfo->pt[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]*sinh(MuonInfo->eta[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]);
            b4P->SetXYZM(mu2px,mu2py,mu2pz,MUON_MASS);
            mu2E = b4P->E();
            Bmu2pt[typesize] = b4P->Pt();
            Bmu2p[typesize] = b4P->P();
            Bmu2eta[typesize] = b4P->Eta();
            Bmu2phi[typesize] = b4P->Phi();
            Bmu2y[typesize] = b4P->Rapidity();

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
            b4P->SetXYZM(BInfo->uj_px[BInfo->rfuj_index[j]],
                    BInfo->uj_py[BInfo->rfuj_index[j]],
                    BInfo->uj_pz[BInfo->rfuj_index[j]],
                    BInfo->uj_mass[BInfo->rfuj_index[j]]);
            Bujpt[typesize] = b4P->Pt();
            Bujeta[typesize] = b4P->PseudoRapidity();
            Bujphi[typesize] = b4P->Phi();
            Bujy[typesize] = b4P->Rapidity();
            Bujlxy[typesize] = ((BInfo->uj_vtxX[BInfo->rfuj_index[j]]-EvtInfo->PVx)*BInfo->uj_px[BInfo->rfuj_index[j]] + (BInfo->uj_vtxY[BInfo->rfuj_index[j]]-EvtInfo->PVy)*BInfo->uj_py[BInfo->rfuj_index[j]])/Bujpt[typesize];
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
                b4P->SetXYZM(BInfo->tktk_px[j],BInfo->tktk_py[j],BInfo->tktk_pz[j],BInfo->tktk_mass[j]);
                Bdoubletpt[typesize] = b4P->Pt();
                Bdoubleteta[typesize] = b4P->PseudoRapidity();
                Bdoubletphi[typesize] = b4P->Phi();
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
                b4P->SetXYZM(BInfo->tktk_px[j],BInfo->tktk_py[j],BInfo->tktk_pz[j],BInfo->tktk_mass[j]);
                Bdoubletpt[typesize] = b4P->Pt();
                Bdoubleteta[typesize] = b4P->PseudoRapidity();
                Bdoubletphi[typesize] = b4P->Phi();
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
        }

        bool signalGen(int Btype, int j, GenInfoBranches *GenInfo)
        {
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
        }

};//}}}

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
        double   Dtrk1Idx[MAX_XB];
        double   Dtrk2Idx[MAX_XB];
        double   Dtrk3Idx[MAX_XB];
        double   Dtrk4Idx[MAX_XB];
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

        void makeDNtuple(int isDchannel[], bool REAL, bool skim, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, TrackInfoBranches *TrackInfo, DInfoBranches *DInfo, GenInfoBranches *GenInfo, TTree* ntD1, TTree* ntD2, TTree* ntD3){
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
        }

        void fillDGenTree(TTree* ntGen, GenInfoBranches *GenInfo, bool gskim=true){
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
        }

        double findMass(int particlePdgId)
        {
            if(TMath::Abs(particlePdgId)==211) return PION_MASS;
            if(TMath::Abs(particlePdgId)==321) return KAON_MASS;
            else
            {
                std::cout<<"ERROR: find particle mass falied >> Particle pdgId: "<<particlePdgId<<std::endl;
                return 0;
            }
        }

        void fillDTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, int j, int typesize, bool REAL, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, TrackInfoBranches *TrackInfo, DInfoBranches *DInfo, GenInfoBranches *GenInfo)
        {
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
            bP->SetXYZ(DInfo->px[j],DInfo->py[j],DInfo->pz[j]);
            bVtx->SetXYZ(DInfo->vtxX[j]-EvtInfo->PVx,
                    DInfo->vtxY[j]-EvtInfo->PVy,
                    DInfo->vtxZ[j]-EvtInfo->PVz);
            b4P->SetXYZM(DInfo->px[j],DInfo->py[j],DInfo->pz[j],DInfo->mass[j]);
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
            Dlxy[typesize] = ((DInfo->vtxX[j]-EvtInfo->PVx)*DInfo->px[j] + (DInfo->vtxY[j]-EvtInfo->PVy)*DInfo->py[j])/DInfo->pt[j];
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
        }//fillDtree

        bool isDsignalGen(int dmesontype, int j, GenInfoBranches *GenInfo)
        {
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
        }
};//}}}

#endif
