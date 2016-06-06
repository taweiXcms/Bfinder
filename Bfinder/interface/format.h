// vim:set ts=4 sw=4 fdm=marker et:
#ifndef _XBFRAMEFORMAT_H_
#define _XBFRAMEFORMAT_H_

//Note, when the array size gett too large, SetBranchAddress will fail, root will abort w/o error msg
#define MAX_XB       20000
#define MAX_MUON     10000
#define MAX_TRACK    6000
#define MAX_GEN      6000
#define MAX_BX       150
#define MAX_Vertices 4000
//
//#define MAX_XB 16384
//#define MAX_MUON 4096
//#define MAX_TRACK 8192
//#define MAX_GEN 4096
//#define MAX_BX 128
//#define MAX_Vertices 4096
//
//#define MAX_XB 8192
//#define MAX_MUON 2048
//#define MAX_TRACK 4096
//#define MAX_GEN 4096
//#define MAX_BX 128
//#define MAX_Vertices 4096
//
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
#define DPLUS_MASS 1.8696
#define DSUBS_MASS 1.9683

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
#define DSTAR_PDGID 413

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
		float	PVx;
		float	PVy;
		float	PVz;
		float	PVxE;
		float	PVyE;
		float	PVzE;
		float	PVnchi2;
		float	PVchi2;
        float  BSx;
        float  BSy;
        float  BSz;
        float  BSxErr;
        float  BSyErr;
        float  BSzErr;
        float  BSdxdz;
        float  BSdydz;
        float  BSdxdzErr;
        float  BSdydzErr;
        float  BSWidthX;
        float  BSWidthXErr;
        float  BSWidthY;
        float  BSWidthYErr;
		//float	PVc2p;
		
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
  		    root->Branch("EvtInfo.PVx"          , &PVx                       , "EvtInfo.PVx/F"			);
  		    root->Branch("EvtInfo.PVy"          , &PVy                       , "EvtInfo.PVy/F"			);
  		    root->Branch("EvtInfo.PVz"          , &PVz                       , "EvtInfo.PVz/F"			);
            root->Branch("EvtInfo.PVxE"         , &PVxE                      , "EvtInfo.PVxE/F"           );
            root->Branch("EvtInfo.PVyE"         , &PVyE                      , "EvtInfo.PVyE/F"           );
            root->Branch("EvtInfo.PVzE"         , &PVzE                      , "EvtInfo.PVzE/F"           );
  		    root->Branch("EvtInfo.PVnchi2"      , &PVnchi2                   , "EvtInfo.PVnchi2/F"		);
  		    root->Branch("EvtInfo.PVchi2"       , &PVchi2                    , "EvtInfo.PVchi2/F"			);
  		    root->Branch("EvtInfo.BSx"          , &BSx                       , "EvtInfo.BSx/F"			);
  		    root->Branch("EvtInfo.BSy"          , &BSy                       , "EvtInfo.BSy/F"			);
  		    root->Branch("EvtInfo.BSz"          , &BSz                       , "EvtInfo.BSz/F"			);
  		    root->Branch("EvtInfo.BSxErr"       , &BSxErr                    , "EvtInfo.BSxErr/F"			);
  		    root->Branch("EvtInfo.BSyErr"       , &BSyErr                    , "EvtInfo.BSyErr/F"			);
  		    root->Branch("EvtInfo.BSzErr"       , &BSzErr                    , "EvtInfo.BSzErr/F"			);
  		    root->Branch("EvtInfo.BSdxdz"       , &BSdxdz                    , "EvtInfo.BSdxdz/F"			);
  		    root->Branch("EvtInfo.BSdydz"       , &BSdydz                    , "EvtInfo.BSdydz/F"			);
  		    root->Branch("EvtInfo.BSdxdzErr"    , &BSdxdzErr                 , "EvtInfo.BSdxdzErr/F"	);
  		    root->Branch("EvtInfo.BSdydzErr"    , &BSdydzErr                 , "EvtInfo.BSdydzErr/F"	);
  		    root->Branch("EvtInfo.BSWidthX"     , &BSWidthX                  , "EvtInfo.BSWidthX/F"		);
  		    root->Branch("EvtInfo.BSWidthXErr"  , &BSWidthXErr               , "EvtInfo.BSWidthXErr/F"	);
  		    root->Branch("EvtInfo.BSWidthY"     , &BSWidthY                  , "EvtInfo.BSWidthY/F"		);
  		    root->Branch("EvtInfo.BSWidthYErr"  , &BSWidthYErr               , "EvtInfo.BSWidthYErr/F"	);
			//root->Branch("EvtInfo.PVc2p"      , &PVc2p                     , "EvtInfo.PVc2p/F"			);//
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
        float 	pt           [ MAX_MUON];
        float	eta          [ MAX_MUON];
        float 	phi          [ MAX_MUON];
        float 	ptErr        [ MAX_MUON];
        float	etaErr       [ MAX_MUON];
        float 	phiErr       [ MAX_MUON];
        bool    isTrackerMuon[ MAX_MUON];
        bool    isGlobalMuon [ MAX_MUON];
        int	    muqual       [ MAX_MUON];
        float  iso_trk      [ MAX_MUON];
        float  iso_ecal     [ MAX_MUON];
        float  iso_hcal     [ MAX_MUON];
        int     type         [ MAX_MUON];
        float  n_matches    [ MAX_MUON];
        bool    TMOneStationTight [MAX_MUON];
        bool    TrackerMuonArbitrated [MAX_MUON];
        bool    isSoftMuon [MAX_MUON];
        int     geninfo_index[ MAX_MUON];
        bool    isNeededMuon[MAX_MUON];//for intermediate Bfinder usage, not stored in output
        bool    BfinderMuID [MAX_MUON];
        bool    SoftMuID [MAX_MUON];

        bool    isStandAloneMuon            [ MAX_MUON];
        int 	StandAloneMuon_charge       [ MAX_MUON];
        float 	StandAloneMuon_pt           [ MAX_MUON];
        float	StandAloneMuon_eta          [ MAX_MUON];
        float 	StandAloneMuon_phi          [ MAX_MUON];
        float	StandAloneMuon_d0           [ MAX_MUON];
        float	StandAloneMuon_dz           [ MAX_MUON];
        float	StandAloneMuon_dzPV         [ MAX_MUON];
        float	StandAloneMuon_dxyPV        [ MAX_MUON];

        bool    outerTrackisNonnull  [MAX_MUON];
        bool    innerTrackisNonnull  [MAX_MUON];
        bool    globalTrackisNonnull [MAX_MUON];
        int     innerTrackQuality    [MAX_MUON];
        float  normchi2     [ MAX_MUON];
        int	    i_striphit   [ MAX_MUON];
        int	    i_pixelhit   [ MAX_MUON];
        int     i_nStripLayer[ MAX_MUON];
        int     i_nPixelLayer[ MAX_MUON];
        float	i_chi2       [ MAX_MUON];
        float	i_ndf        [ MAX_MUON];
        int	    fpbarrelhit  [ MAX_MUON];
        int	    fpendcaphit  [ MAX_MUON];
        float	d0           [ MAX_MUON];
        float	dz           [ MAX_MUON];
        float	dzPV         [ MAX_MUON];
        float	dxyPV        [ MAX_MUON];
        float	g_chi2       [ MAX_MUON];
        float	g_ndf        [ MAX_MUON];
        int	    g_striphit   [ MAX_MUON];
        int	    g_pixelhit   [ MAX_MUON];
        int	    nmuhit       [ MAX_MUON];

        bool    isTriggered  [MAX_MUON];
        int	    MuTrgMatchPathSize;
        std::vector<std::vector<double > > *MuTrgMatchTrgObjE;
        std::vector<std::vector<double > > *MuTrgMatchTrgObjPt;
        std::vector<std::vector<double > > *MuTrgMatchTrgObjEta;
        std::vector<std::vector<double > > *MuTrgMatchTrgObjPhi;

        void regTree(TTree *root, bool detailMode = false){//{{{
            root->Branch("MuonInfo.size"          , &size         , "MuonInfo.size/I"			);
            root->Branch("MuonInfo.index"         , index         , "MuonInfo.index[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.handle_index"  , handle_index  , "MuonInfo.handle_index[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.charge"        , charge        , "MuonInfo.charge[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.pt"            , pt            , "MuonInfo.pt[MuonInfo.size]/F"		);
            root->Branch("MuonInfo.eta"           , eta           , "MuonInfo.eta[MuonInfo.size]/F"	);
            root->Branch("MuonInfo.phi"           , phi           , "MuonInfo.phi[MuonInfo.size]/F"	);
            root->Branch("MuonInfo.ptErr"        , ptErr         , "MuonInfo.ptErr[MuonInfo.size]/F"		);
            root->Branch("MuonInfo.etaErr"        , etaErr        , "MuonInfo.etaErr[MuonInfo.size]/F"	);
            root->Branch("MuonInfo.phiErr"        , phiErr        , "MuonInfo.phiErr[MuonInfo.size]/F"	);
            root->Branch("MuonInfo.isTrackerMuon" , isTrackerMuon , "MuonInfo.isTrackerMuon[MuonInfo.size]/O");
            root->Branch("MuonInfo.isGlobalMuon"  , isGlobalMuon  , "MuonInfo.isGlobalMuon[MuonInfo.size]/O");
            root->Branch("MuonInfo.muqual"        , muqual        , "MuonInfo.muqual[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.type"          , type         ,  "MuonInfo.type[MuonInfo.size]/I"   );
            root->Branch("MuonInfo.n_matches"     , n_matches     , "MuonInfo.n_matches[MuonInfo.size]/F");
            root->Branch("MuonInfo.TMOneStationTight" ,TMOneStationTight, "MuonInfo.TMOneStationTight[MuonInfo.size]/O");
            root->Branch("MuonInfo.TrackerMuonArbitrated" ,TrackerMuonArbitrated, "MuonInfo.TrackerMuonArbitrated[MuonInfo.size]/O");
            root->Branch("MuonInfo.isSoftMuon" ,isSoftMuon, "MuonInfo.isSoftMuon[MuonInfo.size]/O");
            root->Branch("MuonInfo.geninfo_index"    , geninfo_index    , "MuonInfo.geninfo_index[MuonInfo.size]/I");
            root->Branch("MuonInfo.BfinderMuID" ,BfinderMuID, "MuonInfo.BfinderMuID[MuonInfo.size]/O");
            root->Branch("MuonInfo.SoftMuID" ,SoftMuID, "MuonInfo.SoftMuID[MuonInfo.size]/O");

            root->Branch("MuonInfo.innerTrackQuality"    , innerTrackQuality    , "MuonInfo.innerTrackQuality[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.normchi2"      , normchi2      , "MuonInfo.normchi2[MuonInfo.size]/F");
            root->Branch("MuonInfo.i_striphit"    , i_striphit    , "MuonInfo.i_striphit[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.i_pixelhit"    , i_pixelhit    , "MuonInfo.i_pixelhit[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.i_nStripLayer" , i_nStripLayer , "MuonInfo.i_nStripLayer[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.i_nPixelLayer" , i_nPixelLayer , "MuonInfo.i_nPixelLayer[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.i_chi2"        , i_chi2        , "MuonInfo.i_chi2[MuonInfo.size]/F"	);
            root->Branch("MuonInfo.i_ndf"         , i_ndf         , "MuonInfo.i_ndf[MuonInfo.size]/F"	);
            root->Branch("MuonInfo.fpbarrelhit"   , fpbarrelhit   , "MuonInfo.fpbarrelhit[MuonInfo.size]/I");
            root->Branch("MuonInfo.fpendcaphit"   , fpendcaphit   , "MuonInfo.fpendcaphit[MuonInfo.size]/I");
            root->Branch("MuonInfo.d0"            , d0            , "MuonInfo.d0[MuonInfo.size]/F"		);
            root->Branch("MuonInfo.dz"            , dz            , "MuonInfo.dz[MuonInfo.size]/F"		);
            root->Branch("MuonInfo.dzPV"          , dzPV          , "MuonInfo.dzPV[MuonInfo.size]/F"		);
            root->Branch("MuonInfo.dxyPV"         , dxyPV         , "MuonInfo.dxyPV[MuonInfo.size]/F"		);

            root->Branch("MuonInfo.isTriggered" ,isTriggered, "MuonInfo.isTriggered[MuonInfo.size]/O");
            root->Branch("MuonInfo.MuTrgMatchPathSize", &MuTrgMatchPathSize, "MuonInfo.MuTrgMatchPathSize/I");
            root->Branch("MuonInfo.MuTrgMatchTrgObjE",   "std::vector<std::vector<double>>", &MuTrgMatchTrgObjE);
            root->Branch("MuonInfo.MuTrgMatchTrgObjPt",  "std::vector<std::vector<double>>", &MuTrgMatchTrgObjPt);
            root->Branch("MuonInfo.MuTrgMatchTrgObjEta", "std::vector<std::vector<double>>", &MuTrgMatchTrgObjEta);
            root->Branch("MuonInfo.MuTrgMatchTrgObjPhi", "std::vector<std::vector<double>>", &MuTrgMatchTrgObjPhi);

            if(detailMode){
                root->Branch("MuonInfo.iso_trk"       , iso_trk       , "MuonInfo.iso_trk[MuonInfo.size]/F");
                root->Branch("MuonInfo.iso_ecal"      , iso_ecal      , "MuonInfo.iso_ecal[MuonInfo.size]/F");
                root->Branch("MuonInfo.iso_hcal"      , iso_hcal      , "MuonInfo.iso_hcal[MuonInfo.size]/F");
    
                root->Branch("MuonInfo.isStandAloneMuon"             , isStandAloneMuon             , "MuonInfo.isStandAloneMuon[MuonInfo.size]/O"		);
                root->Branch("MuonInfo.StandAloneMuon_charge"        , StandAloneMuon_charge        , "MuonInfo.StandAloneMuon_charge[MuonInfo.size]/I"	);
                root->Branch("MuonInfo.StandAloneMuon_pt"            , StandAloneMuon_pt            , "MuonInfo.StandAloneMuon_pt[MuonInfo.size]/F"		);
                root->Branch("MuonInfo.StandAloneMuon_eta"           , StandAloneMuon_eta           , "MuonInfo.StandAloneMuon_eta[MuonInfo.size]/F"	);
                root->Branch("MuonInfo.StandAloneMuon_phi"           , StandAloneMuon_phi           , "MuonInfo.StandAloneMuon_phi[MuonInfo.size]/F"	);
                root->Branch("MuonInfo.StandAloneMuon_d0"            , StandAloneMuon_d0            , "MuonInfo.StandAloneMuon_d0[MuonInfo.size]/F"		);
                root->Branch("MuonInfo.StandAloneMuon_dz"            , StandAloneMuon_dz            , "MuonInfo.StandAloneMuon_dz[MuonInfo.size]/F"		);
                root->Branch("MuonInfo.StandAloneMuon_dzPV"          , StandAloneMuon_dzPV          , "MuonInfo.StandAloneMuon_dzPV[MuonInfo.size]/F"   );
                root->Branch("MuonInfo.StandAloneMuon_dxyPV"         , StandAloneMuon_dxyPV         , "MuonInfo.StandAloneMuon_dxyPV[MuonInfo.size]/F"	);
    
                root->Branch("MuonInfo.g_chi2"        , g_chi2        , "MuonInfo.g_chi2[MuonInfo.size]/F"	);
                root->Branch("MuonInfo.g_ndf"         , g_ndf         , "MuonInfo.g_ndf[MuonInfo.size]/F"	);
                root->Branch("MuonInfo.g_striphit"    , g_striphit    , "MuonInfo.g_striphit[MuonInfo.size]/I"	);
                root->Branch("MuonInfo.g_pixelhit"    , g_pixelhit    , "MuonInfo.g_pixelhit[MuonInfo.size]/I"	);
                root->Branch("MuonInfo.nmuhit"        , nmuhit        , "MuonInfo.nmuhit[MuonInfo.size]/I"	);

                root->Branch("MuonInfo.outerTrackisNonnull" ,outerTrackisNonnull, "MuonInfo.outerTrackisNonnull[MuonInfo.size]/O");
                root->Branch("MuonInfo.innerTrackisNonnull" ,innerTrackisNonnull, "MuonInfo.innerTrackisNonnull[MuonInfo.size]/O");
                root->Branch("MuonInfo.globalTrackisNonnull" ,globalTrackisNonnull, "MuonInfo.globalTrackisNonnull[MuonInfo.size]/O");
            }
        }//}}}

        void setbranchadd(TTree *root, bool detailMode = false){//{{{
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
            root->SetBranchAddress("MuonInfo.type"          , type          );
            root->SetBranchAddress("MuonInfo.n_matches"     , n_matches);
            root->SetBranchAddress("MuonInfo.TMOneStationTight" , TMOneStationTight);
            root->SetBranchAddress("MuonInfo.TrackerMuonArbitrated" , TrackerMuonArbitrated);
            root->SetBranchAddress("MuonInfo.isSoftMuon" , isSoftMuon);
            root->SetBranchAddress("MuonInfo.geninfo_index"    , geninfo_index);
            root->SetBranchAddress("MuonInfo.BfinderMuID" , BfinderMuID);
            root->SetBranchAddress("MuonInfo.SoftMuID" , SoftMuID);

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

            root->SetBranchAddress("MuonInfo.isTriggered" , isTriggered);
            MuTrgMatchTrgObjE = new  std::vector<std::vector<double > >();
            MuTrgMatchTrgObjPt= new  std::vector<std::vector<double > >();                                                                                                                          
            MuTrgMatchTrgObjEta= new std::vector<std::vector<double > >();
            MuTrgMatchTrgObjPhi= new std::vector<std::vector<double > >();
            root->SetBranchAddress("MuonInfo.MuTrgMatchPathSize", &MuTrgMatchPathSize);
            root->SetBranchAddress("MuonInfo.MuTrgMatchTrgObjE", &MuTrgMatchTrgObjE);
            root->SetBranchAddress("MuonInfo.MuTrgMatchTrgObjPt", &MuTrgMatchTrgObjPt);
            root->SetBranchAddress("MuonInfo.MuTrgMatchTrgObjEta", &MuTrgMatchTrgObjEta);
            root->SetBranchAddress("MuonInfo.MuTrgMatchTrgObjPhi", &MuTrgMatchTrgObjPhi);

            if(detailMode){
                root->SetBranchAddress("MuonInfo.iso_trk"       , iso_trk);
                root->SetBranchAddress("MuonInfo.iso_ecal"      , iso_ecal);
                root->SetBranchAddress("MuonInfo.iso_hcal"      , iso_hcal);
    
                root->SetBranchAddress("MuonInfo.isStandAloneMuon"             , isStandAloneMuon              );
                root->SetBranchAddress("MuonInfo.StandAloneMuon_charge"        , StandAloneMuon_charge         );
                root->SetBranchAddress("MuonInfo.StandAloneMuon_pt"            , StandAloneMuon_pt             );
                root->SetBranchAddress("MuonInfo.StandAloneMuon_eta"           , StandAloneMuon_eta            );
                root->SetBranchAddress("MuonInfo.StandAloneMuon_phi"           , StandAloneMuon_phi            );
                root->SetBranchAddress("MuonInfo.StandAloneMuon_d0"            , StandAloneMuon_d0             );
                root->SetBranchAddress("MuonInfo.StandAloneMuon_dz"            , StandAloneMuon_dz             );
                root->SetBranchAddress("MuonInfo.StandAloneMuon_dzPV"          , StandAloneMuon_dzPV           );
                root->SetBranchAddress("MuonInfo.StandAloneMuon_dxyPV"         , StandAloneMuon_dxyPV          );
    
                root->SetBranchAddress("MuonInfo.g_chi2"        , g_chi2		);
                root->SetBranchAddress("MuonInfo.g_ndf"         , g_ndf		);
                root->SetBranchAddress("MuonInfo.g_striphit"    , g_striphit	);
                root->SetBranchAddress("MuonInfo.g_pixelhit"    , g_pixelhit	);
                root->SetBranchAddress("MuonInfo.nmuhit"        , nmuhit		);

                root->SetBranchAddress("MuonInfo.outerTrackisNonnull"    , outerTrackisNonnull);
                root->SetBranchAddress("MuonInfo.innerTrackisNonnull"    , innerTrackisNonnull);
                root->SetBranchAddress("MuonInfo.globalTrackisNonnull"    , globalTrackisNonnull);
            }
        }//}}}
};//}}}

class TrackInfoBranches{//{{{
    public:
        int     size;
        int     index        [ MAX_TRACK];
        int     handle_index [ MAX_TRACK];
        int     charge       [ MAX_TRACK];
        float  pt           [ MAX_TRACK];
        float  eta          [ MAX_TRACK];
        float  phi          [ MAX_TRACK];
        float  ptErr        [ MAX_TRACK];
        float  etaErr       [ MAX_TRACK];
        float  phiErr       [ MAX_TRACK];
        //float  p            [ MAX_TRACK];
        int     striphit     [ MAX_TRACK];
        int     pixelhit     [ MAX_TRACK];
        int     nStripLayer  [ MAX_TRACK];
        int     nPixelLayer  [ MAX_TRACK];
        int	    fpbarrelhit  [ MAX_TRACK];
        int	    fpendcaphit  [ MAX_TRACK];
        float	chi2         [ MAX_TRACK];
        float	ndf          [ MAX_TRACK];
        float	d0           [ MAX_TRACK];
        float	d0error      [ MAX_TRACK];
        float	dzPV         [ MAX_TRACK];
        float	dxyPV        [ MAX_TRACK];
        int     geninfo_index[ MAX_TRACK];
        int     trackQuality [ MAX_TRACK];
        bool    highPurity   [ MAX_TRACK];
        float   trkMVAVal    [ MAX_TRACK];
        int     trkAlgo      [ MAX_TRACK];
        int   originalTrkAlgo[ MAX_TRACK];

        void regTree(TTree *root, bool detailMode = false){//{{{
            root->Branch("TrackInfo.size"           ,&size		    ,"TrackInfo.size/I"			);
            root->Branch("TrackInfo.index"          ,index          ,"TrackInfo.index[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.handle_index"   ,handle_index   ,"TrackInfo.handle_index[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.charge"  	    ,charge         ,"TrackInfo.charge[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.pt"             ,pt             ,"TrackInfo.pt[TrackInfo.size]/F"	);
            root->Branch("TrackInfo.eta"            ,eta            ,"TrackInfo.eta[TrackInfo.size]/F"	);
            root->Branch("TrackInfo.phi"            ,phi            ,"TrackInfo.phi[TrackInfo.size]/F"	);
            root->Branch("TrackInfo.ptErr"          ,ptErr          ,"TrackInfo.ptErr[TrackInfo.size]/F"	);
            root->Branch("TrackInfo.etaErr"         ,etaErr         ,"TrackInfo.etaErr[TrackInfo.size]/F"	);
            root->Branch("TrackInfo.phiErr"         ,phiErr         ,"TrackInfo.phiErr[TrackInfo.size]/F"	);
            root->Branch("TrackInfo.striphit"	    ,striphit	    ,"TrackInfo.striphit[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.pixelhit"	    ,pixelhit	    ,"TrackInfo.pixelhit[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.nStripLayer"	,nStripLayer	,"TrackInfo.nStripLayer[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.nPixelLayer"	,nPixelLayer	,"TrackInfo.nPixelLayer[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.fpbarrelhit"	,fpbarrelhit	,"TrackInfo.fpbarrelhit[TrackInfo.size]/I");
            root->Branch("TrackInfo.fpendcaphit"	,fpendcaphit	,"TrackInfo.fpendcaphit[TrackInfo.size]/I");
            root->Branch("TrackInfo.chi2"		    ,chi2		    ,"TrackInfo.chi2[TrackInfo.size]/F"	);
            root->Branch("TrackInfo.ndf"		    ,ndf		    ,"TrackInfo.ndf[TrackInfo.size]/F"	);
            root->Branch("TrackInfo.d0"		        ,d0		        ,"TrackInfo.d0[TrackInfo.size]/F"	);
            root->Branch("TrackInfo.d0error"	    ,d0error	    ,"TrackInfo.d0error[TrackInfo.size]/F"	);
            root->Branch("TrackInfo.dzPV"           ,dzPV           ,"TrackInfo.dzPV[TrackInfo.size]/F"		);
            root->Branch("TrackInfo.dxyPV"          ,dxyPV          ,"TrackInfo.dxyPV[TrackInfo.size]/F"		);
            root->Branch("TrackInfo.geninfo_index"  ,geninfo_index  ,"TrackInfo.geninfo_index[TrackInfo.size]/I");
            root->Branch("TrackInfo.trackQuality"   ,trackQuality   ,"TrackInfo.trackQuality[TrackInfo.size]/I");
            root->Branch("TrackInfo.highPurity"     ,highPurity     ,"TrackInfo.highPurity[TrackInfo.size]/O");
            root->Branch("TrackInfo.trkMVAVal"      ,trkMVAVal      ,"TrackInfo.trkMVAVal[TrackInfo.size]/F");
            root->Branch("TrackInfo.trkAlgo"        ,trkAlgo        ,"TrackInfo.trkAlgo[TrackInfo.size]/I");
            root->Branch("TrackInfo.originalTrkAlgo",originalTrkAlgo,"TrackInfo.originalTrkAlgo[TrackInfo.size]/I");

            if(detailMode){
            }
        }//}}}

        void setbranchadd(TTree *root, bool detailMode = false){//{{{
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
            root->SetBranchAddress("TrackInfo.originalTrkAlgo", originalTrkAlgo);

            if(detailMode){
            }
        }//}}}
};//}}}

class BInfoBranches{//{{{
public:
    int	    uj_size;
    int	    uj_index[MAX_XB];
    float  uj_mass[MAX_XB];
    float  uj_pt[MAX_XB];
    float  uj_eta[MAX_XB];
    float  uj_phi[MAX_XB];
    float  uj_px[MAX_XB];
    float  uj_py[MAX_XB];
    float  uj_pz[MAX_XB];
    float	uj_vtxX[MAX_XB];
    float  uj_vtxY[MAX_XB];
    float  uj_vtxZ[MAX_XB];
    float  uj_vtxXErr[MAX_XB];
    float  uj_vtxYErr[MAX_XB];
    float  uj_vtxZErr[MAX_XB];
    float  uj_vtxYXErr[MAX_XB];
    float  uj_vtxZXErr[MAX_XB];
    float  uj_vtxZYErr[MAX_XB];
    float	uj_vtxdof[MAX_XB];
    float	uj_vtxchi2[MAX_XB];
    int     uj_rfmu1_index[MAX_XB];
    int     uj_rfmu2_index[MAX_XB];
    
    float  uj_rfmu1_pt[MAX_XB];
    float  uj_rfmu1_eta[MAX_XB];
    float  uj_rfmu1_phi[MAX_XB];
    float  uj_rfmu2_pt[MAX_XB];
    float  uj_rfmu2_eta[MAX_XB];
    float  uj_rfmu2_phi[MAX_XB];
    
    int	    size;
    int	    index[MAX_XB];
    float	mass[MAX_XB];
    float	unfitted_mass[MAX_XB];
    float	pt[MAX_XB];
    float	eta[MAX_XB];
    float	phi[MAX_XB];
    float	px[MAX_XB];
    float	py[MAX_XB];
    float	pz[MAX_XB];
    float  pxE[MAX_XB];
    float  pyE[MAX_XB];
    float  pzE[MAX_XB];
    float  alpha[MAX_XB];
    float  svpvDistance[MAX_XB];
    float  svpvDisErr[MAX_XB];
    float  svpvDistance_2D[MAX_XB];
    float  svpvDisErr_2D[MAX_XB];
    float  MaxDoca[MAX_XB];
    float  vtxX[MAX_XB];
    float  vtxY[MAX_XB];
    float  vtxZ[MAX_XB];
    float  vtxXErr[MAX_XB];
    float  vtxYErr[MAX_XB];
    float  vtxZErr[MAX_XB];
    float  vtxYXErr[MAX_XB];
    float  vtxZXErr[MAX_XB];
    float  vtxZYErr[MAX_XB];
    float	vtxdof[MAX_XB];
    float	vtxchi2[MAX_XB];
    int     rfuj_index[MAX_XB];
    int     rftk1_index[MAX_XB];
    int     rftk2_index[MAX_XB];
    int     type[MAX_XB];
    
    float  rfmu1_pt[MAX_XB];
    float  rfmu1_eta[MAX_XB];
    float  rfmu1_phi[MAX_XB];
    float  rfmu2_pt[MAX_XB];
    float  rfmu2_eta[MAX_XB];
    float  rfmu2_phi[MAX_XB];
    float  rftk1_pt[MAX_XB];
    float  rftk1_eta[MAX_XB];
    float  rftk1_phi[MAX_XB];
    float  rftk2_pt[MAX_XB];
    float  rftk2_eta[MAX_XB];
    float  rftk2_phi[MAX_XB];
    
    float	tktk_unfitted_mass[MAX_XB];
    float	tktk_mass[MAX_XB];
    float  tktk_pt[MAX_XB];
    float  tktk_eta[MAX_XB];
    float  tktk_phi[MAX_XB];
    float  tktk_px[MAX_XB];
    float  tktk_py[MAX_XB];
    float  tktk_pz[MAX_XB];
    float	tktk_vtxX[MAX_XB];
    float  tktk_vtxY[MAX_XB];
    float  tktk_vtxZ[MAX_XB];
    float  tktk_vtxXErr[MAX_XB];
    float  tktk_vtxYErr[MAX_XB];
    float  tktk_vtxZErr[MAX_XB];
    float  tktk_vtxYXErr[MAX_XB];
    float  tktk_vtxZXErr[MAX_XB];
    float  tktk_vtxZYErr[MAX_XB];
    float	tktk_vtxdof[MAX_XB];
    float	tktk_vtxchi2[MAX_XB];
    float  tktk_rftk1_pt[MAX_XB];
    float  tktk_rftk1_eta[MAX_XB];
    float  tktk_rftk1_phi[MAX_XB];
    float  tktk_rftk2_pt[MAX_XB];
    float  tktk_rftk2_eta[MAX_XB];
    float  tktk_rftk2_phi[MAX_XB];
    
    void regTree(TTree *root, bool detailMode = false){//{{{
        root->Branch("BInfo.uj_size"          , &uj_size       , "BInfo.uj_size/I"			);
        root->Branch("BInfo.uj_index"         , uj_index       , "BInfo.uj_index[BInfo.uj_size]/I"	);
        root->Branch("BInfo.uj_mass"          , uj_mass        , "BInfo.uj_mass[BInfo.uj_size]/F"	);
        root->Branch("BInfo.uj_pt"            , uj_pt          , "BInfo.uj_pt[BInfo.uj_size]/F" );
        root->Branch("BInfo.uj_eta"            , uj_eta          , "BInfo.uj_eta[BInfo.uj_size]/F"  );
        root->Branch("BInfo.uj_phi"            , uj_phi          , "BInfo.uj_phi[BInfo.uj_size]/F"  );
        root->Branch("BInfo.uj_px"            , uj_px          , "BInfo.uj_px[BInfo.uj_size]/F"	);
        root->Branch("BInfo.uj_py"            , uj_py          , "BInfo.uj_py[BInfo.uj_size]/F"	);
        root->Branch("BInfo.uj_pz"            , uj_pz          , "BInfo.uj_pz[BInfo.uj_size]/F"	);
        root->Branch("BInfo.uj_vtxX"          , uj_vtxX        , "BInfo.uj_vtxX[BInfo.uj_size]/F"	);
        root->Branch("BInfo.uj_vtxY"          , uj_vtxY        , "BInfo.uj_vtxY[BInfo.uj_size]/F"	);
        root->Branch("BInfo.uj_vtxZ"          , uj_vtxZ        , "BInfo.uj_vtxZ[BInfo.uj_size]/F"	);
        root->Branch("BInfo.uj_vtxdof"        , uj_vtxdof      , "BInfo.uj_vtxdof[BInfo.uj_size]/F"	);
        root->Branch("BInfo.uj_vtxchi2"       , uj_vtxchi2     , "BInfo.uj_vtxchi2[BInfo.uj_size]/F"	);
        root->Branch("BInfo.uj_rfmu1_index"   , uj_rfmu1_index , "BInfo.uj_rfmu1_index[BInfo.uj_size]/I"	);
        root->Branch("BInfo.uj_rfmu2_index"   , uj_rfmu2_index , "BInfo.uj_rfmu2_index[BInfo.uj_size]/I"	);
        
        root->Branch("BInfo.size"             , &size          , "BInfo.size/I"			);
        root->Branch("BInfo.index"            , index          , "BInfo.index[BInfo.size]/I"		);
        root->Branch("BInfo.mass"             , mass           , "BInfo.mass[BInfo.size]/F"		);
        root->Branch("BInfo.unfitted_mass"    , unfitted_mass  , "BInfo.unfitted_mass[BInfo.size]/F");
        root->Branch("BInfo.pt"               , pt             , "BInfo.pt[BInfo.size]/F"		);
        root->Branch("BInfo.eta"              , eta            , "BInfo.eta[BInfo.size]/F"		);
        root->Branch("BInfo.phi"              , phi            , "BInfo.phi[BInfo.size]/F"		);
        root->Branch("BInfo.px"               , px             , "BInfo.px[BInfo.size]/F"		);
        root->Branch("BInfo.py"               , py             , "BInfo.py[BInfo.size]/F"		);
        root->Branch("BInfo.pz"               , pz             , "BInfo.pz[BInfo.size]/F"		);
        root->Branch("BInfo.pxE"              , pxE            , "BInfo.pxE[BInfo.size]/F"            );
        root->Branch("BInfo.pyE"              , pyE            , "BInfo.pyE[BInfo.size]/F"            );
        root->Branch("BInfo.pzE"              , pzE            , "BInfo.pzE[BInfo.size]/F"            );
        root->Branch("BInfo.alpha"            , alpha          , "BInfo.alpha[BInfo.size]/F"	);
        root->Branch("BInfo.svpvDistance"     , svpvDistance   , "BInfo.svpvDistance[BInfo.size]/F"	);
        root->Branch("BInfo.svpvDisErr"       , svpvDisErr     , "BInfo.svpvDisErr[BInfo.size]/F"	);
        root->Branch("BInfo.svpvDistance_2D"  , svpvDistance_2D, "BInfo.svpvDistance_2D[BInfo.size]/F"	);
        root->Branch("BInfo.svpvDisErr_2D"    , svpvDisErr_2D  , "BInfo.svpvDisErr_2D[BInfo.size]/F"	);
        root->Branch("BInfo.MaxDoca"          , MaxDoca        , "BInfo.MaxDoca[BInfo.size]/F"	);
        root->Branch("BInfo.vtxX"             , vtxX           , "BInfo.vtxX[BInfo.size]/F"		);
        root->Branch("BInfo.vtxY"             , vtxY           , "BInfo.vtxY[BInfo.size]/F"		);
        root->Branch("BInfo.vtxZ"             , vtxZ           , "BInfo.vtxZ[BInfo.size]/F"		);
        root->Branch("BInfo.vtxXErr"          , vtxXErr        , "BInfo.vtxXErr[BInfo.size]/F"          );
        root->Branch("BInfo.vtxYErr"          , vtxYErr        , "BInfo.vtxYErr[BInfo.size]/F"          );
        root->Branch("BInfo.vtxZErr"          , vtxZErr        , "BInfo.vtxZErr[BInfo.size]/F"          );
        root->Branch("BInfo.vtxYXErr"         , vtxYXErr       , "BInfo.vtxYXErr[BInfo.size]/F"          );
        root->Branch("BInfo.vtxZXErr"         , vtxZXErr       , "BInfo.vtxZXErr[BInfo.size]/F"          );
        root->Branch("BInfo.vtxZYErr"         , vtxZYErr       , "BInfo.vtxZYErr[BInfo.size]/F"          );
        root->Branch("BInfo.vtxdof"           , vtxdof         , "BInfo.vtxdof[BInfo.size]/F"		);
        root->Branch("BInfo.vtxchi2"          , vtxchi2        , "BInfo.vtxchi2[BInfo.size]/F"	);
        root->Branch("BInfo.rfuj_index"       , rfuj_index     , "BInfo.rfuj_index[BInfo.size]/I");
        root->Branch("BInfo.rftk1_index"      , rftk1_index    , "BInfo.rftk1_index[BInfo.size]/I");
        root->Branch("BInfo.rftk2_index"      , rftk2_index    , "BInfo.rftk2_index[BInfo.size]/I");
        root->Branch("BInfo.type"             , type           , "BInfo.type[BInfo.size]/I"	);
        
        root->Branch("BInfo.tktk_unfitted_mass" , tktk_unfitted_mass, "BInfo.tktk_unfitted_mass[BInfo.size]/F"     );
        root->Branch("BInfo.tktk_mass"          , tktk_mass        , "BInfo.tktk_mass[BInfo.size]/F"     );
        root->Branch("BInfo.tktk_pt"            , tktk_pt          , "BInfo.tktk_pt[BInfo.size]/F"	);
        root->Branch("BInfo.tktk_eta"            , tktk_eta          , "BInfo.tktk_eta[BInfo.size]/F"	);
        root->Branch("BInfo.tktk_phi"            , tktk_phi          , "BInfo.tktk_phi[BInfo.size]/F"	);
        root->Branch("BInfo.tktk_px"            , tktk_px          , "BInfo.tktk_px[BInfo.size]/F"	);
        root->Branch("BInfo.tktk_py"            , tktk_py          , "BInfo.tktk_py[BInfo.size]/F"	);
        root->Branch("BInfo.tktk_pz"            , tktk_pz          , "BInfo.tktk_pz[BInfo.size]/F"	);

        if(detailMode){
            root->Branch("BInfo.uj_rfmu1_pt"      , uj_rfmu1_pt    , "BInfo.uj_rfmu1_pt[BInfo.uj_size]/F");
            root->Branch("BInfo.uj_rfmu1_eta"     , uj_rfmu1_eta   , "BInfo.uj_rfmu1_eta[BInfo.uj_size]/F");
            root->Branch("BInfo.uj_rfmu1_phi"     , uj_rfmu1_phi   , "BInfo.uj_rfmu1_phi[BInfo.uj_size]/F");
            root->Branch("BInfo.uj_rfmu2_pt"      , uj_rfmu2_pt    , "BInfo.uj_rfmu2_pt[BInfo.uj_size]/F");
            root->Branch("BInfo.uj_rfmu2_eta"     , uj_rfmu2_eta   , "BInfo.uj_rfmu2_eta[BInfo.uj_size]/F");
            root->Branch("BInfo.uj_rfmu2_phi"     , uj_rfmu2_phi   , "BInfo.uj_rfmu2_phi[BInfo.uj_size]/F");

            root->Branch("BInfo.rfmu1_pt"         , rfmu1_pt       , "BInfo.rfmu1_pt[BInfo.size]/F" );
            root->Branch("BInfo.rfmu1_eta"        , rfmu1_eta      , "BInfo.rfmu1_eta[BInfo.size]/F" );
            root->Branch("BInfo.rfmu1_phi"        , rfmu1_phi      , "BInfo.rfmu1_phi[BInfo.size]/F" );
            root->Branch("BInfo.rfmu2_pt"         , rfmu2_pt       , "BInfo.rfmu2_pt[BInfo.size]/F" );
            root->Branch("BInfo.rfmu2_eta"        , rfmu2_eta      , "BInfo.rfmu2_eta[BInfo.size]/F" );
            root->Branch("BInfo.rfmu2_phi"        , rfmu2_phi      , "BInfo.rfmu2_phi[BInfo.size]/F" );
        
            root->Branch("BInfo.rftk1_pt"         , rftk1_pt       , "BInfo.rftk1_pt[BInfo.size]/F"     );
            root->Branch("BInfo.rftk1_eta"        , rftk1_eta      , "BInfo.rftk1_eta[BInfo.size]/F"     );
            root->Branch("BInfo.rftk1_phi"        , rftk1_phi      , "BInfo.rftk1_phi[BInfo.size]/F"     );
            root->Branch("BInfo.rftk2_pt"         , rftk2_pt       , "BInfo.rftk2_pt[BInfo.size]/F"     );
            root->Branch("BInfo.rftk2_eta"        , rftk2_eta      , "BInfo.rftk2_eta[BInfo.size]/F"     );
            root->Branch("BInfo.rftk2_phi"        , rftk2_phi      , "BInfo.rftk2_phi[BInfo.size]/F"     );

            root->Branch("BInfo.uj_vtxXErr"       , uj_vtxXErr     , "BInfo.uj_vtxXErr[BInfo.uj_size]/F"   );
            root->Branch("BInfo.uj_vtxYErr"       , uj_vtxYErr     , "BInfo.uj_vtxYErr[BInfo.uj_size]/F"   );
            root->Branch("BInfo.uj_vtxZErr"       , uj_vtxZErr     , "BInfo.uj_vtxZErr[BInfo.uj_size]/F"   );
            root->Branch("BInfo.uj_vtxYXErr"      , uj_vtxYXErr    , "BInfo.uj_vtxYXErr[BInfo.uj_size]/F"   );
            root->Branch("BInfo.uj_vtxZXErr"      , uj_vtxZXErr    , "BInfo.uj_vtxZXErr[BInfo.uj_size]/F"   );
            root->Branch("BInfo.uj_vtxZYErr"      , uj_vtxZYErr    , "BInfo.uj_vtxZYErr[BInfo.uj_size]/F"   );

            root->Branch("BInfo.tktk_vtxX"          , tktk_vtxX        , "BInfo.tktk_vtxX[BInfo.size]/F"	);
            root->Branch("BInfo.tktk_vtxY"          , tktk_vtxY        , "BInfo.tktk_vtxY[BInfo.size]/F"	);
            root->Branch("BInfo.tktk_vtxZ"          , tktk_vtxZ        , "BInfo.tktk_vtxZ[BInfo.size]/F"	);
            root->Branch("BInfo.tktk_vtxXErr"       , tktk_vtxXErr     , "BInfo.tktk_vtxXErr[BInfo.size]/F"   );
            root->Branch("BInfo.tktk_vtxYErr"       , tktk_vtxYErr     , "BInfo.tktk_vtxYErr[BInfo.size]/F"   );
            root->Branch("BInfo.tktk_vtxZErr"       , tktk_vtxZErr     , "BInfo.tktk_vtxZErr[BInfo.size]/F"   );
            root->Branch("BInfo.tktk_vtxYXErr"      , tktk_vtxYXErr    , "BInfo.tktk_vtxYXErr[BInfo.size]/F"   );
            root->Branch("BInfo.tktk_vtxZXErr"      , tktk_vtxZXErr    , "BInfo.tktk_vtxZXErr[BInfo.size]/F"   );
            root->Branch("BInfo.tktk_vtxZYErr"      , tktk_vtxZYErr    , "BInfo.tktk_vtxZYErr[BInfo.size]/F"   );
            root->Branch("BInfo.tktk_vtxdof"        , tktk_vtxdof      , "BInfo.tktk_vtxdof[BInfo.size]/F"	);
            root->Branch("BInfo.tktk_vtxchi2"       , tktk_vtxchi2     , "BInfo.tktk_vtxchi2[BInfo.size]/F"	);

            root->Branch("BInfo.tktk_rftk1_pt"      , tktk_rftk1_pt    , "BInfo.tktk_rftk1_pt[BInfo.size]/F"     );
            root->Branch("BInfo.tktk_rftk1_eta"     , tktk_rftk1_eta   , "BInfo.tktk_rftk1_eta[BInfo.size]/F"     );
            root->Branch("BInfo.tktk_rftk1_phi"     , tktk_rftk1_phi   , "BInfo.tktk_rftk1_phi[BInfo.size]/F"     );
            root->Branch("BInfo.tktk_rftk2_pt"      , tktk_rftk2_pt    , "BInfo.tktk_rftk2_pt[BInfo.size]/F"     );
            root->Branch("BInfo.tktk_rftk2_eta"     , tktk_rftk2_eta   , "BInfo.tktk_rftk2_eta[BInfo.size]/F"     );
            root->Branch("BInfo.tktk_rftk2_phi"     , tktk_rftk2_phi   , "BInfo.tktk_rftk2_phi[BInfo.size]/F"     );
        }
    }//}}}
    
    void setbranchadd(TTree *root, bool detailMode = false){//{{{
        root->SetBranchAddress("BInfo.uj_size"		    ,&uj_size	);
        root->SetBranchAddress("BInfo.uj_size"		    ,&uj_size	    );
        root->SetBranchAddress("BInfo.uj_index"        ,uj_index   );
        root->SetBranchAddress("BInfo.uj_mass"         ,uj_mass   	);
        root->SetBranchAddress("BInfo.uj_pt"           ,uj_pt       );
        root->SetBranchAddress("BInfo.uj_eta"           ,uj_eta         );
        root->SetBranchAddress("BInfo.uj_phi"           ,uj_phi         );
        root->SetBranchAddress("BInfo.uj_px"           ,uj_px     	);
        root->SetBranchAddress("BInfo.uj_py"           ,uj_py    	);
        root->SetBranchAddress("BInfo.uj_pz"           ,uj_pz   	);
        root->SetBranchAddress("BInfo.uj_vtxX"         ,uj_vtxX    );
        root->SetBranchAddress("BInfo.uj_vtxY"         ,uj_vtxY    );
        root->SetBranchAddress("BInfo.uj_vtxZ"         ,uj_vtxZ    );
        root->SetBranchAddress("BInfo.uj_vtxdof"	   ,uj_vtxdof	);
        root->SetBranchAddress("BInfo.uj_vtxchi2"      ,uj_vtxchi2 );
        root->SetBranchAddress("BInfo.uj_rfmu1_index"  ,uj_rfmu1_index );
        root->SetBranchAddress("BInfo.uj_rfmu2_index"  ,uj_rfmu2_index );
        
        root->SetBranchAddress("BInfo.size"            ,&size        );
        root->SetBranchAddress("BInfo.index"           ,index       	);
        root->SetBranchAddress("BInfo.mass"		    ,mass		);
        root->SetBranchAddress("BInfo.unfitted_mass",unfitted_mass		);
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
        
        root->SetBranchAddress("BInfo.tktk_unfitted_mass",tktk_unfitted_mass    );
        root->SetBranchAddress("BInfo.tktk_mass"         ,tktk_mass    );
        root->SetBranchAddress("BInfo.tktk_pt"           ,tktk_pt     	);
        root->SetBranchAddress("BInfo.tktk_eta"           ,tktk_eta     	);
        root->SetBranchAddress("BInfo.tktk_phi"           ,tktk_phi     	);
        root->SetBranchAddress("BInfo.tktk_px"           ,tktk_px     	);
        root->SetBranchAddress("BInfo.tktk_py"           ,tktk_py    	);
        root->SetBranchAddress("BInfo.tktk_pz"           ,tktk_pz   	);

        if(detailMode){
            root->SetBranchAddress("BInfo.uj_rfmu1_pt"     ,uj_rfmu1_pt     );
            root->SetBranchAddress("BInfo.uj_rfmu1_eta"    ,uj_rfmu1_eta    );
            root->SetBranchAddress("BInfo.uj_rfmu1_phi"    ,uj_rfmu1_phi    );
            root->SetBranchAddress("BInfo.uj_rfmu2_pt"     ,uj_rfmu2_pt     );
            root->SetBranchAddress("BInfo.uj_rfmu2_eta"    ,uj_rfmu2_eta    );
            root->SetBranchAddress("BInfo.uj_rfmu2_phi"    ,uj_rfmu2_phi    );

            root->SetBranchAddress("BInfo.rfmu1_pt"        ,rfmu1_pt 	);
            root->SetBranchAddress("BInfo.rfmu1_eta"       ,rfmu1_eta	);
            root->SetBranchAddress("BInfo.rfmu1_phi"       ,rfmu1_phi   );
            root->SetBranchAddress("BInfo.rfmu2_pt"        ,rfmu2_pt    );
            root->SetBranchAddress("BInfo.rfmu2_eta"       ,rfmu2_eta   );
            root->SetBranchAddress("BInfo.rfmu2_phi"       ,rfmu2_phi   );

            root->SetBranchAddress("BInfo.rftk1_pt"        ,rftk1_pt 	);
            root->SetBranchAddress("BInfo.rftk1_eta"       ,rftk1_eta   );
            root->SetBranchAddress("BInfo.rftk1_phi"       ,rftk1_phi   );
            root->SetBranchAddress("BInfo.rftk2_pt"        ,rftk2_pt    );
            root->SetBranchAddress("BInfo.rftk2_eta"       ,rftk2_eta   );
            root->SetBranchAddress("BInfo.rftk2_phi"       ,rftk2_phi  	);

            root->SetBranchAddress("BInfo.uj_vtxXErr"      ,uj_vtxXErr  );
            root->SetBranchAddress("BInfo.uj_vtxYErr"      ,uj_vtxYErr  );
            root->SetBranchAddress("BInfo.uj_vtxZErr"      ,uj_vtxZErr  );
            root->SetBranchAddress("BInfo.uj_vtxYXErr"     ,uj_vtxYXErr  );
            root->SetBranchAddress("BInfo.uj_vtxZXErr"     ,uj_vtxZXErr  );
            root->SetBranchAddress("BInfo.uj_vtxZYErr"     ,uj_vtxZYErr  );

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

            root->SetBranchAddress("BInfo.tktk_rftk1_pt"     ,tktk_rftk1_pt 	);
            root->SetBranchAddress("BInfo.tktk_rftk1_eta"    ,tktk_rftk1_eta     );
            root->SetBranchAddress("BInfo.tktk_rftk1_phi"    ,tktk_rftk1_phi     );
            root->SetBranchAddress("BInfo.tktk_rftk2_pt"     ,tktk_rftk2_pt    	);
            root->SetBranchAddress("BInfo.tktk_rftk2_eta"    ,tktk_rftk2_eta   	);
            root->SetBranchAddress("BInfo.tktk_rftk2_phi"    ,tktk_rftk2_phi  	);
        }
    }//}}}
};//}}}

class DInfoBranches{//{{{
public:
    int	    size;
    int	    index[MAX_XB];
    int     type[MAX_XB];

    float	tktkRes_mass[MAX_XB];
    float  tktkRes_pt[MAX_XB];
    float  tktkRes_eta[MAX_XB];
    float  tktkRes_phi[MAX_XB];
    float	tktkRes_vtxX[MAX_XB];
    float  tktkRes_vtxY[MAX_XB];
    float  tktkRes_vtxZ[MAX_XB];
    float  tktkRes_vtxXErr[MAX_XB];
    float  tktkRes_vtxYErr[MAX_XB];
    float  tktkRes_vtxZErr[MAX_XB];
    float  tktkRes_vtxYXErr[MAX_XB];
    float  tktkRes_vtxZXErr[MAX_XB];
    float  tktkRes_vtxZYErr[MAX_XB];
    float	tktkRes_vtxdof[MAX_XB];
    float	tktkRes_vtxchi2[MAX_XB];
    float  tktkRes_svpvDistance[MAX_XB];
    float  tktkRes_svpvDisErr[MAX_XB];
	float  tktkRes_alpha[MAX_XB];
    float  tktkRes_rftk1_mass[MAX_XB];
    float  tktkRes_rftk1_pt[MAX_XB];
    float  tktkRes_rftk1_eta[MAX_XB];
    float  tktkRes_rftk1_phi[MAX_XB];
    float  tktkRes_rftk2_mass[MAX_XB];
    float  tktkRes_rftk2_pt[MAX_XB];
    float  tktkRes_rftk2_eta[MAX_XB];
    float  tktkRes_rftk2_phi[MAX_XB];
    float  tktkRes_rftk3_mass[MAX_XB];
    float  tktkRes_rftk3_pt[MAX_XB];
    float  tktkRes_rftk3_eta[MAX_XB];
    float  tktkRes_rftk3_phi[MAX_XB];
    float  tktkRes_rftk4_mass[MAX_XB];
    float  tktkRes_rftk4_pt[MAX_XB];
    float  tktkRes_rftk4_eta[MAX_XB];
    float  tktkRes_rftk4_phi[MAX_XB];
    int     tktkRes_rftk1_index[MAX_XB];
    int     tktkRes_rftk2_index[MAX_XB];
    int     tktkRes_rftk3_index[MAX_XB];
    int     tktkRes_rftk4_index[MAX_XB];

    float	mass[MAX_XB];
    float  pt[MAX_XB];
    float  eta[MAX_XB];
    float  phi[MAX_XB];
    float  px[MAX_XB];
    float  py[MAX_XB];
    float  pz[MAX_XB];
    float  alpha[MAX_XB];
    float  svpvDistance[MAX_XB];
    float  svpvDisErr[MAX_XB];
    float  svpvDistance_2D[MAX_XB];
    float  svpvDisErr_2D[MAX_XB];
    float  MaxDoca[MAX_XB];
    float	vtxX[MAX_XB];
    float  vtxY[MAX_XB];
    float  vtxZ[MAX_XB];
    float  vtxXErr[MAX_XB];
    float  vtxYErr[MAX_XB];
    float  vtxZErr[MAX_XB];
    float  vtxYXErr[MAX_XB];
    float  vtxZXErr[MAX_XB];
    float  vtxZYErr[MAX_XB];
    float	vtxdof[MAX_XB];
    float	vtxchi2[MAX_XB];
    float  rftk1_mass[MAX_XB];
    float  rftk1_pt[MAX_XB];
    float  rftk1_eta[MAX_XB];
    float  rftk1_phi[MAX_XB];
    float  rftk2_mass[MAX_XB];
    float  rftk2_pt[MAX_XB];
    float  rftk2_eta[MAX_XB];
    float  rftk2_phi[MAX_XB];
    float  rftk3_mass[MAX_XB];
    float  rftk3_pt[MAX_XB];
    float  rftk3_eta[MAX_XB];
    float  rftk3_phi[MAX_XB];
    float  rftk4_mass[MAX_XB];
    float  rftk4_pt[MAX_XB];
    float  rftk4_eta[MAX_XB];
    float  rftk4_phi[MAX_XB];
    float  rftk5_mass[MAX_XB];
    float  rftk5_pt[MAX_XB];
    float  rftk5_eta[MAX_XB];
    float  rftk5_phi[MAX_XB];
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
   
    void regTree(TTree *root, bool detailMode = false){//{{{
        root->Branch("DInfo.size"             , &size          , "DInfo.size/I"			);
        root->Branch("DInfo.index"            , index          , "DInfo.index[DInfo.size]/I"		);
        root->Branch("DInfo.type"             , type           , "DInfo.type[DInfo.size]/I"	);
        
        root->Branch("DInfo.tktkRes_mass"          , tktkRes_mass              , "DInfo.tktkRes_mass[DInfo.size]/F"     );
        root->Branch("DInfo.tktkRes_pt"            , tktkRes_pt                , "DInfo.tktkRes_pt[DInfo.size]/F"	);
        root->Branch("DInfo.tktkRes_eta"           , tktkRes_eta               , "DInfo.tktkRes_eta[DInfo.size]/F"	);
        root->Branch("DInfo.tktkRes_phi"           , tktkRes_phi               , "DInfo.tktkRes_phi[DInfo.size]/F"	);
        root->Branch("DInfo.tktkRes_vtxX"          , tktkRes_vtxX              , "DInfo.tktkRes_vtxX[DInfo.size]/F"	);
        root->Branch("DInfo.tktkRes_vtxY"          , tktkRes_vtxY              , "DInfo.tktkRes_vtxY[DInfo.size]/F"	);
        root->Branch("DInfo.tktkRes_vtxZ"          , tktkRes_vtxZ              , "DInfo.tktkRes_vtxZ[DInfo.size]/F"	);
        root->Branch("DInfo.tktkRes_vtxXErr"       , tktkRes_vtxXErr           , "DInfo.tktkRes_vtxXErr[DInfo.size]/F"   );
        root->Branch("DInfo.tktkRes_vtxYErr"       , tktkRes_vtxYErr           , "DInfo.tktkRes_vtxYErr[DInfo.size]/F"   );
        root->Branch("DInfo.tktkRes_vtxZErr"       , tktkRes_vtxZErr           , "DInfo.tktkRes_vtxZErr[DInfo.size]/F"   );
        root->Branch("DInfo.tktkRes_vtxYXErr"      , tktkRes_vtxYXErr          , "DInfo.tktkRes_vtxYXErr[DInfo.size]/F"   );
        root->Branch("DInfo.tktkRes_vtxZXErr"      , tktkRes_vtxZXErr          , "DInfo.tktkRes_vtxZXErr[DInfo.size]/F"   );
        root->Branch("DInfo.tktkRes_vtxZYErr"      , tktkRes_vtxZYErr          , "DInfo.tktkRes_vtxZYErr[DInfo.size]/F"   );
        root->Branch("DInfo.tktkRes_vtxdof"        , tktkRes_vtxdof            , "DInfo.tktkRes_vtxdof[DInfo.size]/F"	);
        root->Branch("DInfo.tktkRes_vtxchi2"       , tktkRes_vtxchi2           , "DInfo.tktkRes_vtxchi2[DInfo.size]/F"	);
        root->Branch("DInfo.tktkRes_svpvDistance"  , tktkRes_svpvDistance      , "DInfo.tktkRes_svpvDistance[DInfo.size]/F"	);
        root->Branch("DInfo.tktkRes_svpvDisErr"    , tktkRes_svpvDisErr        , "DInfo.tktkRes_svpvDisErr[DInfo.size]/F"	);
		root->Branch("DInfo.tktkRes_alpha"         , tktkRes_alpha             , "DInfo.tktkRes_alpha[DInfo.size]/F" );

        root->Branch("DInfo.tktkRes_rftk1_index"   , tktkRes_rftk1_index       , "DInfo.tktkRes_rftk1_index[DInfo.size]/I");
        root->Branch("DInfo.tktkRes_rftk2_index"   , tktkRes_rftk2_index       , "DInfo.tktkRes_rftk2_index[DInfo.size]/I");
        root->Branch("DInfo.tktkRes_rftk3_index"   , tktkRes_rftk3_index       , "DInfo.tktkRes_rftk3_index[DInfo.size]/I");
        root->Branch("DInfo.tktkRes_rftk4_index"   , tktkRes_rftk4_index       , "DInfo.tktkRes_rftk4_index[DInfo.size]/I");

        root->Branch("DInfo.mass"             , mass              , "DInfo.mass[DInfo.size]/F"     );
        root->Branch("DInfo.pt"               , pt                , "DInfo.pt[DInfo.size]/F"	);
        root->Branch("DInfo.eta"              , eta               , "DInfo.eta[DInfo.size]/F"	);
        root->Branch("DInfo.phi"              , phi               , "DInfo.phi[DInfo.size]/F"	);
        root->Branch("DInfo.px"               , px                , "DInfo.px[DInfo.size]/F"	);
        root->Branch("DInfo.py"               , py                , "DInfo.py[DInfo.size]/F"	);
        root->Branch("DInfo.pz"               , pz                , "DInfo.pz[DInfo.size]/F"	);
        root->Branch("DInfo.alpha"            , alpha             , "DInfo.alpha[DInfo.size]/F"	);
        root->Branch("DInfo.svpvDistance"     , svpvDistance      , "DInfo.svpvDistance[DInfo.size]/F"	);
        root->Branch("DInfo.svpvDisErr"       , svpvDisErr        , "DInfo.svpvDisErr[DInfo.size]/F"	);
        root->Branch("DInfo.svpvDistance_2D"  , svpvDistance_2D   , "DInfo.svpvDistance_2D[DInfo.size]/F"	);
        root->Branch("DInfo.svpvDisErr_2D"    , svpvDisErr_2D     , "DInfo.svpvDisErr_2D[DInfo.size]/F"	);
        root->Branch("DInfo.MaxDoca"          , MaxDoca           , "DInfo.MaxDoca[DInfo.size]/F"	);
        root->Branch("DInfo.vtxX"             , vtxX              , "DInfo.vtxX[DInfo.size]/F"	);
        root->Branch("DInfo.vtxY"             , vtxY              , "DInfo.vtxY[DInfo.size]/F"	);
        root->Branch("DInfo.vtxZ"             , vtxZ              , "DInfo.vtxZ[DInfo.size]/F"	);
        root->Branch("DInfo.vtxXErr"          , vtxXErr           , "DInfo.vtxXErr[DInfo.size]/F"   );
        root->Branch("DInfo.vtxYErr"          , vtxYErr           , "DInfo.vtxYErr[DInfo.size]/F"   );
        root->Branch("DInfo.vtxZErr"          , vtxZErr           , "DInfo.vtxZErr[DInfo.size]/F"   );
        root->Branch("DInfo.vtxYXErr"         , vtxYXErr          , "DInfo.vtxYXErr[DInfo.size]/F"   );
        root->Branch("DInfo.vtxZXErr"         , vtxZXErr          , "DInfo.vtxZXErr[DInfo.size]/F"   );
        root->Branch("DInfo.vtxZYErr"         , vtxZYErr          , "DInfo.vtxZYErr[DInfo.size]/F"   );
        root->Branch("DInfo.vtxdof"           , vtxdof            , "DInfo.vtxdof[DInfo.size]/F"	);
        root->Branch("DInfo.vtxchi2"          , vtxchi2           , "DInfo.vtxchi2[DInfo.size]/F"	);

        root->Branch("DInfo.rftk1_MassHypo"   , rftk1_MassHypo       , "DInfo.rftk1_MassHypo[DInfo.size]/I");
        root->Branch("DInfo.rftk2_MassHypo"   , rftk2_MassHypo       , "DInfo.rftk2_MassHypo[DInfo.size]/I");
        root->Branch("DInfo.rftk3_MassHypo"   , rftk3_MassHypo       , "DInfo.rftk3_MassHypo[DInfo.size]/I");
        root->Branch("DInfo.rftk4_MassHypo"   , rftk4_MassHypo       , "DInfo.rftk4_MassHypo[DInfo.size]/I");
        root->Branch("DInfo.rftk5_MassHypo"   , rftk5_MassHypo       , "DInfo.rftk5_MassHypo[DInfo.size]/I");
        root->Branch("DInfo.rftk1_index"      , rftk1_index       , "DInfo.rftk1_index[DInfo.size]/I");
        root->Branch("DInfo.rftk2_index"      , rftk2_index       , "DInfo.rftk2_index[DInfo.size]/I");
        root->Branch("DInfo.rftk3_index"      , rftk3_index       , "DInfo.rftk3_index[DInfo.size]/I");
        root->Branch("DInfo.rftk4_index"      , rftk4_index       , "DInfo.rftk4_index[DInfo.size]/I");
        root->Branch("DInfo.rftk5_index"      , rftk5_index       , "DInfo.rftk5_index[DInfo.size]/I");

        root->Branch("DInfo.rftk1_mass"       , rftk1_mass        , "DInfo.rftk1_mass[DInfo.size]/F"     );
        root->Branch("DInfo.rftk2_mass"       , rftk2_mass        , "DInfo.rftk2_mass[DInfo.size]/F"     );
        root->Branch("DInfo.rftk3_mass"       , rftk3_mass        , "DInfo.rftk3_mass[DInfo.size]/F"     );
        root->Branch("DInfo.rftk4_mass"       , rftk4_mass        , "DInfo.rftk4_mass[DInfo.size]/F"     );
        root->Branch("DInfo.rftk5_mass"       , rftk5_mass        , "DInfo.rftk5_mass[DInfo.size]/F"     );
        root->Branch("DInfo.tktkRes_rftk1_mass"    , tktkRes_rftk1_mass        , "DInfo.tktkRes_rftk1_mass[DInfo.size]/F"     );
        root->Branch("DInfo.tktkRes_rftk2_mass"    , tktkRes_rftk2_mass        , "DInfo.tktkRes_rftk2_mass[DInfo.size]/F"     );
        root->Branch("DInfo.tktkRes_rftk3_mass"    , tktkRes_rftk3_mass        , "DInfo.tktkRes_rftk3_mass[DInfo.size]/F"     );
        root->Branch("DInfo.tktkRes_rftk4_mass"    , tktkRes_rftk4_mass        , "DInfo.tktkRes_rftk4_mass[DInfo.size]/F"     );

        if(detailMode){
            root->Branch("DInfo.tktkRes_rftk1_pt"      , tktkRes_rftk1_pt          , "DInfo.tktkRes_rftk1_pt[DInfo.size]/F"     );
            root->Branch("DInfo.tktkRes_rftk1_eta"     , tktkRes_rftk1_eta         , "DInfo.tktkRes_rftk1_eta[DInfo.size]/F"     );
            root->Branch("DInfo.tktkRes_rftk1_phi"     , tktkRes_rftk1_phi         , "DInfo.tktkRes_rftk1_phi[DInfo.size]/F"     );
            root->Branch("DInfo.tktkRes_rftk2_pt"      , tktkRes_rftk2_pt          , "DInfo.tktkRes_rftk2_pt[DInfo.size]/F"     );
            root->Branch("DInfo.tktkRes_rftk2_eta"     , tktkRes_rftk2_eta         , "DInfo.tktkRes_rftk2_eta[DInfo.size]/F"     );
            root->Branch("DInfo.tktkRes_rftk2_phi"     , tktkRes_rftk2_phi         , "DInfo.tktkRes_rftk2_phi[DInfo.size]/F"     );
            root->Branch("DInfo.tktkRes_rftk3_pt"      , tktkRes_rftk3_pt          , "DInfo.tktkRes_rftk3_pt[DInfo.size]/F"     );
            root->Branch("DInfo.tktkRes_rftk3_eta"     , tktkRes_rftk3_eta         , "DInfo.tktkRes_rftk3_eta[DInfo.size]/F"     );
            root->Branch("DInfo.tktkRes_rftk3_phi"     , tktkRes_rftk3_phi         , "DInfo.tktkRes_rftk3_phi[DInfo.size]/F"     );
            root->Branch("DInfo.tktkRes_rftk4_pt"      , tktkRes_rftk4_pt          , "DInfo.tktkRes_rftk4_pt[DInfo.size]/F"     );
            root->Branch("DInfo.tktkRes_rftk4_eta"     , tktkRes_rftk4_eta         , "DInfo.tktkRes_rftk4_eta[DInfo.size]/F"     );
            root->Branch("DInfo.tktkRes_rftk4_phi"     , tktkRes_rftk4_phi         , "DInfo.tktkRes_rftk4_phi[DInfo.size]/F"     );

            root->Branch("DInfo.rftk1_pt"         , rftk1_pt          , "DInfo.rftk1_pt[DInfo.size]/F"     );
            root->Branch("DInfo.rftk1_eta"        , rftk1_eta         , "DInfo.rftk1_eta[DInfo.size]/F"     );
            root->Branch("DInfo.rftk1_phi"        , rftk1_phi         , "DInfo.rftk1_phi[DInfo.size]/F"     );
            root->Branch("DInfo.rftk2_pt"         , rftk2_pt          , "DInfo.rftk2_pt[DInfo.size]/F"     );
            root->Branch("DInfo.rftk2_eta"        , rftk2_eta         , "DInfo.rftk2_eta[DInfo.size]/F"     );
            root->Branch("DInfo.rftk2_phi"        , rftk2_phi         , "DInfo.rftk2_phi[DInfo.size]/F"     );
            root->Branch("DInfo.rftk3_pt"         , rftk3_pt          , "DInfo.rftk3_pt[DInfo.size]/F"     );
            root->Branch("DInfo.rftk3_eta"        , rftk3_eta         , "DInfo.rftk3_eta[DInfo.size]/F"     );
            root->Branch("DInfo.rftk3_phi"        , rftk3_phi         , "DInfo.rftk3_phi[DInfo.size]/F"     );
            root->Branch("DInfo.rftk4_pt"         , rftk4_pt          , "DInfo.rftk4_pt[DInfo.size]/F"     );
            root->Branch("DInfo.rftk4_eta"        , rftk4_eta         , "DInfo.rftk4_eta[DInfo.size]/F"     );
            root->Branch("DInfo.rftk4_phi"        , rftk4_phi         , "DInfo.rftk4_phi[DInfo.size]/F"     );
            root->Branch("DInfo.rftk5_pt"         , rftk5_pt          , "DInfo.rftk5_pt[DInfo.size]/F"     );
            root->Branch("DInfo.rftk5_eta"        , rftk5_eta         , "DInfo.rftk5_eta[DInfo.size]/F"     );
            root->Branch("DInfo.rftk5_phi"        , rftk5_phi         , "DInfo.rftk5_phi[DInfo.size]/F"     );
        }
    }//}}}
    
    void setbranchadd(TTree *root, bool detailMode = false){//{{{
        root->SetBranchAddress("DInfo.size"            ,&size        );
        root->SetBranchAddress("DInfo.index"           ,index       	);
        root->SetBranchAddress("DInfo.type"            ,type   	);
        
        root->SetBranchAddress("DInfo.tktkRes_mass"            ,tktkRes_mass    );
        root->SetBranchAddress("DInfo.tktkRes_pt"              ,tktkRes_pt     	);
        root->SetBranchAddress("DInfo.tktkRes_eta"             ,tktkRes_eta     	);
        root->SetBranchAddress("DInfo.tktkRes_phi"             ,tktkRes_phi     	);
        root->SetBranchAddress("DInfo.tktkRes_vtxX"            ,tktkRes_vtxX    );
        root->SetBranchAddress("DInfo.tktkRes_vtxY"            ,tktkRes_vtxY    );
        root->SetBranchAddress("DInfo.tktkRes_vtxZ"            ,tktkRes_vtxZ    );
        root->SetBranchAddress("DInfo.tktkRes_vtxdof"          ,tktkRes_vtxdof	);
        root->SetBranchAddress("DInfo.tktkRes_vtxchi2"         ,tktkRes_vtxchi2 );
        root->SetBranchAddress("DInfo.tktkRes_vtxXErr"         ,tktkRes_vtxXErr   );
        root->SetBranchAddress("DInfo.tktkRes_vtxYErr"         ,tktkRes_vtxYErr   );
        root->SetBranchAddress("DInfo.tktkRes_vtxZErr"         ,tktkRes_vtxZErr   );
        root->SetBranchAddress("DInfo.tktkRes_vtxYXErr"        ,tktkRes_vtxYXErr   );
        root->SetBranchAddress("DInfo.tktkRes_vtxZXErr"        ,tktkRes_vtxZXErr   );
        root->SetBranchAddress("DInfo.tktkRes_vtxZYErr"        ,tktkRes_vtxZYErr   );
        root->SetBranchAddress("DInfo.tktkRes_svpvDistance"    ,tktkRes_svpvDistance   	);
        root->SetBranchAddress("DInfo.tktkRes_svpvDisErr"      ,tktkRes_svpvDisErr   	);
		root->SetBranchAddress("DInfo.tktkRes_alpha"           ,tktkRes_alpha    );

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

        root->SetBranchAddress("DInfo.rftk1_mass"      ,rftk1_mass);
        root->SetBranchAddress("DInfo.rftk2_mass"      ,rftk2_mass);
        root->SetBranchAddress("DInfo.rftk3_mass"      ,rftk3_mass);
        root->SetBranchAddress("DInfo.rftk4_mass"      ,rftk4_mass);
        root->SetBranchAddress("DInfo.rftk5_mass"      ,rftk5_mass);
        root->SetBranchAddress("DInfo.tktkRes_rftk1_mass"      ,tktkRes_rftk1_mass  );
        root->SetBranchAddress("DInfo.tktkRes_rftk2_mass"      ,tktkRes_rftk2_mass  );
        root->SetBranchAddress("DInfo.tktkRes_rftk3_mass"      ,tktkRes_rftk3_mass  );
        root->SetBranchAddress("DInfo.tktkRes_rftk4_mass"      ,tktkRes_rftk4_mass  );

        if(detailMode){
            root->SetBranchAddress("DInfo.tktkRes_rftk1_pt"        ,tktkRes_rftk1_pt  );
            root->SetBranchAddress("DInfo.tktkRes_rftk1_eta"       ,tktkRes_rftk1_eta  );
            root->SetBranchAddress("DInfo.tktkRes_rftk1_phi"       ,tktkRes_rftk1_phi  );
            root->SetBranchAddress("DInfo.tktkRes_rftk2_pt"        ,tktkRes_rftk2_pt  );
            root->SetBranchAddress("DInfo.tktkRes_rftk2_eta"       ,tktkRes_rftk2_eta  );
            root->SetBranchAddress("DInfo.tktkRes_rftk2_phi"       ,tktkRes_rftk2_phi  );
            root->SetBranchAddress("DInfo.tktkRes_rftk3_pt"        ,tktkRes_rftk3_pt  );
            root->SetBranchAddress("DInfo.tktkRes_rftk3_eta"       ,tktkRes_rftk3_eta  );
            root->SetBranchAddress("DInfo.tktkRes_rftk3_phi"       ,tktkRes_rftk3_phi  );
            root->SetBranchAddress("DInfo.tktkRes_rftk4_pt"        ,tktkRes_rftk4_pt  );
            root->SetBranchAddress("DInfo.tktkRes_rftk4_eta"       ,tktkRes_rftk4_eta  );
            root->SetBranchAddress("DInfo.tktkRes_rftk4_phi"       ,tktkRes_rftk4_phi  );

            root->SetBranchAddress("DInfo.rftk1_pt"        ,rftk1_pt  );
            root->SetBranchAddress("DInfo.rftk1_eta"       ,rftk1_eta  );
            root->SetBranchAddress("DInfo.rftk1_phi"       ,rftk1_phi  );
            root->SetBranchAddress("DInfo.rftk2_pt"        ,rftk2_pt  );
            root->SetBranchAddress("DInfo.rftk2_eta"       ,rftk2_eta  );
            root->SetBranchAddress("DInfo.rftk2_phi"       ,rftk2_phi  );
            root->SetBranchAddress("DInfo.rftk3_pt"        ,rftk3_pt  );
            root->SetBranchAddress("DInfo.rftk3_eta"       ,rftk3_eta  );
            root->SetBranchAddress("DInfo.rftk3_phi"       ,rftk3_phi  );
            root->SetBranchAddress("DInfo.rftk4_pt"        ,rftk4_pt  );
            root->SetBranchAddress("DInfo.rftk4_eta"       ,rftk4_eta  );
            root->SetBranchAddress("DInfo.rftk4_phi"       ,rftk4_phi  );
            root->SetBranchAddress("DInfo.rftk5_pt"        ,rftk5_pt  );
            root->SetBranchAddress("DInfo.rftk5_eta"       ,rftk5_eta  );
            root->SetBranchAddress("DInfo.rftk5_phi"       ,rftk5_phi  );
        }
    }//}}}
};//}}}

class GenInfoBranches{//{{{
    public:
        int     size;
        int     index       [MAX_GEN];
        int     handle_index[MAX_GEN];
        float  pt          [MAX_GEN];
        float  eta         [MAX_GEN];
        float  phi         [MAX_GEN];
        float  mass        [MAX_GEN];
        int     pdgId       [MAX_GEN];
        int     status      [MAX_GEN];
		int     collisionId [MAX_GEN];//to tell if it is from pythia event or hydjet event
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
            root->Branch("GenInfo.pt"           ,pt             ,"GenInfo.pt[GenInfo.size]/F");
            root->Branch("GenInfo.eta"          ,eta            ,"GenInfo.eta[GenInfo.size]/F");
            root->Branch("GenInfo.phi"          ,phi            ,"GenInfo.phi[GenInfo.size]/F");
            root->Branch("GenInfo.mass"         ,mass           ,"GenInfo.mass[GenInfo.size]/F");
            root->Branch("GenInfo.pdgId"        ,pdgId          ,"GenInfo.pdgId[GenInfo.size]/I");
            root->Branch("GenInfo.status"       ,status         ,"GenInfo.status[GenInfo.size]/I");
            root->Branch("GenInfo.collisionId"  ,collisionId    ,"GenInfo.collisionId[GenInfo.size]/I");
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
            root->SetBranchAddress("GenInfo.collisionId"  ,collisionId    );
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

#endif
