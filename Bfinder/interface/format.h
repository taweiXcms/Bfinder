// vim:set ts=4 sw=4 fdm=marker et:
#ifndef _XBFRAMEFORMAT_H_
#define _XBFRAMEFORMAT_H_

#define MAX_XB 16384
//#define MAX_XB 32768//When the size get too large, SetBranchAddress will fail
#define MAX_MUON 512
#define MAX_TRACK 8192
#define MAX_GEN 4096 
#define MAX_BX 128
#define MAX_Vertices 4096
//#define N_TRIGGER_BOOKINGS 5842

#include <string>
#include <vector>
#include "TTree.h"

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
		//double	PVc2p;
		
		void regTree(TTree *root){//{{{
			root->Branch("EvtInfo.RunNo"      , &RunNo                     , "EvtInfo.RunNo/I"			);
			root->Branch("EvtInfo.EvtNo"      , &EvtNo                     , "EvtInfo.EvtNo/I"			);
			root->Branch("EvtInfo.BxNo"       , &BxNo                      , "EvtInfo.BxNo/I"			);
			root->Branch("EvtInfo.LumiNo"     , &LumiNo                    , "EvtInfo.LumiNo/I"			);
			root->Branch("EvtInfo.Orbit"      , &Orbit                     , "EvtInfo.Orbit/I"			);
			root->Branch("EvtInfo.McFlag"     , &McFlag                    , "EvtInfo.McFlag/O"			);
			root->Branch("EvtInfo.nBX"        , &nBX                       , "EvtInfo.nBX/I" 			);
			root->Branch("EvtInfo.BXPU"       , BXPU                       , "EvtInfo.BXPU[EvtInfo.nBX]/I");
			root->Branch("EvtInfo.nPU"        , nPU                        , "EvtInfo.nPU[EvtInfo.nBX]/I");
			root->Branch("EvtInfo.trueIT"     , trueIT                     , "EvtInfo.trueIT[EvtInfo.nBX]/F");
            //root->Branch("EvtInfo.trgCount"   , &trgCount                  , "EvtInfo.trgCount/I"       );
            //root->Branch("EvtInfo.nTrgBook"   , &nTrgBook                  , "EvtInfo.nTrgBook/I"       );
            //root->Branch("EvtInfo.trgBook"    , trgBook                    , "EvtInfo.trgBook[EvtInfo.nTrgBook]/B");//notice /B
			//root->Branch("EvtInfo.nHLT"       , &nHLT                      , "EvtInfo.nHLT/I"			);
            //root->Branch("EvtInfo.hltBits"    , hltBits                    , "EvtInfo.hltBits[EvtInfo.nHLT]/O");
			//root->Branch("EvtInfo.hltnames" , "std::vector<std::string>" , &hltnames);
			//root->Branch("EvtInfo.hltflag"  , hltflag                    , "EvtInfo.hltflag[EvtInfo.nHLT]/O"	);
			//root->Branch("EvtInfo.nHLTm"    , &nHLTm                     , "EvtInfo.nHLTm/I"			);
			//root->Branch("EvtInfo.hltflagm" , hltflagm                   , "EvtInfo.hltflagm[nHLTm]/O"		);
  		    root->Branch("EvtInfo.PVx"        , &PVx                       , "EvtInfo.PVx/D"			);
  		    root->Branch("EvtInfo.PVy"        , &PVy                       , "EvtInfo.PVy/D"			);
  		    root->Branch("EvtInfo.PVz"        , &PVz                       , "EvtInfo.PVz/D"			);
            root->Branch("EvtInfo.PVxE"       , &PVxE                      , "EvtInfo.PVxE/D"           );
            root->Branch("EvtInfo.PVyE"       , &PVyE                      , "EvtInfo.PVyE/D"           );
            root->Branch("EvtInfo.PVzE"       , &PVzE                      , "EvtInfo.PVzE/D"           );
  		    root->Branch("EvtInfo.PVnchi2"    , &PVnchi2                   , "EvtInfo.PVnchi2/D"		);
  		    root->Branch("EvtInfo.PVchi2"     , &PVchi2                    , "EvtInfo.PVchi2/D"			);
			//root->Branch("EvtInfo.PVc2p"    , &PVc2p                     , "EvtInfo.PVc2p/D"			);//
		}//}}}

    	void setbranchadd(TTree *root){ //{{{
            root->SetBranchAddress("EvtInfo.RunNo"    ,&RunNo	    );
            root->SetBranchAddress("EvtInfo.EvtNo"    ,&EvtNo       );
            root->SetBranchAddress("EvtInfo.BxNo"     ,&BxNo        );
            root->SetBranchAddress("EvtInfo.LumiNo"   ,&LumiNo      );
            root->SetBranchAddress("EvtInfo.Orbit"    ,&Orbit       );
            root->SetBranchAddress("EvtInfo.McFlag"   ,&McFlag      );
			root->SetBranchAddress("EvtInfo.nBX"      ,&nBX         );
			root->SetBranchAddress("EvtInfo.BXPU"     ,BXPU         );
			root->SetBranchAddress("EvtInfo.nPU"      ,nPU          );
			root->SetBranchAddress("EvtInfo.trueIT"   ,trueIT       );
            //root->SetBranchAddress("EvtInfo.trgCount" ,&trgCount    );
            //root->SetBranchAddress("EvtInfo.nTrgBook" ,&nTrgBook    );
            //root->SetBranchAddress("EvtInfo.trgBook"  ,trgBook      );
            //root->SetBranchAddress("EvtInfo.nHLT"     ,&nHLT	    );
            //root->SetBranchAddress("EvtInfo.hltBits"  ,hltBits      );
            //root->SetBranchAddress("EvtInfo.hltnames" ,&hltnames	);
            //root->SetBranchAddress("EvtInfo.hltflag"  ,hltflag	);
            //root->SetBranchAddress("EvtInfo.nHLTm"	  ,&nHLTm	);
            //root->SetBranchAddress("EvtInfo.hltflagm" ,hltflagm	);
            root->SetBranchAddress("EvtInfo.PVx"      ,&PVx		);
            root->SetBranchAddress("EvtInfo.PVy"      ,&PVy		);
            root->SetBranchAddress("EvtInfo.PVz"      ,&PVz		);
            root->SetBranchAddress("EvtInfo.PVxE"     ,&PVxE	);
            root->SetBranchAddress("EvtInfo.PVyE"     ,&PVyE	);
            root->SetBranchAddress("EvtInfo.PVzE"     ,&PVzE	);
            root->SetBranchAddress("EvtInfo.PVnchi2"  ,&PVnchi2	);
            root->SetBranchAddress("EvtInfo.PVchi2"   ,&PVchi2	);
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
        int     isGoodCand   [ MAX_MUON];
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
            root->Branch("MuonInfo.isGoodCand"    , isGoodCand    , "MuonInfo.isGoodCand[MuonInfo.size]/I");
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
            root->SetBranchAddress("MuonInfo.isGoodCand"    , isGoodCand);
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
        int     isGoodCand   [ MAX_TRACK];
        int     geninfo_index[ MAX_TRACK];
        int     trackQuality [ MAX_TRACK];
        bool    highPurity   [ MAX_TRACK];

        void regTree(TTree *root){//{{{
            root->Branch("TrackInfo.size"           ,&size		    ,"TrackInfo.size/I"			);
            root->Branch("TrackInfo.index"          ,index          ,"TrackInfo.index[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.handle_index"   ,handle_index   ,"TrackInfo.handle_index[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.charge"  	    ,charge         ,"TrackInfo.charge[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.pt"             ,pt             ,"TrackInfo.pt[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.eta"            ,eta            ,"TrackInfo.eta[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.phi"            ,phi            ,"TrackInfo.phi[TrackInfo.size]/D"	);
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
            root->Branch("TrackInfo.isGoodCand"     ,isGoodCand     ,"TrackInfo.isGoodCand[TrackInfo.size]/I");
            root->Branch("TrackInfo.geninfo_index"     ,geninfo_index     ,"TrackInfo.geninfo_index[TrackInfo.size]/I");
            root->Branch("TrackInfo.trackQuality"     ,trackQuality     ,"TrackInfo.trackQuality[TrackInfo.size]/I");
            root->Branch("TrackInfo.highPurity"     ,highPurity     ,"TrackInfo.highPurity[TrackInfo.size]/I");
        }//}}}

        void setbranchadd(TTree *root){//{{{
            root->SetBranchAddress("TrackInfo.size"        , &size       );
            root->SetBranchAddress("TrackInfo.index"       , index       );
            root->SetBranchAddress("TrackInfo.handle_index"       , handle_index       );
            root->SetBranchAddress("TrackInfo.charge"      , charge      );
            root->SetBranchAddress("TrackInfo.pt"          , pt          );
            root->SetBranchAddress("TrackInfo.eta"         , eta         );
            root->SetBranchAddress("TrackInfo.phi"         , phi         );
            root->SetBranchAddress("TrackInfo.striphit"    , striphit    );
            root->SetBranchAddress("TrackInfo.pixelhit"    , pixelhit    );
            root->SetBranchAddress("TrackInfo.nStripLayer" , nStripLayer );
            root->SetBranchAddress("TrackInfo.nPixelLayer" , nPixelLayer );
            root->SetBranchAddress("TrackInfo.fpbarrelhit" , fpbarrelhit );
            root->SetBranchAddress("TrackInfo.fpendcaphit" , fpendcaphit );
            root->SetBranchAddress("TrackInfo.chi2"        , chi2        );
            root->SetBranchAddress("TrackInfo.ndf"         , ndf         );
            root->SetBranchAddress("TrackInfo.d0"          , d0          );
            root->SetBranchAddress("TrackInfo.d0error"     , d0error     );
            root->SetBranchAddress("TrackInfo.dzPV"        , dzPV        );
            root->SetBranchAddress("TrackInfo.dxyPV"       , dxyPV       );
            root->SetBranchAddress("TrackInfo.isGoodCand"  , isGoodCand  );
            root->SetBranchAddress("TrackInfo.geninfo_index"  , geninfo_index  );
            root->SetBranchAddress("TrackInfo.trackQuality"  , trackQuality  );
            root->SetBranchAddress("TrackInfo.highPurity"  , highPurity  );
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
    double  uj_vtxXE[MAX_XB];
    double  uj_vtxYE[MAX_XB];
    double  uj_vtxZE[MAX_XB];
    double	uj_vtxdof[MAX_XB];
    double	uj_vtxchi2[MAX_XB];
    int     uj_rfmu1_index[MAX_XB];
    int     uj_rfmu2_index[MAX_XB];
    int     uj_isGoodCand[MAX_XB];
    
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
    double  vtxX[MAX_XB];
    double  vtxY[MAX_XB];
    double  vtxZ[MAX_XB];
    double  vtxXE[MAX_XB];
    double  vtxYE[MAX_XB];
    double  vtxZE[MAX_XB];
    double	vtxdof[MAX_XB];
    double	vtxchi2[MAX_XB];
    int     rfuj_index[MAX_XB];
    int     rftk1_index[MAX_XB];
    int     rftk2_index[MAX_XB];
    int     isGoodCand[MAX_XB];
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
    double  tktk_vtxXE[MAX_XB];
    double  tktk_vtxYE[MAX_XB];
    double  tktk_vtxZE[MAX_XB];
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
        root->Branch("BInfo.uj_vtxXE"         , uj_vtxXE       , "BInfo.uj_vtxXE[BInfo.uj_size]/D"   );
        root->Branch("BInfo.uj_vtxYE"         , uj_vtxYE       , "BInfo.uj_vtxYE[BInfo.uj_size]/D"   );
        root->Branch("BInfo.uj_vtxZE"         , uj_vtxZE       , "BInfo.uj_vtxZE[BInfo.uj_size]/D"   );
        root->Branch("BInfo.uj_vtxdof"        , uj_vtxdof      , "BInfo.uj_vtxdof[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_vtxchi2"       , uj_vtxchi2     , "BInfo.uj_vtxchi2[BInfo.uj_size]/D"	);
        root->Branch("BInfo.uj_rfmu1_index"   , uj_rfmu1_index , "BInfo.uj_rfmu1_index[BInfo.uj_size]/I"	);
        root->Branch("BInfo.uj_rfmu2_index"   , uj_rfmu2_index , "BInfo.uj_rfmu2_index[BInfo.uj_size]/I"	);
        root->Branch("BInfo.uj_isGoodCand"    , uj_isGoodCand  , "BInfo.uj_isGoodCand[BInfo.uj_size]/I"	);
        
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
        root->Branch("BInfo.vtxX"             , vtxX           , "BInfo.vtxX[BInfo.size]/D"		);
        root->Branch("BInfo.vtxY"             , vtxY           , "BInfo.vtxY[BInfo.size]/D"		);
        root->Branch("BInfo.vtxZ"             , vtxZ           , "BInfo.vtxZ[BInfo.size]/D"		);
        root->Branch("BInfo.vtxXE"            , vtxXE          , "BInfo.vtxXE[BInfo.size]/D"          );
        root->Branch("BInfo.vtxYE"            , vtxYE          , "BInfo.vtxYE[BInfo.size]/D"          );
        root->Branch("BInfo.vtxZE"            , vtxZE          , "BInfo.vtxZE[BInfo.size]/D"          );
        root->Branch("BInfo.vtxdof"           , vtxdof         , "BInfo.vtxdof[BInfo.size]/D"		);
        root->Branch("BInfo.vtxchi2"          , vtxchi2        , "BInfo.vtxchi2[BInfo.size]/D"	);
        root->Branch("BInfo.rfuj_index"       , rfuj_index     , "BInfo.rfuj_index[BInfo.size]/I");
        root->Branch("BInfo.rftk1_index"      , rftk1_index    , "BInfo.rftk1_index[BInfo.size]/I");
        root->Branch("BInfo.rftk2_index"      , rftk2_index    , "BInfo.rftk2_index[BInfo.size]/I");
        root->Branch("BInfo.isGoodCand"       , isGoodCand     , "BInfo.isGoodCand[BInfo.size]/I"	);
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
        root->Branch("BInfo.tktk_vtxXE"         , tktk_vtxXE       , "BInfo.tktk_vtxXE[BInfo.size]/D"   );
        root->Branch("BInfo.tktk_vtxYE"         , tktk_vtxYE       , "BInfo.tktk_vtxYE[BInfo.size]/D"   );
        root->Branch("BInfo.tktk_vtxZE"         , tktk_vtxZE       , "BInfo.tktk_vtxZE[BInfo.size]/D"   );
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
        root->SetBranchAddress("BInfo.uj_vtxXE"        ,uj_vtxXE   );
        root->SetBranchAddress("BInfo.uj_vtxYE"        ,uj_vtxYE   );
        root->SetBranchAddress("BInfo.uj_vtxZE"        ,uj_vtxZE   );
        root->SetBranchAddress("BInfo.uj_vtxdof"	    ,uj_vtxdof	);
        root->SetBranchAddress("BInfo.uj_vtxchi2"      ,uj_vtxchi2 );
        root->SetBranchAddress("BInfo.uj_rfmu1_index"  ,uj_rfmu1_index );
        root->SetBranchAddress("BInfo.uj_rfmu2_index"  ,uj_rfmu2_index );
        root->SetBranchAddress("BInfo.uj_isGoodCand"   ,uj_isGoodCand );
        
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
        root->SetBranchAddress("BInfo.vtxX"            ,vtxX       	);
        root->SetBranchAddress("BInfo.vtxY"            ,vtxY      	);
        root->SetBranchAddress("BInfo.vtxZ"            ,vtxZ     	);
        root->SetBranchAddress("BInfo.vtxXE"           ,vtxXE   	);
        root->SetBranchAddress("BInfo.vtxYE"           ,vtxYE        );
        root->SetBranchAddress("BInfo.vtxZE"           ,vtxZE       	);
        root->SetBranchAddress("BInfo.vtxdof"		    ,vtxdof		);
        root->SetBranchAddress("BInfo.vtxchi2"         ,vtxchi2   	);
        root->SetBranchAddress("BInfo.rfuj_index"      ,rfuj_index   	);
        root->SetBranchAddress("BInfo.rftk1_index"     ,rftk1_index   	);
        root->SetBranchAddress("BInfo.rftk2_index"     ,rftk2_index   	);
        root->SetBranchAddress("BInfo.isGoodCand"      ,isGoodCand   	);
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
        root->SetBranchAddress("BInfo.tktk_vtxXE"        ,tktk_vtxXE   );
        root->SetBranchAddress("BInfo.tktk_vtxYE"        ,tktk_vtxYE   );
        root->SetBranchAddress("BInfo.tktk_vtxZE"        ,tktk_vtxZE   );
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
        }//}}}
};//}}}
#endif
