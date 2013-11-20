// vim:set ts=4 sw=4 fdm=marker et:
// Update:
// 2012Mar28    pchen       Add GenInfoBranches,XbInfo.rfuj_index.
//                          Fix Trgresult.
// 2012Nov14    pchen       Kill following branches to reduce size
//                              EvtInfo.PVc2p/F
//                              MuonInfo.p/F
//                              MuonInfo.caloIso/F
//                              TrackInfo.p/F
//                              XbInfo.rfmu1_index/F
//                              XbInfo.rfmu2_index/F
//                              XbInfo.rfmu1_e/F
//                              XbInfo.rfmu2_e/F
//                              XbInfo.rftk1_e/F
//                              XbInfo.rftk2_e/F
//                              XbInfo.jpsi_*
//                          Rename XbInfo.n_uj as XbInfo.uj_size. 
//                          Add MuonInfo.isGoodCand/I
//                              TrackInfo.isGoodCand/I
//                              XbInfo.uj_isGoodCand/I      
//                              XbInfo.isGoodCand/I         
//                          Add option to fit to jpsi/upsilon separately or together.
// 2012Nov28    pchen       Add MuonInfo.i_nStripLayer/I
//                              MuonInfo.i_nPixelLayer/I
// 2013Jan22    pchen       Add MAX_BX 128
//                              EvtInfo.nBX/I
//                              EvtInfo.BXPU[EvtInfo.nBX]/I
//                              EvtInfo.nPU[EvtInfo.nBX]/I
//                              EvtInfo.trueIT[EvtInfo.nBX]/F
//                          Fix GenInfo signal matching
//                              Due to PHOTOS option in evtgen, could have extra photon in (1S)->mumu decay.
#ifndef _XBFRAMEFORMAT_H_
#define _XBFRAMEFORMAT_H_

//#define MAX_XB 2560
#define MAX_XB 8192
#define MAX_MUON 64
#define MAX_TRACK 4096 //default 2048
#define MAX_GEN 4096 //default 2048
#define MAX_BX 128
#define N_TRIGGER_BOOKINGS 5842

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
        int     trgCount;                   //number of successfully triggered HLT path in the booking.
        int     nTrgBook;                   //N_TRIGGER_BOOKING
        char    trgBook[N_TRIGGER_BOOKINGS];//status of booked triggers
		int 	nHLT;                       //# of HLT of the event 
        bool    hltBits[N_TRIGGER_BOOKINGS];//is HLT of the event acceptted?
        int     nBX;
        int     BXPU[MAX_BX];
        int     nPU[MAX_BX];
        float   trueIT[MAX_BX];
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
            root->Branch("EvtInfo.trgCount"   , &trgCount                  , "EvtInfo.trgCount/I"       );
            root->Branch("EvtInfo.nTrgBook"   , &nTrgBook                  , "EvtInfo.nTrgBook/I"       );
            root->Branch("EvtInfo.trgBook"    , trgBook                    , "EvtInfo.trgBook[EvtInfo.nTrgBook]/B");//notice /B
			root->Branch("EvtInfo.nHLT"       , &nHLT                      , "EvtInfo.nHLT/I"			);
            root->Branch("EvtInfo.hltBits"    , hltBits                    , "EvtInfo.hltBits[EvtInfo.nHLT]/O");
			root->Branch("EvtInfo.nBX"        , &nBX                       , "EvtInfo.nBX/I" 			);
			root->Branch("EvtInfo.BXPU"       , BXPU                       , "EvtInfo.BXPU[EvtInfo.nBX]/I");
			root->Branch("EvtInfo.nPU"        , nPU                        , "EvtInfo.nPU[EvtInfo.nBX]/I");
			root->Branch("EvtInfo.trueIT"     , trueIT                     , "EvtInfo.trueIT[EvtInfo.nBX]/F");
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
            root->SetBranchAddress("EvtInfo.trgCount" ,&trgCount    );
            root->SetBranchAddress("EvtInfo.nTrgBook" ,&nTrgBook    );
            root->SetBranchAddress("EvtInfo.trgBook"  ,trgBook      );
            root->SetBranchAddress("EvtInfo.nHLT"     ,&nHLT	    );
            root->SetBranchAddress("EvtInfo.hltBits"  ,hltBits      );
			root->SetBranchAddress("EvtInfo.nBX"      ,&nBX         );
			root->SetBranchAddress("EvtInfo.BXPU"     ,BXPU         );
			root->SetBranchAddress("EvtInfo.nPU"      ,nPU          );
			root->SetBranchAddress("EvtInfo.trueIT"   ,trueIT       );
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

class MuonInfoBranches{//{{{
    public:
        int	    size;
        int     index        [ MAX_MUON];
        int	    handle_index [ MAX_MUON];
        int 	charge       [ MAX_MUON];
        double 	pt           [ MAX_MUON];
        double	eta          [ MAX_MUON];
        double 	phi          [ MAX_MUON];
        //double  p            [ MAX_MUON];
        int	    i_striphit   [ MAX_MUON];
        int	    i_pixelhit   [ MAX_MUON];
        int	    g_striphit   [ MAX_MUON];
        int	    g_pixelhit   [ MAX_MUON];
        int     i_nStripLayer[ MAX_MUON];
        int     i_nPixelLayer[ MAX_MUON];
        double	i_chi2       [ MAX_MUON];
        double	i_ndf        [ MAX_MUON];
        double	g_chi2       [ MAX_MUON];
        double	g_ndf        [ MAX_MUON];
        int	    nmuhit       [ MAX_MUON];
        double	d0           [ MAX_MUON];
        double	dz           [ MAX_MUON];
        double	dzPV         [ MAX_MUON];
        double	dxyPV        [ MAX_MUON];
        int	    fpbarrelhit  [ MAX_MUON];
        int	    fpendcaphit  [ MAX_MUON];
        int	    muqual       [ MAX_MUON];
        double  iso_trk      [ MAX_MUON];
        //double  iso_calo     [ MAX_MUON];
        double  iso_ecal     [ MAX_MUON];
        double  iso_hcal     [ MAX_MUON];
        double  n_matches    [ MAX_MUON];
        int     isGoodCand   [ MAX_MUON];

        void regTree(TTree *root){//{{{
            root->Branch("MuonInfo.size"          , &size         , "MuonInfo.size/I"			);
            root->Branch("MuonInfo.index"         , index         , "MuonInfo.index[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.handle_index"  , handle_index  , "MuonInfo.handle_index[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.charge"        , charge        , "MuonInfo.charge[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.pt"            , pt            , "MuonInfo.pt[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.eta"           , eta           , "MuonInfo.eta[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.phi"           , phi           , "MuonInfo.phi[MuonInfo.size]/D"	);
            //root->Branch("MuonInfo.p"           , p             , "MuonInfo.p[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.i_striphit"    , i_striphit    , "MuonInfo.i_striphit[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.i_pixelhit"    , i_pixelhit    , "MuonInfo.i_pixelhit[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.i_nStripLayer" , i_nStripLayer , "MuonInfo.i_nStripLayer[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.i_nPixelLayer" , i_nPixelLayer , "MuonInfo.i_nPixelLayer[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.g_striphit"    , g_striphit    , "MuonInfo.g_striphit[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.g_pixelhit"    , g_pixelhit    , "MuonInfo.g_pixelhit[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.i_chi2"        , i_chi2        , "MuonInfo.i_chi2[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.i_ndf"         , i_ndf         , "MuonInfo.i_ndf[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.g_chi2"        , g_chi2        , "MuonInfo.g_chi2[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.g_ndf"         , g_ndf         , "MuonInfo.g_ndf[MuonInfo.size]/D"	);
            root->Branch("MuonInfo.nmuhit"        , nmuhit        , "MuonInfo.nmuhit[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.d0"            , d0            , "MuonInfo.d0[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.dz"            , dz            , "MuonInfo.dz[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.dzPV"          , dzPV          , "MuonInfo.dzPV[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.dxyPV"         , dxyPV         , "MuonInfo.dxyPV[MuonInfo.size]/D"		);
            root->Branch("MuonInfo.fpbarrelhit"   , fpbarrelhit   , "MuonInfo.fpbarrelhit[MuonInfo.size]/I");
            root->Branch("MuonInfo.fpendcaphit"   , fpendcaphit   , "MuonInfo.fpendcaphit[MuonInfo.size]/I");
            root->Branch("MuonInfo.muqual"        , muqual        , "MuonInfo.muqual[MuonInfo.size]/I"	);
            root->Branch("MuonInfo.iso_trk"       , iso_trk       , "MuonInfo.iso_trk[MuonInfo.size]/D");
            //root->Branch("MuonInfo.iso_calo"    , iso_calo      , "MuonInfo.iso_calo[MuonInfo.size]/D");
            root->Branch("MuonInfo.iso_ecal"      , iso_ecal      , "MuonInfo.iso_ecal[MuonInfo.size]/D");
            root->Branch("MuonInfo.iso_hcal"      , iso_hcal      , "MuonInfo.iso_hcal[MuonInfo.size]/D");
            root->Branch("MuonInfo.n_matches"     , n_matches     , "MuonInfo.n_matches[MuonInfo.size]/D");
            root->Branch("MuonInfo.isGoodCand"    , isGoodCand    , "MuonInfo.isGoodCand[MuonInfo.size]/I");
        }//}}}

        void setbranchadd(TTree *root){//{{{
            root->SetBranchAddress("MuonInfo.size"          , &size          );
            root->SetBranchAddress("MuonInfo.index"         , index          );
            root->SetBranchAddress("MuonInfo.handle_index"  , handle_index          );
            root->SetBranchAddress("MuonInfo.charge"        , charge         );
            root->SetBranchAddress("MuonInfo.pt"            , pt             );
            root->SetBranchAddress("MuonInfo.eta"           , eta            );
            root->SetBranchAddress("MuonInfo.phi"           , phi            );
            //root->SetBranchAddress("MuonInfo.p"           , p              );
            root->SetBranchAddress("MuonInfo.i_striphit"    , i_striphit	);
            root->SetBranchAddress("MuonInfo.i_pixelhit"    , i_pixelhit	);
            root->SetBranchAddress("MuonInfo.i_nStripLayer" , i_nStripLayer);
            root->SetBranchAddress("MuonInfo.i_nPixelLayer" , i_nPixelLayer);
            root->SetBranchAddress("MuonInfo.g_striphit"    , g_striphit	);
            root->SetBranchAddress("MuonInfo.i_pixelhit"    , g_pixelhit	);
            root->SetBranchAddress("MuonInfo.i_chi2"        , i_chi2		);
            root->SetBranchAddress("MuonInfo.i_ndf"         , i_ndf		);
            root->SetBranchAddress("MuonInfo.g_chi2"        , g_chi2		);
            root->SetBranchAddress("MuonInfo.g_ndf"         , g_ndf		);
            root->SetBranchAddress("MuonInfo.nmuhit"        , nmuhit		);
            root->SetBranchAddress("MuonInfo.d0"            , d0             );
            root->SetBranchAddress("MuonInfo.dz"            , dz             );
            root->SetBranchAddress("MuonInfo.dzPV"          , dzPV             );
            root->SetBranchAddress("MuonInfo.dxyPV"         , dxyPV             );
            root->SetBranchAddress("MuonInfo.fpbarrelhit"   , fpbarrelhit    );
            root->SetBranchAddress("MuonInfo.fpendcaphit"   , fpendcaphit    );
            root->SetBranchAddress("MuonInfo.muqual"        , muqual         );
            root->SetBranchAddress("MuonInfo.iso_trk"       , iso_trk);
            //root->SetBranchAddress("MuonInfo.iso_calo"    , iso_calo);
            root->SetBranchAddress("MuonInfo.iso_ecal"      , iso_ecal);
            root->SetBranchAddress("MuonInfo.iso_hcal"      , iso_hcal);
            root->SetBranchAddress("MuonInfo.n_matches"     , n_matches);
            root->SetBranchAddress("MuonInfo.isGoodCand"    , isGoodCand);

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
        int	    fpbarrelhit  [ MAX_TRACK];
        int	    fpendcaphit  [ MAX_TRACK];
        double	chi2         [ MAX_TRACK];
        double	ndf          [ MAX_TRACK];
        double	d0           [ MAX_TRACK];
        double	d0error      [ MAX_TRACK];
        double	dzPV         [ MAX_TRACK];
        double	dxyPV        [ MAX_TRACK];
        int     isGoodCand   [ MAX_TRACK];

        void regTree(TTree *root){//{{{
            root->Branch("TrackInfo.size"           ,&size		    ,"TrackInfo.size/I"			);
            root->Branch("TrackInfo.index"          ,index          ,"TrackInfo.index[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.handle_index"   ,handle_index   ,"TrackInfo.handle_index[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.charge"  	    ,charge         ,"TrackInfo.charge[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.pt"             ,pt             ,"TrackInfo.pt[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.eta"            ,eta            ,"TrackInfo.eta[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.phi"            ,phi            ,"TrackInfo.phi[TrackInfo.size]/D"	);
            //root->Branch("TrackInfo.p"              ,p              ,"TrackInfo.p[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.striphit"	    ,striphit	    ,"TrackInfo.striphit[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.pixelhit"	    ,pixelhit	    ,"TrackInfo.pixelhit[TrackInfo.size]/I"	);
            root->Branch("TrackInfo.fpbarrelhit"	,fpbarrelhit	,"TrackInfo.fpbarrelhit[TrackInfo.size]/I");
            root->Branch("TrackInfo.fpendcaphit"	,fpendcaphit	,"TrackInfo.fpendcaphit[TrackInfo.size]/I");
            root->Branch("TrackInfo.chi2"		    ,chi2		    ,"TrackInfo.chi2[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.ndf"		    ,ndf		    ,"TrackInfo.ndf[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.d0"		        ,d0		        ,"TrackInfo.d0[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.d0error"	    ,d0error	    ,"TrackInfo.d0error[TrackInfo.size]/D"	);
            root->Branch("TrackInfo.dzPV"           ,dzPV           ,"TrackInfo.dzPV[TrackInfo.size]/D"		);
            root->Branch("TrackInfo.dxyPV"          ,dxyPV          ,"TrackInfo.dxyPV[TrackInfo.size]/D"		);
            root->Branch("TrackInfo.isGoodCand"     ,isGoodCand     ,"TrackInfo.isGoodCand[TrackInfo.size]/I");
        }//}}}

        void setbranchadd(TTree *root){//{{{
            root->SetBranchAddress("TrackInfo.size"        , &size       );
            root->SetBranchAddress("TrackInfo.index"       , index       );
            root->SetBranchAddress("TrackInfo.handle_index"       , handle_index       );
            root->SetBranchAddress("TrackInfo.charge"      , charge      );
            root->SetBranchAddress("TrackInfo.pt"          , pt          );
            root->SetBranchAddress("TrackInfo.eta"         , eta         );
            root->SetBranchAddress("TrackInfo.phi"         , phi         );
            //root->SetBranchAddress("TrackInfo.p"           , p           );
            root->SetBranchAddress("TrackInfo.striphit"    , striphit    );
            root->SetBranchAddress("TrackInfo.pixelhit"    , pixelhit    );
            root->SetBranchAddress("TrackInfo.fpbarrelhit" , fpbarrelhit );
            root->SetBranchAddress("TrackInfo.fpendcaphit" , fpendcaphit );
            root->SetBranchAddress("TrackInfo.chi2"        , chi2        );
            root->SetBranchAddress("TrackInfo.ndf"         , ndf         );
            root->SetBranchAddress("TrackInfo.d0"          , d0          );
            root->SetBranchAddress("TrackInfo.d0error"     , d0error     );
            root->SetBranchAddress("TrackInfo.dzPV"        , dzPV        );
            root->SetBranchAddress("TrackInfo.dxyPV"       , dxyPV       );
            root->SetBranchAddress("TrackInfo.isGoodCand"  , isGoodCand  );
        }//}}}
};//}}}

class XbInfoBranches{//{{{
    public:
        int	    uj_size;
        int     n_uj;//Compatible with uj_size
        int	    uj_index[MAX_XB];
        double  uj_mass[MAX_XB];
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
        //int     rfmu1_index[MAX_XB];
        //int     rfmu2_index[MAX_XB];
        int     rftk1_index[MAX_XB];
        int     rftk2_index[MAX_XB];
        int     isGoodCand[MAX_XB];

        double  rfmu1_px[MAX_XB];
        double  rfmu1_py[MAX_XB];
        double  rfmu1_pz[MAX_XB];
        //double  rfmu1_e[MAX_XB];
        double  rfmu2_px[MAX_XB];
        double  rfmu2_py[MAX_XB];
        double  rfmu2_pz[MAX_XB];
        //double  rfmu2_e[MAX_XB];
        double  rftk1_px[MAX_XB];
        double  rftk1_py[MAX_XB];
        double  rftk1_pz[MAX_XB];
        //double  rftk1_e[MAX_XB];
        double  rftk2_px[MAX_XB];
        double  rftk2_py[MAX_XB];
        double  rftk2_pz[MAX_XB];
        //double  rftk2_e[MAX_XB];//38
    
        ////Jpsi_fit
        //int	    jpsi_fitstatus[MAX_XB];
        //double	jpsi_mass[MAX_XB];
        //double	jpsi_px[MAX_XB];
        //double	jpsi_py[MAX_XB];
        //double	jpsi_pz[MAX_XB];
        //double    jpsi_vtxX[MAX_XB];
        //double    jpsi_vtxY[MAX_XB];
        //double    jpsi_vtxZ[MAX_XB];
        //double	jpsi_vtxdof[MAX_XB];
        //double	jpsi_vtxchi2[MAX_XB];//10


        void regTree(TTree *root){//{{{
            root->Branch("XbInfo.uj_size"          , &uj_size       , "XbInfo.uj_size/I"			);
            root->Branch("XbInfo.uj_index"         , uj_index       , "XbInfo.uj_index[XbInfo.uj_size]/I"	);
            root->Branch("XbInfo.uj_mass"          , uj_mass        , "XbInfo.uj_mass[XbInfo.uj_size]/D"	);
            root->Branch("XbInfo.uj_px"            , uj_px          , "XbInfo.uj_px[XbInfo.uj_size]/D"	);
            root->Branch("XbInfo.uj_py"            , uj_py          , "XbInfo.uj_py[XbInfo.uj_size]/D"	);
            root->Branch("XbInfo.uj_pz"            , uj_pz          , "XbInfo.uj_pz[XbInfo.uj_size]/D"	);
            root->Branch("XbInfo.uj_vtxX"          , uj_vtxX        , "XbInfo.uj_vtxX[XbInfo.uj_size]/D"	);
            root->Branch("XbInfo.uj_vtxY"          , uj_vtxY        , "XbInfo.uj_vtxY[XbInfo.uj_size]/D"	);
            root->Branch("XbInfo.uj_vtxZ"          , uj_vtxZ        , "XbInfo.uj_vtxZ[XbInfo.uj_size]/D"	);
            root->Branch("XbInfo.uj_vtxXE"         , uj_vtxXE       , "XbInfo.uj_vtxXE[XbInfo.uj_size]/D"   );
            root->Branch("XbInfo.uj_vtxYE"         , uj_vtxYE       , "XbInfo.uj_vtxYE[XbInfo.uj_size]/D"   );
            root->Branch("XbInfo.uj_vtxZE"         , uj_vtxZE       , "XbInfo.uj_vtxZE[XbInfo.uj_size]/D"   );
            root->Branch("XbInfo.uj_vtxdof"        , uj_vtxdof      , "XbInfo.uj_vtxdof[XbInfo.uj_size]/D"	);
            root->Branch("XbInfo.uj_vtxchi2"       , uj_vtxchi2     , "XbInfo.uj_vtxchi2[XbInfo.uj_size]/D"	);
            root->Branch("XbInfo.uj_rfmu1_index"   , uj_rfmu1_index , "XbInfo.uj_rfmu1_index[XbInfo.uj_size]/I"	);
            root->Branch("XbInfo.uj_rfmu2_index"   , uj_rfmu2_index , "XbInfo.uj_rfmu2_index[XbInfo.uj_size]/I"	);
            root->Branch("XbInfo.uj_isGoodCand"    , uj_isGoodCand  , "XbInfo.uj_isGoodCand[XbInfo.uj_size]/I"	);

            root->Branch("XbInfo.uj_rfmu1_px"      , uj_rfmu1_px    , "XbInfo.uj_rfmu1_px[XbInfo.uj_size]/D");
            root->Branch("XbInfo.uj_rfmu1_py"      , uj_rfmu1_py    , "XbInfo.uj_rfmu1_py[XbInfo.uj_size]/D");
            root->Branch("XbInfo.uj_rfmu1_pz"      , uj_rfmu1_pz    , "XbInfo.uj_rfmu1_pz[XbInfo.uj_size]/D");
            root->Branch("XbInfo.uj_rfmu2_px"      , uj_rfmu2_px    , "XbInfo.uj_rfmu2_px[XbInfo.uj_size]/D");
            root->Branch("XbInfo.uj_rfmu2_py"      , uj_rfmu2_py    , "XbInfo.uj_rfmu2_py[XbInfo.uj_size]/D");
            root->Branch("XbInfo.uj_rfmu2_pz"      , uj_rfmu2_pz    , "XbInfo.uj_rfmu2_pz[XbInfo.uj_size]/D");

            root->Branch("XbInfo.size"             , &size          , "XbInfo.size/I"			);
            root->Branch("XbInfo.index"            , index          , "XbInfo.index[XbInfo.size]/I"		);
            root->Branch("XbInfo.mass"             , mass           , "XbInfo.mass[XbInfo.size]/D"		);
            root->Branch("XbInfo.px"               , px             , "XbInfo.px[XbInfo.size]/D"		);
            root->Branch("XbInfo.py"               , py             , "XbInfo.py[XbInfo.size]/D"		);
            root->Branch("XbInfo.pz"               , pz             , "XbInfo.pz[XbInfo.size]/D"		);
            root->Branch("XbInfo.pxE"              , pxE            , "XbInfo.pxE[XbInfo.size]/D"            );
            root->Branch("XbInfo.pyE"              , pyE            , "XbInfo.pyE[XbInfo.size]/D"            );
            root->Branch("XbInfo.pzE"              , pzE            , "XbInfo.pzE[XbInfo.size]/D"            );
            root->Branch("XbInfo.vtxX"             , vtxX           , "XbInfo.vtxX[XbInfo.size]/D"		);
            root->Branch("XbInfo.vtxY"             , vtxY           , "XbInfo.vtxY[XbInfo.size]/D"		);
            root->Branch("XbInfo.vtxZ"             , vtxZ           , "XbInfo.vtxZ[XbInfo.size]/D"		);
            root->Branch("XbInfo.vtxXE"            , vtxXE          , "XbInfo.vtxXE[XbInfo.size]/D"          );
            root->Branch("XbInfo.vtxYE"            , vtxYE          , "XbInfo.vtxYE[XbInfo.size]/D"          );
            root->Branch("XbInfo.vtxZE"            , vtxZE          , "XbInfo.vtxZE[XbInfo.size]/D"          );
            root->Branch("XbInfo.vtxdof"           , vtxdof         , "XbInfo.vtxdof[XbInfo.size]/D"		);
            root->Branch("XbInfo.vtxchi2"          , vtxchi2        , "XbInfo.vtxchi2[XbInfo.size]/D"	);
            root->Branch("XbInfo.rfuj_index"       , rfuj_index     , "XbInfo.rfuj_index[XbInfo.size]/I");
            //root->Branch("XbInfo.rfmu1_index"    , rfmu1_index    , "XbInfo.rfmu1_index[XbInfo.size]/I");
            //root->Branch("XbInfo.rfmu2_index"    , rfmu2_index    , "XbInfo.rfmu2_index[XbInfo.size]/I");
            root->Branch("XbInfo.rftk1_index"      , rftk1_index    , "XbInfo.rftk1_index[XbInfo.size]/I");
            root->Branch("XbInfo.rftk2_index"      , rftk2_index    , "XbInfo.rftk2_index[XbInfo.size]/I");
            root->Branch("XbInfo.isGoodCand"       , isGoodCand     , "XbInfo.isGoodCand[XbInfo.size]/I"	);

            root->Branch("XbInfo.rfmu1_px"         , rfmu1_px       , "XbInfo.rfmu1_px[XbInfo.size]/D"	);
            root->Branch("XbInfo.rfmu1_py"         , rfmu1_py       , "XbInfo.rfmu1_py[XbInfo.size]/D"	);
            root->Branch("XbInfo.rfmu1_pz"         , rfmu1_pz       , "XbInfo.rfmu1_pz[XbInfo.size]/D"	);
            //root->Branch("XbInfo.rfmu1_e"        , rfmu1_e        , "XbInfo.rfmu1_e[XbInfo.size]/D" 	);
            root->Branch("XbInfo.rfmu2_px"         , rfmu2_px       , "XbInfo.rfmu2_px[XbInfo.size]/D"	);
            root->Branch("XbInfo.rfmu2_py"         , rfmu2_py       , "XbInfo.rfmu2_py[XbInfo.size]/D"	);
            root->Branch("XbInfo.rfmu2_pz"         , rfmu2_pz       , "XbInfo.rfmu2_pz[XbInfo.size]/D"	);
            //root->Branch("XbInfo.rfmu2_e"        , rfmu2_e        , "XbInfo.rfmu2_e[XbInfo.size]/D" 	);
            root->Branch("XbInfo.rftk1_px"         , rftk1_px       , "XbInfo.rftk1_px[XbInfo.size]/D"     );
            root->Branch("XbInfo.rftk1_py"         , rftk1_py       , "XbInfo.rftk1_py[XbInfo.size]/D"     );
            root->Branch("XbInfo.rftk1_pz"         , rftk1_pz       , "XbInfo.rftk1_pz[XbInfo.size]/D"     );
            //root->Branch("XbInfo.rftk1_e"        , rftk1_e        , "XbInfo.rftk1_e[XbInfo.size]/D"      );
            root->Branch("XbInfo.rftk2_px"         , rftk2_px       , "XbInfo.rftk2_px[XbInfo.size]/D"     );
            root->Branch("XbInfo.rftk2_py"         , rftk2_py       , "XbInfo.rftk2_py[XbInfo.size]/D"     );
            root->Branch("XbInfo.rftk2_pz"         , rftk2_pz       , "XbInfo.rftk2_pz[XbInfo.size]/D"     );
            //root->Branch("XbInfo.rftk2_e"        , rftk2_e        , "XbInfo.rftk2_e[XbInfo.size]/D"      );

            //root->Branch("XbInfo.jpsi_fitstatus" , jpsi_fitstatus , "XbInfo.jpsi_fitstatus[XbInfo.size]/I");
            //root->Branch("XbInfo.jpsi_mass"      , jpsi_mass      , "XbInfo.jpsi_mass[XbInfo.size]/D"		);
            //root->Branch("XbInfo.jpsi_px"        , jpsi_px        , "XbInfo.jpsi_px[XbInfo.size]/D"		);
            //root->Branch("XbInfo.jpsi_py"        , jpsi_py        , "XbInfo.jpsi_py[XbInfo.size]/D"		);
            //root->Branch("XbInfo.jpsi_pz"        , jpsi_pz        , "XbInfo.jpsi_pz[XbInfo.size]/D"		);
            //root->Branch("XbInfo.jpsi_vtxX"      , jpsi_vtxX      , "XbInfo.jpsi_vtxX[XbInfo.size]/D"		);
            //root->Branch("XbInfo.jpsi_vtxY"      , jpsi_vtxY      , "XbInfo.jpsi_vtxY[XbInfo.size]/D"		);
            //root->Branch("XbInfo.jpsi_vtxZ"      , jpsi_vtxZ      , "XbInfo.jpsi_vtxZ[XbInfo.size]/D"		);
            //root->Branch("XbInfo.jpsi_vtxdof"    , jpsi_vtxdof    , "XbInfo.jpsi_vtxdof[XbInfo.size]/D"		);
            //root->Branch("XbInfo.jpsi_vtxchi2"   , jpsi_vtxchi2   , "XbInfo.jpsi_vtxchi2[XbInfo.size]/D"	);
        }//}}}

        void setbranchadd(TTree *root){//{{{
            root->SetBranchAddress("XbInfo.uj_size"		    ,&uj_size	);   
            root->SetBranchAddress("XbInfo.uj_size"		    ,&n_uj	    );   
            root->SetBranchAddress("XbInfo.uj_index"        ,uj_index   ); 
            root->SetBranchAddress("XbInfo.uj_mass"         ,uj_mass   	); 
            root->SetBranchAddress("XbInfo.uj_px"           ,uj_px     	);
            root->SetBranchAddress("XbInfo.uj_py"           ,uj_py    	);
            root->SetBranchAddress("XbInfo.uj_pz"           ,uj_pz   	);
            root->SetBranchAddress("XbInfo.uj_vtxX"         ,uj_vtxX    );  
            root->SetBranchAddress("XbInfo.uj_vtxY"         ,uj_vtxY    ); 
            root->SetBranchAddress("XbInfo.uj_vtxZ"         ,uj_vtxZ    );
            root->SetBranchAddress("XbInfo.uj_vtxXE"        ,uj_vtxXE   );  
            root->SetBranchAddress("XbInfo.uj_vtxYE"        ,uj_vtxYE   ); 
            root->SetBranchAddress("XbInfo.uj_vtxZE"        ,uj_vtxZE   );  
            root->SetBranchAddress("XbInfo.uj_vtxdof"	    ,uj_vtxdof	);
            root->SetBranchAddress("XbInfo.uj_vtxchi2"      ,uj_vtxchi2 );  
            root->SetBranchAddress("XbInfo.uj_rfmu1_index"  ,uj_rfmu1_index );  
            root->SetBranchAddress("XbInfo.uj_rfmu2_index"  ,uj_rfmu2_index );  
            root->SetBranchAddress("XbInfo.uj_isGoodCand"   ,uj_isGoodCand );  
            //root->SetBranchAddress("XbInfo.uj_rfmu1_px"     ,uj_rfmu1_px);  
            //root->SetBranchAddress("XbInfo.uj_rfmu1_py"     ,uj_rfmu1_py);  
            //root->SetBranchAddress("XbInfo.uj_rfmu1_pz"     ,uj_rfmu1_pz);  
            //root->SetBranchAddress("XbInfo.uj_rfmu2_px"     ,uj_rfmu2_px); 
            //root->SetBranchAddress("XbInfo.uj_rfmu2_py"     ,uj_rfmu2_py); 
            //root->SetBranchAddress("XbInfo.uj_rfmu2_pz"     ,uj_rfmu2_pz);

            root->SetBranchAddress("XbInfo.size"            ,&size        );
            root->SetBranchAddress("XbInfo.index"           ,index       	);
            root->SetBranchAddress("XbInfo.mass"		    ,mass		);
            root->SetBranchAddress("XbInfo.px"              ,px         	);
            root->SetBranchAddress("XbInfo.py"              ,py        	);
            root->SetBranchAddress("XbInfo.pz"              ,pz           );  
            root->SetBranchAddress("XbInfo.pxE"             ,pxE          ); 
            root->SetBranchAddress("XbInfo.pyE"             ,pyE          );
            root->SetBranchAddress("XbInfo.pzE"             ,pzE         	);
            root->SetBranchAddress("XbInfo.vtxX"            ,vtxX       	);
            root->SetBranchAddress("XbInfo.vtxY"            ,vtxY      	);
            root->SetBranchAddress("XbInfo.vtxZ"            ,vtxZ     	);
            root->SetBranchAddress("XbInfo.vtxXE"           ,vtxXE   	);
            root->SetBranchAddress("XbInfo.vtxYE"           ,vtxYE        );
            root->SetBranchAddress("XbInfo.vtxZE"           ,vtxZE       	);
            root->SetBranchAddress("XbInfo.vtxdof"		    ,vtxdof		);	
            root->SetBranchAddress("XbInfo.vtxchi2"         ,vtxchi2   	);
            root->SetBranchAddress("XbInfo.rfuj_index"      ,rfuj_index   	);
            //root->SetBranchAddress("XbInfo.rfmu1_index"     ,rfmu1_index   	);
            //root->SetBranchAddress("XbInfo.rfmu2_index"     ,rfmu2_index   	);
            root->SetBranchAddress("XbInfo.rftk1_index"     ,rftk1_index   	);
            root->SetBranchAddress("XbInfo.rftk2_index"     ,rftk2_index   	);
            root->SetBranchAddress("XbInfo.isGoodCand"      ,isGoodCand   	);
            root->SetBranchAddress("XbInfo.rfmu1_px"        ,rfmu1_px 	);
            root->SetBranchAddress("XbInfo.rfmu1_py"        ,rfmu1_py	);
            root->SetBranchAddress("XbInfo.rfmu1_pz"        ,rfmu1_pz     );
            //root->SetBranchAddress("XbInfo.rfmu1_e"         ,rfmu1_e      );
            root->SetBranchAddress("XbInfo.rfmu2_px"        ,rfmu2_px     );
            root->SetBranchAddress("XbInfo.rfmu2_py"        ,rfmu2_py    	);
            root->SetBranchAddress("XbInfo.rfmu2_pz"        ,rfmu2_pz   	);
            //root->SetBranchAddress("XbInfo.rfmu2_e"         ,rfmu2_e   	);
            root->SetBranchAddress("XbInfo.rftk1_px"        ,rftk1_px 	);
            root->SetBranchAddress("XbInfo.rftk1_py"        ,rftk1_py     );
            root->SetBranchAddress("XbInfo.rftk1_pz"        ,rftk1_pz     );
            //root->SetBranchAddress("XbInfo.rftk1_e"         ,rftk1_e      );
            root->SetBranchAddress("XbInfo.rftk2_px"        ,rftk2_px    	);
            root->SetBranchAddress("XbInfo.rftk2_py"        ,rftk2_py   	);
            root->SetBranchAddress("XbInfo.rftk2_pz"        ,rftk2_pz  	);
            //root->SetBranchAddress("XbInfo.rftk2_e"         ,rftk2_e  	);
            

            //root->SetBranchAddress("XbInfo.jpsi_fitstatus"  ,jpsi_fitstatus);
            //root->SetBranchAddress("XbInfo.jpsi_mass"		  ,jpsi_mass		);
            //root->SetBranchAddress("XbInfo.jpsi_px"         ,jpsi_px         	);
            //root->SetBranchAddress("XbInfo.jpsi_py"         ,jpsi_py        	);
            //root->SetBranchAddress("XbInfo.jpsi_pz"         ,jpsi_pz           );  
            //root->SetBranchAddress("XbInfo.jpsi_vtxX"       ,jpsi_vtxX       	);
            //root->SetBranchAddress("XbInfo.jpsi_vtxY"       ,jpsi_vtxY      	);
            //root->SetBranchAddress("XbInfo.jpsi_vtxZ"       ,jpsi_vtxZ     	);
            //root->SetBranchAddress("XbInfo.jpsi_vtxdof"	  ,jpsi_vtxdof		);	
            //root->SetBranchAddress("XbInfo.jpsi_vtxchi2"    ,jpsi_vtxchi2   	);
        }//}}}
};//}}}

class BInfoBranches{//{{{
public:
    int	    uj_size;
    int     n_uj;
    int	    uj_index[MAX_XB];
    double  uj_mass[MAX_XB];
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
    
    void regTree(TTree *root){//{{{
        root->Branch("BInfo.uj_size"          , &uj_size       , "BInfo.uj_size/I"			);
        root->Branch("BInfo.uj_index"         , uj_index       , "BInfo.uj_index[BInfo.uj_size]/I"	);
        root->Branch("BInfo.uj_mass"          , uj_mass        , "BInfo.uj_mass[BInfo.uj_size]/D"	);
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
        
    }//}}}
    
    void setbranchadd(TTree *root){//{{{
        root->SetBranchAddress("BInfo.uj_size"		    ,&uj_size	);
        root->SetBranchAddress("BInfo.uj_size"		    ,&n_uj	    );
        root->SetBranchAddress("BInfo.uj_index"        ,uj_index   );
        root->SetBranchAddress("BInfo.uj_mass"         ,uj_mass   	);
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
        int     mhmu1_index;
        int     mhmu2_index;
        int     mhtk1_index;
        int     mhtk2_index;
        double  mhujMass;
        double  mhxbMass;

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
            root->Branch("Geninfo.mhmu1_index"  ,&mhmu1_index   ,"GenInfo.mhmu1_index/I");
            root->Branch("Geninfo.mhmu2_index"  ,&mhmu2_index   ,"GenInfo.mhmu2_index/I");
            root->Branch("Geninfo.mhtk1_index"  ,&mhtk1_index   ,"GenInfo.mhtk1_index/I");
            root->Branch("Geninfo.mhtk2_index"  ,&mhtk2_index   ,"GenInfo.mhtk2_index/I");
            root->Branch("GenInfo.mhujMass"     ,&mhujMass      ,"GenInfo.mhujMass/D");
            root->Branch("GenInfo.mhxbMass"     ,&mhxbMass      ,"GenInfo.mhxbMass/D");
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
            root->SetBranchAddress("Geninfo.mhmu1_index"  ,&mhmu1_index   );
            root->SetBranchAddress("Geninfo.mhmu2_index"  ,&mhmu2_index   );
            root->SetBranchAddress("Geninfo.mhtk1_index"  ,&mhtk1_index   );
            root->SetBranchAddress("Geninfo.mhtk2_index"  ,&mhtk2_index   );
            root->SetBranchAddress("GenInfo.mhujMass"     ,&mhujMass      );
            root->SetBranchAddress("GenInfo.mhxbMass"     ,&mhxbMass      );
        }//}}}
};//}}}
#endif
