// vim:set ts=4 sw=4 fdm=marker et:
// Ntuplt creator for B meson related analysis.
// Maintain and contact: ta-wei wang
// Email: "tawei@mit.edu" or "ta-wei.wang@cern.ch"
#include "Bfinder/Bfinder/interface/format.h"
#include "Bfinder/Bfinder/interface/Bntuple.h"
#include "Bfinder/Bfinder/interface/utilities.h"
//
// class declaration
//

class Bfinder : public edm::EDAnalyzer
{//{{{
    public:
        explicit Bfinder(const edm::ParameterSet&);
        ~Bfinder();
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 
    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
 
        virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        virtual void endRun(edm::Run const&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        
        virtual void BranchOut2MuTk(
            BInfoBranches &BInfo,
            std::vector<pat::GenericParticle> input_tracks,
            reco::Vertex thePrimaryV,
            std::vector<bool> isNeededTrack,
            TLorentzVector v4_mu1,
            TLorentzVector v4_mu2,
            reco::TransientTrack muonPTT,
            reco::TransientTrack muonMTT,
            std::vector<int> &B_counter,
            float *mass_window,
            float MuMu_MASS,
            float Tk_MASS,
            int channel_number
        );
        virtual void BranchOut2MuX_XtoTkTk(
            BInfoBranches &BInfo,
            std::vector<pat::GenericParticle> input_tracks,
            reco::Vertex thePrimaryV,
            std::vector<bool> isNeededTrack,
            TLorentzVector v4_mu1,
            TLorentzVector v4_mu2,
            reco::TransientTrack muonPTT,
            reco::TransientTrack muonMTT,
            std::vector<int> &B_counter,
            float *mass_window,
            float MuMu_MASS,
            float TkTk_MASS,
            float TkTk_window,
            float Tk1_MASS,
            float Tk2_MASS,     
            int channel_number,
            int fit_option
        );
 
        // ----------member data ---------------------------
        edm::ESHandle<MagneticField> bField;
        edm::ParameterSet theConfig;
        //std::vector<std::string> TriggersForMatching_;
        std::vector<int> Bchannel_;
        std::vector<std::string> MuonTriggerMatchingPath_;
        //edm::InputTag hltLabel_;
        edm::InputTag genLabel_;
        edm::InputTag muonLabel_;
        edm::InputTag trackLabel_;
        edm::InputTag puInfoLabel_;
        edm::InputTag bsLabel_;
        edm::InputTag pvLabel_;
        double tkPtCut_;
        double tkEtaCut_;
        double jpsiPtCut_;
        double bPtCut_;
        double bEtaCut_;
        double VtxChiProbCut_;
        double svpvDistanceCut_;
        double MaxDocaCut_;
        double alphaCut_;
        bool RunOnMC_;
        bool doTkPreCut_;
        bool doMuPreCut_;
        bool doBntupleSkim_;
        std::string MVAMapLabel_;

        edm::Service<TFileService> fs;
        TTree *root;
        EvtInfoBranches     EvtInfo;
        VtxInfoBranches     VtxInfo;
        MuonInfoBranches    MuonInfo;
        TrackInfoBranches   TrackInfo;
        BInfoBranches       BInfo;
        GenInfoBranches     GenInfo;
        CommonFuncts        Functs;
        BntupleBranches     *Bntuple = new BntupleBranches;
        TTree* nt0;
        TTree* nt1;
        TTree* nt2;
        TTree* nt3;
        TTree* nt5;
        TTree* nt6;
        TTree* ntGen;

        //histograms
        TH1F *MuonCutLevel;
        TH1F *TrackCutLevel;
        TH1F *XbujCutLevel;
        //How many channel
        static int const Nchannel = 20;
        std::vector<TH1F*> XbMassCutLevel;

};//}}}

void Bfinder::beginJob()
{//{{{
    root = fs->make<TTree>("root","root");
    nt0   = fs->make<TTree>("ntKp","");     Bntuple->buildBranch(nt0);
    nt1   = fs->make<TTree>("ntpi","");     Bntuple->buildBranch(nt1);
    nt2   = fs->make<TTree>("ntKs","");     Bntuple->buildBranch(nt2);
    nt3   = fs->make<TTree>("ntKstar","");  Bntuple->buildBranch(nt3);
    nt5   = fs->make<TTree>("ntphi","");    Bntuple->buildBranch(nt5);
    nt6   = fs->make<TTree>("ntmix","");    Bntuple->buildBranch(nt6);
    ntGen = fs->make<TTree>("ntGen","");    Bntuple->buildGenBranch(ntGen);
    EvtInfo.regTree(root);
    VtxInfo.regTree(root);
    MuonInfo.regTree(root);
    TrackInfo.regTree(root);
    BInfo.regTree(root);
    GenInfo.regTree(root);
}//}}}

Bfinder::Bfinder(const edm::ParameterSet& iConfig):theConfig(iConfig)
{//{{{
    //now do what ever initialization is needed
//  TriggersForMatching_= iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching");
    Bchannel_= iConfig.getParameter<std::vector<int> >("Bchannel");
    MuonTriggerMatchingPath_ = iConfig.getParameter<std::vector<std::string> >("MuonTriggerMatchingPath");
    genLabel_           = iConfig.getParameter<edm::InputTag>("GenLabel");
    trackLabel_         = iConfig.getParameter<edm::InputTag>("TrackLabel");
    muonLabel_          = iConfig.getParameter<edm::InputTag>("MuonLabel");
    //hltLabel_           = iConfig.getParameter<edm::InputTag>("HLTLabel");
    puInfoLabel_        = iConfig.getParameter<edm::InputTag>("PUInfoLabel");
    bsLabel_        = iConfig.getParameter<edm::InputTag>("BSLabel");
    pvLabel_        = iConfig.getParameter<edm::InputTag>("PVLabel");

    tkPtCut_ = iConfig.getParameter<double>("tkPtCut");
    tkEtaCut_ = iConfig.getParameter<double>("tkEtaCut");
    jpsiPtCut_ = iConfig.getParameter<double>("jpsiPtCut");
    bPtCut_ = iConfig.getParameter<double>("bPtCut");
    bEtaCut_ = iConfig.getParameter<double>("bEtaCut");
    VtxChiProbCut_ = iConfig.getParameter<double>("VtxChiProbCut");
    svpvDistanceCut_ = iConfig.getParameter<double>("svpvDistanceCut");
    MaxDocaCut_ = iConfig.getParameter<double>("MaxDocaCut");
    alphaCut_ = iConfig.getParameter<double>("alphaCut");
    RunOnMC_ = iConfig.getParameter<bool>("RunOnMC");
    doTkPreCut_ = iConfig.getParameter<bool>("doTkPreCut");
    doMuPreCut_ = iConfig.getParameter<bool>("doMuPreCut");
    doBntupleSkim_ = iConfig.getParameter<bool>("doBntupleSkim");
    MVAMapLabel_  = iConfig.getParameter<std::string>("MVAMapLabel");

    MuonCutLevel        = fs->make<TH1F>("MuonCutLevel"     , "MuonCutLevel"    , 10, 0, 10);
    TrackCutLevel       = fs->make<TH1F>("TrackCutLevel"    , "TrackCutLevel"   , 10, 0, 10);
    XbujCutLevel        = fs->make<TH1F>("XbujCutLevel"     , "XbujCutLevel"    , 10, 0, 10);
    for(unsigned int i = 0; i < Bchannel_.size(); i++){
        TH1F* XbMassCutLevel_temp      = fs->make<TH1F>(TString::Format("XbMassCutLevel_i")   ,TString::Format("XbMassCutLevel_i")  , 10, 0, 10);
        XbMassCutLevel.push_back(XbMassCutLevel_temp);
    }
}//}}}

Bfinder::~Bfinder()
{//{{{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}//}}}

//
// member functions
//

// ------------ method called for each event  ------------
void Bfinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //std::cout << "*************************\nReconstructing event number: " << iEvent.id() << "\n";
    using namespace edm;
    using namespace reco;
    //ESHandle<MagneticField> bField;
    iSetup.get<IdealMagneticFieldRecord>().get(bField);

    // Change used muon and track collections
    edm::Handle< std::vector<pat::Muon> > muons;
    iEvent.getByLabel(muonLabel_,muons);
    edm::Handle< std::vector<pat::GenericParticle> > tks;
    iEvent.getByLabel(trackLabel_, tks);

    //CLEAN all memory
    memset(&EvtInfo     ,0x00,sizeof(EvtInfo)   );
    memset(&VtxInfo     ,0x00,sizeof(VtxInfo)   );
    memset(&MuonInfo    ,0x00,sizeof(MuonInfo)  );
    memset(&TrackInfo   ,0x00,sizeof(TrackInfo) );
    memset(&BInfo       ,0x00,sizeof(BInfo)    );
    memset(&GenInfo     ,0x00,sizeof(GenInfo)   );

    // EvtInfo section{{{
    EvtInfo.RunNo   = iEvent.id().run();
    EvtInfo.EvtNo   = iEvent.id().event();
    //std::cout<<"(EvtInfo.EvtNo)"<<EvtInfo.EvtNo<<std::endl;
    EvtInfo.BxNo    = iEvent.bunchCrossing();
    EvtInfo.LumiNo  = iEvent.luminosityBlock();
    EvtInfo.Orbit   = iEvent.orbitNumber();
    EvtInfo.McFlag  = !iEvent.isRealData();
    //EvtInfo.hltnames->clear();
    //EvtInfo.nTrgBook= N_TRIGGER_BOOKINGS;

    //Using HI HLT analysis now
    /*
    //HLT{{{
    edm::Handle<TriggerResults> TrgResultsHandle; //catch triggerresults
    bool with_TriggerResults = iEvent.getByLabel(hltLabel_,TrgResultsHandle);
    if(!with_TriggerResults){//
        std::cout << "Sorry there is no TriggerResult in the file" << std::endl;
    }else{
        //get the names of the triggers
        const edm::TriggerNames &TrgNames = iEvent.triggerNames(*TrgResultsHandle);
        EvtInfo.trgCount = 0;
        for(int i=0; i< N_TRIGGER_BOOKINGS; i++){
            unsigned int TrgIndex = TrgNames.triggerIndex(TriggerBooking[i]);
            if (TrgIndex == TrgNames.size()) {
                EvtInfo.trgBook[i] = -4; // The trigger path is not known in this event.
            }else if ( !TrgResultsHandle->wasrun( TrgIndex ) ) {
                EvtInfo.trgBook[i] = -3; // The trigger path was not included in this event.
            }else if ( !TrgResultsHandle->accept( TrgIndex ) ) {
                EvtInfo.trgBook[i] = -2; // The trigger path was not accepted in this event.
            }else if (  TrgResultsHandle->error ( TrgIndex ) ) {
                EvtInfo.trgBook[i] = -1; // The trigger path has an error in this event.
            }else {
                EvtInfo.trgBook[i] = +1; // It's triggered.
                EvtInfo.trgCount++; 
            }
        }
        EvtInfo.nHLT = TrgNames.size();
        for(unsigned int i=0; i<TrgNames.size(); i++){
            EvtInfo.hltBits[i] = (TrgResultsHandle->accept(i) == true) ? 1:0;
        }
    }//end(!with_TriggerResults)}}}
    */

    // Handle primary vertex properties
    Vertex thePrimaryV;
    math::XYZPoint RefVtx;
    //get beamspot information
    Vertex theBeamSpotV;
    reco::BeamSpot beamSpot;
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    //iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
    iEvent.getByLabel(bsLabel_, beamSpotHandle);
    if (beamSpotHandle.isValid()){
        beamSpot = *beamSpotHandle;
        theBeamSpotV = Vertex(beamSpot.position(), beamSpot.covariance3D());
        EvtInfo.BSx             = beamSpot.x0();
        EvtInfo.BSy             = beamSpot.y0();
        EvtInfo.BSz             = beamSpot.z0();
        EvtInfo.BSxErr          = beamSpot.x0Error();
        EvtInfo.BSyErr          = beamSpot.y0Error();
        EvtInfo.BSzErr          = beamSpot.z0Error();
        EvtInfo.BSdxdz          = beamSpot.dxdz();
        EvtInfo.BSdydz          = beamSpot.dydz();
        EvtInfo.BSdxdzErr       = beamSpot.dxdzError();
        EvtInfo.BSdydzErr       = beamSpot.dydzError();
        EvtInfo.BSWidthX        = beamSpot.BeamWidthX();
        EvtInfo.BSWidthXErr     = beamSpot.BeamWidthXError();
        EvtInfo.BSWidthY        = beamSpot.BeamWidthY();
        EvtInfo.BSWidthYErr     = beamSpot.BeamWidthYError();
    }else{
        std::cout<< "No beam spot available from EventSetup \n";
    }

    //get vertex informationa
    edm::Handle<reco::VertexCollection> VertexHandle;
    iEvent.getByLabel(pvLabel_, VertexHandle);

    /*  
    if (!VertexHandle.failedToGet() && VertexHandle->size()>0){
        //int nVtxTrks = 0;//outdated PV definition
        double max_tkSt = 0;
        for(std::vector<reco::Vertex>::const_iterator it_vtx = VertexHandle->begin(); it_vtx != VertexHandle->end(); it_vtx++){
            if (!it_vtx->isValid()) continue;
            //find primary vertex with largest St
            double tkSt = 0;
            for(std::vector<reco::TrackBaseRef>::const_iterator it_tk = it_vtx->tracks_begin();
                it_tk != it_vtx->tracks_end(); it_tk++){
                tkSt += it_tk->get()->pt();
            }
            if (tkSt > max_tkSt){
                max_tkSt = tkSt;
                thePrimaryV = Vertex(*it_vtx);
            }
        }
    }else{ 
        thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }
    RefVtx = thePrimaryV.position();
    */

    double PVBS_Pt_Max = -100.;
    reco::Vertex PVtx_BS;
    if( VertexHandle.isValid() && !VertexHandle.failedToGet() && VertexHandle->size() > 0) {
        //const vector<reco::Vertex> VerticesBS = *VertexHandle;
        for(std::vector<reco::Vertex>::const_iterator it_vtx = VertexHandle->begin();it_vtx != VertexHandle->end(); it_vtx++ ) {
        if (VtxInfo.Size>=MAX_Vertices) {
            std::cout << "PVBS " << VtxInfo.Size << std::endl;
            fprintf(stderr,"ERROR: number of  Vertices exceeds the size of array.\n");
            break;//exit(0);
        }
        VtxInfo.isValid[VtxInfo.Size] = it_vtx->isValid();
        VtxInfo.isFake[VtxInfo.Size] = it_vtx->isFake();
        VtxInfo.Ndof[VtxInfo.Size] = it_vtx->ndof();
        VtxInfo.NormalizedChi2[VtxInfo.Size] = it_vtx->normalizedChi2();
        VtxInfo.x[VtxInfo.Size] = it_vtx->x(); 
        VtxInfo.y[VtxInfo.Size] = it_vtx->y();
        VtxInfo.z[VtxInfo.Size] = it_vtx->z();
        VtxInfo.Pt_Sum[VtxInfo.Size] = 0.;
        VtxInfo.Pt_Sum2[VtxInfo.Size] = 0.;
        //if its hiSelectedVertex, then there will be only one vertex and will have no associated tracks
        if(int(VertexHandle->end()-VertexHandle->begin())==1){
            thePrimaryV = *it_vtx;
            VtxInfo.Size++;
            break;
        }
        for (reco::Vertex::trackRef_iterator it = it_vtx->tracks_begin(); it != it_vtx->tracks_end(); it++) {
           VtxInfo.Pt_Sum[VtxInfo.Size] += (*it)->pt();
           VtxInfo.Pt_Sum2[VtxInfo.Size] += ((*it)->pt() * (*it)->pt());
        }
        if( VtxInfo.Pt_Sum[VtxInfo.Size] >= PVBS_Pt_Max ){
            PVBS_Pt_Max = VtxInfo.Pt_Sum[VtxInfo.Size];
            thePrimaryV = *it_vtx;
        }            
        VtxInfo.Size++;
        }
    }else{ 
        thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }
    RefVtx = thePrimaryV.position();

    EvtInfo.PVx     = thePrimaryV.position().x();
    EvtInfo.PVy     = thePrimaryV.position().y();
    EvtInfo.PVz     = thePrimaryV.position().z();
    EvtInfo.PVxE    = thePrimaryV.xError();
    EvtInfo.PVyE    = thePrimaryV.yError();
    EvtInfo.PVzE    = thePrimaryV.zError();
    EvtInfo.PVnchi2 = thePrimaryV.normalizedChi2();
    EvtInfo.PVchi2  = thePrimaryV.chi2();

    // get pile-up information
    if (!iEvent.isRealData() && RunOnMC_){
        edm::Handle<std::vector< PileupSummaryInfo > >  PUHandle;
        iEvent.getByLabel(puInfoLabel_, PUHandle);
        std::vector<PileupSummaryInfo>::const_iterator PVI;
        for(PVI = PUHandle->begin(); PVI != PUHandle->end(); ++PVI) {
            EvtInfo.nPU[EvtInfo.nBX]   = PVI->getPU_NumInteractions();
            EvtInfo.BXPU[EvtInfo.nBX]  = PVI->getBunchCrossing();
            EvtInfo.trueIT[EvtInfo.nBX]= PVI->getTrueNumInteractions();
            EvtInfo.nBX += 1;
        }
    }else{
        EvtInfo.nBX = 0;
    }

    //}}}
    //printf("-----*****DEBUG:End of EvtInfo.\n");

    // Double check size=0.
    MuonInfo.size   = 0;
    TrackInfo.size  = 0;
    BInfo.uj_size   = 0;
    BInfo.size      = 0;
    GenInfo.size    = 0;
    
    std::vector<int> B_counter;
    for(unsigned int i = 0; i < Bchannel_.size(); i++){
        B_counter.push_back(0);
    }

    //Branches of type std::vector pointer must reserve memory before using
    MuonInfo.MuTrgMatchTrgObjE = new std::vector<std::vector<double>>();
    MuonInfo.MuTrgMatchTrgObjPt = new std::vector<std::vector<double>>();
    MuonInfo.MuTrgMatchTrgObjEta = new std::vector<std::vector<double>>();
    MuonInfo.MuTrgMatchTrgObjPhi = new std::vector<std::vector<double>>();

    std::vector<pat::Muon>              input_muons;
    std::vector<pat::GenericParticle>   input_tracks;
    input_muons = *muons;
    input_tracks = *tks;
    try{
        const reco::GenParticle* genMuonPtr[MAX_MUON];
        memset(genMuonPtr,0x00,MAX_MUON);
        const reco::GenParticle* genTrackPtr[MAX_TRACK];
        memset(genTrackPtr,0x00,MAX_GEN);
        //standard check for validity of input data
        if (input_muons.size() == 0){
            std::cout << "There's no muon : " << iEvent.id() << std::endl;
        }else{
            std::cout << "Got " << input_muons.size() << " muons / ";
            if (input_tracks.size() == 0){
                std::cout << "There's no track: " << iEvent.id() << std::endl;
            }else{
                std::cout << "Got " << input_tracks.size() << " tracks" << std::endl;
                if (input_tracks.size() > 0 && input_muons.size() > 1){

                    //MuonInfo section{{{
                    int mu_hindex = -1;
                    for(std::vector<pat::Muon>::const_iterator mu_it=input_muons.begin();
                        mu_it != input_muons.end() ; mu_it++){
                        mu_hindex = int(mu_it - input_muons.begin());
                        if(MuonInfo.size >= MAX_MUON){
                            fprintf(stderr,"ERROR: number of muons exceeds the size of array.\n");
                            break;//exit(0);
                        }
 
                        //Muon cut level
                        MuonCutLevel->Fill(0);
                        if (!(mu_it->isTrackerMuon() || mu_it->isGlobalMuon())) ;
                        else {
                            MuonCutLevel->Fill(1);
                            if (!muon::isGoodMuon(*mu_it,muon::TMOneStationTight)) ;
                            else {
                                MuonCutLevel->Fill(2);
                                if(!mu_it->innerTrack().isNonnull()) ;
                                else {
                                    MuonCutLevel->Fill(3);
                                    if (  fabs(mu_it->innerTrack()->dxy(thePrimaryV.position())) >= 3.        || 
                                          fabs(mu_it->innerTrack()->dz(thePrimaryV.position()))  >= 30.       
                                       ) ;
                                    else {
                                        MuonCutLevel->Fill(4);
                                        if (mu_it->innerTrack()->hitPattern().pixelLayersWithMeasurement()<1    ||
                                            mu_it->innerTrack()->normalizedChi2()>1.8
                                           ) ;
                                        else {
                                            MuonCutLevel->Fill(5);
                                            if (mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement()<6) ;
                                            else 
                                                MuonCutLevel->Fill(6);
                                             }   
                                         }
                                     }
                                 }
                            }

                        //Muon Id flag
                        MuonInfo.BfinderMuID[MuonInfo.size] = false;
                        if(mu_it->innerTrack().isNonnull()){
                            if( (mu_it->isTrackerMuon() || mu_it->isGlobalMuon()) 
                            && (muon::isGoodMuon(*mu_it,muon::TMOneStationTight)) 
                            && fabs(mu_it->innerTrack()->dxy(thePrimaryV.position())) < 3.
                            && fabs(mu_it->innerTrack()->dz(thePrimaryV.position()))  < 30.
                            && mu_it->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 
                            && mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
                            && mu_it->innerTrack()->normalizedChi2() <= 1.8
                            )
                                MuonInfo.BfinderMuID[MuonInfo.size] = true;
                        }
                        
                        MuonInfo.SoftMuID[MuonInfo.size] = false;
                        if(mu_it->innerTrack().isNonnull()){
                            if( (muon::isGoodMuon(*mu_it,muon::TMOneStationTight)) 
                            && mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
                            && mu_it->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 
                            && mu_it->innerTrack()->quality(reco::TrackBase::highPurity)
                            && fabs(mu_it->innerTrack()->dxy(thePrimaryV.position())) < 0.3
                            && fabs(mu_it->innerTrack()->dz(thePrimaryV.position()))  < 20.
                            )
                                MuonInfo.SoftMuID[MuonInfo.size] = true;
                        }

                        //outdated selections
                        //if(!(mu_it->innerTrack().isNonnull()*mu_it->globalTrack().isNonnull())) {continue;}
                        //if (!(mu_it->isGlobalMuon()*mu_it->track().isAvailable()*mu_it->globalTrack().isAvailable())) continue;
                        //if (mu_it->p()>200 || mu_it->pt()>200)                  continue;
                        //if (mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement()<6  &&
                        //    mu_it->innerTrack()->hitPattern().numberOfValidStripHits()<11
                        //   ) continue;

                        //Get Muon HLT Trigger matching
                        MuonInfo.MuTrgMatchPathSize = MuonTriggerMatchingPath_.size();
                        std::vector<double> trgobjE;
                        std::vector<double> trgobjPt;
                        std::vector<double> trgobjEta;
                        std::vector<double> trgobjPhi;
                        MuonInfo.isTriggered[MuonInfo.size] = false;
                        for(int _m = 0; _m < MuonInfo.MuTrgMatchPathSize; _m++){
                            pat::TriggerObjectStandAloneCollection match = mu_it->triggerObjectMatchesByPath(MuonTriggerMatchingPath_[_m].c_str());
                            //pat::TriggerObjectStandAloneCollection match = mu_it->triggerObjectMatchesByPath(MuonTriggerMatchingPath_[_m].c_str(), true, false);
                            if (match.empty()) {
                                trgobjE.push_back(-999.);
                                trgobjPt.push_back(-999.);
                                trgobjEta.push_back(-999.);
                                trgobjPhi.push_back(-999.);
                                //std::cout << "Muon didn't match Trigger Object" << std::endl;
                            } else {
                                trgobjE.push_back(match[0].energy());
                                trgobjPt.push_back(match[0].pt());
                                trgobjEta.push_back(match[0].eta());
                                trgobjPhi.push_back(match[0].phi());
                                //std::cout << "Propagation succeeeded; eta = " << match[0].eta() << ", phi = " << match[0].phi() << std::endl;
                                MuonInfo.isTriggered[MuonInfo.size] = true;
                            }
                        }
                        MuonInfo.MuTrgMatchTrgObjE->push_back(trgobjE);
                        MuonInfo.MuTrgMatchTrgObjPt->push_back(trgobjPt);
                        MuonInfo.MuTrgMatchTrgObjEta->push_back(trgobjEta);
                        MuonInfo.MuTrgMatchTrgObjPhi->push_back(trgobjPhi);
                        
                        //Muon general info.
                        MuonInfo.index          [MuonInfo.size] = MuonInfo.size;
                        MuonInfo.handle_index   [MuonInfo.size] = mu_hindex;
                        MuonInfo.charge         [MuonInfo.size] = mu_it->charge();
                        MuonInfo.pt             [MuonInfo.size] = mu_it->pt();
                        MuonInfo.eta            [MuonInfo.size] = mu_it->eta();
                        MuonInfo.phi            [MuonInfo.size] = mu_it->phi();
                        MuonInfo.isTrackerMuon  [MuonInfo.size] = mu_it->isTrackerMuon();
                        MuonInfo.isGlobalMuon   [MuonInfo.size] = mu_it->isGlobalMuon();
                        MuonInfo.iso_trk        [MuonInfo.size] = mu_it->trackIso();//R<0.3
                        MuonInfo.iso_ecal       [MuonInfo.size] = mu_it->ecalIso();
                        MuonInfo.iso_hcal       [MuonInfo.size] = mu_it->hcalIso();
                        MuonInfo.type           [MuonInfo.size] = mu_it->type();//CaloMuon = 1<<4  GlobalMuon = 1<<1  PFMuon = 1<<5  StandAloneMuon = 1<<3  TrackerMuon = 1<<2
                        MuonInfo.n_matches      [MuonInfo.size] = mu_it->numberOfMatches();//only in chamber
                        MuonInfo.geninfo_index  [MuonInfo.size] = -1;//initialize for later use
                        MuonInfo.TMOneStationTight[MuonInfo.size] = muon::isGoodMuon(*mu_it,muon::TMOneStationTight);//For Muon ID for convenience
                        MuonInfo.TrackerMuonArbitrated[MuonInfo.size] = muon::isGoodMuon(*mu_it,muon::TrackerMuonArbitrated);//For Muon ID for convenience
                        if (!iEvent.isRealData() && RunOnMC_) genMuonPtr [MuonInfo.size] = mu_it->genParticle();

                        //Muon standalone info.
                        MuonInfo.isStandAloneMuon[MuonInfo.size] = false;
                        if(mu_it->isStandAloneMuon()){
                            MuonInfo.isStandAloneMuon[MuonInfo.size] = true;
                            reco::TrackRef tkref;
                            tkref = mu_it->standAloneMuon();
                            const reco::Track &trk = *tkref;
                            MuonInfo.StandAloneMuon_charge         [MuonInfo.size] = trk.charge();
                            MuonInfo.StandAloneMuon_pt             [MuonInfo.size] = trk.pt();
                            MuonInfo.StandAloneMuon_eta            [MuonInfo.size] = trk.eta();
                            MuonInfo.StandAloneMuon_phi            [MuonInfo.size] = trk.phi();
                            MuonInfo.StandAloneMuon_d0             [MuonInfo.size] = trk.d0();
                            MuonInfo.StandAloneMuon_dz             [MuonInfo.size] = trk.dz();
                            MuonInfo.StandAloneMuon_dzPV           [MuonInfo.size] = trk.dz(RefVtx);
                            MuonInfo.StandAloneMuon_dxyPV          [MuonInfo.size] = trk.dxy(RefVtx);
                            //std::cout<<"sta pt: "<<trk.pt()<<std::endl; std::cout<<"pt: "<<mu_it->pt()<<std::endl;
                        }

                        MuonInfo.outerTrackisNonnull[MuonInfo.size] = mu_it->outerTrack().isNonnull();
                        MuonInfo.innerTrackisNonnull[MuonInfo.size] = mu_it->innerTrack().isNonnull();
                        //Muon inner track info.
                        
                        if(mu_it->innerTrack().isNonnull()){
                            //Muon inner track track quality
                            //enum TrackQuality { undefQuality = -1, loose = 0, tight = 1, highPurity = 2, confirmed = 3, goodIterative = 4, looseSetWithPV = 5, highPuritySetWithPV = 6, qualitySize = 7}
                            for(int tq = 0; tq < reco::TrackBase::qualitySize; tq++){
                                if (mu_it->innerTrack()->quality(static_cast<reco::TrackBase::TrackQuality>(tq))) MuonInfo.innerTrackQuality[MuonInfo.size] += 1 << (tq);
                                //std::cout<<"type: "<<mu_it->innerTrack()->quality(static_cast<reco::TrackBase::TrackQuality>(tq))<<std::endl;
                            }
                            MuonInfo.normchi2                [MuonInfo.size] = mu_it->innerTrack()->normalizedChi2();
                            MuonInfo.i_striphit              [MuonInfo.size] = mu_it->innerTrack()->hitPattern().numberOfValidStripHits();
                            MuonInfo.i_pixelhit              [MuonInfo.size] = mu_it->innerTrack()->hitPattern().numberOfValidPixelHits();
                            MuonInfo.i_nStripLayer           [MuonInfo.size] = mu_it->innerTrack()->hitPattern().stripLayersWithMeasurement();
                            MuonInfo.i_nPixelLayer           [MuonInfo.size] = mu_it->innerTrack()->hitPattern().pixelLayersWithMeasurement();
                            MuonInfo.i_chi2                  [MuonInfo.size] = mu_it->innerTrack()->chi2();
                            MuonInfo.i_ndf                   [MuonInfo.size] = mu_it->innerTrack()->ndof();
                            MuonInfo.fpbarrelhit             [MuonInfo.size] = mu_it->innerTrack()->hitPattern().hasValidHitInFirstPixelBarrel();
                            MuonInfo.fpendcaphit             [MuonInfo.size] = mu_it->innerTrack()->hitPattern().hasValidHitInFirstPixelEndcap();
                            MuonInfo.ptErr                   [MuonInfo.size] = mu_it->track()->ptError();
                            MuonInfo.etaErr                  [MuonInfo.size] = mu_it->track()->etaError();
                            MuonInfo.phiErr                  [MuonInfo.size] = mu_it->track()->phiError();
                            MuonInfo.d0                      [MuonInfo.size] = mu_it->track()->d0();
                            MuonInfo.dz                      [MuonInfo.size] = mu_it->track()->dz();
                            MuonInfo.dzPV                    [MuonInfo.size] = mu_it->track()->dz(RefVtx);//==mu_it->innerTrack()->dxy(thePrimaryV.position());
                            MuonInfo.dxyPV                   [MuonInfo.size] = mu_it->track()->dxy(RefVtx);//==mu_it->innerTrack()->dz(thePrimaryV.position());
                            //mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement() == MuonInfo.i_nStripLayer + MuonInfo.i_nPixelLayer
                        }
                        //Muon global track info.
                        MuonInfo.globalTrackisNonnull[MuonInfo.size] = mu_it->globalTrack().isNonnull();
                        if(mu_it->isGlobalMuon()){
                            MuonInfo.g_striphit [MuonInfo.size] = mu_it->globalTrack()->hitPattern().numberOfValidStripHits();
                            MuonInfo.g_pixelhit [MuonInfo.size] = mu_it->globalTrack()->hitPattern().numberOfValidPixelHits();
                            MuonInfo.g_chi2     [MuonInfo.size] = mu_it->globalTrack()->chi2();
                            MuonInfo.g_ndf      [MuonInfo.size] = mu_it->globalTrack()->ndof();
                            MuonInfo.nmuhit     [MuonInfo.size] = mu_it->globalTrack()->hitPattern().numberOfValidMuonHits();
                        }else{
                            MuonInfo.g_striphit [MuonInfo.size] = -1;
                            MuonInfo.g_pixelhit [MuonInfo.size] = -1;
                            MuonInfo.g_chi2     [MuonInfo.size] = -1;
                            MuonInfo.g_ndf      [MuonInfo.size] = -1;
                            MuonInfo.nmuhit     [MuonInfo.size] = -1;
                        }
                        //Muon quality
                        //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMuonAnalysis
                        int qm = 0;
                        for(int qi=1; qi!= 24; ++qi){
                            if (muon::isGoodMuon(*mu_it, muon::SelectionType(qi))){
                                qm += 1 << qi;
                            }
                        }
                        MuonInfo.muqual         [MuonInfo.size] = qm;   

                        // Basic Muon selection for fitting
                        //Can not be just CaloMuon or empty type
                        if((MuonInfo.type[MuonInfo.size]|(1<<4))==(1<<4)){
                            MuonInfo.isNeededMuon[MuonInfo.size] = false;
                        }
                        else if(doMuPreCut_ &&
                                (   !muon::isGoodMuon(*mu_it,muon::TMOneStationTight)
                                //||  some other cut
                                )
                        ){
                            MuonInfo.isNeededMuon[MuonInfo.size] = false;
                        }
                        else MuonInfo.isNeededMuon[MuonInfo.size] = true;

                        MuonInfo.size++;
                    }//end of MuonInfo}}}
                    //printf("-----*****DEBUG:End of MuonInfo.\n");

                    //Preselect tracks{{{
                    std::vector<bool> isNeededTrack;// Are the tracks redundant?
                    int PassedTrk = 0;
                    for(std::vector<pat::GenericParticle>::const_iterator tk_it=input_tracks.begin();
                        tk_it != input_tracks.end(); tk_it++){
                        if(PassedTrk >= MAX_TRACK){
                            fprintf(stderr,"ERROR: number of tracks exceeds the size of array.\n");
                            break;
                        }
                        isNeededTrack.push_back(false);
                        TrackCutLevel->Fill(0);//number of all tracks
                        bool isMuonTrack = false; //remove muon track
                        for(std::vector<pat::Muon>::iterator it=input_muons.begin() ; 
                            it != input_muons.end(); it++){
                            if (!it->track().isNonnull())                   continue;
                            if((it->type()|(1<<4))==(1<<4)) continue;//Don't clean track w.r.t. calo muon 
                            if (fabs(tk_it->pt() -it->track()->pt() )<0.00001 &&
                                fabs(tk_it->eta()-it->track()->eta())<0.00001 &&
                                fabs(tk_it->phi()-it->track()->phi())<0.00001 ){
                                    isMuonTrack = true;
                                    break;
                            }
                        }
                        if (isMuonTrack)                                    continue;
                        TrackCutLevel->Fill(1);//number of non muon tracks
                        if (tk_it->pt()<tkPtCut_)                           continue;
                        TrackCutLevel->Fill(2);
                        if (fabs(tk_it->eta())>tkEtaCut_)                   continue;
                        TrackCutLevel->Fill(3);
                        //if (fabs(tk_it->eta()) > 2.5)                       continue;
                        TrackCutLevel->Fill(4);
                        if(doTkPreCut_){
                            if( !(tk_it->track()->quality(reco::TrackBase::highPurity))) continue;
                            //outdated selections
                            //if (tk_it->track()->normalizedChi2()>5)             continue;
                            //if (tk_it->p()>200 || tk_it->pt()>200)              continue;
                            //if (tk_it->track()->hitPattern().numberOfValidStripHits()<10)continue;
                            //if (tk_it->track()->hitPattern().numberOfValidPixelHits()<2) continue;
                            TrackCutLevel->Fill(5);
                        }
                        isNeededTrack[tk_it-input_tracks.begin()] = true;
                        PassedTrk++;
                    }//end of track preselection}}}
                    //printf("-----*****DEBUG:End of track preselection.\n");

                    // BInfo section{{{
                    int mu1_index = -1;
                    int mu1_hindex = -1;
                    bool gogogo = false;
                    for(std::vector<pat::Muon>::const_iterator mu_it1=input_muons.begin();
                        mu_it1 != input_muons.end(); mu_it1++){
                        //check if muon track is non null
                        if(!mu_it1->track().isNonnull()) continue;
                        //Check if it in MuonInfo and isNeedeMuon
                        mu1_hindex = int(mu_it1 - input_muons.begin());
                        gogogo = false;
                        for(int i=0; i < MuonInfo.size; i++){
                            if (mu1_hindex == MuonInfo.handle_index[i] && MuonInfo.isNeededMuon[i]){
                                gogogo = true;
                                break;
                            }
                        }
                        if (!gogogo) continue;
                        //Get the corrisponding index in MuonInfo
                        mu1_index ++;
                        if (mu_it1->charge()<0) continue;
                        
                        int mu2_index = -1;
                        int mu2_hindex = -1; 
                        for(std::vector<pat::Muon>::const_iterator mu_it2=input_muons.begin();
                            mu_it2 != input_muons.end(); mu_it2++){
                            //check if muon track is non null
                            if(!mu_it2->track().isNonnull()) continue;
                            //Check if it in MuonInfo and isNeedeMuon
                            mu2_hindex = int(mu_it2 - input_muons.begin()); 
                            gogogo = false;
                            for(int j=0; j < MuonInfo.size; j++){
                                if(mu2_hindex == MuonInfo.handle_index[j] && MuonInfo.isNeededMuon[j]){
                                    gogogo = true;
                                    break;
                                }
                            }
                            if (!gogogo) continue;
                            mu2_index ++;   
                            if (mu_it2->charge()>0) continue;
                            XbujCutLevel->Fill(0);
                            
                            TLorentzVector v4_mu1,v4_mu2;
                            v4_mu1.SetPtEtaPhiM(mu_it1->pt(),mu_it1->eta(),mu_it1->phi(),MUON_MASS);
                            v4_mu2.SetPtEtaPhiM(mu_it2->pt(),mu_it2->eta(),mu_it2->phi(),MUON_MASS);
                            if (fabs((v4_mu1+v4_mu2).Mag()-JPSI_MASS)>0.4) continue;
                            if((v4_mu1+v4_mu2).Pt()<jpsiPtCut_)continue;
    
                            //Fit 2 muon
                            reco::TransientTrack muonPTT(mu_it1->track(), &(*bField) );
                            reco::TransientTrack muonMTT(mu_it2->track(), &(*bField) );
                            if(!muonPTT.isValid()) continue;
                            if(!muonMTT.isValid()) continue;
                            XbujCutLevel->Fill(1);
    
                            const reco::Muon* rmu1 = dynamic_cast<const reco::Muon * >(mu_it1->originalObject());
                            const reco::Muon* rmu2 = dynamic_cast<const reco::Muon * >(mu_it2->originalObject());
                            if(muon::overlap(*rmu1, *rmu2)) continue;
                            XbujCutLevel->Fill(2);
    
                            KinematicParticleFactoryFromTransientTrack pFactory;
                            ParticleMass muon_mass = MUON_MASS; //pdg mass
                            float muon_sigma = muon_mass*1.e-6;
                            float chi = 0.;
                            float ndf = 0.;
                            std::vector<RefCountedKinematicParticle> muonParticles;
                            muonParticles.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
                            muonParticles.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));
    
                            KinematicParticleVertexFitter   fitter;   
                            RefCountedKinematicTree         ujVFT;
                            ujVFT = fitter.fit(muonParticles); 
                            if (!ujVFT->isValid()) continue;
                            XbujCutLevel->Fill(3); 

                            ujVFT->movePointerToTheTop();
    
                            RefCountedKinematicParticle ujVFP       = ujVFT->currentParticle();
                            RefCountedKinematicVertex   ujVFPvtx    = ujVFT->currentDecayVertex();
                            ujVFT->movePointerToTheFirstChild();
                            KinematicParameters         ujmu1KP     = ujVFT->currentParticle()->currentState().kinematicParameters();
                            ujVFT->movePointerToTheNextChild();
                            KinematicParameters         ujmu2KP     = ujVFT->currentParticle()->currentState().kinematicParameters();
                            double chi2_prob_uj = TMath::Prob(ujVFPvtx->chiSquared(), ujVFPvtx->degreesOfFreedom());
                            if(chi2_prob_uj < VtxChiProbCut_) continue;
                            XbujCutLevel->Fill(4);

                            if (fabs(ujVFP->currentState().mass()-JPSI_MASS)>0.3) continue;

                            TLorentzVector uj_4vec,uj_mu1_4vec,uj_mu2_4vec;
                            uj_4vec.SetPxPyPzE(ujVFP->currentState().kinematicParameters().momentum().x(),
                                               ujVFP->currentState().kinematicParameters().momentum().y(),
                                               ujVFP->currentState().kinematicParameters().momentum().z(),
                                               ujVFP->currentState().kinematicParameters().energy());
                            uj_mu1_4vec.SetPxPyPzE( ujmu1KP.momentum().x(),
                                                    ujmu1KP.momentum().y(),
                                                    ujmu1KP.momentum().z(),
                                                    ujmu1KP.energy());
                            uj_mu2_4vec.SetPxPyPzE( ujmu2KP.momentum().x(),
                                                    ujmu2KP.momentum().y(),
                                                    ujmu2KP.momentum().z(),
                                                    ujmu2KP.energy());
                            //uj_4vec.Print();

                            BInfo.uj_index         [BInfo.uj_size]= BInfo.uj_size;
                            BInfo.uj_mass          [BInfo.uj_size]= uj_4vec.Mag();
                            BInfo.uj_pt            [BInfo.uj_size]= uj_4vec.Pt();
                            BInfo.uj_eta            [BInfo.uj_size]= uj_4vec.Eta();
                            BInfo.uj_phi            [BInfo.uj_size]= uj_4vec.Phi();
                            BInfo.uj_px            [BInfo.uj_size]= uj_4vec.Px();
                            BInfo.uj_py            [BInfo.uj_size]= uj_4vec.Py();
                            BInfo.uj_pz            [BInfo.uj_size]= uj_4vec.Pz();
                            BInfo.uj_vtxX          [BInfo.uj_size]= ujVFPvtx->position().x();
                            BInfo.uj_vtxY          [BInfo.uj_size]= ujVFPvtx->position().y();
                            BInfo.uj_vtxZ          [BInfo.uj_size]= ujVFPvtx->position().z();
                            BInfo.uj_vtxXErr       [BInfo.uj_size]= ujVFPvtx->error().cxx();
                            BInfo.uj_vtxYErr       [BInfo.uj_size]= ujVFPvtx->error().cyy();
                            BInfo.uj_vtxZErr       [BInfo.uj_size]= ujVFPvtx->error().czz();
                            BInfo.uj_vtxYXErr      [BInfo.uj_size]= ujVFPvtx->error().cyx();
                            BInfo.uj_vtxZXErr      [BInfo.uj_size]= ujVFPvtx->error().czx();
                            BInfo.uj_vtxZYErr      [BInfo.uj_size]= ujVFPvtx->error().czy();
                            BInfo.uj_vtxdof        [BInfo.uj_size]= ujVFPvtx->degreesOfFreedom();
                            BInfo.uj_vtxchi2       [BInfo.uj_size]= ujVFPvtx->chiSquared();
                            BInfo.uj_rfmu1_index   [BInfo.uj_size]= mu1_hindex;
                            BInfo.uj_rfmu2_index   [BInfo.uj_size]= mu2_hindex;
                            BInfo.uj_rfmu1_px      [BInfo.uj_size]= uj_mu1_4vec.Px();
                            BInfo.uj_rfmu1_py      [BInfo.uj_size]= uj_mu1_4vec.Py();
                            BInfo.uj_rfmu1_pz      [BInfo.uj_size]= uj_mu1_4vec.Pz();
                            BInfo.uj_rfmu2_px      [BInfo.uj_size]= uj_mu2_4vec.Px();
                            BInfo.uj_rfmu2_py      [BInfo.uj_size]= uj_mu2_4vec.Py();
                            BInfo.uj_rfmu2_pz      [BInfo.uj_size]= uj_mu2_4vec.Pz();

                            BInfo.uj_size++;
                            muonParticles.clear();

                            //////////////////////////////////////////////////////////////////////////
                            // RECONSTRUCTION: J/psi + K
                            //////////////////////////////////////////////////////////////////////////
                            //float mass_window[2] = {4.3, 6.4};
                            float mass_window[2] = {5., 6.};
                            if(Bchannel_[0] == 1){
                                BranchOut2MuTk(
                                    BInfo,
                                    input_tracks,
                                    thePrimaryV,
                                    isNeededTrack,
                                    v4_mu1,
                                    v4_mu2,
                                    muonPTT,
                                    muonMTT,
                                    B_counter,
                                    mass_window,
                                    JPSI_MASS,
                                    KAON_MASS,
                                    1
                                );
                            }                            
                            //////////////////////////////////////////////////////////////////////////
                            // RECONSTRUCTION: J/psi + Pi
                            //////////////////////////////////////////////////////////////////////////
                            if(Bchannel_[1] == 1){
                                BranchOut2MuTk(
                                    BInfo,
                                    input_tracks,
                                    thePrimaryV,
                                    isNeededTrack,
                                    v4_mu1,
                                    v4_mu2,
                                    muonPTT,
                                    muonMTT,
                                    B_counter,
                                    mass_window,
                                    JPSI_MASS,
                                    PION_MASS,
                                    2
                                );
                            }                            

                            //////////////////////////////////////////////////////////////////////////
                            // RECONSTRUCTION: J/psi + Ks
                            //////////////////////////////////////////////////////////////////////////

                            float TkTk_window = 0;
                            TkTk_window = 0.3;
                            if(Bchannel_[2] == 1){
                                BranchOut2MuX_XtoTkTk(
                                    BInfo,
                                    input_tracks,
                                    thePrimaryV,
                                    isNeededTrack,
                                    v4_mu1,
                                    v4_mu2,
                                    muonPTT,
                                    muonMTT,
                                    B_counter,
                                    mass_window,
                                    JPSI_MASS,
                                    KSHORT_MASS,
                                    TkTk_window,
                                    PION_MASS,        
                                    PION_MASS,
                                    3,
                                    1
                                );
                            }                            
                            
                            //////////////////////////////////////////////////////////////////////////
                            // RECONSTRUCTION: J/psi + K* (K+, Pi-)
                            //////////////////////////////////////////////////////////////////////////

                            //TkTk_window = 0.4;
                            TkTk_window = 0.25;
                            if(Bchannel_[3] == 1){
                                BranchOut2MuX_XtoTkTk(
                                    BInfo,
                                    input_tracks,
                                    thePrimaryV,
                                    isNeededTrack,
                                    v4_mu1,
                                    v4_mu2,
                                    muonPTT,
                                    muonMTT,
                                    B_counter,
                                    mass_window,
                                    JPSI_MASS,
                                    KSTAR_MASS,
                                    TkTk_window,
                                    KAON_MASS,        
                                    PION_MASS,
                                    4,
                                    0
                                );
                            }
                            
                            //////////////////////////////////////////////////////////////////////////
                            // RECONSTRUCTION: J/psi + K* (K-, Pi+)
                            //////////////////////////////////////////////////////////////////////////

                            //TkTk_window = 0.4;
                            TkTk_window = 0.25;
                            if(Bchannel_[4] == 1){
                                BranchOut2MuX_XtoTkTk(
                                    BInfo,
                                    input_tracks,
                                    thePrimaryV,
                                    isNeededTrack,
                                    v4_mu1,
                                    v4_mu2,
                                    muonPTT,
                                    muonMTT,
                                    B_counter,
                                    mass_window,
                                    JPSI_MASS,
                                    KSTAR_MASS,
                                    TkTk_window,
                                    PION_MASS,        
                                    KAON_MASS,
                                    5,
                                    0
                                );
                            }
                            
                            //////////////////////////////////////////////////////////////////////////
                            // RECONSTRUCTION: J/psi + phi
                            //////////////////////////////////////////////////////////////////////////
                            
                            TkTk_window = 0.1;
                            if(Bchannel_[5] == 1){
                                BranchOut2MuX_XtoTkTk(
                                    BInfo,
                                    input_tracks,
                                    thePrimaryV,
                                    isNeededTrack,
                                    v4_mu1,
                                    v4_mu2,
                                    muonPTT,
                                    muonMTT,
                                    B_counter,
                                    mass_window,
                                    JPSI_MASS,
                                    PHI_MASS,
                                    TkTk_window,
                                    KAON_MASS,        
                                    KAON_MASS,
                                    6,
                                    0
                                );
                            }

                            //////////////////////////////////////////////////////////////////////////
                            // RECONSTRUCTION: J/psi + pi pi <= psi', X(3872), Bs->J/psi f0
                            //////////////////////////////////////////////////////////////////////////
                            //mass_window[0] = 3;
                            //mass_window[1] = 6.4;
                            //TkTk_window = 1.6;
                            mass_window[0] = 3.4;
                            mass_window[1] = 4.2;
                            TkTk_window = 0;
                            if(Bchannel_[6] == 1){
                                BranchOut2MuX_XtoTkTk(
                                    BInfo,
                                    input_tracks,
                                    thePrimaryV,
                                    isNeededTrack,
                                    v4_mu1,
                                    v4_mu2,
                                    muonPTT,
                                    muonMTT,
                                    B_counter,
                                    mass_window,
                                    JPSI_MASS,
                                    -1,
                                    TkTk_window,
                                    PION_MASS,        
                                    PION_MASS,
                                    7,
                                    0
                                );
                            }
                            
                        }//Mu2
                    }//Mu1
                    printf("B_counter: ");
                    for(unsigned int i = 0; i < Bchannel_.size(); i++){
                        printf("%d/", B_counter[i]);
                    }
                    printf("\n");//}}}
                    //printf("-----*****DEBUG:End of BInfo.\n");

                    // TrackInfo section {{{
                    Handle<edm::ValueMap<float> > mvaoutput;
                    iEvent.getByLabel(MVAMapLabel_, "MVAVals", mvaoutput);
                    for(std::vector<pat::GenericParticle>::const_iterator tk_it=input_tracks.begin();
                        tk_it != input_tracks.end() ; tk_it++){
                        int tk_hindex = int(tk_it - input_tracks.begin());
                        if(tk_hindex>=int(isNeededTrack.size())) break;
                        if (isNeededTrack[tk_hindex]==false) continue;

                        //Create list of relative xb candidates for later filling
                        std::vector<int> listOfRelativeXbCands1;//1~nXb
                        std::vector<int> listOfRelativeXbCands2;//1~nXb
                        for(int iXb=0; iXb < BInfo.size; iXb++){
                            if(BInfo.rftk1_index[iXb] == tk_hindex){
                                listOfRelativeXbCands1.push_back(iXb+1);
                            }if(BInfo.rftk2_index[iXb] == tk_hindex){
                                listOfRelativeXbCands2.push_back(iXb+1);
                            }
                        }
                        //if (listOfRelativeXbCands1.size() == 0 && listOfRelativeXbCands2.size() == 0) continue;//drop unused tracks
                        
                        TrackInfo.index          [TrackInfo.size] = TrackInfo.size;
                        TrackInfo.handle_index   [TrackInfo.size] = tk_hindex;
                        TrackInfo.charge         [TrackInfo.size] = tk_it->charge();
                        TrackInfo.pt             [TrackInfo.size] = tk_it->pt();
                        TrackInfo.eta            [TrackInfo.size] = tk_it->eta();
                        TrackInfo.phi            [TrackInfo.size] = tk_it->phi();
                        TrackInfo.ptErr          [TrackInfo.size] = tk_it->track()->ptError();
                        TrackInfo.etaErr         [TrackInfo.size] = tk_it->track()->etaError();
                        TrackInfo.phiErr         [TrackInfo.size] = tk_it->track()->phiError();
                        //TrackInfo.p              [TrackInfo.size] = tk_it->p();
                        TrackInfo.striphit       [TrackInfo.size] = tk_it->track()->hitPattern().numberOfValidStripHits();
                        TrackInfo.pixelhit       [TrackInfo.size] = tk_it->track()->hitPattern().numberOfValidPixelHits();
                        TrackInfo.nStripLayer    [TrackInfo.size] = tk_it->track()->hitPattern().stripLayersWithMeasurement();
                        TrackInfo.nPixelLayer    [TrackInfo.size] = tk_it->track()->hitPattern().pixelLayersWithMeasurement();
                        TrackInfo.fpbarrelhit    [TrackInfo.size] = tk_it->track()->hitPattern().hasValidHitInFirstPixelBarrel();
                        TrackInfo.fpendcaphit    [TrackInfo.size] = tk_it->track()->hitPattern().hasValidHitInFirstPixelEndcap();
                        TrackInfo.chi2           [TrackInfo.size] = tk_it->track()->chi2();
                        TrackInfo.ndf            [TrackInfo.size] = tk_it->track()->ndof();
                        TrackInfo.d0             [TrackInfo.size] = tk_it->track()->d0();
                        TrackInfo.d0error        [TrackInfo.size] = tk_it->track()->d0Error();
                        TrackInfo.dzPV           [TrackInfo.size] = tk_it->track()->dz(RefVtx);
                        TrackInfo.dxyPV          [TrackInfo.size] = tk_it->track()->dxy(RefVtx);
                        TrackInfo.highPurity     [TrackInfo.size] = tk_it->track()->quality(reco::TrackBase::highPurity);
                        TrackInfo.geninfo_index  [TrackInfo.size] = -1;//initialize for later use
                        TrackInfo.trkMVAVal      [TrackInfo.size] = (*mvaoutput)[tk_it->track()];
                        TrackInfo.trkAlgo        [TrackInfo.size] = tk_it->track()->algo();


                        //https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_5_patch1/DataFormats/TrackReco/interface/TrackBase.h#L149
                        if(tk_it->track().isNonnull()){
                            for(int tq = 0; tq < reco::TrackBase::qualitySize; tq++){
                            if (tk_it->track()->quality(static_cast<reco::TrackBase::TrackQuality>(tq))) TrackInfo.trackQuality[TrackInfo.size] += 1 << (tq);
                        }}

                        if (!iEvent.isRealData() && RunOnMC_)
                            genTrackPtr [TrackInfo.size] = tk_it->genParticle();

                        //Fill correct track index and track quality to correspond Xb candidate
                        for(unsigned int iCands=0; iCands < listOfRelativeXbCands1.size(); iCands++){
                                BInfo.rftk1_index[listOfRelativeXbCands1[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeXbCands2.size(); iCands++){
                                BInfo.rftk2_index[listOfRelativeXbCands2[iCands]-1] = TrackInfo.size;
                        }

                        TrackInfo.size++;
                    }//end of TrackInfo}}}
                    //printf("-----*****DEBUG:End of TrackInfo.\n");
                }//has nTracks>1
            }//if no Tracks
        }//if no Muons

        // GenInfo section{{{
        if (!iEvent.isRealData() && RunOnMC_){
        //if (RunOnMC_){
        //if (1){
            //edm::Handle< std::vector<reco::GenParticle> > gens;
            edm::Handle<reco::GenParticleCollection> gens;
            iEvent.getByLabel(genLabel_, gens);

            std::vector<const reco::Candidate *> sel_cands;
            //deprecated
            /*
            std::vector<const reco::Candidate *> cands;
            for(std::vector<reco::GenParticle>::const_iterator it_gen = gens->begin();
                it_gen != gens->end(); it_gen++ ){
                cands.push_back(&*it_gen);
            }
            */

            for(std::vector<reco::GenParticle>::const_iterator it_gen=gens->begin();
                it_gen != gens->end(); it_gen++){
                if (it_gen->status() > 2 && it_gen->status() != 8) continue;//only status 1, 2, 8(simulated)
                if(GenInfo.size >= MAX_GEN){
                    fprintf(stderr,"ERROR: number of gens exceeds the size of array.\n");
                    break;;
                }

                /*
                if (
                    //(abs(it_gen->pdgId()) == 111 && it_gen->status() == 2) ||//pi 0
                    (abs(it_gen->pdgId()) == 211 && it_gen->status() == 2) ||//pi +-
                    //(abs(it_gen->pdgId()) == 311 && it_gen->status() == 2) ||//K0
                    (abs(it_gen->pdgId()) == 321 && it_gen->status() == 2) //K+-
                ) continue;//only status=1 pi+- and K+-
                */

                bool isGenSignal = false;
                //save target intermediat state particle
                if (
                    abs(int(it_gen->pdgId()/100) % 100) == 4  ||//c menson
                    abs(int(it_gen->pdgId()/100) % 100) == 5  ||//b menson
                    abs(it_gen->pdgId()) == 511 ||//B_0
                    abs(it_gen->pdgId()) == 521 ||//B_+-
                    abs(it_gen->pdgId()) == 531 ||//B_s
                    abs(it_gen->pdgId()) == 130 ||//KL
                    //abs(it_gen->pdgId()) == 311 ||//K0
                    //abs(it_gen->pdgId()) == 321 ||//K+
                    //abs(it_gen->pdgId()) == 310 ||//KS
                    //abs(it_gen->pdgId()) == 313 ||//K*0(892)
                    //abs(it_gen->pdgId()) == 323 ||//K*+-(892)
                    //abs(it_gen->pdgId()) == 333 ||//phi(1020)
                    it_gen->pdgId() == 443      ||//Jpsi
                    it_gen->pdgId() == 100443   ||//Psi(2S)
                    it_gen->pdgId() == 553      ||//Upsilon
                    it_gen->pdgId() == 100553     //Upsilon(2S)
                   ) isGenSignal = true;//b, c, s mesons

                if (abs(it_gen->pdgId()) == 13) isGenSignal = true;//all mu

                if (
                    abs(int(it_gen->pdgId()/100) % 100) == 3  ||//s menson
                    abs(it_gen->pdgId()) == 111 || //pi0
                    abs(it_gen->pdgId()) == 211 //pi+
                    ){
                    reco::GenParticle _deRef = (*it_gen);
                    reco::Candidate* Myself = dynamic_cast<reco::Candidate*>(&_deRef);
                    //std::cout<<Myself->pdgId()<<"-----------"<<std::endl;
                    isGenSignal = (Functs.GetAncestor(Myself, 5) | Functs.GetAncestor(Myself, 4));
                }//all pi and K from b or c meson
                if (!isGenSignal) continue;

                /*deprecated
                int iMo1 = -1,  iMo2 = -1,  iDa1 = -1,  iDa2 = -1;
                for(std::vector<const reco::Candidate *>::iterator iCands = cands.begin();
                    iCands != cands.end(); iCands++){
                    if (it_gen->numberOfMothers() >= 2){
                        if (it_gen->mother(0) == *iCands)
                            iMo1 = iCands - cands.begin();
                        if (it_gen->mother(1) == *iCands)
                            iMo2 = iCands - cands.begin();
                    }else if(it_gen->numberOfMothers() == 1){
                        if (it_gen->mother(0) == *iCands)
                            iMo1 = iCands - cands.begin();
                    }
                    if (it_gen->numberOfDaughters() >= 2){
                        if (it_gen->daughter(0) == *iCands)
                            iDa1 = iCands - cands.begin();
                        else if (it_gen->daughter(1) == *iCands)
                            iDa2 = iCands - cands.begin();
                    }else if(it_gen->numberOfDaughters() == 1){
                        if (it_gen->daughter(0) == *iCands)
                            iDa1 = iCands - cands.begin();
                    }
                }
                */
                //Find all other particle in TrackInfo
                //printf("-----*****DEBUG:Start of matching.\n");
                //if (abs(it_gen->pdgId()) == 13){
                    //printf("-----*****DEBUG:Entered muon matching block.\n");
                    for(int muonIdx = 0; muonIdx < MuonInfo.size; muonIdx++){
                        // match by pat::Muon
                        if (genMuonPtr[muonIdx] == 0) continue;
                        if (it_gen->p4() == genMuonPtr[muonIdx]->p4()){
                            MuonInfo.geninfo_index[muonIdx] = GenInfo.size;
                            //printf("-----*****DEBUG:[Mu]Tar.Pt /Ref.Pt = %9f/%9f\n",
                            //    MuonInfo.pt [muonIdx],it_gen->pt ());
                            break;
                        }
                    }
                //}
                //else{
                    //printf("-----*****DEBUG:Entered pion matching block.\n");
                    for(int trackIdx = 0; trackIdx < TrackInfo.size; trackIdx++){
                        //match by pat::GenericParticle
                        if (genTrackPtr[trackIdx] == 0 ) continue;
                        if (it_gen->p4() == genTrackPtr[trackIdx]->p4()){
                            TrackInfo.geninfo_index[trackIdx] = GenInfo.size;
                            break;
                        }
                    }
                //}

                GenInfo.index[GenInfo.size]         = GenInfo.size;
                GenInfo.handle_index[GenInfo.size]  = it_gen-gens->begin();
                GenInfo.pt[GenInfo.size]            = it_gen->pt();
                GenInfo.eta[GenInfo.size]           = it_gen->eta();
                GenInfo.phi[GenInfo.size]           = it_gen->phi();
                GenInfo.mass[GenInfo.size]          = it_gen->mass();
                GenInfo.pdgId[GenInfo.size]         = it_gen->pdgId();
                GenInfo.status[GenInfo.size]        = it_gen->status();
                GenInfo.nMo[GenInfo.size]           = it_gen->numberOfMothers();
                GenInfo.nDa[GenInfo.size]           = it_gen->numberOfDaughters();
                //GenInfo.mo1[GenInfo.size]           = iMo1;//To be matched later.
                //GenInfo.mo2[GenInfo.size]           = iMo2;
                //GenInfo.da1[GenInfo.size]           = iDa1;
                //GenInfo.da2[GenInfo.size]           = iDa2;
                GenInfo.mo1[GenInfo.size]           = -1;//To be matched later.
                GenInfo.mo2[GenInfo.size]           = -1;
                GenInfo.da1[GenInfo.size]           = -1;
                GenInfo.da2[GenInfo.size]           = -1;
                GenInfo.da3[GenInfo.size]           = -1;
                GenInfo.da4[GenInfo.size]           = -1;
                GenInfo.size++;
                sel_cands.push_back(&*it_gen);
            }
            //printf("-----*****DEBUG:End of gens loop.\n");

            int geninfo_idx = 0;
            for(std::vector<const reco::Candidate *>::iterator sCands = sel_cands.begin();
                sCands != sel_cands.end(); sCands++){
                geninfo_idx = int(sCands-sel_cands.begin());
                for(int nGenMo = 0; nGenMo < std::min(2,int((*sCands)->numberOfMothers())); nGenMo++){
                //if((*sCands)->numberOfMothers()==1){
                    for(std::vector<const reco::Candidate *>::iterator mCands = sel_cands.begin();
                    mCands != sel_cands.end(); mCands++){
                        if((*sCands)->mother(nGenMo) == *mCands){
                        //if((*sCands)->mother(0) == *mCands){
                            if(nGenMo == 0) GenInfo.mo1[geninfo_idx] = int(mCands-sel_cands.begin());
                            if(nGenMo == 1) GenInfo.mo2[geninfo_idx] = int(mCands-sel_cands.begin());
                        }
                    }
                }
                for(int nGenDa = 0; nGenDa < std::min(4,int((*sCands)->numberOfDaughters())); nGenDa++){
                    for(std::vector<const reco::Candidate *>::iterator mCands = sel_cands.begin();
                    mCands != sel_cands.end(); mCands++){
                        if((*sCands)->daughter(nGenDa) == *mCands){
                            if(nGenDa == 0) GenInfo.da1[geninfo_idx] = int(mCands-sel_cands.begin());
                            if(nGenDa == 1) GenInfo.da2[geninfo_idx] = int(mCands-sel_cands.begin());
                            if(nGenDa == 2) GenInfo.da3[geninfo_idx] = int(mCands-sel_cands.begin());
                            if(nGenDa == 3) GenInfo.da4[geninfo_idx] = int(mCands-sel_cands.begin());
                        }
                    }
                }
            }
            /*deprecated
            //Pass handle_index to igen
            for(int igen = 0; igen < GenInfo.size; igen++){
                int iMo1 = GenInfo.mo1[igen];
                int iMo2 = GenInfo.mo2[igen];
                int iDa1 = GenInfo.da1[igen];
                int iDa2 = GenInfo.da2[igen];
                for(int k = 0; k < GenInfo.size; k++){
                    if (iMo1 == GenInfo.handle_index[k])
                        GenInfo.mo1[igen] = k;
                    else if (iMo2 == GenInfo.handle_index[k])
                        GenInfo.mo2[igen] = k;
                    else if (iDa1 == GenInfo.handle_index[k])
                        GenInfo.da1[igen] = k;
                    else if (iDa2 == GenInfo.handle_index[k])
                        GenInfo.da2[igen] = k;
                }
                //In case that GEN particles are omitted from GenInfo
                //handle_index couldn't be the same as igen
                //since the very first proton pair has status 3.
                if (iMo1 == GenInfo.mo1[igen])
                    GenInfo.mo1[igen] = -1;
                if (iMo2 == GenInfo.mo2[igen])
                    GenInfo.mo2[igen] = -1;
                if (iDa1 == GenInfo.da1[igen])
                    GenInfo.da1[igen] = -1;
                if (iDa2 == GenInfo.da2[igen])
                   GenInfo.da2[igen] = -1;
            }
            //printf("-----*****DEBUG:End of IndexToIgen\n");
            */
        }//isRealData}}}
        //printf("-----*****DEBUG:End of GenInfo.\n");
        //std::cout<<"Start to fill!\n";

    }//try
    catch (std::exception & err){
            std::cout  << "Exception during event number: " << iEvent.id()
                << "\n" << err.what() << "\n";
    }//catch 
    root->Fill();
    //std::cout<<"filled!\n";
    
    //Made a Bntuple on the fly
    int ifchannel[7];
    ifchannel[0] = 1; //jpsi+Kp
    ifchannel[1] = 1; //jpsi+pi
    ifchannel[2] = 1; //jpsi+Ks(pi+,pi-)
    ifchannel[3] = 1; //jpsi+K*(K+,pi-)
    ifchannel[4] = 1; //jpsi+K*(K-,pi+)
    ifchannel[5] = 1; //jpsi+phi(K+,K-)
    ifchannel[6] = 1; //jpsi+pi pi <= psi', X(3872), Bs->J/psi f0
    bool REAL = ((!iEvent.isRealData() && RunOnMC_) ? false:true);
    Bntuple->makeNtuple(ifchannel, REAL, &EvtInfo, &VtxInfo, &MuonInfo, &TrackInfo, &BInfo, &GenInfo, nt0, nt1, nt2, nt3, nt5, nt6);
    if(!REAL) Bntuple->fillGenTree(ntGen, &GenInfo);
}

// ------------ method called once each job just after ending the event loop  ------------{{{
void Bfinder::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void Bfinder::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void Bfinder::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void Bfinder::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void Bfinder::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Bfinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}//}}}

//BranchOut2MuTk{{{
void Bfinder::BranchOut2MuTk(
    BInfoBranches &BInfo, 
    std::vector<pat::GenericParticle> input_tracks, 
    reco::Vertex thePrimaryV,
    std::vector<bool> isNeededTrack,
    TLorentzVector v4_mu1, 
    TLorentzVector v4_mu2,
    reco::TransientTrack muonPTT,
    reco::TransientTrack muonMTT,
    std::vector<int> &B_counter,
    float *mass_window,
    float MuMu_MASS,
    float Tk_MASS,
    int channel_number

){
  if(channel_number > (int)Bchannel_.size()){ printf("Exceeding defined # of channel, exit"); return;}
  float chi = 0.;
  float ndf = 0.;
  int tk1_hindex = -1;
  KinematicParticleFactoryFromTransientTrack pFactory;
  ParticleMass muon_mass = MUON_MASS; //pdg mass
  float muon_sigma = Functs.getParticleSigma(muon_mass);

  for(std::vector<pat::GenericParticle>::const_iterator tk_it1=input_tracks.begin();
      tk_it1 != input_tracks.end() ; tk_it1++){
      if (BInfo.size >= MAX_XB) break;
      tk1_hindex = int(tk_it1 - input_tracks.begin());
      if(tk1_hindex>=int(isNeededTrack.size())) break;
      if (!isNeededTrack[tk1_hindex]) continue;
      if (abs(tk_it1->charge()) != 1) continue;
      
      TLorentzVector v4_tk1;
      v4_tk1.SetPtEtaPhiM(tk_it1->pt(),tk_it1->eta(),tk_it1->phi(),KAON_MASS);
  
      //if ((v4_mu1+v4_mu2+v4_tk1).Mag()<mass_window[0]-0.2 || (v4_mu1+v4_mu2+v4_tk1).Mag()>mass_window[1]+0.2) continue;
      if ((v4_mu1+v4_mu2+v4_tk1).Mag()<mass_window[0] || (v4_mu1+v4_mu2+v4_tk1).Mag()>mass_window[1]) continue;
      XbMassCutLevel[channel_number-1]->Fill(0);
      if((v4_mu1+v4_mu2+v4_tk1).Pt()<bPtCut_)continue;
      XbMassCutLevel[channel_number-1]->Fill(1);
      
      reco::TransientTrack kaonTT(tk_it1->track(), &(*bField) );
      if (!kaonTT.isValid()) continue;
      XbMassCutLevel[channel_number-1]->Fill(2);
      
      ParticleMass kaon_mass = Tk_MASS;
      float kaon_sigma = Functs.getParticleSigma(kaon_mass);
      
      std::vector<RefCountedKinematicParticle> Xb_candidate;
      Xb_candidate.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
      Xb_candidate.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));
      Xb_candidate.push_back(pFactory.particle(kaonTT,kaon_mass,chi,ndf,kaon_sigma));
      RefCountedKinematicTree xbVFT;

      double MaximumDoca = Functs.getMaxDoca(Xb_candidate);
      if (MaximumDoca > MaxDocaCut_) continue;
      XbMassCutLevel[channel_number-1]->Fill(3);
      
      ParticleMass uj_mass = MuMu_MASS;
      MultiTrackKinematicConstraint *uj_c = new TwoTrackMassKinematicConstraint(uj_mass);
      KinematicConstrainedVertexFitter kcvFitter;
      xbVFT = kcvFitter.fit(Xb_candidate, uj_c);
      if (!xbVFT->isValid()) continue;
      XbMassCutLevel[channel_number-1]->Fill(4);

      xbVFT->movePointerToTheTop();
      RefCountedKinematicParticle     xbVFP       = xbVFT->currentParticle();
      RefCountedKinematicVertex       xbVFPvtx    = xbVFT->currentDecayVertex();
      if (!xbVFPvtx->vertexIsValid()) continue;
      XbMassCutLevel[channel_number-1]->Fill(5);

      std::vector<RefCountedKinematicParticle> xCands  = xbVFT->finalStateParticles();
      
      double chi2_prob = TMath::Prob(xbVFPvtx->chiSquared(),xbVFPvtx->degreesOfFreedom());
      if (chi2_prob < VtxChiProbCut_) continue;
      XbMassCutLevel[channel_number-1]->Fill(6);
      
      //if (xbVFP->currentState().mass()<mass_window[0] || xbVFP->currentState().mass()>mass_window[1]) continue;
      
      TLorentzVector xb_4vec,xb_mu1_4vec,xb_mu2_4vec,xb_tk1_4vec,xb_tk2_4vec;
      xb_4vec.SetPxPyPzE(xbVFP->currentState().kinematicParameters().momentum().x(),
                         xbVFP->currentState().kinematicParameters().momentum().y(),
                         xbVFP->currentState().kinematicParameters().momentum().z(),
                         xbVFP->currentState().kinematicParameters().energy());
      xb_mu1_4vec.SetPxPyPzE(xCands[0]->currentState().kinematicParameters().momentum().x(),
                             xCands[0]->currentState().kinematicParameters().momentum().y(),
                             xCands[0]->currentState().kinematicParameters().momentum().z(),
                             xCands[0]->currentState().kinematicParameters().energy());
      xb_mu2_4vec.SetPxPyPzE(xCands[1]->currentState().kinematicParameters().momentum().x(),
                             xCands[1]->currentState().kinematicParameters().momentum().y(),
                             xCands[1]->currentState().kinematicParameters().momentum().z(),
                             xCands[1]->currentState().kinematicParameters().energy());
      xb_tk1_4vec.SetPxPyPzE(xCands[2]->currentState().kinematicParameters().momentum().x(),
                             xCands[2]->currentState().kinematicParameters().momentum().y(),
                             xCands[2]->currentState().kinematicParameters().momentum().z(),
                             xCands[2]->currentState().kinematicParameters().energy());
      
      BInfo.index[BInfo.size]   = BInfo.size;
      BInfo.mass[BInfo.size]    = xb_4vec.Mag();
      BInfo.pt[BInfo.size]    = xb_4vec.Pt();
      BInfo.eta[BInfo.size]    = xb_4vec.Eta();
      BInfo.phi[BInfo.size]    = xb_4vec.Phi();
      BInfo.px[BInfo.size]      = xb_4vec.Px();
      BInfo.py[BInfo.size]      = xb_4vec.Py();
      BInfo.pz[BInfo.size]      = xb_4vec.Pz();
      BInfo.pxE[BInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(3,3));
      BInfo.pyE[BInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(4,4));
      BInfo.pzE[BInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(5,5));
      BInfo.MaxDoca[BInfo.size]         = MaximumDoca;

      VertexDistance3D a3d;
      //https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_0/RecoVertex/VertexTools/src/VertexDistance3D.cc
      BInfo.svpvDistance[BInfo.size] = a3d.distance(thePrimaryV,xbVFPvtx->vertexState()).value();
      BInfo.svpvDisErr[BInfo.size] = a3d.distance(thePrimaryV,xbVFPvtx->vertexState()).error();
      if( (BInfo.svpvDistance[BInfo.size]/BInfo.svpvDisErr[BInfo.size]) < svpvDistanceCut_) continue;
      XbMassCutLevel[channel_number-1]->Fill(7);

      reco::Vertex::Point vp1(thePrimaryV.position().x(), thePrimaryV.position().y(), 0.);
      reco::Vertex::Point vp2(xbVFPvtx->vertexState().position().x(), xbVFPvtx->vertexState().position().y(), 0.);
      ROOT::Math::SVector<double, 6> sv1(thePrimaryV.covariance(0,0), thePrimaryV.covariance(0,1), thePrimaryV.covariance(1,1), 0., 0., 0.);
      ROOT::Math::SVector<double, 6> sv2(xbVFPvtx->vertexState().error().cxx(), xbVFPvtx->vertexState().error().cyx(), xbVFPvtx->vertexState().error().cyy(), 0., 0., 0.);
      reco::Vertex::Error ve1(sv1);
      reco::Vertex::Error ve2(sv2);
      reco::Vertex v1(vp1, ve1);
      reco::Vertex v2(vp2, ve2);
      BInfo.svpvDistance_2D[BInfo.size] = a3d.distance(v1, v2).value();
      BInfo.svpvDisErr_2D[BInfo.size] = a3d.distance(v1, v2).error();

      BInfo.vtxX[BInfo.size]    = xbVFPvtx->position().x();
      BInfo.vtxY[BInfo.size]    = xbVFPvtx->position().y();
      BInfo.vtxZ[BInfo.size]    = xbVFPvtx->position().z();
      BInfo.vtxXErr[BInfo.size]   = xbVFPvtx->error().cxx();
      BInfo.vtxYErr[BInfo.size]   = xbVFPvtx->error().cyy();
      BInfo.vtxZErr[BInfo.size]   = xbVFPvtx->error().czz();
      BInfo.vtxYXErr[BInfo.size]  = xbVFPvtx->error().cyx();
      BInfo.vtxZXErr[BInfo.size]  = xbVFPvtx->error().czx();
      BInfo.vtxZYErr[BInfo.size]  = xbVFPvtx->error().czy();
      BInfo.vtxdof[BInfo.size]  = xbVFPvtx->degreesOfFreedom();
      BInfo.vtxchi2[BInfo.size] = xbVFPvtx->chiSquared();

      TVector3 svpvVec;
      svpvVec.SetXYZ(BInfo.vtxX[BInfo.size]-EvtInfo.PVx, BInfo.vtxY[BInfo.size]-EvtInfo.PVy, BInfo.vtxZ[BInfo.size]-EvtInfo.PVz);
      TVector3 dVec;
      dVec.SetXYZ(BInfo.px[BInfo.size], BInfo.py[BInfo.size], BInfo.pz[BInfo.size]);
      BInfo.alpha[BInfo.size] = svpvVec.Angle(dVec);
      if( BInfo.alpha[BInfo.size] > alphaCut_) continue;
      XbMassCutLevel[channel_number-1]->Fill(8);
      
      BInfo.rftk1_index[BInfo.size] = -2;
      BInfo.rftk2_index[BInfo.size] = -2;
      BInfo.rfuj_index[BInfo.size]  = BInfo.uj_size-1;
      BInfo.rftk1_index[BInfo.size] = tk1_hindex;
      BInfo.rftk2_index[BInfo.size] = tk1_hindex;
      
      BInfo.rfmu1_px[BInfo.size]=xb_mu1_4vec.Px();
      BInfo.rfmu1_py[BInfo.size]=xb_mu1_4vec.Py();
      BInfo.rfmu1_pz[BInfo.size]=xb_mu1_4vec.Pz();
      BInfo.rfmu2_px[BInfo.size]=xb_mu2_4vec.Px();
      BInfo.rfmu2_py[BInfo.size]=xb_mu2_4vec.Py();
      BInfo.rfmu2_pz[BInfo.size]=xb_mu2_4vec.Pz();
      BInfo.rftk1_px[BInfo.size]=xb_tk1_4vec.Px();
      BInfo.rftk1_py[BInfo.size]=xb_tk1_4vec.Py();
      BInfo.rftk1_pz[BInfo.size]=xb_tk1_4vec.Pz();
      BInfo.rftk2_px[BInfo.size]=-999.;
      BInfo.rftk2_py[BInfo.size]=-999.;
      BInfo.rftk2_pz[BInfo.size]=-999.;
      
      BInfo.type[BInfo.size] = channel_number;
      B_counter[channel_number-1]++;
      
      Xb_candidate.clear();
      xCands.clear();
      BInfo.size++;
  }//Tk1
}
//}}}

//BranchOut2MuX{{{
void Bfinder::BranchOut2MuX_XtoTkTk(
    BInfoBranches &BInfo, 
    std::vector<pat::GenericParticle> input_tracks, 
    reco::Vertex thePrimaryV,
    std::vector<bool> isNeededTrack,
    TLorentzVector v4_mu1, 
    TLorentzVector v4_mu2,
    reco::TransientTrack muonPTT,
    reco::TransientTrack muonMTT,
    std::vector<int> &B_counter,
    float *mass_window,
    float MuMu_MASS,
    float TkTk_MASS,
    float TkTk_window,
    float Tk1_MASS,
    float Tk2_MASS,
    int channel_number,
    int fit_option
){
    if(channel_number > (int)Bchannel_.size()){ printf("Exceeding defined # of channel, exit"); return;}
    float chi = 0.;
    float ndf = 0.;
    KinematicParticleFactoryFromTransientTrack pFactory;
    ParticleMass muon_mass = MUON_MASS; //pdg mass
    float muon_sigma = Functs.getParticleSigma(muon_mass);
    int tk1_hindex = -1;
    int tk2_hindex = -1;

    for(std::vector<pat::GenericParticle>::const_iterator tk_it1=input_tracks.begin();
        tk_it1 != input_tracks.end() ; tk_it1++){
        tk1_hindex = int(tk_it1 - input_tracks.begin());
        if(tk1_hindex>=int(isNeededTrack.size())) break;
        if (!isNeededTrack[tk1_hindex]) continue;
        if (tk_it1->charge()<0) continue;
        
        for(std::vector<pat::GenericParticle>::const_iterator tk_it2=input_tracks.begin();
            tk_it2 != input_tracks.end() ; tk_it2++){
            if (BInfo.size >= MAX_XB) break;
            if(tk2_hindex>=int(isNeededTrack.size())) break;
            tk2_hindex = int(tk_it2 - input_tracks.begin());
            if (!isNeededTrack[tk2_hindex]) continue;
            if (tk_it2->charge()>0) continue;
            
            TLorentzVector v4_tk1,v4_tk2;
            v4_tk1.SetPtEtaPhiM(tk_it1->pt(),tk_it1->eta(),tk_it1->phi(),Tk1_MASS);
            v4_tk2.SetPtEtaPhiM(tk_it2->pt(),tk_it2->eta(),tk_it2->phi(),Tk2_MASS);
            if(TkTk_MASS > 0) {if (fabs((v4_tk1+v4_tk2).Mag()-TkTk_MASS)>TkTk_window) continue;}
            //else {if (fabs((v4_tk1+v4_tk2).Mag())>TkTk_window) continue;}//if no tktk mass constrain, require it to be at least < some mass value
            XbMassCutLevel[channel_number-1]->Fill(0);
            
            //if ((v4_mu1+v4_mu2+v4_tk1+v4_tk2).Mag()<mass_window[0]-0.2 || (v4_mu1+v4_mu2+v4_tk1+v4_tk2).Mag()>mass_window[1]+0.2) continue;
            if ((v4_mu1+v4_mu2+v4_tk1+v4_tk2).Mag()<mass_window[0] || (v4_mu1+v4_mu2+v4_tk1+v4_tk2).Mag()>mass_window[1]) continue;
            XbMassCutLevel[channel_number-1]->Fill(1);
            if((v4_mu1+v4_mu2+v4_tk1+v4_tk2).Pt()<bPtCut_)continue;
            XbMassCutLevel[channel_number-1]->Fill(2);
            
            reco::TransientTrack tk1PTT(tk_it1->track(), &(*bField) );
            reco::TransientTrack tk2MTT(tk_it2->track(), &(*bField) );
            if (!tk1PTT.isValid()) continue;
            if (!tk2MTT.isValid()) continue;
            XbMassCutLevel[channel_number-1]->Fill(3);
            
            ParticleMass tk1_mass = Tk1_MASS;
            float tk1_sigma = Functs.getParticleSigma(tk1_mass);
            ParticleMass tk2_mass = Tk2_MASS;
            float tk2_sigma = Functs.getParticleSigma(tk2_mass);

            //doing tktk fit
            std::vector<RefCountedKinematicParticle> tktk_candidate;
            tktk_candidate.push_back(pFactory.particle(tk1PTT,tk1_mass,chi,ndf,tk1_sigma));
            tktk_candidate.push_back(pFactory.particle(tk2MTT,tk2_mass,chi,ndf,tk2_sigma));
            
            KinematicParticleVertexFitter   tktk_fitter;
            RefCountedKinematicTree         tktk_VFT;
            tktk_VFT = tktk_fitter.fit(tktk_candidate);
            if(TkTk_MASS > 0 && !tktk_VFT->isValid()) continue;
            XbMassCutLevel[channel_number-1]->Fill(4);
            
            tktk_VFT->movePointerToTheTop();
            RefCountedKinematicParticle tktk_VFP   = tktk_VFT->currentParticle();
            RefCountedKinematicVertex   tktk_VFPvtx = tktk_VFT->currentDecayVertex();
            if (TkTk_MASS > 0 && !tktk_VFPvtx->vertexIsValid()) continue;
            XbMassCutLevel[channel_number-1]->Fill(5);

            double chi2_prob_tktk = TMath::Prob(tktk_VFPvtx->chiSquared(),tktk_VFPvtx->degreesOfFreedom());
            if(TkTk_MASS > 0 && chi2_prob_tktk < VtxChiProbCut_) continue;
            XbMassCutLevel[channel_number-1]->Fill(6);

            std::vector<RefCountedKinematicParticle> Xb_candidate;
            Xb_candidate.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
            Xb_candidate.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));
            if(fit_option == 0){
                Xb_candidate.push_back(pFactory.particle(tk1PTT,tk1_mass,chi,ndf,tk1_sigma));
                Xb_candidate.push_back(pFactory.particle(tk2MTT,tk2_mass,chi,ndf,tk2_sigma));
            }
            else if(fit_option == 1){
                VirtualKinematicParticleFactory vFactory;
                float tktkchi = tktk_VFPvtx->chiSquared();
                float tktkndf = tktk_VFPvtx->degreesOfFreedom();
                Xb_candidate.push_back(vFactory.particle(tktk_VFP->currentState(),tktkchi,tktkndf,tktk_VFP));
            }
            RefCountedKinematicTree xbVFT;

            double MaximumDoca = Functs.getMaxDoca(Xb_candidate);
            if (MaximumDoca > MaxDocaCut_) continue;
            XbMassCutLevel[channel_number-1]->Fill(7);
            
            ParticleMass uj_mass = MuMu_MASS;
            MultiTrackKinematicConstraint *uj_c = new  TwoTrackMassKinematicConstraint(uj_mass);
            KinematicConstrainedVertexFitter kcvFitter;
            xbVFT = kcvFitter.fit(Xb_candidate, uj_c);
            if (!xbVFT->isValid()) continue;
            XbMassCutLevel[channel_number-1]->Fill(8);

            xbVFT->movePointerToTheTop();
            RefCountedKinematicParticle     xbVFP       = xbVFT->currentParticle();
            RefCountedKinematicVertex       xbVFPvtx    = xbVFT->currentDecayVertex();
            if (!xbVFPvtx->vertexIsValid()) continue;
            XbMassCutLevel[channel_number-1]->Fill(9);
            
            double chi2_prob = TMath::Prob(xbVFPvtx->chiSquared(),xbVFPvtx->degreesOfFreedom());
            if (chi2_prob < VtxChiProbCut_) continue;
            XbMassCutLevel[channel_number-1]->Fill(10);
            
            //Cut out a mass window
            //if (xbVFP->currentState().mass()<mass_window[0]|| xbVFP->currentState().mass()>mass_window[1]) continue;
            //
            std::vector<RefCountedKinematicParticle> tktkCands  = tktk_VFT->finalStateParticles();
            std::vector<RefCountedKinematicParticle> xCands  = xbVFT->finalStateParticles();
            
            TLorentzVector xb_4vec,xb_mu1_4vec,xb_mu2_4vec,tktk_4vec,xb_tk1_4vec,xb_tk2_4vec,tktk_tk1_4vec, tktk_tk2_4vec;
            xb_4vec.SetPxPyPzE(xbVFP->currentState().kinematicParameters().momentum().x(),
                               xbVFP->currentState().kinematicParameters().momentum().y(),
                               xbVFP->currentState().kinematicParameters().momentum().z(),
                               xbVFP->currentState().kinematicParameters().energy());
            tktk_4vec.SetPxPyPzE(tktk_VFP->currentState().kinematicParameters().momentum().x(),
                                 tktk_VFP->currentState().kinematicParameters().momentum().y(),
                                 tktk_VFP->currentState().kinematicParameters().momentum().z(),
                                 tktk_VFP->currentState().kinematicParameters().energy());
            
            xb_mu1_4vec.SetPxPyPzE(xCands[0]->currentState().kinematicParameters().momentum().x(),
                                   xCands[0]->currentState().kinematicParameters().momentum().y(),
                                   xCands[0]->currentState().kinematicParameters().momentum().z(),
                                   xCands[0]->currentState().kinematicParameters().energy());
            xb_mu2_4vec.SetPxPyPzE(xCands[1]->currentState().kinematicParameters().momentum().x(),
                                   xCands[1]->currentState().kinematicParameters().momentum().y(),
                                   xCands[1]->currentState().kinematicParameters().momentum().z(),
                                   xCands[1]->currentState().kinematicParameters().energy());
            tktk_tk1_4vec.SetPxPyPzE(tktkCands[0]->currentState().kinematicParameters().momentum().x(),
                                   tktkCands[0]->currentState().kinematicParameters().momentum().y(),
                                   tktkCands[0]->currentState().kinematicParameters().momentum().z(),
                                   tktkCands[0]->currentState().kinematicParameters().energy());
            tktk_tk2_4vec.SetPxPyPzE(tktkCands[1]->currentState().kinematicParameters().momentum().x(),
                                   tktkCands[1]->currentState().kinematicParameters().momentum().y(),
                                   tktkCands[1]->currentState().kinematicParameters().momentum().z(),
                                   tktkCands[1]->currentState().kinematicParameters().energy());
            if(fit_option == 0){
                xb_tk1_4vec.SetPxPyPzE(xCands[2]->currentState().kinematicParameters().momentum().x(),
                                       xCands[2]->currentState().kinematicParameters().momentum().y(),
                                       xCands[2]->currentState().kinematicParameters().momentum().z(),
                                       xCands[2]->currentState().kinematicParameters().energy());
                xb_tk2_4vec.SetPxPyPzE(xCands[3]->currentState().kinematicParameters().momentum().x(),
                                       xCands[3]->currentState().kinematicParameters().momentum().y(),
                                       xCands[3]->currentState().kinematicParameters().momentum().z(),
                                       xCands[3]->currentState().kinematicParameters().energy());
            }
            
            if(fit_option == 1){
                xb_tk1_4vec.SetPxPyPzE(xCands[2]->currentState().kinematicParameters().momentum().x(),
                                       xCands[2]->currentState().kinematicParameters().momentum().y(),
                                       xCands[2]->currentState().kinematicParameters().momentum().z(),
                                       xCands[2]->currentState().kinematicParameters().energy());
            }
            
            BInfo.index[BInfo.size]   = BInfo.size;
            BInfo.mass[BInfo.size]    = xb_4vec.Mag();
            BInfo.pt[BInfo.size]    = xb_4vec.Pt();
            BInfo.eta[BInfo.size]    = xb_4vec.Eta();
            BInfo.phi[BInfo.size]    = xb_4vec.Phi();
            BInfo.px[BInfo.size]      = xb_4vec.Px();
            BInfo.py[BInfo.size]      = xb_4vec.Py();
            BInfo.pz[BInfo.size]      = xb_4vec.Pz();
            BInfo.pxE[BInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(3,3));
            BInfo.pyE[BInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(4,4));
            BInfo.pzE[BInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(5,5));
            BInfo.MaxDoca[BInfo.size]         = MaximumDoca;

            VertexDistance3D a3d;
            //https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_0/RecoVertex/VertexTools/src/VertexDistance3D.cc
            BInfo.svpvDistance[BInfo.size] = a3d.distance(thePrimaryV,xbVFPvtx->vertexState()).value();
            BInfo.svpvDisErr[BInfo.size] = a3d.distance(thePrimaryV,xbVFPvtx->vertexState()).error();
            if( (BInfo.svpvDistance[BInfo.size]/BInfo.svpvDisErr[BInfo.size]) < svpvDistanceCut_) continue;
            XbMassCutLevel[channel_number-1]->Fill(11);
           
      
            reco::Vertex::Point vp1(thePrimaryV.position().x(), thePrimaryV.position().y(), 0.);
            reco::Vertex::Point vp2(xbVFPvtx->vertexState().position().x(), xbVFPvtx->vertexState().position().y(), 0.);
            ROOT::Math::SVector<double, 6> sv1(thePrimaryV.covariance(0,0), thePrimaryV.covariance(0,1), thePrimaryV.covariance(1,1), 0., 0., 0.);
            ROOT::Math::SVector<double, 6> sv2(xbVFPvtx->vertexState().error().cxx(), xbVFPvtx->vertexState().error().cyx(), xbVFPvtx->vertexState().error().cyy(), 0., 0., 0.);
            reco::Vertex::Error ve1(sv1);
            reco::Vertex::Error ve2(sv2);
            reco::Vertex v1(vp1, ve1);
            reco::Vertex v2(vp2, ve2);
            BInfo.svpvDistance_2D[BInfo.size] = a3d.distance(v1, v2).value();
            BInfo.svpvDisErr_2D[BInfo.size] = a3d.distance(v1, v2).error();

            BInfo.vtxX[BInfo.size]    = xbVFPvtx->position().x();
            BInfo.vtxY[BInfo.size]    = xbVFPvtx->position().y();
            BInfo.vtxZ[BInfo.size]    = xbVFPvtx->position().z();
            BInfo.vtxXErr[BInfo.size]        = xbVFPvtx->error().cxx();
            BInfo.vtxYErr[BInfo.size]        = xbVFPvtx->error().cyy();
            BInfo.vtxZErr[BInfo.size]        = xbVFPvtx->error().czz();
            BInfo.vtxYXErr[BInfo.size]     = xbVFPvtx->error().cyx();
            BInfo.vtxZXErr[BInfo.size]     = xbVFPvtx->error().czx();
            BInfo.vtxZYErr[BInfo.size]     = xbVFPvtx->error().czy();
            BInfo.vtxdof[BInfo.size]  = xbVFPvtx->degreesOfFreedom();
            BInfo.vtxchi2[BInfo.size] = xbVFPvtx->chiSquared();

            TVector3 svpvVec;
            svpvVec.SetXYZ(BInfo.vtxX[BInfo.size]-EvtInfo.PVx, BInfo.vtxY[BInfo.size]-EvtInfo.PVy, BInfo.vtxZ[BInfo.size]-EvtInfo.PVz);
            TVector3 dVec;
            dVec.SetXYZ(BInfo.px[BInfo.size], BInfo.py[BInfo.size], BInfo.pz[BInfo.size]);
            BInfo.alpha[BInfo.size] = svpvVec.Angle(dVec);
            if( BInfo.alpha[BInfo.size] > alphaCut_) continue;
            XbMassCutLevel[channel_number-1]->Fill(12);
            
            BInfo.rftk1_index[BInfo.size] = -2;
            BInfo.rftk2_index[BInfo.size] = -2;
            BInfo.rfuj_index[BInfo.size]  = BInfo.uj_size-1;
            BInfo.rftk1_index[BInfo.size] = tk1_hindex;
            BInfo.rftk2_index[BInfo.size] = tk2_hindex;
            
            //tktk fit info
            BInfo.tktk_mass[BInfo.size]    = tktk_4vec.Mag();
            BInfo.tktk_pt[BInfo.size]      = tktk_4vec.Pt();
            BInfo.tktk_eta[BInfo.size]     = tktk_4vec.Eta();
            BInfo.tktk_phi[BInfo.size]     = tktk_4vec.Phi();
            BInfo.tktk_px[BInfo.size]      = tktk_4vec.Px();
            BInfo.tktk_py[BInfo.size]      = tktk_4vec.Py();
            BInfo.tktk_pz[BInfo.size]      = tktk_4vec.Pz();
            BInfo.tktk_vtxX[BInfo.size]    = tktk_VFPvtx->position().x();
            BInfo.tktk_vtxY[BInfo.size]    = tktk_VFPvtx->position().y();
            BInfo.tktk_vtxZ[BInfo.size]    = tktk_VFPvtx->position().z();
            BInfo.tktk_vtxXErr[BInfo.size] = tktk_VFPvtx->error().cxx();
            BInfo.tktk_vtxYErr[BInfo.size] = tktk_VFPvtx->error().cyy();
            BInfo.tktk_vtxZErr[BInfo.size] = tktk_VFPvtx->error().czz();
            BInfo.tktk_vtxYXErr[BInfo.size]= tktk_VFPvtx->error().cyx();
            BInfo.tktk_vtxZXErr[BInfo.size]= tktk_VFPvtx->error().czx();
            BInfo.tktk_vtxZYErr[BInfo.size]= tktk_VFPvtx->error().czy();
            BInfo.tktk_vtxdof[BInfo.size]  = tktk_VFPvtx->degreesOfFreedom();
            BInfo.tktk_vtxchi2[BInfo.size] = tktk_VFPvtx->chiSquared();
            BInfo.tktk_rftk1_px[BInfo.size]=tktk_tk1_4vec.Px();
            BInfo.tktk_rftk1_py[BInfo.size]=tktk_tk1_4vec.Py();
            BInfo.tktk_rftk1_pz[BInfo.size]=tktk_tk1_4vec.Pz();
            BInfo.tktk_rftk2_px[BInfo.size]=tktk_tk2_4vec.Px();
            BInfo.tktk_rftk2_py[BInfo.size]=tktk_tk2_4vec.Py();
            BInfo.tktk_rftk2_pz[BInfo.size]=tktk_tk2_4vec.Pz();
            
            BInfo.rfmu1_px[BInfo.size]=xb_mu1_4vec.Px();
            BInfo.rfmu1_py[BInfo.size]=xb_mu1_4vec.Py();
            BInfo.rfmu1_pz[BInfo.size]=xb_mu1_4vec.Pz();
            BInfo.rfmu2_px[BInfo.size]=xb_mu2_4vec.Px();
            BInfo.rfmu2_py[BInfo.size]=xb_mu2_4vec.Py();
            BInfo.rfmu2_pz[BInfo.size]=xb_mu2_4vec.Pz();
            //If option == 1, this momentum is the tktk virtual particle p.
            BInfo.rftk1_px[BInfo.size]=xb_tk1_4vec.Px();
            BInfo.rftk1_py[BInfo.size]=xb_tk1_4vec.Py();
            BInfo.rftk1_pz[BInfo.size]=xb_tk1_4vec.Pz();
            if(fit_option == 0){
                BInfo.rftk2_px[BInfo.size]=xb_tk2_4vec.Px();
                BInfo.rftk2_py[BInfo.size]=xb_tk2_4vec.Py();
                BInfo.rftk2_pz[BInfo.size]=xb_tk2_4vec.Pz();
            }
            else if(fit_option == 1){
                BInfo.rftk2_px[BInfo.size]=-999;
                BInfo.rftk2_py[BInfo.size]=-999;
                BInfo.rftk2_pz[BInfo.size]=-999;
            }
            
            BInfo.type[BInfo.size] = channel_number;
            B_counter[channel_number-1]++;
            
            Xb_candidate.clear();
            xCands.clear();
            BInfo.size++;
        }//Tk2
    }//Tk1
}
//}}}

//define this as a plug-in
DEFINE_FWK_MODULE(Bfinder);
