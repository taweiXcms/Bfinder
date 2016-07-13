// vim:set ts=4 sw=4 fdm=marker et:
// Ntuplt creator for D meson related analysis.
// Maintain and contact: ta-wei wang
// Email: "tawei@mit.edu" or "ta-wei.wang@cern.ch"
#include "Bfinder/Bfinder/interface/format.h"
#include "Bfinder/Bfinder/interface/Dntuple.h"
#include "Bfinder/Bfinder/interface/utilities.h"
//
// class declaration
//

class Dfinder : public edm::EDAnalyzer
{//{{{
    public:
        explicit Dfinder(const edm::ParameterSet&);
        ~Dfinder();
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 
    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
 
        virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        virtual void endRun(edm::Run const&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        
        virtual std::vector< std::vector< std::pair<float, int> > > GetPermu(std::vector< std::pair<float, int> > InVec);
        virtual std::vector< std::vector< std::pair<float, int> > > DelDuplicate(std::vector< std::vector< std::pair<float, int> > > InVec);
        
		virtual void BranchOutNTk(
            DInfoBranches &DInfo, 
            std::vector<pat::GenericParticle> input_tracks, 
            reco::Vertex thePrimaryV,
            std::vector<int> isNeededTrackIdx,
            std::vector<int> &D_counter,
            float *mass_window,
            std::vector< std::pair<float, int> > TkMassCharge,
            double tktkRes_mass,
            double tktkRes_mass_window,
            bool doConstrainFit,
            bool SequentialFit,
            int Dchannel_number,
			int TkCombinationMethod = -1
        );

        virtual void TkCombinationPermutation(
            std::vector<pat::GenericParticle> input_tracks, 
            std::vector<int> isNeededTrackIdx,
            float *mass_window,
            std::vector< std::pair<float, int> > TkMassCharge,
            double tktkRes_mass,
            double tktkRes_mass_window,
        	std::vector< std::vector<int> > &selectedTkhidxSet,
			int Dchannel_number
        );

        virtual void TkCombinationResFast(
            std::vector<pat::GenericParticle> input_tracks, 
            std::vector<int> isNeededTrackIdx,
            float *mass_window,
            std::vector< std::pair<float, int> > TkMassCharge,
            double tktkRes_mass,
            double tktkRes_mass_window,
        	std::vector< std::vector<int> > &selectedTkhidxSet,
			int Dchannel_number
        );

        // ----------member data ---------------------------
        edm::ESHandle<MagneticField> bField;
        edm::ParameterSet theConfig;

        bool detailMode_;
        bool dropUnusedTracks_;
        std::vector<int> Dchannel_;
        //edm::InputTag hltLabel_;
        edm::EDGetTokenT<reco::GenParticleCollection> genLabel_;
        edm::EDGetTokenT<std::vector<pat::GenericParticle>> trackLabel_;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfoLabel_;
        edm::EDGetTokenT<reco::BeamSpot> bsLabel_;
        edm::EDGetTokenT<reco::VertexCollection> pvLabel_;
        double tkPtCut_;
        double tkEtaCut_;
        std::vector<double> dPtCut_;
        std::vector<double> dEtaCut_;
        std::vector<double> VtxChiProbCut_;
        std::vector<double> dCutSeparating_PtVal_;
        std::vector<double> tktkRes_svpvDistanceCut_lowptD_;
        std::vector<double> tktkRes_svpvDistanceCut_highptD_;
		std::vector<double> tktkRes_alphaCut_;
        std::vector<double> svpvDistanceCut_lowptD_;
        std::vector<double> svpvDistanceCut_highptD_;
        std::vector<double> MaxDocaCut_;
        std::vector<double> alphaCut_;
        bool RunOnMC_;
        bool doTkPreCut_;
        bool makeDntuple_;
        bool doDntupleSkim_;
        //edm::EDGetTokenT<edm::ValueMap<float>> MVAMapLabel_;


        edm::Service<TFileService> fs;
        TTree *root;
        EvtInfoBranches     EvtInfo;
        VtxInfoBranches     VtxInfo;
        TrackInfoBranches   TrackInfo;
        DInfoBranches       DInfo;
        GenInfoBranches     GenInfo;
        CommonFuncts        Functs;
        DntupleBranches     *Dntuple = new DntupleBranches;
        TTree* ntD1; 
        TTree* ntD2;
        TTree* ntD3; 
        TTree* ntD4; 
        TTree* ntD5; 
        TTree* ntD6; 
        TTree* ntGen;

        //histograms
        TH1F *TrackCutLevel;
        //How many channel
        static int const Nchannel = 20;
        std::vector<TH1F*> DMassCutLevel;
        
};//}}}

void Dfinder::beginJob()
{//{{{
    root = fs->make<TTree>("root","root");
    ntD1 = fs->make<TTree>("ntDkpi","");           Dntuple->buildDBranch(ntD1);
    ntD2 = fs->make<TTree>("ntDkpipi","");         Dntuple->buildDBranch(ntD2);
    ntD3 = fs->make<TTree>("ntDkpipipi","");       Dntuple->buildDBranch(ntD3);
    ntD4 = fs->make<TTree>("ntDPhikkpi","");       Dntuple->buildDBranch(ntD4);
    ntD5 = fs->make<TTree>("ntDD0kpipi","");       Dntuple->buildDBranch(ntD5);
    ntD6 = fs->make<TTree>("ntDD0kpipipipi","");   Dntuple->buildDBranch(ntD6);
    ntGen = fs->make<TTree>("ntGen","");           Dntuple->buildGenBranch(ntGen);
    EvtInfo.regTree(root);
    VtxInfo.regTree(root);
    TrackInfo.regTree(root, detailMode_);
    DInfo.regTree(root, detailMode_);
    GenInfo.regTree(root);
}//}}}

Dfinder::Dfinder(const edm::ParameterSet& iConfig):theConfig(iConfig)
{//{{{
    //now do what ever initialization is needed
    detailMode_ = iConfig.getParameter<bool>("detailMode");
    dropUnusedTracks_ = iConfig.getParameter<bool>("dropUnusedTracks");

    Dchannel_ = iConfig.getParameter<std::vector<int> >("Dchannel");
    genLabel_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenLabel"));
    trackLabel_ = consumes<std::vector<pat::GenericParticle>>(iConfig.getParameter<edm::InputTag>("TrackLabel"));
    //hltLabel_           = iConfig.getParameter<edm::InputTag>("HLTLabel");
    puInfoLabel_ = consumes<std::vector< PileupSummaryInfo > >(iConfig.getParameter<edm::InputTag>("PUInfoLabel"));
    bsLabel_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BSLabel"));
    pvLabel_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PVLabel"));

    tkPtCut_ = iConfig.getParameter<double>("tkPtCut");
    tkEtaCut_ = iConfig.getParameter<double>("tkEtaCut");
    dPtCut_ = iConfig.getParameter<std::vector<double> >("dPtCut");
    dEtaCut_ = iConfig.getParameter<std::vector<double> >("dEtaCut");
    VtxChiProbCut_ = iConfig.getParameter<std::vector<double> >("VtxChiProbCut");
    dCutSeparating_PtVal_ = iConfig.getParameter<std::vector<double> >("dCutSeparating_PtVal");
    tktkRes_svpvDistanceCut_lowptD_ = iConfig.getParameter<std::vector<double> >("tktkRes_svpvDistanceCut_lowptD");
    tktkRes_svpvDistanceCut_highptD_ = iConfig.getParameter<std::vector<double> >("tktkRes_svpvDistanceCut_highptD");
	tktkRes_alphaCut_ = iConfig.getParameter<std::vector<double> >("tktkRes_alphaCut");
    svpvDistanceCut_lowptD_ = iConfig.getParameter<std::vector<double> >("svpvDistanceCut_lowptD");
    svpvDistanceCut_highptD_ = iConfig.getParameter<std::vector<double> >("svpvDistanceCut_highptD");
    MaxDocaCut_ = iConfig.getParameter<std::vector<double> >("MaxDocaCut");
    alphaCut_ = iConfig.getParameter<std::vector<double> >("alphaCut");
    RunOnMC_ = iConfig.getParameter<bool>("RunOnMC");
    doTkPreCut_ = iConfig.getParameter<bool>("doTkPreCut");
    makeDntuple_ = iConfig.getParameter<bool>("makeDntuple");
    doDntupleSkim_ = iConfig.getParameter<bool>("doDntupleSkim");
    //MVAMapLabel_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("MVAMapLabel"));

    TrackCutLevel       = fs->make<TH1F>("TrackCutLevel"    , "TrackCutLevel"   , 10, 0, 10);
    for(unsigned int i = 0; i < Dchannel_.size(); i++){
        TH1F* DMassCutLevel_temp      = fs->make<TH1F>(TString::Format("DMassCutLevel_i")   ,TString::Format("DMassCutLevel_i")  , 10, 0, 10);
        DMassCutLevel.push_back(DMassCutLevel_temp);
    }
}//}}}

Dfinder::~Dfinder()
{//{{{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}//}}}

//
// member functions
//

// ------------ method called for each event  ------------
void Dfinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //checking input parameter size
    if( (Dchannel_.size() != dPtCut_.size()) || (dPtCut_.size() != dEtaCut_.size()) || (dEtaCut_.size() != VtxChiProbCut_.size()) || (VtxChiProbCut_.size() != tktkRes_svpvDistanceCut_lowptD_.size()) || (tktkRes_svpvDistanceCut_lowptD_.size() != tktkRes_svpvDistanceCut_highptD_.size()) || (tktkRes_svpvDistanceCut_highptD_.size() != svpvDistanceCut_lowptD_.size()) || (svpvDistanceCut_lowptD_.size() != svpvDistanceCut_highptD_.size()) || (svpvDistanceCut_highptD_.size() != MaxDocaCut_.size()) || (MaxDocaCut_.size() != alphaCut_.size()) || (alphaCut_.size() != tktkRes_alphaCut_.size())){
        std::cout<<"Unmatched input parameter vector size, EXIT"<<std::endl;
        return;
    }

    //std::cout << "*************************\nReconstructing event number: " << iEvent.id() << "\n";
    using namespace edm;
    using namespace reco;
    //ESHandle<MagneticField> bField;
    iSetup.get<IdealMagneticFieldRecord>().get(bField);

    // Change used muon and track collections
    edm::Handle< std::vector<pat::GenericParticle> > tks;
    iEvent.getByToken(trackLabel_, tks);

    //CLEAN all memory
    memset(&EvtInfo     ,0x00,sizeof(EvtInfo)   );
    memset(&VtxInfo     ,0x00,sizeof(VtxInfo)   );
    memset(&TrackInfo   ,0x00,sizeof(TrackInfo) );
    memset(&DInfo       ,0x00,sizeof(DInfo)    );
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
    iEvent.getByToken(bsLabel_, beamSpotHandle);
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
    iEvent.getByToken(pvLabel_, VertexHandle);

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
        iEvent.getByToken(puInfoLabel_, PUHandle);
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
    TrackInfo.size  = 0;
    DInfo.size      = 0;
    GenInfo.size    = 0;
    
    std::vector<int> D_counter;
    for(unsigned int i = 0; i < Dchannel_.size(); i++){
        D_counter.push_back(0);
    }

    std::vector<pat::GenericParticle>   input_tracks;
    input_tracks = *tks;
    try{
        const reco::GenParticle* genMuonPtr[MAX_MUON];
        memset(genMuonPtr,0x00,MAX_MUON);
        const reco::GenParticle* genTrackPtr[MAX_TRACK];
        memset(genTrackPtr,0x00,MAX_GEN);
        //standard check for validity of input data
        if (0){
        //if (input_muons.size() == 0){
            std::cout << "There's no muon : " << iEvent.id() << std::endl;
        }else{
            //std::cout << "Got " << input_muons.size() << " muons / ";
            if (input_tracks.size() == 0){
                std::cout << "There's no track: " << iEvent.id() << std::endl;
            }else{
                std::cout << "Got " << input_tracks.size() << " tracks" << std::endl;
                if (input_tracks.size() > 0){

                    //Preselect tracks{{{
                    std::vector<bool> isNeededTrack;// Are the tracks redundant?
                    std::vector<int> isNeededTrackIdx;
                    int PassedTrk = 0;
                    for(std::vector<pat::GenericParticle>::const_iterator tk_it=input_tracks.begin();
                        tk_it != input_tracks.end(); tk_it++){
                        if(PassedTrk >= MAX_TRACK){
                            fprintf(stderr,"ERROR: number of tracks exceeds the size of array.\n");
                            break;
                        }
                        isNeededTrack.push_back(false);
                        TrackCutLevel->Fill(0);//number of all tracks
                        TrackCutLevel->Fill(1);//
                        if (tk_it->pt()<tkPtCut_)                           continue;
                        TrackCutLevel->Fill(2);
                        if (fabs(tk_it->eta())>tkEtaCut_)                   continue;
                        TrackCutLevel->Fill(3);
                        //if (fabs(tk_it->eta()) > 2.5)                       continue;
                        TrackCutLevel->Fill(4);
                        if(doTkPreCut_){
                            if( !(tk_it->track()->quality(reco::TrackBase::highPurity))) continue;
                            //d0 analysis cuts
                            //if(tk_it->track()->hitPattern().numberOfValidHits() < 12) continue;
                            //if(tk_it->track()->ptError()/tk_it->track()->pt() > 0.075) continue;
                            //if(tk_it->track()->normalizedChi2()/tk_it->track()->hitPattern().trackerLayersWithMeasurement() > 0.25) continue;
                            //outdated selections
                            //if (tk_it->track()->normalizedChi2()>5)             continue;
                            //if (tk_it->p()>200 || tk_it->pt()>200)              continue;
                            //if (tk_it->track()->hitPattern().numberOfValidStripHits()<10)continue;
                            //if (tk_it->track()->hitPattern().numberOfValidPixelHits()<2) continue;
                            TrackCutLevel->Fill(5);
                        }
                        isNeededTrack[tk_it-input_tracks.begin()] = true;
                        isNeededTrackIdx.push_back(tk_it-input_tracks.begin());
                        PassedTrk++;
                    }//end of track preselection}}}
                    //printf("-----*****DEBUG:End of track preselection.\n");
                    std::cout<<"PassedTrk: "<<PassedTrk<<std::endl;
                    
                    // DInfo section{{{
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+pi-
                    //////////////////////////////////////////////////////////////////////////
                    float d0_mass_window[2] = {D0_MASS-0.2,D0_MASS+0.2};
                    
                    if(Dchannel_[0] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, d0_mass_window, PermuVec[i], -1, -1, false, false, 1, 0);
                        }
                        //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, d0_mass_window, InVec, -1, -1, false, false, 1, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K-pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[1] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(-KAON_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, d0_mass_window, PermuVec[i], -1, -1, false, false, 2, 0);
                        }
                        //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, d0_mass_window, InVec, -1, -1, false, false, 2, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K-pi+pi+
                    //////////////////////////////////////////////////////////////////////////
					float dplus_mass_window[2] = {DPLUS_MASS-0.2,DPLUS_MASS+0.2};
                    if(Dchannel_[2] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(-KAON_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(PION_MASS, 0);
                        std::pair<float, int> tk3 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dplus_mass_window, PermuVec[i], -1, -1, false, false, 3, 0);
                        }
                        //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dplus_mass_window, InVec, -1, -1, false, false, 3, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+pi-pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[3] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(-PION_MASS, 0);
                        std::pair<float, int> tk3 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dplus_mass_window, PermuVec[i], -1, -1, false, false, 4, 0);
                        }
                        //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dplus_mass_window, InVec, -1, -1, false, false, 4, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K-pi-pi+pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[4] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(-KAON_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(-PION_MASS, 0);
                        std::pair<float, int> tk3 = std::make_pair(PION_MASS, 0);
                        std::pair<float, int> tk4 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        InVec.push_back(tk4);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, d0_mass_window, PermuVec[i], -1, -1, false, false, 5, 0);
                        }
                        //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, d0_mass_window, InVec, -1, -1, false, false, 5, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+pi+pi-pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[5] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(PION_MASS, 0);
                        std::pair<float, int> tk3 = std::make_pair(-PION_MASS, 0);
                        std::pair<float, int> tk4 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        InVec.push_back(tk4);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, d0_mass_window, PermuVec[i], -1, -1, false, false, 6, 0);
                        }
                        //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, d0_mass_window, InVec, -1, -1, false, false, 6, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+K-(Phi)pi+
                    //////////////////////////////////////////////////////////////////////////
					float dsubs_mass_window[2] = {DSUBS_MASS-0.2,DSUBS_MASS+0.2};
                    if(Dchannel_[6] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(-KAON_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                        //    Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dsubs_mass_window, PermuVec[i], PHI_MASS, 0.1, false, false, 7, 0);
                        //}
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dsubs_mass_window, InVec, PHI_MASS, 0.1, false, false, 7, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+K-(Phi)pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[7] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(-KAON_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                        //    Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dsubs_mass_window, PermuVec[i], PHI_MASS, 0.1, false, false, 8, 0);
                        //}
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dsubs_mass_window, InVec, PHI_MASS, 0.1, false, false, 8, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: D0(K-pi+)pi+
                    //////////////////////////////////////////////////////////////////////////
                    float dstar_mass_window[2] = {DSTAR_MASS-0.2,DSTAR_MASS+0.2};
                    if(Dchannel_[8] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(-KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(PION_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                        //    Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dstar_mass_window, PermuVec[i], D0_MASS, 0.1, false, true, 9, 0);
                        //}
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dstar_mass_window, InVec, D0_MASS, 0.1, false, true, 9, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: D0bar(K+pi-)pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[9] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(-PION_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                        //    Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dstar_mass_window, PermuVec[i], D0_MASS, 0.1, false, true, 10, 0);
                        //}
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dstar_mass_window, InVec, D0_MASS, 0.1, false, true, 10, 1);
                    }

                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: D0(K-pi-pi+pi+)pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[10] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(-KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(-PION_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(PION_MASS, 1);
                        std::pair<float, int> tk4 = std::make_pair(PION_MASS, 1);
                        std::pair<float, int> tk5 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        InVec.push_back(tk4);
                        InVec.push_back(tk5);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                        //    Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dstar_mass_window, PermuVec[i], D0_MASS, 0.1, false, true, 11, 0);
                        //}
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dstar_mass_window, InVec, D0_MASS, 0.1, false, true, 11, 1);
                    }

                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: D0bar(K+pi+pi-pi-)pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[11] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(PION_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(-PION_MASS, 1);
                        std::pair<float, int> tk4 = std::make_pair(-PION_MASS, 1);
                        std::pair<float, int> tk5 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        InVec.push_back(tk4);
                        InVec.push_back(tk5);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                        //    Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dstar_mass_window, PermuVec[i], D0_MASS, 0.1, false, true, 12, 0);
                        //}
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dstar_mass_window, InVec, D0_MASS, 0.1, false, true, 12, 1);
                    }

                    printf("D_counter: ");
                    for(unsigned int i = 0; i < Dchannel_.size(); i++){
                        printf("%d/", D_counter[i]);
                    }
                    printf("\n");//}}}
                    //printf("-----*****DEBUG:End of DInfo.\n");

                    // TrackInfo section {{{
                    edm::Handle<edm::ValueMap<float> > mvaoutput;
                   // iEvent.getByToken(MVAMapLabel_, mvaoutput);
                    for(std::vector<pat::GenericParticle>::const_iterator tk_it=input_tracks.begin();
                        tk_it != input_tracks.end() ; tk_it++){
                        int tk_hindex = int(tk_it - input_tracks.begin());
                        if(tk_hindex>=int(isNeededTrack.size())) break;
                        if (isNeededTrack[tk_hindex]==false) continue;

                        //Create list of relative xb candidates for later filling
                        std::vector<int> listOfRelativeDCand1;//1~nXb
                        std::vector<int> listOfRelativeDCand2;//1~nXb
                        std::vector<int> listOfRelativeDCand3;//1~nXb
                        std::vector<int> listOfRelativeDCand4;//1~nXb
                        std::vector<int> listOfRelativeDCand5;//1~nXb
                        std::vector<int> listOfRelativeDResCand1;//1~nXb
                        std::vector<int> listOfRelativeDResCand2;//1~nXb
                        std::vector<int> listOfRelativeDResCand3;//1~nXb
                        std::vector<int> listOfRelativeDResCand4;//1~nXb
                        for(int d=0; d < DInfo.size; d++){
                            if(DInfo.rftk1_index[d] == tk_hindex){
                                listOfRelativeDCand1.push_back(d+1);
                            }
                            if(DInfo.rftk2_index[d] == tk_hindex){
                                listOfRelativeDCand2.push_back(d+1);
                            }
                            if(DInfo.rftk3_index[d] == tk_hindex){
                                listOfRelativeDCand3.push_back(d+1);
                            }
                            if(DInfo.rftk4_index[d] == tk_hindex){
                                listOfRelativeDCand4.push_back(d+1);
                            }
                            if(DInfo.rftk5_index[d] == tk_hindex){
                                listOfRelativeDCand5.push_back(d+1);
                            }

                            if(DInfo.tktkRes_rftk1_index[d] == tk_hindex){
                                listOfRelativeDResCand1.push_back(d+1);
                            }
                            if(DInfo.tktkRes_rftk2_index[d] == tk_hindex){
                                listOfRelativeDResCand2.push_back(d+1);
                            }
                            if(DInfo.tktkRes_rftk3_index[d] == tk_hindex){
                                listOfRelativeDResCand3.push_back(d+1);
                            }
                            if(DInfo.tktkRes_rftk4_index[d] == tk_hindex){
                                listOfRelativeDResCand4.push_back(d+1);
                            }
                        }
                        
                        if(dropUnusedTracks_ && listOfRelativeDCand1.size() == 0 && listOfRelativeDCand2.size() == 0 && listOfRelativeDCand3.size() == 0 && listOfRelativeDCand4.size() == 0 && listOfRelativeDCand5.size() == 0 && listOfRelativeDResCand1.size() == 0 && listOfRelativeDResCand2.size() == 0 && listOfRelativeDResCand3.size() == 0 && listOfRelativeDResCand4.size() == 0) continue;//drop unused tracks

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
                        //TrackInfo.trkMVAVal      [TrackInfo.size] = (*mvaoutput)[tk_it->track()];
                        TrackInfo.trkAlgo        [TrackInfo.size] = tk_it->track()->algo();
                        TrackInfo.originalTrkAlgo[TrackInfo.size] = tk_it->track()->originalAlgo();

                        if(tk_it->track().isNonnull()){
                            for(int tq = 0; tq < reco::TrackBase::qualitySize; tq++){
                            if (tk_it->track()->quality(static_cast<reco::TrackBase::TrackQuality>(tq))) TrackInfo.trackQuality[TrackInfo.size] += 1 << (tq);
                        }}

                        if (!iEvent.isRealData() && RunOnMC_)
                            genTrackPtr [TrackInfo.size] = tk_it->genParticle();

                        // Fill the same list for DInfo
                        for(unsigned int iCands=0; iCands < listOfRelativeDCand1.size(); iCands++){
                            DInfo.rftk1_index[listOfRelativeDCand1[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDCand2.size(); iCands++){
                            DInfo.rftk2_index[listOfRelativeDCand2[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDCand3.size(); iCands++){
                            DInfo.rftk3_index[listOfRelativeDCand3[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDCand4.size(); iCands++){
                            DInfo.rftk4_index[listOfRelativeDCand4[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDCand5.size(); iCands++){
                            DInfo.rftk5_index[listOfRelativeDCand5[iCands]-1] = TrackInfo.size;
                        }

                        for(unsigned int iCands=0; iCands < listOfRelativeDResCand1.size(); iCands++){
                            DInfo.tktkRes_rftk1_index[listOfRelativeDResCand1[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDResCand2.size(); iCands++){
                            DInfo.tktkRes_rftk2_index[listOfRelativeDResCand2[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDResCand3.size(); iCands++){
                            DInfo.tktkRes_rftk3_index[listOfRelativeDResCand3[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDResCand4.size(); iCands++){
                            DInfo.tktkRes_rftk4_index[listOfRelativeDResCand4[iCands]-1] = TrackInfo.size;
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
            iEvent.getByToken(genLabel_, gens);

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
				
				if( !isGenSignal && abs(it_gen->pdgId()) > 400 ){ //should be OK to require the PID > 400
					reco::GenParticle _deRef = (*it_gen);
					reco::Candidate* Myself = dynamic_cast<reco::Candidate*>(&_deRef);
					isGenSignal = Functs.GetDescendant(Myself, 4);
				}//other particles (with pid > 400) have D meson in descendant chain

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
                //printf("-----*****DEBUG:Entered pion matching block.\n");
                for(int trackIdx = 0; trackIdx < TrackInfo.size; trackIdx++){
                    //match by pat::GenericParticle
                    if (genTrackPtr[trackIdx] == 0 ) continue;
                    if (it_gen->p4() == genTrackPtr[trackIdx]->p4()){
                        TrackInfo.geninfo_index[trackIdx] = GenInfo.size;
                        break;
                    }
                }

                GenInfo.index[GenInfo.size]         = GenInfo.size;
                GenInfo.handle_index[GenInfo.size]  = it_gen-gens->begin();
                GenInfo.pt[GenInfo.size]            = it_gen->pt();
                GenInfo.eta[GenInfo.size]           = it_gen->eta();
                GenInfo.phi[GenInfo.size]           = it_gen->phi();
                GenInfo.mass[GenInfo.size]          = it_gen->mass();
                GenInfo.pdgId[GenInfo.size]         = it_gen->pdgId();
                GenInfo.status[GenInfo.size]        = it_gen->status();
                GenInfo.collisionId[GenInfo.size]   = it_gen->collisionId();
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
 
    //Made a Dntuple on the fly   
    if(makeDntuple_){
        int isDchannel[12];
        isDchannel[0] = 1; //k+pi-
        isDchannel[1] = 1; //k-pi+
        isDchannel[2] = 1; //k-pi+pi+
        isDchannel[3] = 1; //k+pi-pi-
        isDchannel[4] = 1; //k-pi-pi+pi+
        isDchannel[5] = 1; //k+pi+pi-pi-
        isDchannel[6] = 1; 
        isDchannel[7] = 1; 
        isDchannel[8] = 1; 
        isDchannel[9] = 1; 
        isDchannel[10] = 1; 
        isDchannel[11] = 1;
        bool REAL = ((!iEvent.isRealData() && RunOnMC_) ? false:true);
        Dntuple->makeDNtuple(isDchannel, REAL, doDntupleSkim_, &EvtInfo, &VtxInfo, &TrackInfo, &DInfo, &GenInfo, ntD1, ntD2, ntD3, ntD4, ntD5, ntD6);
        if(!REAL) Dntuple->fillDGenTree(ntGen, &GenInfo);
    }

}

// ------------ method called once each job just after ending the event loop  ------------{{{
void Dfinder::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void Dfinder::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void Dfinder::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void Dfinder::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void Dfinder::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Dfinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}//}}}

//Functions{{{
std::vector< std::vector< std::pair<float, int> > > Dfinder::GetPermu(std::vector< std::pair<float, int> > InVec){
    if(InVec.size() == 1){
        std::vector< std::vector< std::pair<float, int> > > OneEntryVec;
        OneEntryVec.push_back(InVec);
        return OneEntryVec;
    }
    std::vector< std::vector< std::pair<float, int> > > NPermu;
    for(unsigned int i = 0; i < InVec.size(); i++){
        std::vector< std::pair<float, int> > copy;
        copy = InVec;
        copy.erase(copy.begin()+i);
        std::vector< std::vector< std::pair<float, int> > > Nminus1Permu;
        Nminus1Permu = GetPermu(copy);
        for(unsigned int j = 0; j < Nminus1Permu.size(); j++){
            Nminus1Permu[j].push_back(InVec[i]);
        }
        NPermu.insert(NPermu.end(), Nminus1Permu.begin(), Nminus1Permu.end());
    }
    return NPermu;
}

std::vector< std::vector< std::pair<float, int> > > Dfinder::DelDuplicate(std::vector< std::vector< std::pair<float, int> > > InVec){
    std::vector< std::vector< std::pair<float, int> > > CleanedVec;
    for(unsigned int i = 0; i < InVec.size(); i++){
        bool IsDuplicate = false;
        for(unsigned int j = 0; j < CleanedVec.size(); j++){
            bool ADuplicate = true;
            for(unsigned int k = 0; k < InVec[i].size(); k++){
                if(InVec[i][k] != CleanedVec[j][k]) {
                    ADuplicate = false;
                    break;
                }
            }
            if(ADuplicate) {
                IsDuplicate = ADuplicate;
                break;
            }
        }
        if(!IsDuplicate) CleanedVec.push_back(InVec[i]);
    }
    return CleanedVec;
}

void Dfinder::TkCombinationPermutation(
    std::vector<pat::GenericParticle> input_tracks, 
    std::vector<int> isNeededTrackIdx,
    float *mass_window,
    std::vector< std::pair<float, int> > TkMassCharge,
    double tktkRes_mass,
    double tktkRes_mass_window,
	std::vector< std::vector<int> > &selectedTkhidxSet,
	int Dchannel_number
){
    int tk1_hindex = -1;
    int tk2_hindex = -1;
    int tk3_hindex = -1;
    int tk4_hindex = -1;
    int tk5_hindex = -1;

    std::vector<int> selectedTkhidx;
    TLorentzVector v4_D, v4_Res;//unfitted D and Res
    TLorentzVector v4_tk1, v4_tk2, v4_tk3, v4_tk4, v4_tk5;
    TLorentzVector v4_Restk1, v4_Restk2, v4_Restk3, v4_Restk4, v4_Restk5;
    //for(std::vector<pat::GenericParticle>::const_iterator tk_it1=input_tracks.begin();
    //        tk_it1 != input_tracks.end() ; tk_it1++){
    for(int tk1idx = 0; tk1idx < (int)isNeededTrackIdx.size(); tk1idx++){
        v4_D.Clear(); v4_Res.Clear();
        //tk1_hindex = int(tk_it1 - input_tracks.begin());
        //if(tk1_hindex>=int(isNeededTrack.size())) break;
        //if(!isNeededTrack[tk1_hindex]) continue;
        //if(tk_it1->charge()*TkMassCharge[0].first<0) continue;
        tk1_hindex = isNeededTrackIdx[tk1idx];
        if(input_tracks[tk1_hindex].charge()*TkMassCharge[0].first<0) continue;
        v4_tk1.SetPtEtaPhiM(input_tracks[tk1_hindex].pt(),input_tracks[tk1_hindex].eta(),input_tracks[tk1_hindex].phi(),fabs(TkMassCharge[0].first));
        v4_D = v4_tk1;
        if(TkMassCharge[0].second==1) v4_Restk1.SetPtEtaPhiM(input_tracks[tk1_hindex].pt(),input_tracks[tk1_hindex].eta(),input_tracks[tk1_hindex].phi(),fabs(TkMassCharge[0].first));
        else v4_Restk1.SetPtEtaPhiM(0,0,0,0);
        v4_Res = v4_Restk1;

        //for(std::vector<pat::GenericParticle>::const_iterator tk_it2=tk_it1+1;
        //        tk_it2 != input_tracks.end() ; tk_it2++){
        for(int tk2idx = tk1idx+1; tk2idx < (int)isNeededTrackIdx.size(); tk2idx++){
            //tk2_hindex = int(tk_it2 - input_tracks.begin());
            //if(tk2_hindex>=int(isNeededTrack.size())) break;
            //if(!isNeededTrack[tk2_hindex]) continue;
            //if(tk_it2->charge()*TkMassCharge[1].first<0) continue;
            tk2_hindex = isNeededTrackIdx[tk2idx];
            if(input_tracks[tk2_hindex].charge()*TkMassCharge[1].first<0) continue;
            if(tk2_hindex==tk1_hindex) continue;
            v4_tk2.SetPtEtaPhiM(input_tracks[tk2_hindex].pt(),input_tracks[tk2_hindex].eta(),input_tracks[tk2_hindex].phi(),fabs(TkMassCharge[1].first));
            v4_D = v4_tk1 + v4_tk2;
            if(TkMassCharge[1].second==1) v4_Restk2.SetPtEtaPhiM(input_tracks[tk2_hindex].pt(),input_tracks[tk2_hindex].eta(),input_tracks[tk2_hindex].phi(),fabs(TkMassCharge[1].first));
            else v4_Restk2.SetPtEtaPhiM(0,0,0,0);
            v4_Res = v4_Restk1 + v4_Restk2;
            if(TkMassCharge.size()==2){
                //cut mass window before fit
                if(tktkRes_mass > 0) {if (fabs(v4_Res.Mag()-tktkRes_mass) > tktkRes_mass_window) continue;}
                DMassCutLevel[Dchannel_number-1]->Fill(0);
                if(v4_D.Mag()<mass_window[0] || v4_D.Mag()>mass_window[1]) continue;
                if(v4_D.Pt() < dPtCut_[Dchannel_number-1])continue;
                DMassCutLevel[Dchannel_number-1]->Fill(1);
                //if(fabs(v4_D.Eta()) > dEtaCut_[Dchannel_number-1])continue;
                DMassCutLevel[Dchannel_number-1]->Fill(2);
                selectedTkhidx.push_back(tk1_hindex);
                selectedTkhidx.push_back(tk2_hindex);
                selectedTkhidxSet.push_back(selectedTkhidx);
                selectedTkhidx.clear();
                continue;
            }
            //for(std::vector<pat::GenericParticle>::const_iterator tk_it3=tk_it2+1;
            //        tk_it3 != input_tracks.end() ; tk_it3++){
            for(int tk3idx = tk2idx+1; tk3idx < (int)isNeededTrackIdx.size(); tk3idx++){
                //tk3_hindex = int(tk_it3 - input_tracks.begin());
                //if(tk3_hindex>=int(isNeededTrack.size())) break;
                //if(!isNeededTrack[tk3_hindex]) continue;
                //if(tk_it3->charge()*TkMassCharge[2].first<0) continue;
                tk3_hindex = isNeededTrackIdx[tk3idx];
                if(input_tracks[tk3_hindex].charge()*TkMassCharge[2].first<0) continue;
                if(tk3_hindex==tk1_hindex) continue;
                if(tk3_hindex==tk2_hindex) continue;
                v4_tk3.SetPtEtaPhiM(input_tracks[tk3_hindex].pt(),input_tracks[tk3_hindex].eta(),input_tracks[tk3_hindex].phi(),fabs(TkMassCharge[2].first));
                v4_D = v4_tk1 + v4_tk2 + v4_tk3;
                if(TkMassCharge[2].second==1) v4_Restk3.SetPtEtaPhiM(input_tracks[tk3_hindex].pt(),input_tracks[tk3_hindex].eta(),input_tracks[tk3_hindex].phi(),fabs(TkMassCharge[2].first));
                else v4_Restk3.SetPtEtaPhiM(0,0,0,0);
                v4_Res = v4_Restk1 + v4_Restk2 + v4_Restk3;
                if(TkMassCharge.size()==3){
                    if(tktkRes_mass > 0) {if (fabs(v4_Res.Mag()-tktkRes_mass) > tktkRes_mass_window) continue;}
                    DMassCutLevel[Dchannel_number-1]->Fill(0);
                    if(v4_D.Mag()<mass_window[0] || v4_D.Mag()>mass_window[1]) continue;
                    if(v4_D.Pt() < dPtCut_[Dchannel_number-1])continue;
                    DMassCutLevel[Dchannel_number-1]->Fill(1);
                    //if(fabs(v4_D.Eta()) > dEtaCut_[Dchannel_number-1])continue;
                    DMassCutLevel[Dchannel_number-1]->Fill(2);
                    selectedTkhidx.push_back(tk1_hindex);
                    selectedTkhidx.push_back(tk2_hindex);
                    selectedTkhidx.push_back(tk3_hindex);
                    selectedTkhidxSet.push_back(selectedTkhidx);
                    selectedTkhidx.clear();
                    continue;
                }
                //for(std::vector<pat::GenericParticle>::const_iterator tk_it4=tk_it3+1;
                //        tk_it4 != input_tracks.end() ; tk_it4++){
                for(int tk4idx = tk3idx+1; tk4idx < (int)isNeededTrackIdx.size(); tk4idx++){
                    //tk4_hindex = int(tk_it4 - input_tracks.begin());
                    //if(tk4_hindex>=int(isNeededTrack.size())) break;
                    //if(!isNeededTrack[tk4_hindex]) continue;
                    //if(tk_it4->charge()*TkMassCharge[3].first<0) continue;
                    tk4_hindex = isNeededTrackIdx[tk4idx];
                    if(input_tracks[tk4_hindex].charge()*TkMassCharge[3].first<0) continue;
                    if(tk4_hindex==tk1_hindex) continue;
                    if(tk4_hindex==tk2_hindex) continue;
                    if(tk4_hindex==tk3_hindex) continue;
                    v4_tk4.SetPtEtaPhiM(input_tracks[tk4_hindex].pt(),input_tracks[tk4_hindex].eta(),input_tracks[tk4_hindex].phi(),fabs(TkMassCharge[3].first));
                    v4_D = v4_tk1 + v4_tk2 + v4_tk3 + v4_tk4;
                    if(TkMassCharge[3].second==1) v4_Restk4.SetPtEtaPhiM(input_tracks[tk4_hindex].pt(),input_tracks[tk4_hindex].eta(),input_tracks[tk4_hindex].phi(),fabs(TkMassCharge[3].first));
                    else v4_Restk4.SetPtEtaPhiM(0,0,0,0);
                    v4_Res = v4_Restk1 + v4_Restk2 + v4_Restk3 + v4_Restk4;
                    if(TkMassCharge.size()==4){
                        if(tktkRes_mass > 0) {if (fabs(v4_Res.Mag()-tktkRes_mass) > tktkRes_mass_window) continue;}
                        DMassCutLevel[Dchannel_number-1]->Fill(0);
                        if(v4_D.Mag()<mass_window[0] || v4_D.Mag()>mass_window[1]) continue;
                        if(v4_D.Pt() < dPtCut_[Dchannel_number-1])continue;
                        DMassCutLevel[Dchannel_number-1]->Fill(1);
                        //if(fabs(v4_D.Eta()) > dEtaCut_[Dchannel_number-1])continue;
                        DMassCutLevel[Dchannel_number-1]->Fill(2);
                        selectedTkhidx.push_back(tk1_hindex);
                        selectedTkhidx.push_back(tk2_hindex);
                        selectedTkhidx.push_back(tk3_hindex);
                        selectedTkhidx.push_back(tk4_hindex);
                        selectedTkhidxSet.push_back(selectedTkhidx);
                        selectedTkhidx.clear();
                        continue;
                    }
                    //for(std::vector<pat::GenericParticle>::const_iterator tk_it5=tk_it4+1;
                    //        tk_it5 != input_tracks.end() ; tk_it5++){
                    for(int tk5idx = tk4idx+1; tk5idx < (int)isNeededTrackIdx.size(); tk5idx++){
                        //tk5_hindex = int(tk_it5 - input_tracks.begin());
                        //if(tk5_hindex>=int(isNeededTrack.size())) break;
                        //if(!isNeededTrack[tk5_hindex]) continue;
                        //if(tk_it5->charge()*TkMassCharge[4].first<0) continue;
                        tk5_hindex = isNeededTrackIdx[tk5idx];
                        if(input_tracks[tk5_hindex].charge()*TkMassCharge[4].first<0) continue;
                        if(tk5_hindex==tk1_hindex) continue;
                        if(tk5_hindex==tk2_hindex) continue;
                        if(tk5_hindex==tk3_hindex) continue;
                        if(tk5_hindex==tk4_hindex) continue;
                        v4_tk5.SetPtEtaPhiM(input_tracks[tk5_hindex].pt(),input_tracks[tk5_hindex].eta(),input_tracks[tk5_hindex].phi(),fabs(TkMassCharge[4].first));
                        v4_D = v4_tk1 + v4_tk2 + v4_tk3 + v4_tk4 + v4_tk5;
                        if(TkMassCharge[4].second==1) v4_Restk5.SetPtEtaPhiM(input_tracks[tk5_hindex].pt(),input_tracks[tk5_hindex].eta(),input_tracks[tk5_hindex].phi(),fabs(TkMassCharge[4].first));
                        else v4_Restk5.SetPtEtaPhiM(0,0,0,0);
                        v4_Res = v4_Restk1 + v4_Restk2 + v4_Restk3 + v4_Restk4 + v4_Restk5;
                        if(TkMassCharge.size()==5){
                            if(tktkRes_mass > 0) {if (fabs(v4_Res.Mag()-tktkRes_mass) > tktkRes_mass_window) continue;}
                            DMassCutLevel[Dchannel_number-1]->Fill(0);
                            if(v4_D.Mag()<mass_window[0] || v4_D.Mag()>mass_window[1]) continue;
                            if(v4_D.Pt() < dPtCut_[Dchannel_number-1])continue;
                            DMassCutLevel[Dchannel_number-1]->Fill(1);
                            //if(fabs(v4_D.Eta()) > dEtaCut_[Dchannel_number-1])continue;
                            DMassCutLevel[Dchannel_number-1]->Fill(2);
                            selectedTkhidx.push_back(tk1_hindex);
                            selectedTkhidx.push_back(tk2_hindex);
                            selectedTkhidx.push_back(tk3_hindex);
                            selectedTkhidx.push_back(tk4_hindex);
                            selectedTkhidx.push_back(tk5_hindex);
                            selectedTkhidxSet.push_back(selectedTkhidx);
                            selectedTkhidx.clear();
                            continue;
                        }
                    }
                }
            }
        }
    }
    //std::cout<<"TkCombinationPermutation, selectedTkhidxSet.size: "<<selectedTkhidxSet.size()<<std::endl;
	return;
}

void Dfinder::TkCombinationResFast(
    std::vector<pat::GenericParticle> input_tracks, 
    std::vector<int> isNeededTrackIdx,
    float *mass_window,
    std::vector< std::pair<float, int> > TkMassCharge,
    double tktkRes_mass,
    double tktkRes_mass_window,
	std::vector< std::vector<int> > &selectedTkhidxSet,
	int Dchannel_number
){
    int tk1_hindex = -1;
    int tk2_hindex = -1;
    int tk3_hindex = -1;
    int tk4_hindex = -1;
    int tk5_hindex = -1;

    std::vector<int> selectedTkhidx;
    TLorentzVector v4_D, v4_Res;//unfitted D and Res
    TLorentzVector v4_tk1, v4_tk2, v4_tk3, v4_tk4, v4_tk5;
    TLorentzVector v4_Restk1, v4_Restk2, v4_Restk3, v4_Restk4, v4_Restk5;
	for(int tk1idx = 0; tk1idx < (int)isNeededTrackIdx.size(); tk1idx++){
        v4_D.Clear(); v4_Res.Clear();
        tk1_hindex = isNeededTrackIdx[tk1idx];
        if(input_tracks[tk1_hindex].charge()*TkMassCharge[0].first<0) continue;
        v4_tk1.SetPtEtaPhiM(input_tracks[tk1_hindex].pt(),input_tracks[tk1_hindex].eta(),input_tracks[tk1_hindex].phi(),fabs(TkMassCharge[0].first));
        v4_D = v4_tk1;
        if(TkMassCharge[0].second==1) v4_Restk1.SetPtEtaPhiM(input_tracks[tk1_hindex].pt(),input_tracks[tk1_hindex].eta(),input_tracks[tk1_hindex].phi(),fabs(TkMassCharge[0].first));
        else v4_Restk1.SetPtEtaPhiM(0,0,0,0);
        v4_Res = v4_Restk1;

        for(int tk2idx = ((TkMassCharge[1]==TkMassCharge[0]) ? tk1idx+1:0); tk2idx < (int)isNeededTrackIdx.size(); tk2idx++){
            tk2_hindex = isNeededTrackIdx[tk2idx];
            if(input_tracks[tk2_hindex].charge()*TkMassCharge[1].first<0) continue;
            if(tk2_hindex==tk1_hindex) continue;
            v4_tk2.SetPtEtaPhiM(input_tracks[tk2_hindex].pt(),input_tracks[tk2_hindex].eta(),input_tracks[tk2_hindex].phi(),fabs(TkMassCharge[1].first));
            v4_D = v4_tk1 + v4_tk2;
            if(TkMassCharge[1].second==1) v4_Restk2.SetPtEtaPhiM(input_tracks[tk2_hindex].pt(),input_tracks[tk2_hindex].eta(),input_tracks[tk2_hindex].phi(),fabs(TkMassCharge[1].first));
            else v4_Restk2.SetPtEtaPhiM(0,0,0,0);
            v4_Res = v4_Restk1 + v4_Restk2;
            if(TkMassCharge.size()>2)
            if(TkMassCharge[1].second == 1 && TkMassCharge[2].second == 0 ){
                if(tktkRes_mass > 0) {if (fabs(v4_Res.Mag()-tktkRes_mass) > tktkRes_mass_window) continue;}
                DMassCutLevel[Dchannel_number-1]->Fill(0);
            }
            if(TkMassCharge.size()==2){
                //cut mass window before fit
                if(v4_D.Mag()<mass_window[0] || v4_D.Mag()>mass_window[1]) continue;
                if(v4_D.Pt() < dPtCut_[Dchannel_number-1])continue;
                DMassCutLevel[Dchannel_number-1]->Fill(1);
                //if(fabs(v4_D.Eta()) > dEtaCut_[Dchannel_number-1])continue;
                DMassCutLevel[Dchannel_number-1]->Fill(2);
                selectedTkhidx.push_back(tk1_hindex);
                selectedTkhidx.push_back(tk2_hindex);
                selectedTkhidxSet.push_back(selectedTkhidx);
                selectedTkhidx.clear();
                continue;
            }
            for(int tk3idx = ((TkMassCharge[2]==TkMassCharge[1]) ? tk2idx+1:0); tk3idx < (int)isNeededTrackIdx.size(); tk3idx++){
                tk3_hindex = isNeededTrackIdx[tk3idx];
                if(input_tracks[tk3_hindex].charge()*TkMassCharge[2].first<0) continue;
                if(tk3_hindex==tk1_hindex) continue;
                if(tk3_hindex==tk2_hindex) continue;
                v4_tk3.SetPtEtaPhiM(input_tracks[tk3_hindex].pt(),input_tracks[tk3_hindex].eta(),input_tracks[tk3_hindex].phi(),fabs(TkMassCharge[2].first));
                v4_D = v4_tk1 + v4_tk2 + v4_tk3;
                if(TkMassCharge[2].second==1) v4_Restk3.SetPtEtaPhiM(input_tracks[tk3_hindex].pt(),input_tracks[tk3_hindex].eta(),input_tracks[tk3_hindex].phi(),fabs(TkMassCharge[2].first));
                else v4_Restk3.SetPtEtaPhiM(0,0,0,0);
                v4_Res = v4_Restk1 + v4_Restk2 + v4_Restk3;
                if(TkMassCharge.size()>3)
                if(TkMassCharge[2].second == 1 && TkMassCharge[3].second == 0 ){
                    if(tktkRes_mass > 0) {if (fabs(v4_Res.Mag()-tktkRes_mass) > tktkRes_mass_window) continue;}
                    DMassCutLevel[Dchannel_number-1]->Fill(0);
                }
                if(TkMassCharge.size()==3){
                    if(v4_D.Mag()<mass_window[0] || v4_D.Mag()>mass_window[1]) continue;
                    if(v4_D.Pt() < dPtCut_[Dchannel_number-1])continue;
                    DMassCutLevel[Dchannel_number-1]->Fill(1);
                    //if(fabs(v4_D.Eta()) > dEtaCut_[Dchannel_number-1])continue;
                    DMassCutLevel[Dchannel_number-1]->Fill(2);
                    selectedTkhidx.push_back(tk1_hindex);
                    selectedTkhidx.push_back(tk2_hindex);
                    selectedTkhidx.push_back(tk3_hindex);
                    selectedTkhidxSet.push_back(selectedTkhidx);
                    selectedTkhidx.clear();
                    continue;
                }
                for(int tk4idx = ((TkMassCharge[3]==TkMassCharge[2]) ? tk3idx+1:0); tk4idx < (int)isNeededTrackIdx.size(); tk4idx++){
                    tk4_hindex = isNeededTrackIdx[tk4idx];
                    if(input_tracks[tk4_hindex].charge()*TkMassCharge[3].first<0) continue;
                    if(tk4_hindex==tk1_hindex) continue;
                    if(tk4_hindex==tk2_hindex) continue;
                    if(tk4_hindex==tk3_hindex) continue;
                    v4_tk4.SetPtEtaPhiM(input_tracks[tk4_hindex].pt(),input_tracks[tk4_hindex].eta(),input_tracks[tk4_hindex].phi(),fabs(TkMassCharge[3].first));
                    v4_D = v4_tk1 + v4_tk2 + v4_tk3 + v4_tk4;
                    if(TkMassCharge[3].second==1) v4_Restk4.SetPtEtaPhiM(input_tracks[tk4_hindex].pt(),input_tracks[tk4_hindex].eta(),input_tracks[tk4_hindex].phi(),fabs(TkMassCharge[3].first));
                    else v4_Restk4.SetPtEtaPhiM(0,0,0,0);
                    v4_Res = v4_Restk1 + v4_Restk2 + v4_Restk3 + v4_Restk4;
                    if(TkMassCharge.size()>4)
                    if(TkMassCharge[3].second == 1 && TkMassCharge[4].second == 0 ){
                        if(tktkRes_mass > 0) {if (fabs(v4_Res.Mag()-tktkRes_mass) > tktkRes_mass_window) continue;}
                        DMassCutLevel[Dchannel_number-1]->Fill(0);
                    }
                    if(TkMassCharge.size()==4){
                        if(v4_D.Mag()<mass_window[0] || v4_D.Mag()>mass_window[1]) continue;
                        if(v4_D.Pt() < dPtCut_[Dchannel_number-1])continue;
                        DMassCutLevel[Dchannel_number-1]->Fill(1);
                        //if(fabs(v4_D.Eta()) > dEtaCut_[Dchannel_number-1])continue;
                        DMassCutLevel[Dchannel_number-1]->Fill(2);
                        selectedTkhidx.push_back(tk1_hindex);
                        selectedTkhidx.push_back(tk2_hindex);
                        selectedTkhidx.push_back(tk3_hindex);
                        selectedTkhidx.push_back(tk4_hindex);
                        selectedTkhidxSet.push_back(selectedTkhidx);
                        selectedTkhidx.clear();
                        continue;
                    }
                    for(int tk5idx = ((TkMassCharge[4]==TkMassCharge[3]) ? tk4idx+1:0); tk5idx < (int)isNeededTrackIdx.size(); tk5idx++){
                        tk5_hindex = isNeededTrackIdx[tk5idx];
                        if(input_tracks[tk5_hindex].charge()*TkMassCharge[4].first<0) continue;
                        if(tk5_hindex==tk1_hindex) continue;
                        if(tk5_hindex==tk2_hindex) continue;
                        if(tk5_hindex==tk3_hindex) continue;
                        if(tk5_hindex==tk4_hindex) continue;
                        v4_tk5.SetPtEtaPhiM(input_tracks[tk5_hindex].pt(),input_tracks[tk5_hindex].eta(),input_tracks[tk5_hindex].phi(),fabs(TkMassCharge[4].first));
                        v4_D = v4_tk1 + v4_tk2 + v4_tk3 + v4_tk4 + v4_tk5;
                        if(TkMassCharge[4].second==1) v4_Restk5.SetPtEtaPhiM(input_tracks[tk5_hindex].pt(),input_tracks[tk5_hindex].eta(),input_tracks[tk5_hindex].phi(),fabs(TkMassCharge[4].first));
                        else v4_Restk5.SetPtEtaPhiM(0,0,0,0);
                        v4_Res = v4_Restk1 + v4_Restk2 + v4_Restk3 + v4_Restk4 + v4_Restk5;
                        if(TkMassCharge.size()==5){
                            if(tktkRes_mass > 0) {if (fabs(v4_Res.Mag()-tktkRes_mass) > tktkRes_mass_window) continue;}
                            DMassCutLevel[Dchannel_number-1]->Fill(0);
                            if(v4_D.Mag()<mass_window[0] || v4_D.Mag()>mass_window[1]) continue;
                            if(v4_D.Pt() < dPtCut_[Dchannel_number-1])continue;
                            DMassCutLevel[Dchannel_number-1]->Fill(1);
                            //if(fabs(v4_D.Eta()) > dEtaCut_[Dchannel_number-1])continue;
                            DMassCutLevel[Dchannel_number-1]->Fill(2);
                            selectedTkhidx.push_back(tk1_hindex);
                            selectedTkhidx.push_back(tk2_hindex);
                            selectedTkhidx.push_back(tk3_hindex);
                            selectedTkhidx.push_back(tk4_hindex);
                            selectedTkhidx.push_back(tk5_hindex);
                            selectedTkhidxSet.push_back(selectedTkhidx);
                            selectedTkhidx.clear();
                            continue;
                        }
                    }
                }
            }
        }
    }
    //std::cout<<"TkCombinationResFast, selectedTkhidxSet.size: "<<selectedTkhidxSet.size()<<std::endl;
	return;
}
//}}}

//BranchOutNTk{{{
void Dfinder::BranchOutNTk(//input 2~4 tracks
    DInfoBranches &DInfo, 
    std::vector<pat::GenericParticle> input_tracks, 
    reco::Vertex thePrimaryV,
    std::vector<int> isNeededTrackIdx,
    std::vector<int> &D_counter,
    float *mass_window,
    std::vector< std::pair<float, int> > TkMassCharge,
    double tktkRes_mass,
    double tktkRes_mass_window,
    bool doConstrainFit,
    bool SequentialFit,
    int Dchannel_number,
	int TkCombinationMethod
){
    if(Dchannel_number > (int)Dchannel_.size()){ printf("Exceeding defined # of channel, exit"); return;}

	std::vector< std::vector<int> > selectedTkhidxSet;
	if( TkCombinationMethod == 0 ) 
        TkCombinationPermutation( input_tracks, isNeededTrackIdx, mass_window, TkMassCharge, tktkRes_mass, tktkRes_mass_window, selectedTkhidxSet, Dchannel_number);
    else if( TkCombinationMethod == 1 )
        TkCombinationResFast( input_tracks, isNeededTrackIdx, mass_window, TkMassCharge, tktkRes_mass, tktkRes_mass_window, selectedTkhidxSet, Dchannel_number);
    else
		{ printf("unknown method of track combination, exit"); return;}
//    std::cout<<"selectedTkhidxSet.size: "<<selectedTkhidxSet.size()<<std::endl;
	
	float chi = 0.;
	float ndf = 0.;

    //particle factory: produce transient tracks
    KinematicParticleFactoryFromTransientTrack pFactory;
    VirtualKinematicParticleFactory vFactory;
    //fitter for D
    KinematicParticleVertexFitter   tktk_fitter;
    RefCountedKinematicTree         tktk_VFT;
    RefCountedKinematicParticle     tktk_VFP;
    RefCountedKinematicVertex       tktk_VFPvtx;
    //constrain fit fitter
    KinematicConstrainedVertexFitter kcv_tktk_fitter;
    //fitter for Res
    KinematicParticleVertexFitter   tktkRes_fitter;
    RefCountedKinematicTree         tktkRes_VFT;
    RefCountedKinematicParticle     tktkRes_VFP;
    RefCountedKinematicVertex       tktkRes_VFPvtx;

    TLorentzVector v4_tk;
    std::vector<TLorentzVector> tktk_4vecs;//fitted tks
    TLorentzVector tktk_4vec;//fitted D
    std::vector<TLorentzVector> tktkRes_4vecs;//fitted Res tks
    TLorentzVector tktkRes_4vec;//fitted Res
    std::vector<RefCountedKinematicParticle> tktk_candidate;//input tracks to D fitter
    std::vector<RefCountedKinematicParticle> tktkRes_candidate;//input tracks to Res fitter
    std::vector<RefCountedKinematicParticle> tktkCands;//output tracks from D fitter
    std::vector<RefCountedKinematicParticle> tktkResCands;//output tracks from Res fitter

    for(int i = 0; i < int(selectedTkhidxSet.size()); i++){
        if (DInfo.size >= MAX_XB) break;

        //clear before using
        v4_tk.Clear();
        tktk_4vecs.clear();
        tktk_4vec.Clear();
        tktkRes_4vecs.clear();
        tktkRes_4vec.Clear();
        tktk_candidate.clear();
        tktkRes_candidate.clear();
        tktkCands.clear();
        tktkResCands.clear();

        //push back the Res tracks as first tracks
        ParticleMass tk_mass;
        //keep track of the push_back track index
        std::vector<int> pushbackTrkIdx;
        std::vector<int> pushbackResTrkIdx;
        float tk_sigma;
        for(int p = 0; p < int(selectedTkhidxSet[0].size()); p++){        
            if(TkMassCharge[p].second==0) continue;
            reco::TransientTrack tkTT(input_tracks[selectedTkhidxSet[i][p]].track(), &(*bField) );
            if (!tkTT.isValid()) continue;
            tk_mass = fabs(TkMassCharge[p].first);
            tk_sigma = Functs.getParticleSigma(tk_mass);
            if(!SequentialFit){
                tktk_candidate.push_back(pFactory.particle(tkTT,tk_mass,chi,ndf,tk_sigma));
                pushbackTrkIdx.push_back(selectedTkhidxSet[i][p]);
            }
            if(tktkRes_mass>0){
                tktkRes_candidate.push_back(pFactory.particle(tkTT,tk_mass,chi,ndf,tk_sigma));
                pushbackResTrkIdx.push_back(selectedTkhidxSet[i][p]);
            }
        }

        //if there is tktk Res, check its fit validity also
        if(tktkRes_mass>0){
            tktkRes_VFT = tktkRes_fitter.fit(tktkRes_candidate);
            if(!tktkRes_VFT->isValid()) continue;
            DMassCutLevel[Dchannel_number-1]->Fill(3);
            tktkRes_VFP   = tktkRes_VFT->currentParticle();
            tktkRes_VFPvtx = tktkRes_VFT->currentDecayVertex();
            double chi2_prob_tktkRes = TMath::Prob(tktkRes_VFPvtx->chiSquared(),tktkRes_VFPvtx->degreesOfFreedom());

            if(chi2_prob_tktkRes < VtxChiProbCut_[Dchannel_number-1]) continue;
            DMassCutLevel[Dchannel_number-1]->Fill(4);

            if(SequentialFit){
                float tktkchi = tktkRes_VFPvtx->chiSquared();
                float tktkndf = tktkRes_VFPvtx->degreesOfFreedom();
                tktk_candidate.push_back(vFactory.particle(tktkRes_VFP->currentState(),tktkchi,tktkndf,tktkRes_VFP));
                //tktk_candidate.push_back(tktkRes_VFP);
                pushbackTrkIdx.push_back(-1);//means its a resonance particle
            }
        }
        //push back the other tracks
        for(int p = 0; p < int(selectedTkhidxSet[0].size()); p++){        
            if(TkMassCharge[p].second==1) continue;
            reco::TransientTrack tkTT(input_tracks[selectedTkhidxSet[i][p]].track(), &(*bField) );
            if (!tkTT.isValid()) continue;
            tk_mass = fabs(TkMassCharge[p].first);
            tk_sigma = Functs.getParticleSigma(tk_mass);
            tktk_candidate.push_back(pFactory.particle(tkTT,tk_mass,chi,ndf,tk_sigma));
            pushbackTrkIdx.push_back(selectedTkhidxSet[i][p]);
        }
        DMassCutLevel[Dchannel_number-1]->Fill(5);

        double MaximumDoca = Functs.getMaxDoca(tktk_candidate);
        if (MaximumDoca > MaxDocaCut_[Dchannel_number-1]) continue;
        DMassCutLevel[Dchannel_number-1]->Fill(6);

        if(tktkRes_mass>0){
            if(doConstrainFit){
                ParticleMass tktkResMass = tktkRes_mass;
                MultiTrackKinematicConstraint *tktkResConstraint = new TwoTrackMassKinematicConstraint(tktkResMass);
                tktk_VFT = kcv_tktk_fitter.fit(tktk_candidate, tktkResConstraint);
            }
            else tktk_VFT = tktk_fitter.fit(tktk_candidate);
        }
        else{
            tktk_VFT = tktk_fitter.fit(tktk_candidate);
        }

        if(!tktk_VFT->isValid()) continue;
        DMassCutLevel[Dchannel_number-1]->Fill(7);

        tktk_VFT->movePointerToTheTop();
        tktk_VFP   = tktk_VFT->currentParticle();
        tktk_VFPvtx = tktk_VFT->currentDecayVertex();
        if (!tktk_VFPvtx->vertexIsValid()) continue;
        DMassCutLevel[Dchannel_number-1]->Fill(8);

        double chi2_prob_tktk = TMath::Prob(tktk_VFPvtx->chiSquared(),tktk_VFPvtx->degreesOfFreedom());
        if(chi2_prob_tktk < VtxChiProbCut_[Dchannel_number-1]) continue;
        DMassCutLevel[Dchannel_number-1]->Fill(9);

        tktkCands  = tktk_VFT->finalStateParticles();

        //cut mass window after fit
        //if (tktk_VFP->currentState().mass()<mass_window[0] || tktk_VFP->currentState().mass()>mass_window[1]) continue;
        DMassCutLevel[Dchannel_number-1]->Fill(10);

        for(unsigned int k = 0; k < tktkCands.size(); k++){
            v4_tk.SetPxPyPzE(tktkCands[k]->currentState().kinematicParameters().momentum().x(),
                                     tktkCands[k]->currentState().kinematicParameters().momentum().y(),
                                     tktkCands[k]->currentState().kinematicParameters().momentum().z(),
                                     tktkCands[k]->currentState().kinematicParameters().energy());
            tktk_4vecs.push_back(v4_tk);
            v4_tk.Clear();
        }

        tktk_4vec.SetPxPyPzE(tktk_VFP->currentState().kinematicParameters().momentum().x(),
                tktk_VFP->currentState().kinematicParameters().momentum().y(),
                tktk_VFP->currentState().kinematicParameters().momentum().z(),
                tktk_VFP->currentState().kinematicParameters().energy());

        //tktkRes fit info
        if(tktkRes_mass>0){
            tktkResCands = tktkRes_VFT->finalStateParticles();

            for(unsigned int k = 0; k < tktkResCands.size(); k++){
                v4_tk.SetPxPyPzE(tktkResCands[k]->currentState().kinematicParameters().momentum().x(),
                                        tktkResCands[k]->currentState().kinematicParameters().momentum().y(),
                                        tktkResCands[k]->currentState().kinematicParameters().momentum().z(),
                                        tktkResCands[k]->currentState().kinematicParameters().energy());
                tktkRes_4vecs.push_back(v4_tk);
                v4_tk.Clear();
            }

            tktkRes_4vec.SetPxPyPzE(tktkRes_VFP->currentState().kinematicParameters().momentum().x(),
                    tktkRes_VFP->currentState().kinematicParameters().momentum().y(),
                    tktkRes_VFP->currentState().kinematicParameters().momentum().z(),
                    tktkRes_VFP->currentState().kinematicParameters().energy());

            DInfo.tktkRes_mass[DInfo.size]            = tktkRes_4vec.Mag();
            DInfo.tktkRes_pt[DInfo.size]              = tktkRes_4vec.Pt();
            DInfo.tktkRes_eta[DInfo.size]             = tktkRes_4vec.Eta();
            DInfo.tktkRes_phi[DInfo.size]             = tktkRes_4vec.Phi();
            DInfo.tktkRes_vtxX[DInfo.size]            = tktkRes_VFPvtx->position().x();
            DInfo.tktkRes_vtxY[DInfo.size]            = tktkRes_VFPvtx->position().y();
            DInfo.tktkRes_vtxZ[DInfo.size]            = tktkRes_VFPvtx->position().z();
            DInfo.tktkRes_vtxXErr[DInfo.size]         = tktkRes_VFPvtx->error().cxx();
            DInfo.tktkRes_vtxYErr[DInfo.size]         = tktkRes_VFPvtx->error().cyy();
            DInfo.tktkRes_vtxZErr[DInfo.size]         = tktkRes_VFPvtx->error().czz();
            DInfo.tktkRes_vtxYXErr[DInfo.size]        = tktkRes_VFPvtx->error().cyx();
            DInfo.tktkRes_vtxZXErr[DInfo.size]        = tktkRes_VFPvtx->error().czx();
            DInfo.tktkRes_vtxZYErr[DInfo.size]        = tktkRes_VFPvtx->error().czy();
            DInfo.tktkRes_vtxdof[DInfo.size]          = tktkRes_VFPvtx->degreesOfFreedom();
            DInfo.tktkRes_vtxchi2[DInfo.size]         = tktkRes_VFPvtx->chiSquared();
        
            TVector3 Res_svpvVec;
            Res_svpvVec.SetXYZ(DInfo.tktkRes_vtxX[DInfo.size]-EvtInfo.PVx, DInfo.tktkRes_vtxY[DInfo.size]-EvtInfo.PVy, DInfo.tktkRes_vtxZ[DInfo.size]-EvtInfo.PVz);
            TVector3 Res_dVec;
            Res_dVec.SetPtEtaPhi(DInfo.tktkRes_pt[DInfo.size], DInfo.tktkRes_eta[DInfo.size], DInfo.tktkRes_phi[DInfo.size]);
            DInfo.tktkRes_alpha[DInfo.size] = Res_svpvVec.Angle(Res_dVec);
            if( DInfo.tktkRes_alpha[DInfo.size] > tktkRes_alphaCut_[Dchannel_number-1]) continue;

            VertexDistance3D Res_a3d;
            DInfo.tktkRes_svpvDistance[DInfo.size] = Res_a3d.distance(thePrimaryV,tktkRes_VFPvtx->vertexState()).value();
            DInfo.tktkRes_svpvDisErr[DInfo.size] = Res_a3d.distance(thePrimaryV,tktkRes_VFPvtx->vertexState()).error();
            if(DInfo.tktkRes_pt[DInfo.size] <= dCutSeparating_PtVal_[Dchannel_number-1] && (DInfo.tktkRes_svpvDistance[DInfo.size]/DInfo.tktkRes_svpvDisErr[DInfo.size]) < tktkRes_svpvDistanceCut_lowptD_[Dchannel_number-1]) continue;
            else if( DInfo.tktkRes_pt[DInfo.size] > dCutSeparating_PtVal_[Dchannel_number-1] && (DInfo.tktkRes_svpvDistance[DInfo.size]/DInfo.tktkRes_svpvDisErr[DInfo.size]) < tktkRes_svpvDistanceCut_highptD_[Dchannel_number-1]) continue;
            DMassCutLevel[Dchannel_number-1]->Fill(11);

            //index initialization to -2
            DInfo.tktkRes_rftk1_index[DInfo.size]     = -2;
            DInfo.tktkRes_rftk2_index[DInfo.size]     = -2;
            DInfo.tktkRes_rftk3_index[DInfo.size]     = -2;
            DInfo.tktkRes_rftk4_index[DInfo.size]     = -2;

            DInfo.tktkRes_rftk1_mass[DInfo.size]      = tktkRes_4vecs[0].Mag();
            DInfo.tktkRes_rftk1_pt[DInfo.size]        = tktkRes_4vecs[0].Pt();
            DInfo.tktkRes_rftk1_eta[DInfo.size]       = tktkRes_4vecs[0].Eta();
            DInfo.tktkRes_rftk1_phi[DInfo.size]       = tktkRes_4vecs[0].Phi();
            DInfo.tktkRes_rftk1_index[DInfo.size]     = pushbackResTrkIdx[0];
            DInfo.tktkRes_rftk2_mass[DInfo.size]      = tktkRes_4vecs[1].Mag();
            DInfo.tktkRes_rftk2_pt[DInfo.size]        = tktkRes_4vecs[1].Pt();
            DInfo.tktkRes_rftk2_eta[DInfo.size]       = tktkRes_4vecs[1].Eta();
            DInfo.tktkRes_rftk2_phi[DInfo.size]       = tktkRes_4vecs[1].Phi();
            DInfo.tktkRes_rftk2_index[DInfo.size]     = pushbackResTrkIdx[1];
            if(tktkResCands.size()>2){
                DInfo.tktkRes_rftk3_mass[DInfo.size]      = tktkRes_4vecs[2].Mag();
                DInfo.tktkRes_rftk3_pt[DInfo.size]        = tktkRes_4vecs[2].Pt();
                DInfo.tktkRes_rftk3_eta[DInfo.size]       = tktkRes_4vecs[2].Eta();
                DInfo.tktkRes_rftk3_phi[DInfo.size]       = tktkRes_4vecs[2].Phi();
                DInfo.tktkRes_rftk3_index[DInfo.size]     = pushbackResTrkIdx[2];
            }
            if(tktkResCands.size()>3){
                DInfo.tktkRes_rftk4_mass[DInfo.size]      = tktkRes_4vecs[3].Mag();
                DInfo.tktkRes_rftk4_pt[DInfo.size]        = tktkRes_4vecs[3].Pt();
                DInfo.tktkRes_rftk4_eta[DInfo.size]       = tktkRes_4vecs[3].Eta();
                DInfo.tktkRes_rftk4_phi[DInfo.size]       = tktkRes_4vecs[3].Phi();
                DInfo.tktkRes_rftk4_index[DInfo.size]     = pushbackResTrkIdx[3];
            }
        }

        //fit info
        DInfo.index[DInfo.size]           = DInfo.size;
        DInfo.mass[DInfo.size]            = tktk_4vec.Mag();
        DInfo.pt[DInfo.size]              = tktk_4vec.Pt();
        DInfo.eta[DInfo.size]             = tktk_4vec.Eta();
        DInfo.phi[DInfo.size]             = tktk_4vec.Phi();
        DInfo.px[DInfo.size]              = tktk_4vec.Px();
        DInfo.py[DInfo.size]              = tktk_4vec.Py();
        DInfo.pz[DInfo.size]              = tktk_4vec.Pz();
        DInfo.MaxDoca[DInfo.size]         = MaximumDoca;

        VertexDistance3D a3d;
        //https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_0/RecoVertex/VertexTools/src/VertexDistance3D.cc
        DInfo.svpvDistance[DInfo.size] = a3d.distance(thePrimaryV,tktk_VFPvtx->vertexState()).value();
        DInfo.svpvDisErr[DInfo.size] = a3d.distance(thePrimaryV,tktk_VFPvtx->vertexState()).error();
        if( DInfo.pt[DInfo.size] <= dCutSeparating_PtVal_[Dchannel_number-1] && (DInfo.svpvDistance[DInfo.size]/DInfo.svpvDisErr[DInfo.size]) < svpvDistanceCut_lowptD_[Dchannel_number-1]) continue;
        else if( DInfo.pt[DInfo.size] > dCutSeparating_PtVal_[Dchannel_number-1] && (DInfo.svpvDistance[DInfo.size]/DInfo.svpvDisErr[DInfo.size]) < svpvDistanceCut_highptD_[Dchannel_number-1]) continue;
        DMassCutLevel[Dchannel_number-1]->Fill(12);

        reco::Vertex::Point vp1(thePrimaryV.position().x(), thePrimaryV.position().y(), 0.);
        reco::Vertex::Point vp2(tktk_VFPvtx->vertexState().position().x(), tktk_VFPvtx->vertexState().position().y(), 0.);
        ROOT::Math::SVector<double, 6> sv1(thePrimaryV.covariance(0,0), thePrimaryV.covariance(0,1), thePrimaryV.covariance(1,1), 0., 0., 0.);
        ROOT::Math::SVector<double, 6> sv2(tktk_VFPvtx->vertexState().error().cxx(), tktk_VFPvtx->vertexState().error().cyx(), tktk_VFPvtx->vertexState().error().cyy(), 0., 0., 0.);
        reco::Vertex::Error ve1(sv1);
        reco::Vertex::Error ve2(sv2);
        reco::Vertex v1(vp1, ve1);
        reco::Vertex v2(vp2, ve2);
        DInfo.svpvDistance_2D[DInfo.size] = a3d.distance(v1, v2).value();
        DInfo.svpvDisErr_2D[DInfo.size] = a3d.distance(v1, v2).error();

        DInfo.vtxX[DInfo.size]            = tktk_VFPvtx->position().x();
        DInfo.vtxY[DInfo.size]            = tktk_VFPvtx->position().y();
        DInfo.vtxZ[DInfo.size]            = tktk_VFPvtx->position().z();
        DInfo.vtxXErr[DInfo.size]         = tktk_VFPvtx->error().cxx();
        DInfo.vtxYErr[DInfo.size]         = tktk_VFPvtx->error().cyy();
        DInfo.vtxZErr[DInfo.size]         = tktk_VFPvtx->error().czz();
        DInfo.vtxYXErr[DInfo.size]        = tktk_VFPvtx->error().cyx();
        DInfo.vtxZXErr[DInfo.size]        = tktk_VFPvtx->error().czx();
        DInfo.vtxZYErr[DInfo.size]        = tktk_VFPvtx->error().czy();
        DInfo.vtxdof[DInfo.size]          = tktk_VFPvtx->degreesOfFreedom();
        DInfo.vtxchi2[DInfo.size]         = tktk_VFPvtx->chiSquared();
        
        TVector3 svpvVec;
        svpvVec.SetXYZ(DInfo.vtxX[DInfo.size]-EvtInfo.PVx, DInfo.vtxY[DInfo.size]-EvtInfo.PVy, DInfo.vtxZ[DInfo.size]-EvtInfo.PVz);
        TVector3 dVec;
        dVec.SetXYZ(DInfo.px[DInfo.size], DInfo.py[DInfo.size], DInfo.pz[DInfo.size]);
        DInfo.alpha[DInfo.size] = svpvVec.Angle(dVec);
        if( DInfo.alpha[DInfo.size] > alphaCut_[Dchannel_number-1]) continue;
        DMassCutLevel[Dchannel_number-1]->Fill(13);

        DInfo.rftk1_mass[DInfo.size]      = tktk_4vecs[0].Mag();
        DInfo.rftk1_pt[DInfo.size]        = tktk_4vecs[0].Pt();
        DInfo.rftk1_eta[DInfo.size]       = tktk_4vecs[0].Eta();
        DInfo.rftk1_phi[DInfo.size]       = tktk_4vecs[0].Phi();
        DInfo.rftk2_mass[DInfo.size]      = tktk_4vecs[1].Mag();
        DInfo.rftk2_pt[DInfo.size]        = tktk_4vecs[1].Pt();
        DInfo.rftk2_eta[DInfo.size]       = tktk_4vecs[1].Eta();
        DInfo.rftk2_phi[DInfo.size]       = tktk_4vecs[1].Phi();

        //Index initialize to -2
        DInfo.rftk1_index[DInfo.size] = -2;
        DInfo.rftk2_index[DInfo.size] = -2;
        DInfo.rftk3_index[DInfo.size] = -2;
        DInfo.rftk4_index[DInfo.size] = -2;
        DInfo.rftk5_index[DInfo.size] = -2;

        DInfo.rftk1_index[DInfo.size]     = pushbackTrkIdx[0];
        DInfo.rftk2_index[DInfo.size]     = pushbackTrkIdx[1];

        //document the mass hypothesis
        if( fabs(tktk_4vecs[0].Mag()-PION_MASS) < fabs(tktk_4vecs[0].Mag()-KAON_MASS) ) DInfo.rftk1_MassHypo[DInfo.size] = 211;
        else DInfo.rftk1_MassHypo[DInfo.size] = 321;
        //If its a Res particle, save it as D0
        if( DInfo.rftk1_index[DInfo.size] == -1) DInfo.rftk1_MassHypo[DInfo.size] = 421;

        if( fabs(tktk_4vecs[1].Mag()-PION_MASS) < fabs(tktk_4vecs[1].Mag()-KAON_MASS) ) DInfo.rftk2_MassHypo[DInfo.size] = 211;
        else DInfo.rftk2_MassHypo[DInfo.size] = 321;

        if(tktkCands.size()>2){
            DInfo.rftk3_mass[DInfo.size]  = tktk_4vecs[2].Mag();
            DInfo.rftk3_pt[DInfo.size]    = tktk_4vecs[2].Pt();
            DInfo.rftk3_eta[DInfo.size]   = tktk_4vecs[2].Eta();
            DInfo.rftk3_phi[DInfo.size]   = tktk_4vecs[2].Phi();
            DInfo.rftk3_index[DInfo.size] = pushbackTrkIdx[2];
            if( fabs(tktk_4vecs[2].Mag()-PION_MASS) < fabs(tktk_4vecs[2].Mag()-KAON_MASS) ) DInfo.rftk3_MassHypo[DInfo.size] = 211;
            else DInfo.rftk3_MassHypo[DInfo.size] = 321;
        }
        if(tktkCands.size()>3){
            DInfo.rftk4_mass[DInfo.size]  = tktk_4vecs[3].Mag();
            DInfo.rftk4_pt[DInfo.size]    = tktk_4vecs[3].Pt();
            DInfo.rftk4_eta[DInfo.size]   = tktk_4vecs[3].Eta();
            DInfo.rftk4_phi[DInfo.size]   = tktk_4vecs[3].Phi();
            DInfo.rftk4_index[DInfo.size] = pushbackTrkIdx[3];
            if( fabs(tktk_4vecs[3].Mag()-PION_MASS) < fabs(tktk_4vecs[3].Mag()-KAON_MASS) ) DInfo.rftk4_MassHypo[DInfo.size] = 211;
            else DInfo.rftk4_MassHypo[DInfo.size] = 321;
        }
        if(tktkCands.size()>4){
            DInfo.rftk5_mass[DInfo.size]  = tktk_4vecs[4].Mag();
            DInfo.rftk5_pt[DInfo.size]    = tktk_4vecs[4].Pt();
            DInfo.rftk5_eta[DInfo.size]   = tktk_4vecs[4].Eta();
            DInfo.rftk5_phi[DInfo.size]   = tktk_4vecs[4].Phi();
            DInfo.rftk5_index[DInfo.size] = pushbackTrkIdx[4];
            if( fabs(tktk_4vecs[4].Mag()-PION_MASS) < fabs(tktk_4vecs[4].Mag()-KAON_MASS) ) DInfo.rftk5_MassHypo[DInfo.size] = 211;
            else DInfo.rftk5_MassHypo[DInfo.size] = 321;
        }

        DInfo.type[DInfo.size] = Dchannel_number;
        D_counter[Dchannel_number-1]++;

        tktk_candidate.clear();
        tktkCands.clear();
        DInfo.size++;
    }
}
//}}}

//define this as a plug-in
DEFINE_FWK_MODULE(Dfinder);
