// vim:set ts=4 sw=4 fdm=marker et:
// Ntuplt creator for B meson related analysis.
// Maintain and contact: ta-wei wang
// Email: "tawei@mit.edu" or "ta-wei.wang@cern.ch"
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"

        //message logger
#include "FWCore/MessageLogger/interface/MessageLogger.h"
        //magnetic field
#include "MagneticField/Engine/interface/MagneticField.h"
        //vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/Ref.h"
        //For Kalman vertex fitters
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
        //KalmanTrimmedVertexFinder
#include "RecoVertex/TrimmedKalmanVertexFinder/interface/KalmanTrimmedVertexFinder.h"
        //ROOT
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
        //others
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <iostream>
        //Bfinder
#include "Bfinder/Bfinder/interface/format.h"
//#include "Bfinder/Bfinder/interface/TriggerBooking.h"

#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TBranch.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <math.h>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <list>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#define MUON_MASS   0.10565837
#define PION_MASS   0.13957018
#define KAON_MASS   0.493677
#define KSHORT_MASS 0.497614
#define KSTAR_MASS  0.89594
#define PHI_MASS    1.019455
#define JPSI_MASS   3.096916
#define PSI2S_MASS  3.686109

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
        
        virtual bool GetAncestor(const reco::Candidate* p, int PDGprefix);
        virtual void BranchOut2MuTk(
            BInfoBranches &BInfo,
            std::vector<pat::GenericParticle> input_tracks,
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
        virtual void BranchOutNTk(
            DInfoBranches &DInfo, 
            std::vector<pat::GenericParticle> input_tracks, 
            std::vector<bool> isNeededTrack,
            std::vector<int> &D_counter,
            float *mass_window,
            std::vector<double> TkMass,
            std::vector<int> TkCharge,
            double tktkRes_mass,
            double tktkRes_mass_window,
            double *ResIndex,
            bool doConstrainFit,
            int Dchannel_number
        );
 
        // ----------member data ---------------------------
        edm::ESHandle<MagneticField> bField;
        edm::ParameterSet theConfig;
//      std::vector<std::string> TriggersForMatching_;
        std::vector<int> Bchannel_;
        std::vector<int> Dchannel_;
        std::vector<std::string> MuonTriggerMatchingPath_;
//        edm::InputTag hltLabel_;
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
        double dPtCut_;
        double bEtaCut_;
        double dEtaCut_;
        bool RunOnMC_;
        bool doTkPreCut_;
        bool doMuPreCut_;

        edm::Service<TFileService> fs;
        TTree *root;
        EvtInfoBranches     EvtInfo;
        VtxInfoBranches     VtxInfo;
        MuonInfoBranches    MuonInfo;
        TrackInfoBranches   TrackInfo;
        BInfoBranches       BInfo;
        DInfoBranches       DInfo;
        GenInfoBranches     GenInfo;

        //histograms
        TH1F *MuonCutLevel;
        TH1F *TrackCutLevel;
        TH1F *XbujCutLevel;
        //How many channel
        static int const Nchannel = 20;
        std::vector<TH1F*> XbMassCutLevel;
        std::vector<TH1F*> DMassCutLevel;
        
};//}}}

void Bfinder::beginJob()
{//{{{
    root = fs->make<TTree>("root","root");
    EvtInfo.regTree(root);
    VtxInfo.regTree(root);
    MuonInfo.regTree(root);
    TrackInfo.regTree(root);
    BInfo.regTree(root);
    DInfo.regTree(root);
    GenInfo.regTree(root);
}//}}}

Bfinder::Bfinder(const edm::ParameterSet& iConfig):theConfig(iConfig)
{//{{{
    //now do what ever initialization is needed
//  TriggersForMatching_= iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching");
    Bchannel_= iConfig.getParameter<std::vector<int> >("Bchannel");
    Dchannel_= iConfig.getParameter<std::vector<int> >("Dchannel");
    MuonTriggerMatchingPath_ = iConfig.getParameter<std::vector<std::string> >("MuonTriggerMatchingPath");
    genLabel_           = iConfig.getParameter<edm::InputTag>("GenLabel");
    trackLabel_         = iConfig.getParameter<edm::InputTag>("TrackLabel");
    muonLabel_          = iConfig.getParameter<edm::InputTag>("MuonLabel");
//    hltLabel_           = iConfig.getParameter<edm::InputTag>("HLTLabel");
    puInfoLabel_        = iConfig.getParameter<edm::InputTag>("PUInfoLabel");
    bsLabel_        = iConfig.getParameter<edm::InputTag>("BSLabel");
    pvLabel_        = iConfig.getParameter<edm::InputTag>("PVLabel");

    doTkPreCut_ = iConfig.getParameter<bool>("doTkPreCut");
    doMuPreCut_ = iConfig.getParameter<bool>("doMuPreCut");
    tkPtCut_ = iConfig.getParameter<double>("tkPtCut");
    tkEtaCut_ = iConfig.getParameter<double>("tkEtaCut");
    jpsiPtCut_ = iConfig.getParameter<double>("jpsiPtCut");
    bPtCut_ = iConfig.getParameter<double>("bPtCut");
    dPtCut_ = iConfig.getParameter<double>("dPtCut");
    bEtaCut_ = iConfig.getParameter<double>("bEtaCut");
    dEtaCut_ = iConfig.getParameter<double>("dEtaCut");
    RunOnMC_ = iConfig.getParameter<bool>("RunOnMC");

    MuonCutLevel        = fs->make<TH1F>("MuonCutLevel"     , "MuonCutLevel"    , 10, 0, 10);
    TrackCutLevel       = fs->make<TH1F>("TrackCutLevel"    , "TrackCutLevel"   , 10, 0, 10);
    XbujCutLevel        = fs->make<TH1F>("XbujCutLevel"     , "XbujCutLevel"    , 10, 0, 10);
    for(unsigned int i = 0; i < Bchannel_.size(); i++){
        TH1F* XbMassCutLevel_temp      = fs->make<TH1F>(TString::Format("XbMassCutLevel_i")   ,TString::Format("XbMassCutLevel_i")  , 10, 0, 10);
        XbMassCutLevel.push_back(XbMassCutLevel_temp);
    }
    for(unsigned int i = 0; i < Dchannel_.size(); i++){
        TH1F* DMassCutLevel_temp      = fs->make<TH1F>(TString::Format("DMassCutLevel_i")   ,TString::Format("DMassCutLevel_i")  , 10, 0, 10);
        DMassCutLevel.push_back(DMassCutLevel_temp);
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
    iEvent.getByLabel(bsLabel_, beamSpotHandle);
    if (beamSpotHandle.isValid()){
        beamSpot = *beamSpotHandle;
        theBeamSpotV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }else{
        std::cout<< "No beam spot available from EventSetup \n";
    }

    //get vertex informationa
    edm::Handle<reco::VertexCollection> VertexHandle;
    //iEvent.getByLabel("offlinePrimaryVertexHandle", VertexHandle);
    //iEvent.getByLabel("offlinePrimaryVerticesWithBS", VertexHandle);
    iEvent.getByLabel(pvLabel_, VertexHandle);

    /*  
    if (!VertexHandle.failedToGet() && VertexHandle->size()>0){
        //int nVtxTrks = 0;//outdated PV definition
        double max_tkSt = 0;
        for(std::vector<reco::Vertex>::const_iterator it_vtx = VertexHandle->begin(); it_vtx != VertexHandle->end(); it_vtx++){
            if (!it_vtx->isValid()) continue;
            //find primary primary vertex with largest St
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
    DInfo.size      = 0;
    GenInfo.size    = 0;
    
    std::vector<int> B_counter;
    for(unsigned int i = 0; i < Bchannel_.size(); i++){
        B_counter.push_back(0);
    }
    std::vector<int> D_counter;
    for(unsigned int i = 0; i < Dchannel_.size(); i++){
        D_counter.push_back(0);
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
        if (0){
        //if (input_muons.size() == 0){
            std::cout << "There's no muon : " << iEvent.id() << std::endl;
        }else{
            std::cout << "Got " << input_muons.size() << " muons / ";
            if (0){
            //if (input_tracks.size() == 0){
                std::cout << "There's no track: " << iEvent.id() << std::endl;
            }else{
                std::cout << "Got " << input_tracks.size() << " tracks" << std::endl;
                //if (input_tracks.size() > 0 && input_muons.size() > 1){
                if (input_tracks.size() > 0){

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

                        //1.every digit after this must be 1 since only binary of the form 000111 will have: (b&b+1)==0 
                        //0.pass all cuts
                        MuonInfo.isGoodCand[MuonInfo.size] += 1 << 1;
                        if (((MuonInfo.isGoodCand[MuonInfo.size]>>1)&((MuonInfo.isGoodCand[MuonInfo.size]>>1)+1)) == 0)
                            MuonInfo.isGoodCand[MuonInfo.size] += 1;

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
                            if (!(tk_it->track()->qualityByName("highPurity")))        continue;
                            TrackCutLevel->Fill(5);
                            //outdated selections
                            //if (tk_it->track()->normalizedChi2()>5)             continue;
                            //if (tk_it->p()>200 || tk_it->pt()>200)              continue;
                            //if (tk_it->track()->hitPattern().numberOfValidStripHits()<10)continue;
                            //if (tk_it->track()->hitPattern().numberOfValidPixelHits()<2) continue;
                        }
                        isNeededTrack[tk_it-input_tracks.begin()] = true;
                        PassedTrk++;
                    }//end of track preselection}}}
                    //printf("-----*****DEBUG:End of track preselection.\n");
                    
                    // DInfo section{{{
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+pi-
                    //////////////////////////////////////////////////////////////////////////
                    float mass_window[2] = {1.7,2.0};
                    double ResIndex[2] = {0, 1};
                    std::vector<double> TkMass;
                    std::vector<int> TkCharge;
                    if(Dchannel_[0] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 1);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 1);
                        TkMass.clear();
                        TkCharge.clear();
                    }

                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K-pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[1] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 2);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 2);
                        TkMass.clear();
                        TkCharge.clear();
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K-pi+pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[2] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 3);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 3);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 3);
                        TkMass.clear();
                        TkCharge.clear();
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+pi-pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[3] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 4);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 4);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 4);
                        TkMass.clear();
                        TkCharge.clear();
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K-pi-pi+pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[4] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+pi+pi-pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[5] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+K-(Phi)pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[6] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        ResIndex[0] = 0; ResIndex[1] = 1;
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 7);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 7);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        ResIndex[0] = 0; ResIndex[1] = 2;
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 7);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 7);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        ResIndex[0] = 1; ResIndex[1] = 2;
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 7);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 7);
                        TkMass.clear();
                        TkCharge.clear();
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+K-(Phi)pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[7] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        ResIndex[0] = 0; ResIndex[1] = 1;
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 8);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 8);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        ResIndex[0] = 0; ResIndex[1] = 2;
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 8);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 8);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        ResIndex[0] = 1; ResIndex[1] = 2;
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 8);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Bfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 8);
                        TkMass.clear();
                        TkCharge.clear();
                    }
                    //}}}
                    //printf("-----*****DEBUG:End of DInfo.\n");

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
                            if(chi2_prob_uj < 0.01) continue;
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
                            BInfo.uj_vtxXE         [BInfo.uj_size]= sqrt(ujVFPvtx->error().cxx());
                            BInfo.uj_vtxYE         [BInfo.uj_size]= sqrt(ujVFPvtx->error().cyy());
                            BInfo.uj_vtxZE         [BInfo.uj_size]= sqrt(ujVFPvtx->error().czz());
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

                            //1.good muon candidates?
                            //0.passed all selections
                            BInfo.uj_isGoodCand[BInfo.uj_size] += 1 << 2;
                            if (MuonInfo.isGoodCand[mu1_hindex]%2 == 1 && MuonInfo.isGoodCand[mu2_hindex]%2 == 1)
                                BInfo.uj_isGoodCand[BInfo.uj_size] += 1 << 1;
                            if (((BInfo.uj_isGoodCand[BInfo.uj_size]>>1)&((BInfo.uj_isGoodCand[BInfo.uj_size]>>1)+1))==0)
                                BInfo.uj_isGoodCand[BInfo.uj_size] += 1;

                            BInfo.uj_size++;
                            muonParticles.clear();

                            //////////////////////////////////////////////////////////////////////////
                            // RECONSTRUCTION: J/psi + K
                            //////////////////////////////////////////////////////////////////////////
                            float mass_window[2] = {4.3, 6.4};
                            if(Bchannel_[0] == 1){
                                BranchOut2MuTk(
                                    BInfo,
                                    input_tracks,
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

                            if(Bchannel_[2] == 1){
                                BranchOut2MuX_XtoTkTk(
                                    BInfo,
                                    input_tracks,
                                    isNeededTrack,
                                    v4_mu1,
                                    v4_mu2,
                                    muonPTT,
                                    muonMTT,
                                    B_counter,
                                    mass_window,
                                    JPSI_MASS,
                                    KSHORT_MASS,
                                    0.3,
                                    PION_MASS,        
                                    PION_MASS,
                                    3,
                                    1
                                );
                            }                            
                            
                            //////////////////////////////////////////////////////////////////////////
                            // RECONSTRUCTION: J/psi + K* (K+, Pi-)
                            //////////////////////////////////////////////////////////////////////////

                            if(Bchannel_[3] == 1){
                                BranchOut2MuX_XtoTkTk(
                                    BInfo,
                                    input_tracks,
                                    isNeededTrack,
                                    v4_mu1,
                                    v4_mu2,
                                    muonPTT,
                                    muonMTT,
                                    B_counter,
                                    mass_window,
                                    JPSI_MASS,
                                    KSTAR_MASS,
                                    0.4,
                                    KAON_MASS,        
                                    PION_MASS,
                                    4,
                                    0
                                );
                            }
                            
                            //////////////////////////////////////////////////////////////////////////
                            // RECONSTRUCTION: J/psi + K* (K-, Pi+)
                            //////////////////////////////////////////////////////////////////////////

                            if(Bchannel_[4] == 1){
                                BranchOut2MuX_XtoTkTk(
                                    BInfo,
                                    input_tracks,
                                    isNeededTrack,
                                    v4_mu1,
                                    v4_mu2,
                                    muonPTT,
                                    muonMTT,
                                    B_counter,
                                    mass_window,
                                    JPSI_MASS,
                                    KSTAR_MASS,
                                    0.4,
                                    PION_MASS,        
                                    KAON_MASS,
                                    5,
                                    0
                                );
                            }
                            
                            //////////////////////////////////////////////////////////////////////////
                            // RECONSTRUCTION: J/psi + phi
                            //////////////////////////////////////////////////////////////////////////
                            
                            if(Bchannel_[5] == 1){
                                BranchOut2MuX_XtoTkTk(
                                    BInfo,
                                    input_tracks,
                                    isNeededTrack,
                                    v4_mu1,
                                    v4_mu2,
                                    muonPTT,
                                    muonMTT,
                                    B_counter,
                                    mass_window,
                                    JPSI_MASS,
                                    PHI_MASS,
                                    0.1,
                                    KAON_MASS,        
                                    KAON_MASS,
                                    6,
                                    0
                                );
                            }

                            //////////////////////////////////////////////////////////////////////////
                            // RECONSTRUCTION: J/psi + pi pi <= psi', X(3872), Bs->J/psi f0
                            //////////////////////////////////////////////////////////////////////////
                            mass_window[0] = 3;
                            mass_window[1] = 6.4;
                            if(Bchannel_[6] == 1){
                                BranchOut2MuX_XtoTkTk(
                                    BInfo,
                                    input_tracks,
                                    isNeededTrack,
                                    v4_mu1,
                                    v4_mu2,
                                    muonPTT,
                                    muonMTT,
                                    B_counter,
                                    mass_window,
                                    JPSI_MASS,
                                    -1,
                                    1.6,
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
                    printf("\n");
                    printf("D_counter: ");
                    for(unsigned int i = 0; i < Dchannel_.size(); i++){
                        printf("%d/", D_counter[i]);
                    }
                    printf("\n");//}}}
                    //printf("-----*****DEBUG:End of BInfo.\n");

                    // TrackInfo section {{{
                    for(std::vector<pat::GenericParticle>::const_iterator tk_it=input_tracks.begin();
                        tk_it != input_tracks.end() ; tk_it++){
                        int tk_hindex = int(tk_it - input_tracks.begin());
                        if(tk_hindex>=int(isNeededTrack.size())) break;
                        if (isNeededTrack[tk_hindex]==false) continue;

                        //Create list of relative xb candidates for later filling
                        std::vector<int> listOfRelativeXbCands;//1~nXb
                        std::vector<int> listOfRelativeDCand1;//1~nXb
                        std::vector<int> listOfRelativeDCand2;//1~nXb
                        std::vector<int> listOfRelativeDCand3;//1~nXb
                        std::vector<int> listOfRelativeDCand4;//1~nXb
                        for(int iXb=0; iXb < BInfo.size; iXb++){
                            if(BInfo.rftk1_index[iXb] == -tk_hindex-1){
                                listOfRelativeXbCands.push_back(iXb+1);
                            //}else if(BInfo.rftk2_index[iXb] == -tk_hindex-1){
                            }if(BInfo.rftk2_index[iXb] == -tk_hindex-1){
                                listOfRelativeXbCands.push_back(-iXb-1);
                            }
                        }
                        for(int d=0; d < DInfo.size; d++){
                            if(DInfo.rftk1_index[d] == -tk_hindex-1){
                                listOfRelativeDCand1.push_back(d+1);
                            }
                            if(DInfo.rftk2_index[d] == -tk_hindex-1){
                                listOfRelativeDCand2.push_back(d+1);
                            }
                            if(DInfo.rftk3_index[d] == -tk_hindex-1){
                                listOfRelativeDCand3.push_back(d+1);
                            }
                            if(DInfo.rftk4_index[d] == -tk_hindex-1){
                                listOfRelativeDCand4.push_back(d+1);
                            }
                        }
                        //if (listOfRelativeXbCands.size() == 0) continue;//drop unused tracks

                        TrackInfo.index          [TrackInfo.size] = TrackInfo.size;
                        TrackInfo.handle_index   [TrackInfo.size] = tk_hindex;
                        TrackInfo.charge         [TrackInfo.size] = tk_it->charge();
                        TrackInfo.pt             [TrackInfo.size] = tk_it->pt();
                        TrackInfo.eta            [TrackInfo.size] = tk_it->eta();
                        TrackInfo.phi            [TrackInfo.size] = tk_it->phi();
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
                        TrackInfo.highPurity     [TrackInfo.size] = tk_it->track()->qualityByName("highPurity");
                        TrackInfo.geninfo_index  [TrackInfo.size] = -1;//initialize for later use
                        if(tk_it->track().isNonnull()){
                            for(int tq = 0; tq < reco::TrackBase::qualitySize; tq++){
                            if (tk_it->track()->quality(static_cast<reco::TrackBase::TrackQuality>(tq))) TrackInfo.trackQuality[TrackInfo.size] += 1 << (tq);
                        }}

                        if (!iEvent.isRealData() && RunOnMC_)
                            genTrackPtr [TrackInfo.size] = tk_it->genParticle();

                        //0.passed all selections
                        TrackInfo.isGoodCand     [TrackInfo.size]+= 1 << 1;
                        if (((TrackInfo.isGoodCand[TrackInfo.size]>>1)&((TrackInfo.isGoodCand[TrackInfo.size]>>1)+1))==0)
                            TrackInfo.isGoodCand     [TrackInfo.size]+= 1;

                        //Fill correct track index and track quality to correspond Xb candidate
                        for(unsigned int iCands=0; iCands < listOfRelativeXbCands.size(); iCands++){
                            if (listOfRelativeXbCands[iCands]>0){
                                BInfo.rftk1_index[listOfRelativeXbCands[iCands]-1] = TrackInfo.size;
                                //if ((TrackInfo.isGoodCand[TrackInfo.size]&1) == 1){
                                //    BInfo.isGoodCand[listOfRelativeXbCands[iCands]-1] += 1 << 1;
                                //    if (BInfo.rftk2_index[listOfRelativeXbCands[iCands]-1]>=0){
                                //        BInfo.isGoodCand[listOfRelativeXbCands[iCands]-1] += ((BInfo.isGoodCand[listOfRelativeXbCands[iCands]-1]>>2)&1);
                                //    }
                                //    if (BInfo.rftk2_index[listOfRelativeXbCands[iCands]-1]>=0){
                                //        BInfo.isGoodCand[listOfRelativeXbCands[iCands]-1] += ((BInfo.isGoodCand[listOfRelativeXbCands[iCands]-1]>>2)&1)<<2;
                                //    }else{
                                //        BInfo.isGoodCand[listOfRelativeXbCands[iCands]-1] += 1 << 2;
                                //    }
                                //}
                            }else{
                                BInfo.rftk2_index[-listOfRelativeXbCands[iCands]-1] = TrackInfo.size;
                                //if ((TrackInfo.isGoodCand[TrackInfo.size]&1) == 1){
                                //    BInfo.isGoodCand[-listOfRelativeXbCands[iCands]-1] += 1 << 2;
                                //    if (BInfo.rftk1_index[-listOfRelativeXbCands[iCands]-1]>=0){
                                //        BInfo.isGoodCand[-listOfRelativeXbCands[iCands]-1] += ((BInfo.isGoodCand[-listOfRelativeXbCands[iCands]-1]>>1)&1);
                                //    }
                                //    if (BInfo.rftk1_index[-listOfRelativeXbCands[iCands]-1]>=0){
                                //        BInfo.isGoodCand[-listOfRelativeXbCands[iCands]-1] += ((BInfo.isGoodCand[-listOfRelativeXbCands[iCands]-1]>>2)&1)<<2;
                                //    }else{
                                //        BInfo.isGoodCand[-listOfRelativeXbCands[iCands]-1] += 1 << 2;
                                //    }
                                //}
                            }
                        }
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
                    isGenSignal = (GetAncestor(Myself, 5) | GetAncestor(Myself, 4));
                }//all pi and K from b or c meson

                /*
                if (abs(it_gen->pdgId()) == 13                              &&
                    it_gen->numberOfMothers() == 1                          &&
                    (it_gen->mother()->pdgId() == 443 ||
                     it_gen->mother()->pdgId() == 553 )                     &&
                    it_gen->mother()->numberOfDaughters() >= 2              &&
                    it_gen->mother()->numberOfMothers() == 1                &&
                    (abs(it_gen->mother()->mother()->pdgId()) == 511 ||
                     abs(it_gen->mother()->mother()->pdgId()) == 521 ||
                     abs(it_gen->mother()->mother()->pdgId()) == 531 )
                    //(it_gen->mother()->mother()->pdgId() == 100553 ||
                    // it_gen->mother()->mother()->pdgId() == 100443 )        &&
                    //it_gen->mother()->mother()->numberOfDaughters() == 3    &&
                    //abs(it_gen->mother()->mother()->daughter(1)->pdgId()) == 211 
                   ) isGenSignal = true;//signal mu

                if ((abs(it_gen->pdgId()) == 111 || 
                     abs(it_gen->pdgId()) == 211 || 
                     abs(it_gen->pdgId()) == 310 || 
                     abs(it_gen->pdgId()) == 311 || 
                     abs(it_gen->pdgId()) == 321) &&
                    it_gen->numberOfMothers() == 1 ){
                    //Mother is B
                    if (abs(it_gen->mother()->pdgId()) == 511 || abs(it_gen->mother()->pdgId()) == 521 || abs(it_gen->mother()->pdgId()) == 531){
                    //Mom's 1st dau is Upsilon/Jpsi
                        if (it_gen->mother()->daughter(0)->pdgId() == 553 || it_gen->mother()->daughter(0)->pdgId() == 443) {
                            //Mom's 1st dau's dau are muons
                            if (it_gen->mother()->daughter(0)->numberOfDaughters()>=2){
                                if (abs(it_gen->mother()->daughter(0)->daughter(0)->pdgId()) == 13)
                                  isGenSignal = true;     
                            }
                        }
                    }
                    //A level up
                    if(it_gen->mother()->numberOfMothers() == 1){
                        if(abs(it_gen->mother()->mother()->pdgId()) == 511 || abs(it_gen->mother()->mother()->pdgId()) == 521 || abs(it_gen->mother()->mother()->pdgId()) == 531
                        //Specific for Ks from K0
                        || (abs(it_gen->pdgId()) == 310 && abs(it_gen->mother()->mother()->pdgId()) == 311)
                        ){
                            if(it_gen->mother()->mother()->daughter(0)->pdgId() == 553 ||  it_gen->mother()->mother()->daughter(0)->pdgId() == 443){
                                if (it_gen->mother()->mother()->daughter(0)->numberOfDaughters()>=2) {
                                    if (abs(it_gen->mother()->mother()->daughter(0)->daughter(0)->pdgId()) == 13){
                                    isGenSignal = true;
                                    }
                                }
                            }
                        }
                    }
                }
                //signal pion/kaon

                //Specific for pi from Ks (K0->Ks->pi)                           
                if (abs(it_gen->pdgId()) == 211 && it_gen->numberOfMothers() == 1){
                    if(it_gen->mother()->numberOfMothers() == 1 && abs(it_gen->mother()->pdgId()) == 310){
                        if(it_gen->mother()->mother()->numberOfMothers() == 1 && abs(it_gen->mother()->mother()->pdgId()) == 311){
                            if(abs(it_gen->mother()->mother()->mother()->pdgId()) == 511){
                                if(it_gen->mother()->mother()->mother()->daughter(0)->pdgId() == 553 ||  it_gen->mother()->mother()->mother()->daughter(0)->pdgId() == 443){
                                    if (it_gen->mother()->mother()->mother()->daughter(0)->numberOfDaughters()>=2) {
                                        if (abs(it_gen->mother()->mother()->mother()->daughter(0)->daughter(0)->pdgId()) == 13){
                                            isGenSignal = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                */
                /*
                if ((it_gen->pdgId() == 443 || it_gen->pdgId() == 553)      &&
                    (it_gen->mother()->pdgId() == 100443 ||
                     it_gen->mother()->pdgId() == 100553 )                  &&
                    it_gen->mother()->numberOfDaughters() == 3              &&
                    abs(it_gen->mother()->daughter(1)->pdgId()) == 211      &&
                    it_gen->numberOfDaughters() >= 2                        &&
                    abs(it_gen->daughter(0)->pdgId()) == 13    
                   ) isGenSignal = true;//signal uj

                if ((it_gen->pdgId()==100443 || it_gen->pdgId()==100553)    &&
                    it_gen->numberOfDaughters() == 3                        &&
                    it_gen->daughter(0)->numberOfDaughters() >= 2           &&
                    abs(it_gen->daughter(0)->daughter(0)->pdgId()) == 13    &&
                    abs(it_gen->daughter(1)->pdgId()) == 211 
                   ) isGenSignal = true;//signal xb
                */
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
                for(int nGenDa = 0; nGenDa < std::min(2,int((*sCands)->numberOfDaughters())); nGenDa++){
                    for(std::vector<const reco::Candidate *>::iterator mCands = sel_cands.begin();
                    mCands != sel_cands.end(); mCands++){
                        if((*sCands)->daughter(nGenDa) == *mCands){
                            if(nGenDa == 0) GenInfo.da1[geninfo_idx] = int(mCands-sel_cands.begin());
                            if(nGenDa == 1) GenInfo.da2[geninfo_idx] = int(mCands-sel_cands.begin());
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
    std::cout << "BInfo.size=" << BInfo.size << std::endl;
    root->Fill();
    //std::cout<<"filled!\n";
}

// ------------ method called once each job just after ending the event loop  ------------
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
}

bool Bfinder::GetAncestor(const reco::Candidate* p, int PDGprefix)
{
    if(p->numberOfMothers()==0) return false;
    else{
        const reco::Candidate* MyMom = p->mother(0);
        int mpid = abs(MyMom->pdgId());
        if(abs(int(mpid/100) % 100) == PDGprefix) return true;
        else return GetAncestor(MyMom, PDGprefix);
    }
}

//{{{
void Bfinder::BranchOut2MuTk(
    BInfoBranches &BInfo, 
    std::vector<pat::GenericParticle> input_tracks, 
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
  float muon_sigma = muon_mass*1.e-6;

  for(std::vector<pat::GenericParticle>::const_iterator tk_it1=input_tracks.begin();
      tk_it1 != input_tracks.end() ; tk_it1++){
      tk1_hindex = int(tk_it1 - input_tracks.begin());
      if(tk1_hindex>=int(isNeededTrack.size())) break;
      if (!isNeededTrack[tk1_hindex]) continue;
      if (abs(tk_it1->charge()) != 1) continue;
      
      TLorentzVector v4_tk1;
      v4_tk1.SetPtEtaPhiM(tk_it1->pt(),tk_it1->eta(),tk_it1->phi(),KAON_MASS);
  
      if ((v4_mu1+v4_mu2+v4_tk1).Mag()<mass_window[0]-0.2 || (v4_mu1+v4_mu2+v4_tk1).Mag()>mass_window[1]+0.2) continue;
      XbMassCutLevel[channel_number-1]->Fill(0);
      if((v4_mu1+v4_mu2+v4_tk1).Pt()<bPtCut_)continue;
      XbMassCutLevel[channel_number-1]->Fill(1);
      
      reco::TransientTrack kaonTT(tk_it1->track(), &(*bField) );
      if (!kaonTT.isValid()) continue;
      
      ParticleMass kaon_mass = Tk_MASS;
      float kaon_sigma = kaon_mass*1.e-6;
      
      if (BInfo.size >= MAX_XB) continue;
      
      std::vector<RefCountedKinematicParticle> Xb_candidate;
      Xb_candidate.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
      Xb_candidate.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));
      Xb_candidate.push_back(pFactory.particle(kaonTT,kaon_mass,chi,ndf,kaon_sigma));
      RefCountedKinematicTree xbVFT;
      
      ParticleMass uj_mass = MuMu_MASS;
      MultiTrackKinematicConstraint *uj_c = new TwoTrackMassKinematicConstraint(uj_mass);
      KinematicConstrainedVertexFitter kcvFitter;
      xbVFT = kcvFitter.fit(Xb_candidate, uj_c);
      
      if (!xbVFT->isValid()) continue;
      XbMassCutLevel[channel_number-1]->Fill(2);

      xbVFT->movePointerToTheTop();
      RefCountedKinematicParticle     xbVFP       = xbVFT->currentParticle();
      RefCountedKinematicVertex       xbVFPvtx    = xbVFT->currentDecayVertex();
      if (!xbVFPvtx->vertexIsValid()) continue;
      XbMassCutLevel[channel_number-1]->Fill(3);

      std::vector<RefCountedKinematicParticle> xCands  = xbVFT->finalStateParticles();
      
      double chi2_prob = TMath::Prob(xbVFPvtx->chiSquared(),xbVFPvtx->degreesOfFreedom());
      if (chi2_prob < 0.01) continue;
      XbMassCutLevel[channel_number-1]->Fill(4);
      
      if (xbVFP->currentState().mass()<mass_window[0] || xbVFP->currentState().mass()>mass_window[1]) continue;
      XbMassCutLevel[channel_number-1]->Fill(5);
      
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
      BInfo.vtxX[BInfo.size]    = xbVFPvtx->position().x();
      BInfo.vtxY[BInfo.size]    = xbVFPvtx->position().y();
      BInfo.vtxZ[BInfo.size]    = xbVFPvtx->position().z();
      BInfo.vtxXE[BInfo.size]   = sqrt(xbVFPvtx->error().cxx());
      BInfo.vtxYE[BInfo.size]   = sqrt(xbVFPvtx->error().cyy());
      BInfo.vtxZE[BInfo.size]   = sqrt(xbVFPvtx->error().czz());
      BInfo.vtxdof[BInfo.size]  = xbVFPvtx->degreesOfFreedom();
      BInfo.vtxchi2[BInfo.size] = xbVFPvtx->chiSquared();
      
      BInfo.rfuj_index[BInfo.size]  = BInfo.uj_size-1;
      BInfo.rftk1_index[BInfo.size] = -tk1_hindex-1;
      BInfo.rftk2_index[BInfo.size] = -tk1_hindex-1;
      
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
//{{{
void Bfinder::BranchOut2MuX_XtoTkTk(
    BInfoBranches &BInfo, 
    std::vector<pat::GenericParticle> input_tracks, 
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
    float muon_sigma = muon_mass*1.e-6;
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
            if(tk2_hindex>=int(isNeededTrack.size())) break;
            tk2_hindex = int(tk_it2 - input_tracks.begin());
            if (!isNeededTrack[tk2_hindex]) continue;
            if (tk_it2->charge()>0) continue;
            
            TLorentzVector v4_tk1,v4_tk2;
            v4_tk1.SetPtEtaPhiM(tk_it1->pt(),tk_it1->eta(),tk_it1->phi(),Tk1_MASS);
            v4_tk2.SetPtEtaPhiM(tk_it2->pt(),tk_it2->eta(),tk_it2->phi(),Tk2_MASS);
            if(TkTk_MASS > 0) {if (fabs((v4_tk1+v4_tk2).Mag()-TkTk_MASS)>TkTk_window) continue;}
            else {if (fabs((v4_tk1+v4_tk2).Mag())>TkTk_window) continue;}//if no tktk mass constrain, require it to be at least < some window
            
            if ((v4_mu1+v4_mu2+v4_tk1+v4_tk2).Mag()<mass_window[0]-0.2 || (v4_mu1+v4_mu2+v4_tk1+v4_tk2).Mag()>mass_window[1]+0.2) continue;
            XbMassCutLevel[channel_number-1]->Fill(0);
            if((v4_mu1+v4_mu2+v4_tk1+v4_tk2).Pt()<bPtCut_)continue;
            XbMassCutLevel[channel_number-1]->Fill(1);
            
            reco::TransientTrack tk1PTT(tk_it1->track(), &(*bField) );
            reco::TransientTrack tk2MTT(tk_it2->track(), &(*bField) );
            if (!tk1PTT.isValid()) continue;
            if (!tk2MTT.isValid()) continue;
            
            ParticleMass tk1_mass = Tk1_MASS;
            float tk1_sigma = tk1_mass*1.e-6;
            ParticleMass tk2_mass = Tk2_MASS;
            float tk2_sigma = tk2_mass*1.e-6;
            
            if (BInfo.size >= MAX_XB) continue;
            //doing tktk fit
            std::vector<RefCountedKinematicParticle> tktk_candidate;
            tktk_candidate.push_back(pFactory.particle(tk1PTT,tk1_mass,chi,ndf,tk1_sigma));
            tktk_candidate.push_back(pFactory.particle(tk2MTT,tk2_mass,chi,ndf,tk2_sigma));
            
            KinematicParticleVertexFitter   tktk_fitter;
            RefCountedKinematicTree         tktk_VFT;
            tktk_VFT = tktk_fitter.fit(tktk_candidate);
            if(!tktk_VFT->isValid()) continue;
            XbMassCutLevel[channel_number-1]->Fill(2);
            
            tktk_VFT->movePointerToTheTop();
            RefCountedKinematicParticle tktk_VFP   = tktk_VFT->currentParticle();
            RefCountedKinematicVertex   tktk_VFPvtx = tktk_VFT->currentDecayVertex();
            double chi2_prob_tktk = TMath::Prob(tktk_VFPvtx->chiSquared(),
                                                tktk_VFPvtx->degreesOfFreedom());
            if(chi2_prob_tktk < 0.01) continue;
            XbMassCutLevel[channel_number-1]->Fill(3);

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
            
            ParticleMass uj_mass = MuMu_MASS;
            MultiTrackKinematicConstraint *uj_c = new  TwoTrackMassKinematicConstraint(uj_mass);
            KinematicConstrainedVertexFitter kcvFitter;
            xbVFT = kcvFitter.fit(Xb_candidate, uj_c);
            
            if (!xbVFT->isValid()) continue;
            XbMassCutLevel[channel_number-1]->Fill(4);

            xbVFT->movePointerToTheTop();
            RefCountedKinematicParticle     xbVFP       = xbVFT->currentParticle();
            RefCountedKinematicVertex       xbVFPvtx    = xbVFT->currentDecayVertex();
            if (!xbVFPvtx->vertexIsValid()) continue;
            XbMassCutLevel[channel_number-1]->Fill(5);
            
            std::vector<RefCountedKinematicParticle> tktkCands  = tktk_VFT->finalStateParticles();
            std::vector<RefCountedKinematicParticle> xCands  = xbVFT->finalStateParticles();
            
            double chi2_prob = TMath::Prob(xbVFPvtx->chiSquared(),xbVFPvtx->degreesOfFreedom());
            if (chi2_prob < 0.01) continue;
            XbMassCutLevel[channel_number-1]->Fill(6);
            
            //Cut out a mass window
            if (xbVFP->currentState().mass()<mass_window[0]|| xbVFP->currentState().mass()>mass_window[1]) continue;
            XbMassCutLevel[channel_number-1]->Fill(7);
            
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
            BInfo.vtxX[BInfo.size]    = xbVFPvtx->position().x();
            BInfo.vtxY[BInfo.size]    = xbVFPvtx->position().y();
            BInfo.vtxZ[BInfo.size]    = xbVFPvtx->position().z();
            BInfo.vtxXE[BInfo.size]   = sqrt(xbVFPvtx->error().cxx());
            BInfo.vtxYE[BInfo.size]   = sqrt(xbVFPvtx->error().cyy());
            BInfo.vtxZE[BInfo.size]   = sqrt(xbVFPvtx->error().czz());
            BInfo.vtxdof[BInfo.size]  = xbVFPvtx->degreesOfFreedom();
            BInfo.vtxchi2[BInfo.size] = xbVFPvtx->chiSquared();
            
            BInfo.rfuj_index[BInfo.size]  = BInfo.uj_size-1;
            BInfo.rftk1_index[BInfo.size] = -tk1_hindex-1;
            BInfo.rftk2_index[BInfo.size] = -tk2_hindex-1;
            
            //tktk fit info
            BInfo.tktk_mass[BInfo.size]    = tktk_4vec.Mag();
            BInfo.tktk_pt[BInfo.size]      = tktk_4vec.Pt();
            BInfo.tktk_eta[BInfo.size]      = tktk_4vec.Eta();
            BInfo.tktk_phi[BInfo.size]      = tktk_4vec.Phi();
            BInfo.tktk_px[BInfo.size]      = tktk_4vec.Px();
            BInfo.tktk_py[BInfo.size]      = tktk_4vec.Py();
            BInfo.tktk_pz[BInfo.size]      = tktk_4vec.Pz();
            BInfo.tktk_vtxX[BInfo.size]    = tktk_VFPvtx->position().x();
            BInfo.tktk_vtxY[BInfo.size]    = tktk_VFPvtx->position().y();
            BInfo.tktk_vtxZ[BInfo.size]    = tktk_VFPvtx->position().z();
            BInfo.tktk_vtxXE[BInfo.size]   = sqrt(tktk_VFPvtx->error().cxx());
            BInfo.tktk_vtxYE[BInfo.size]   = sqrt(tktk_VFPvtx->error().cyy());
            BInfo.tktk_vtxZE[BInfo.size]   = sqrt(tktk_VFPvtx->error().czz());
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
//{{{
void Bfinder::BranchOutNTk(//input 2~4 tracks
    DInfoBranches &DInfo, 
    std::vector<pat::GenericParticle> input_tracks, 
    std::vector<bool> isNeededTrack,
    std::vector<int> &D_counter,
    float *mass_window,
    std::vector<double> TkMass,
    std::vector<int> TkCharge,
    double tktkRes_mass,
    double tktkRes_mass_window,
    double *ResIndex,
    bool doConstrainFit,
    int Dchannel_number
){
    if(Dchannel_number > (int)Dchannel_.size()){ printf("Exceeding defined # of channel, exit"); return;}
    float chi = 0.;
    float ndf = 0.;
    KinematicParticleFactoryFromTransientTrack pFactory;
    //    ParticleMass muon_mass = MUON_MASS; //pdg mass
    //    float muon_sigma = muon_mass*1.e-6;

    int tk1_hindex = -1;
    int tk2_hindex = -1;
    int tk3_hindex = -1;
    int tk4_hindex = -1;
    std::vector< std::vector<double> > selectedTkhidxSet;
    for(std::vector<pat::GenericParticle>::const_iterator tk_it1=input_tracks.begin();
            tk_it1 != input_tracks.end() ; tk_it1++){
        std::vector<double> selectedTkhidx;
        tk1_hindex = int(tk_it1 - input_tracks.begin());
        if(tk1_hindex>=int(isNeededTrack.size())) break;
        if(!isNeededTrack[tk1_hindex]) continue;
        if(tk_it1->charge()!=TkCharge[0]) continue;

        //for(std::vector<pat::GenericParticle>::const_iterator tk_it2=input_tracks.begin();
        for(std::vector<pat::GenericParticle>::const_iterator tk_it2=tk_it1+1;
                tk_it2 != input_tracks.end() ; tk_it2++){
            tk2_hindex = int(tk_it2 - input_tracks.begin());
            if(tk2_hindex>=int(isNeededTrack.size())) break;
            if(!isNeededTrack[tk2_hindex]) continue;
            if(tk_it2->charge()!=TkCharge[1]) continue;
            if(tk2_hindex==tk1_hindex) continue;
            if(TkMass.size()==2){
                selectedTkhidx.push_back(tk1_hindex);
                selectedTkhidx.push_back(tk2_hindex);
                selectedTkhidxSet.push_back(selectedTkhidx);
                selectedTkhidx.clear();
                continue;
            }

            //for(std::vector<pat::GenericParticle>::const_iterator tk_it3=input_tracks.begin();
            for(std::vector<pat::GenericParticle>::const_iterator tk_it3=tk_it2+1;
                    tk_it3 != input_tracks.end() ; tk_it3++){
                tk3_hindex = int(tk_it3 - input_tracks.begin());
                if(tk3_hindex>=int(isNeededTrack.size())) break;
                if(!isNeededTrack[tk3_hindex]) continue;
                if(tk_it3->charge()!=TkCharge[2]) continue;
                if(tk3_hindex==tk1_hindex) continue;
                if(tk3_hindex==tk2_hindex) continue;
                if(TkMass.size()==3){
                    selectedTkhidx.push_back(tk1_hindex);
                    selectedTkhidx.push_back(tk2_hindex);
                    selectedTkhidx.push_back(tk3_hindex);
                    selectedTkhidxSet.push_back(selectedTkhidx);
                    selectedTkhidx.clear();
                    continue;
                }

                //for(std::vector<pat::GenericParticle>::const_iterator tk_it4=input_tracks.begin();
                for(std::vector<pat::GenericParticle>::const_iterator tk_it4=tk_it3+1;
                        tk_it4 != input_tracks.end() ; tk_it4++){
                    tk4_hindex = int(tk_it4 - input_tracks.begin());
                    if(tk4_hindex>=int(isNeededTrack.size())) break;
                    if(!isNeededTrack[tk4_hindex]) continue;
                    if(tk_it4->charge()!=TkCharge[3]) continue;
                    if(tk4_hindex==tk1_hindex) continue;
                    if(tk4_hindex==tk2_hindex) continue;
                    if(tk4_hindex==tk3_hindex) continue;
                    if(TkMass.size()==4){
                        selectedTkhidx.push_back(tk1_hindex);
                        selectedTkhidx.push_back(tk2_hindex);
                        selectedTkhidx.push_back(tk3_hindex);
                        selectedTkhidx.push_back(tk4_hindex);
                        selectedTkhidxSet.push_back(selectedTkhidx);
                        selectedTkhidx.clear();
                        continue;
                    }
                }
            }
        }
    }

    TLorentzVector v4_tk1,v4_tk2,v4_tk3,v4_tk4, v4_D;
    TLorentzVector v4_tkRes1,v4_tkRes2;
    //std::cout<<"selectedTkhidxSet.size: "<<selectedTkhidxSet.size()<<std::endl;
    for(int i = 0; i < int(selectedTkhidxSet.size()); i++){
        //check mass before fit
        v4_tk1.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i][0]].pt(),input_tracks[selectedTkhidxSet[i][0]].eta(),input_tracks[selectedTkhidxSet[i][0]].phi(),TkMass[0]);
        v4_tk2.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i][1]].pt(),input_tracks[selectedTkhidxSet[i][1]].eta(),input_tracks[selectedTkhidxSet[i][1]].phi(),TkMass[1]);
        v4_D = v4_tk1 + v4_tk2;
        if(TkMass.size()>2) {
            v4_tk3.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i][2]].pt(),input_tracks[selectedTkhidxSet[i][2]].eta(),input_tracks[selectedTkhidxSet[i][2]].phi(),TkMass[2]);
            v4_D = v4_tk1 + v4_tk2 + v4_tk3;
        }
        if(TkMass.size()>3) {
            v4_tk4.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i][3]].pt(),input_tracks[selectedTkhidxSet[i][3]].eta(),input_tracks[selectedTkhidxSet[i][3]].phi(),TkMass[3]);
            v4_D = v4_tk1 + v4_tk2 + v4_tk3 + v4_tk4;
        }
        if(v4_D.Mag()<mass_window[0]-0.05 || v4_D.Mag()>mass_window[1]+0.05) continue;

        //if there's tktk Res, also check tktk res mass window
        if(tktkRes_mass > 0) {
            v4_tkRes1.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i][ResIndex[0]]].pt(),input_tracks[selectedTkhidxSet[i][ResIndex[0]]].eta(),input_tracks[selectedTkhidxSet[i][ResIndex[0]]].phi(),TkMass[ResIndex[0]]);
            v4_tkRes2.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i][ResIndex[1]]].pt(),input_tracks[selectedTkhidxSet[i][ResIndex[1]]].eta(),input_tracks[selectedTkhidxSet[i][ResIndex[1]]].phi(),TkMass[ResIndex[1]]);        
            if (fabs((v4_tkRes1+v4_tkRes2).Mag()-tktkRes_mass) > tktkRes_mass_window) continue;
        }
        DMassCutLevel[Dchannel_number-1]->Fill(0);

        if(v4_D.Pt() < dPtCut_)continue;
        DMassCutLevel[Dchannel_number-1]->Fill(1);

        //if(fabs(v4_D.Eta()) > dEtaCut_)continue;
        DMassCutLevel[Dchannel_number-1]->Fill(2);

        if (DInfo.size >= MAX_XB) break;
        std::vector<RefCountedKinematicParticle> tktk_candidate;

        //push back the Res tracks as first tow tracks
        reco::TransientTrack tkTT1(input_tracks[selectedTkhidxSet[i][ResIndex[0]]].track(), &(*bField) );
        if (!tkTT1.isValid()) continue;
        ParticleMass tk_mass = TkMass[ResIndex[0]];
        float tk_sigma = tk_mass*1.e-6;
        tktk_candidate.push_back(pFactory.particle(tkTT1,tk_mass,chi,ndf,tk_sigma));

        reco::TransientTrack tkTT2(input_tracks[selectedTkhidxSet[i][ResIndex[1]]].track(), &(*bField) );
        if (!tkTT2.isValid()) continue;
        tk_mass = TkMass[ResIndex[1]];
        tk_sigma = tk_mass*1.e-6;
        tktk_candidate.push_back(pFactory.particle(tkTT2,tk_mass,chi,ndf,tk_sigma));
        
        //push back the other tracks
        for(int p = 0; p < int(selectedTkhidxSet[0].size()); p++){        
            if(p == ResIndex[0] || p == ResIndex[1]) continue;
            reco::TransientTrack tkTT(input_tracks[selectedTkhidxSet[i][p]].track(), &(*bField) );
            if (!tkTT.isValid()) continue;
            tk_mass = TkMass[p];
            tk_sigma = tk_mass*1.e-6;
            tktk_candidate.push_back(pFactory.particle(tkTT,tk_mass,chi,ndf,tk_sigma));
        }
        DMassCutLevel[Dchannel_number-1]->Fill(3);

        //doing tktk fit
        KinematicParticleVertexFitter   tktk_fitter;
        RefCountedKinematicTree         tktk_VFT;
        KinematicConstrainedVertexFitter kcv_tktk_fitter;
        if(tktkRes_mass<=0){
            tktk_VFT = tktk_fitter.fit(tktk_candidate);
        }
            
        //if these's a tktk Res, check its fit validity also
        KinematicParticleVertexFitter   tktkRes_fitter;
        RefCountedKinematicTree         tktkRes_VFT;
        RefCountedKinematicParticle tktkRes_VFP;
        RefCountedKinematicVertex   tktkRes_VFPvtx;
        std::vector<RefCountedKinematicParticle> tktkResCands;
        if(tktkRes_mass>0){
            reco::TransientTrack tk1Res(input_tracks[selectedTkhidxSet[i][ResIndex[0]]].track(), &(*bField) );
            reco::TransientTrack tk2Res(input_tracks[selectedTkhidxSet[i][ResIndex[1]]].track(), &(*bField) );
            ParticleMass tk1_mass = TkMass[ResIndex[0]];
            float tk1_sigma = tk1_mass*1.e-6;
            ParticleMass tk2_mass = TkMass[ResIndex[1]];
            float tk2_sigma = tk2_mass*1.e-6;
            std::vector<RefCountedKinematicParticle> tktkRes_candidate;
            tktkRes_candidate.push_back(pFactory.particle(tk1Res,tk1_mass,chi,ndf,tk1_sigma));
            tktkRes_candidate.push_back(pFactory.particle(tk2Res,tk2_mass,chi,ndf,tk2_sigma));
            tktkRes_VFT = tktkRes_fitter.fit(tktkRes_candidate);
            if(!tktkRes_VFT->isValid()) continue;
            DMassCutLevel[Dchannel_number-1]->Fill(4);
            tktkRes_VFP   = tktkRes_VFT->currentParticle();
            tktkRes_VFPvtx = tktkRes_VFT->currentDecayVertex();
            double chi2_prob_tktkRes = TMath::Prob(tktkRes_VFPvtx->chiSquared(),tktkRes_VFPvtx->degreesOfFreedom());
            //if(chi2_prob_tktkRes < 0.01) continue;
            if(chi2_prob_tktkRes < 0.05) continue;
            DMassCutLevel[Dchannel_number-1]->Fill(5);

            if(doConstrainFit){
                ParticleMass tktkResMass = tktkRes_mass;
                MultiTrackKinematicConstraint *tktkResConstraint = new TwoTrackMassKinematicConstraint(tktkResMass);
                tktk_VFT = kcv_tktk_fitter.fit(tktk_candidate, tktkResConstraint);
            }
            else tktk_VFT = tktk_fitter.fit(tktk_candidate);
        }

        if(!tktk_VFT->isValid()) continue;
        DMassCutLevel[Dchannel_number-1]->Fill(6);

        tktk_VFT->movePointerToTheTop();
        RefCountedKinematicParticle tktk_VFP   = tktk_VFT->currentParticle();
        RefCountedKinematicVertex   tktk_VFPvtx = tktk_VFT->currentDecayVertex();
        if (!tktk_VFPvtx->vertexIsValid()) continue;
        DMassCutLevel[Dchannel_number-1]->Fill(7);

        double chi2_prob_tktk = TMath::Prob(tktk_VFPvtx->chiSquared(),tktk_VFPvtx->degreesOfFreedom());
        //if(chi2_prob_tktk < 0.01) continue;
        if(chi2_prob_tktk < 0.05) continue;
        DMassCutLevel[Dchannel_number-1]->Fill(8);

        std::vector<RefCountedKinematicParticle> tktkCands  = tktk_VFT->finalStateParticles();

        //Cut out a mass window
        if (tktk_VFP->currentState().mass()<mass_window[0] || tktk_VFP->currentState().mass()>mass_window[1]) continue;
        DMassCutLevel[Dchannel_number-1]->Fill(9);

        TLorentzVector tktk_4vec,tktk_tk1_4vec, tktk_tk2_4vec, tktk_tk3_4vec, tktk_tk4_4vec;
        tktk_4vec.SetPxPyPzE(tktk_VFP->currentState().kinematicParameters().momentum().x(),
                tktk_VFP->currentState().kinematicParameters().momentum().y(),
                tktk_VFP->currentState().kinematicParameters().momentum().z(),
                tktk_VFP->currentState().kinematicParameters().energy());

        tktk_tk1_4vec.SetPxPyPzE(tktkCands[0]->currentState().kinematicParameters().momentum().x(),
                tktkCands[0]->currentState().kinematicParameters().momentum().y(),
                tktkCands[0]->currentState().kinematicParameters().momentum().z(),
                tktkCands[0]->currentState().kinematicParameters().energy());
        tktk_tk2_4vec.SetPxPyPzE(tktkCands[1]->currentState().kinematicParameters().momentum().x(),
                tktkCands[1]->currentState().kinematicParameters().momentum().y(),
                tktkCands[1]->currentState().kinematicParameters().momentum().z(),
                tktkCands[1]->currentState().kinematicParameters().energy());
        if(TkMass.size()>2){
            tktk_tk3_4vec.SetPxPyPzE(tktkCands[2]->currentState().kinematicParameters().momentum().x(),
                    tktkCands[2]->currentState().kinematicParameters().momentum().y(),
                    tktkCands[2]->currentState().kinematicParameters().momentum().z(),
                    tktkCands[2]->currentState().kinematicParameters().energy());
        }
        if(TkMass.size()>3){
            tktk_tk4_4vec.SetPxPyPzE(tktkCands[3]->currentState().kinematicParameters().momentum().x(),
                    tktkCands[3]->currentState().kinematicParameters().momentum().y(),
                    tktkCands[3]->currentState().kinematicParameters().momentum().z(),
                    tktkCands[3]->currentState().kinematicParameters().energy());
        }

        //tktkRes fit info
        if(tktkRes_mass>0){
            tktkResCands  = tktkRes_VFT->finalStateParticles();
            TLorentzVector tktkRes_tk1_4vec, tktkRes_tk2_4vec,  tktkRes_4vec;
            tktkRes_tk1_4vec.SetPxPyPzE(tktkResCands[0]->currentState().kinematicParameters().momentum().x(),
                    tktkResCands[0]->currentState().kinematicParameters().momentum().y(),
                    tktkResCands[0]->currentState().kinematicParameters().momentum().z(),
                    tktkResCands[0]->currentState().kinematicParameters().energy());
            tktkRes_tk2_4vec.SetPxPyPzE(tktkResCands[1]->currentState().kinematicParameters().momentum().x(),
                    tktkCands[1]->currentState().kinematicParameters().momentum().y(),
                    tktkCands[1]->currentState().kinematicParameters().momentum().z(),
                    tktkCands[1]->currentState().kinematicParameters().energy());
            tktkRes_4vec.SetPxPyPzE(tktkRes_VFP->currentState().kinematicParameters().momentum().x(),
                    tktkRes_VFP->currentState().kinematicParameters().momentum().y(),
                    tktkRes_VFP->currentState().kinematicParameters().momentum().z(),
                    tktkRes_VFP->currentState().kinematicParameters().energy());
            DInfo.tktkRes_mass[DInfo.size]            = tktkRes_4vec.Mag();
            DInfo.tktkRes_pt[DInfo.size]              = tktkRes_4vec.Pt();
            DInfo.tktkRes_eta[DInfo.size]             = tktkRes_4vec.Eta();
            DInfo.tktkRes_phi[DInfo.size]             = tktkRes_4vec.Phi();
            DInfo.tktkRes_px[DInfo.size]              = tktkRes_4vec.Px();
            DInfo.tktkRes_py[DInfo.size]              = tktkRes_4vec.Py();
            DInfo.tktkRes_pz[DInfo.size]              = tktkRes_4vec.Pz();
            DInfo.tktkRes_vtxX[DInfo.size]            = tktkRes_VFPvtx->position().x();
            DInfo.tktkRes_vtxY[DInfo.size]            = tktkRes_VFPvtx->position().y();
            DInfo.tktkRes_vtxZ[DInfo.size]            = tktkRes_VFPvtx->position().z();
            DInfo.tktkRes_vtxXE[DInfo.size]           = sqrt(tktkRes_VFPvtx->error().cxx());
            DInfo.tktkRes_vtxYE[DInfo.size]           = sqrt(tktkRes_VFPvtx->error().cyy());
            DInfo.tktkRes_vtxZE[DInfo.size]           = sqrt(tktkRes_VFPvtx->error().czz());
            DInfo.tktkRes_vtxdof[DInfo.size]          = tktkRes_VFPvtx->degreesOfFreedom();
            DInfo.tktkRes_vtxchi2[DInfo.size]         = tktkRes_VFPvtx->chiSquared();
    
            DInfo.tktkRes_rftk1_pt[DInfo.size]        = tktkRes_tk1_4vec.Pt();
            DInfo.tktkRes_rftk1_eta[DInfo.size]       = tktkRes_tk1_4vec.Eta();
            DInfo.tktkRes_rftk1_phi[DInfo.size]       = tktkRes_tk1_4vec.Phi();
            DInfo.tktkRes_rftk2_pt[DInfo.size]        = tktkRes_tk2_4vec.Pt();
            DInfo.tktkRes_rftk2_eta[DInfo.size]       = tktkRes_tk2_4vec.Eta();
            DInfo.tktkRes_rftk2_phi[DInfo.size]       = tktkRes_tk2_4vec.Phi();
        }
        //fit info
        DInfo.index[DInfo.size]           = DInfo.size;
        DInfo.b4fit_mass[DInfo.size]      = v4_D.Mag();
        DInfo.b4fit_pt[DInfo.size]        = v4_D.Pt();
        DInfo.b4fit_eta[DInfo.size]       = v4_D.Eta();
        DInfo.b4fit_phi[DInfo.size]       = v4_D.Phi();
        DInfo.mass[DInfo.size]            = tktk_4vec.Mag();
        DInfo.pt[DInfo.size]              = tktk_4vec.Pt();
        DInfo.eta[DInfo.size]             = tktk_4vec.Eta();
        DInfo.phi[DInfo.size]             = tktk_4vec.Phi();
        DInfo.px[DInfo.size]              = tktk_4vec.Px();
        DInfo.py[DInfo.size]              = tktk_4vec.Py();
        DInfo.pz[DInfo.size]              = tktk_4vec.Pz();
        DInfo.vtxX[DInfo.size]            = tktk_VFPvtx->position().x();
        DInfo.vtxY[DInfo.size]            = tktk_VFPvtx->position().y();
        DInfo.vtxZ[DInfo.size]            = tktk_VFPvtx->position().z();
        DInfo.vtxXE[DInfo.size]           = sqrt(tktk_VFPvtx->error().cxx());
        DInfo.vtxYE[DInfo.size]           = sqrt(tktk_VFPvtx->error().cyy());
        DInfo.vtxZE[DInfo.size]           = sqrt(tktk_VFPvtx->error().czz());
        DInfo.vtxdof[DInfo.size]          = tktk_VFPvtx->degreesOfFreedom();
        DInfo.vtxchi2[DInfo.size]         = tktk_VFPvtx->chiSquared();

        DInfo.rftk1_px[DInfo.size]        = tktk_tk1_4vec.Px();
        DInfo.rftk1_py[DInfo.size]        = tktk_tk1_4vec.Py();
        DInfo.rftk1_pz[DInfo.size]        = tktk_tk1_4vec.Pz();
        DInfo.rftk2_px[DInfo.size]        = tktk_tk2_4vec.Px();
        DInfo.rftk2_py[DInfo.size]        = tktk_tk2_4vec.Py();
        DInfo.rftk2_pz[DInfo.size]        = tktk_tk2_4vec.Pz();
        DInfo.rftk1_pt[DInfo.size]        = tktk_tk1_4vec.Pt();
        DInfo.rftk1_eta[DInfo.size]       = tktk_tk1_4vec.Eta();
        DInfo.rftk1_phi[DInfo.size]       = tktk_tk1_4vec.Phi();
        DInfo.rftk2_pt[DInfo.size]        = tktk_tk2_4vec.Pt();
        DInfo.rftk2_eta[DInfo.size]       = tktk_tk2_4vec.Eta();
        DInfo.rftk2_phi[DInfo.size]       = tktk_tk2_4vec.Phi();
        DInfo.rftk1_index[DInfo.size]     = -selectedTkhidxSet[i][0]-1;
        DInfo.rftk2_index[DInfo.size]     = -selectedTkhidxSet[i][1]-1;
        if(TkMass.size()>2){
            DInfo.rftk3_px[DInfo.size]    = tktk_tk3_4vec.Px();
            DInfo.rftk3_py[DInfo.size]    = tktk_tk3_4vec.Py();
            DInfo.rftk3_pz[DInfo.size]    = tktk_tk3_4vec.Pz();
            DInfo.rftk3_pt[DInfo.size]    = tktk_tk3_4vec.Pt();
            DInfo.rftk3_eta[DInfo.size]   = tktk_tk3_4vec.Eta();
            DInfo.rftk3_phi[DInfo.size]   = tktk_tk3_4vec.Phi();
            DInfo.rftk3_index[DInfo.size] = -selectedTkhidxSet[i][2]-1;
        }
        if(TkMass.size()>3){
            DInfo.rftk4_px[DInfo.size]    = tktk_tk4_4vec.Px();
            DInfo.rftk4_py[DInfo.size]    = tktk_tk4_4vec.Py();
            DInfo.rftk4_pz[DInfo.size]    = tktk_tk4_4vec.Pz();
            DInfo.rftk4_pt[DInfo.size]    = tktk_tk4_4vec.Pt();
            DInfo.rftk4_eta[DInfo.size]   = tktk_tk4_4vec.Eta();
            DInfo.rftk4_phi[DInfo.size]   = tktk_tk4_4vec.Phi();
            DInfo.rftk4_index[DInfo.size] = -selectedTkhidxSet[i][3]-1;
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
DEFINE_FWK_MODULE(Bfinder);
