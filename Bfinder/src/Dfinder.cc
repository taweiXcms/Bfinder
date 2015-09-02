// vim:set ts=4 sw=4 fdm=marker et:
// Ntuplt creator for D meson related analysis.
// Maintain and contact: ta-wei wang
// Email: "tawei@mit.edu" or "ta-wei.wang@cern.ch"
#include <memory>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <list>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/TrimmedKalmanVertexFinder/interface/KalmanTrimmedVertexFinder.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Bfinder/Bfinder/interface/format.h"
//#include "Bfinder/Bfinder/interface/TriggerBooking.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>

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
        
        virtual bool GetAncestor(const reco::Candidate* p, int PDGprefix);
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
        std::vector<int> Dchannel_;
//        edm::InputTag hltLabel_;
        edm::InputTag genLabel_;
        edm::InputTag trackLabel_;
        edm::InputTag puInfoLabel_;
        edm::InputTag bsLabel_;
        edm::InputTag pvLabel_;
        double tkPtCut_;
        double tkEtaCut_;
        double dPtCut_;
        double dEtaCut_;
        double VtxChiProbCut_;
        bool RunOnMC_;
        bool doTkPreCut_;

        edm::Service<TFileService> fs;
        TTree *root;
        EvtInfoBranches     EvtInfo;
        VtxInfoBranches     VtxInfo;
        TrackInfoBranches   TrackInfo;
        DInfoBranches       DInfo;
        GenInfoBranches     GenInfo;

        //histograms
        TH1F *TrackCutLevel;
        //How many channel
        static int const Nchannel = 20;
        std::vector<TH1F*> DMassCutLevel;
        
};//}}}

void Dfinder::beginJob()
{//{{{
    root = fs->make<TTree>("root","root");
    EvtInfo.regTree(root);
    VtxInfo.regTree(root);
    TrackInfo.regTree(root);
    DInfo.regTree(root);
    GenInfo.regTree(root);
}//}}}

Dfinder::Dfinder(const edm::ParameterSet& iConfig):theConfig(iConfig)
{//{{{
    //now do what ever initialization is needed
    Dchannel_= iConfig.getParameter<std::vector<int> >("Dchannel");
    genLabel_           = iConfig.getParameter<edm::InputTag>("GenLabel");
    trackLabel_         = iConfig.getParameter<edm::InputTag>("TrackLabel");
//    hltLabel_           = iConfig.getParameter<edm::InputTag>("HLTLabel");
    puInfoLabel_        = iConfig.getParameter<edm::InputTag>("PUInfoLabel");
    bsLabel_        = iConfig.getParameter<edm::InputTag>("BSLabel");
    pvLabel_        = iConfig.getParameter<edm::InputTag>("PVLabel");

    doTkPreCut_ = iConfig.getParameter<bool>("doTkPreCut");
    tkPtCut_ = iConfig.getParameter<double>("tkPtCut");
    tkEtaCut_ = iConfig.getParameter<double>("tkEtaCut");
    dPtCut_ = iConfig.getParameter<double>("dPtCut");
    dEtaCut_ = iConfig.getParameter<double>("dEtaCut");
    VtxChiProbCut_ = iConfig.getParameter<double>("VtxChiProbCut");
    RunOnMC_ = iConfig.getParameter<bool>("RunOnMC");

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
    //std::cout << "*************************\nReconstructing event number: " << iEvent.id() << "\n";
    using namespace edm;
    using namespace reco;
    //ESHandle<MagneticField> bField;
    iSetup.get<IdealMagneticFieldRecord>().get(bField);

    // Change used muon and track collections
    edm::Handle< std::vector<pat::GenericParticle> > tks;
    iEvent.getByLabel(trackLabel_, tks);

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
    iEvent.getByLabel(bsLabel_, beamSpotHandle);
    if (beamSpotHandle.isValid()){
        beamSpot = *beamSpotHandle;
        theBeamSpotV = Vertex(beamSpot.position(), beamSpot.covariance3D());
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
                    
                    // DInfo section{{{
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+pi-
                    //////////////////////////////////////////////////////////////////////////
                    float mass_window[2] = {1.7,2.1};
                    double ResIndex[2] = {0, 1};
                    std::vector<double> TkMass;
                    std::vector<int> TkCharge;
                    if(Dchannel_[0] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 1);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 1);
                        TkMass.clear();
                        TkCharge.clear();
                    }

                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K-pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[1] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 2);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 2);
                        TkMass.clear();
                        TkCharge.clear();
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K-pi+pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[2] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 3);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 3);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 3);
                        TkMass.clear();
                        TkCharge.clear();
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+pi-pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[3] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 4);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 4);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 4);
                        TkMass.clear();
                        TkCharge.clear();
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K-pi-pi+pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[4] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 5);
                        TkMass.clear();
                        TkCharge.clear();
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+pi+pi-pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[5] == 1){
                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, -1, -1, ResIndex, false, 6);
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
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 7);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 7);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        ResIndex[0] = 0; ResIndex[1] = 2;
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 7);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 7);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        ResIndex[0] = 1; ResIndex[1] = 2;
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 7);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 7);
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
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 8);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 8);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(1); TkCharge.push_back(-1); TkCharge.push_back(-1);  
                        ResIndex[0] = 0; ResIndex[1] = 2;
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 8);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(KAON_MASS); TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 8);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(1); TkCharge.push_back(-1);  
                        ResIndex[0] = 1; ResIndex[1] = 2;
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 8);
                        TkMass.clear();
                        TkCharge.clear();

                        TkMass.push_back(PION_MASS); TkMass.push_back(KAON_MASS); TkMass.push_back(KAON_MASS);
                        TkCharge.push_back(-1); TkCharge.push_back(-1); TkCharge.push_back(1);  
                        Dfinder::BranchOutNTk( DInfo, input_tracks, isNeededTrack, D_counter, mass_window, TkMass, TkCharge, PHI_MASS, 0.1, ResIndex, false, 8);
                        TkMass.clear();
                        TkCharge.clear();
                    }
                    printf("D_counter: ");
                    for(unsigned int i = 0; i < Dchannel_.size(); i++){
                        printf("%d/", D_counter[i]);
                    }
                    printf("\n");//}}}
                    //}}}
                    //printf("-----*****DEBUG:End of DInfo.\n");

                    // TrackInfo section {{{
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
                        //
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
}

// ------------ method called once each job just after ending the event loop  ------------
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
}

bool Dfinder::GetAncestor(const reco::Candidate* p, int PDGprefix)
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
void Dfinder::BranchOutNTk(//input 2~4 tracks
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
            if(chi2_prob_tktkRes < VtxChiProbCut_) continue;
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
        if(chi2_prob_tktk < VtxChiProbCut_) continue;
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
        //document the mass hypothesis
        if( fabs(tktk_tk1_4vec.Mag()-PION_MASS) < fabs(tktk_tk1_4vec.Mag()-KAON_MASS) ) DInfo.rftk1_MassHypo[DInfo.size] = 211;
        else DInfo.rftk1_MassHypo[DInfo.size] = 321;
        if( fabs(tktk_tk2_4vec.Mag()-PION_MASS) < fabs(tktk_tk2_4vec.Mag()-KAON_MASS) ) DInfo.rftk2_MassHypo[DInfo.size] = 211;
        else DInfo.rftk2_MassHypo[DInfo.size] = 321;

        if(TkMass.size()>2){
            DInfo.rftk3_px[DInfo.size]    = tktk_tk3_4vec.Px();
            DInfo.rftk3_py[DInfo.size]    = tktk_tk3_4vec.Py();
            DInfo.rftk3_pz[DInfo.size]    = tktk_tk3_4vec.Pz();
            DInfo.rftk3_pt[DInfo.size]    = tktk_tk3_4vec.Pt();
            DInfo.rftk3_eta[DInfo.size]   = tktk_tk3_4vec.Eta();
            DInfo.rftk3_phi[DInfo.size]   = tktk_tk3_4vec.Phi();
            DInfo.rftk3_index[DInfo.size] = -selectedTkhidxSet[i][2]-1;
            if( fabs(tktk_tk3_4vec.Mag()-PION_MASS) < fabs(tktk_tk3_4vec.Mag()-KAON_MASS) ) DInfo.rftk3_MassHypo[DInfo.size] = 211;
            else DInfo.rftk3_MassHypo[DInfo.size] = 321;
        }
        if(TkMass.size()>3){
            DInfo.rftk4_px[DInfo.size]    = tktk_tk4_4vec.Px();
            DInfo.rftk4_py[DInfo.size]    = tktk_tk4_4vec.Py();
            DInfo.rftk4_pz[DInfo.size]    = tktk_tk4_4vec.Pz();
            DInfo.rftk4_pt[DInfo.size]    = tktk_tk4_4vec.Pt();
            DInfo.rftk4_eta[DInfo.size]   = tktk_tk4_4vec.Eta();
            DInfo.rftk4_phi[DInfo.size]   = tktk_tk4_4vec.Phi();
            DInfo.rftk4_index[DInfo.size] = -selectedTkhidxSet[i][3]-1;
            if( fabs(tktk_tk4_4vec.Mag()-PION_MASS) < fabs(tktk_tk4_4vec.Mag()-KAON_MASS) ) DInfo.rftk4_MassHypo[DInfo.size] = 211;
            else DInfo.rftk4_MassHypo[DInfo.size] = 321;
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
