// vim:set ts=4 sw=4 fdm=marker et:
#ifndef _UTILITIES_H_
#define _UTILITIES_H_

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
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"//calculate trajectory distance

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

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
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"//proper covariance error calculation
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
#include "DataFormats/TrackReco/interface/DeDxData.h"

class CommonFuncts{//{{{
    public:
        void test(){
    }   
        
    bool GetAncestor(const reco::Candidate* p, int PDGprefix)
    {
        if(p->numberOfMothers()==0) return false;
        else{
            const reco::Candidate* MyMom = p->mother(0);
            int mpid = abs(MyMom->pdgId());
            if(abs(int(mpid/100) % 100) == PDGprefix) return true;
            else return GetAncestor(MyMom, PDGprefix);
        }
    }

	bool GetDescendant(const reco::Candidate* p, int PDGprefix)
	{
		bool Saveit = false;
		if( p->numberOfDaughters() == 0 ) return false; //not necessary, keep it to make clear
		const reco::Candidate * daughterparticle = NULL;
		for( unsigned int idau = 0; idau < p->numberOfDaughters(); idau++ )
		{
			daughterparticle = p->daughter(idau);
			if( abs( int(daughterparticle->pdgId()/100) % 100 ) == PDGprefix ) return true;
			if( daughterparticle == p )   continue; //protection
			Saveit = GetDescendant(daughterparticle, PDGprefix);
			if( Saveit )  return true;
		}
		return Saveit;
	}

    float getParticleSigma(double mass)
    {
        if(mass == ELECTRON_MASS)
            return 0.013E-9f;
        else if(mass == MUON_MASS)
            return 4E-9f;
        else if(mass == PION_MASS)
            return 3.5E-7f;
        else if(mass == KAON_MASS)
            return 1.6E-5f;
        else if(mass == PROTON_MASS)
            return 8E-8f;
        else
            return 1E-6;
    }

    double getMaxDoca(std::vector<RefCountedKinematicParticle> &kinParticles)
    {
        double maxDoca = -1.0;
        TwoTrackMinimumDistance md;
        std::vector<RefCountedKinematicParticle>::iterator in_it, out_it;
        for (out_it = kinParticles.begin(); out_it != kinParticles.end(); ++out_it) {
            for (in_it = out_it + 1; in_it != kinParticles.end(); ++in_it) {
                md.calculate((*out_it)->currentState().freeTrajectoryState(),(*in_it)->currentState().freeTrajectoryState());
                if (md.distance() > maxDoca)
                    maxDoca = md.distance();
            }
        }
        return maxDoca;
    }
};//}}}

#endif
