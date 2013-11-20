#ifndef _Xbuse_H_
#define _Xbuse_H_

#include "TLorentzVector.h"
#include "TTree.h"
#include <math.h>
#include <string>
#include <vector>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "Bfinder/Bfinder/interface/format.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"


//a namespace used for xb frame 
namespace xb{

using namespace std;
//define particle masses
#define PI_MASS 0.13957
#define MU_MASS 0.10566

//function to construct TLorentzVector with data (most common-seen format can use)
template< typename T >
void FillLVec(const T* particle, TLorentzVector& vec, double mass)
{
        vec.SetPxPyPzE(particle->px(), particle->py(), particle->pz(), sqrt(particle->p()*particle->p()+mass*mass));
}

//function to calculate deltaR between 2 particles
template< typename T1 , typename T2>
double dR(const T1& p1, const T2& p2)
{
        return sqrt((p1.eta()-p2.eta())*(p1.eta()-p2.eta()) + (p1.phi()-p2.phi())*(p1.phi()-p2.phi()));
}

//vector of TLorentzVector
//typedef std::vector<TLorentzVector> TLVCollection;


}




#endif
