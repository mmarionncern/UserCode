// -*- C++ -*-
//
// Package:    AdvLeptonFilter
// Class:      AdvLeptonFilter
// 
/**\class AdvLeptonFilter AdvLeptonFilter.h MMarionneau/AdvLeptonFilter/src/AdvLeptonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthieu Pierre Marionneau,8 R-019,+41227675765,
//         Created:  Tue Aug  7 11:56:39 CEST 2012
// $Id: AdvLeptonFilter.h,v 1.1 2012/08/07 10:51:51 mmarionn Exp $
//
//



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "PhysicsTools/CandUtils/interface/CandCombiner.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/METReco/interface/PFMET.h" 
#include "DataFormats/METReco/interface/PFMETFwd.h" 

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TVector3.h"

//#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#ifndef __CombUtils
#include "MMarionneau/AdvLeptonFilter/utils/CombUtils.cc"
#endif

//
// class declaration
//

using namespace edm;
using namespace std;
using namespace cms;

typedef map<string, int> idMap;
typedef map<string, int>::iterator idMIter;

class AdvLeptonFilter : public edm::EDFilter {
   public:
      explicit AdvLeptonFilter(const edm::ParameterSet&);
      ~AdvLeptonFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
  InputTag LeptonCollection_[3];
  InputTag ZllCollection_[2];
  InputTag PhotonCollection_;
  InputTag PFJetCollection_;

  InputTag MetCollection_;

  int NFailed_,Nevt_;
  
  //counters
  int NSingleEl;
  int NSingleMu;
  int NSingleTau;
  int NDoubleEl;
  int NDoubleMu;
  int NDoubleTau;
  int NXElMu;
  int NXElTau;
  int NXMuTau;


  double ePt_, mPt_,tPt_,pPt_;
  vector<double> dePt_,dmPt_,dpPt_,dtPt_;
  vector<double> XemPt_,XetPt_,XmtPt_;

  string muID(std::vector<pat::Muon>::const_iterator mu);
  string elID(std::vector<pat::Electron>::const_iterator it);

  string tauID(const pat::TauRef tau);
  string photonID(std::vector<pat::Photon>::const_iterator ph);
  string jetID(const pat::JetRef jet);

  int findPhEtaBin( float eta );

  void  parser(vector<string> sels);
  void SetCateg(string T, int pT, string Id);
  bool acceptEvent();

  vector<string> sels;

  vector< idMap > _selMaps;
  vector< idMap > _selections;
  map<string, int>::const_iterator iter;

  vector<pair<string,int> > counter;

  vector<int> thrs;
  
  edm::Handle< std::vector<pat::Electron> > e_h;
  edm::Handle< pat::MuonCollection  >      m_h;
  edm::Handle< pat::TauCollection  >      t_h;
  edm::Handle< pat::PhotonCollection  >      p_h;
  edm::Handle< pat::JetCollection  >      jets_h;
  edm::Handle< std::vector<pat::MET> > met_h;


  edm::Handle<reco::BeamSpot> beamSpot_h_;
  edm::Handle<reco::ConversionCollection> conversions_h_;
  edm::Handle<reco::VertexCollection> vertices_h_;
  edm::Handle<double> rho_h_;
  edm::Handle< reco::PFCandidateCollection > pfcs_h_;

  //photon isolator
  PFIsolationEstimator isolator;

  vector<float> phEACh;
  vector<float> phEANe;
  vector<float> phEAPh;
  
  //Combination
  CombUtils combU;

  void MergeSelections(vector<vector<idMap> > inSels);
  idMap MergeSels( idMap s1, idMap s2);

};
