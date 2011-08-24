
// Package:    CleanMETProducer
// Class:      CleanMETProducer
// 
/**\class CleanMETProducer CleanMETProducer.hh MMarionneau/CleanMETProducer/src/CleanMETProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthieu Pierre Marionneau,27 2-005,+41227673174,
//         Created:  Tue Aug 23 09:36:19 CEST 2011
// $Id: CleanMETProducer.h,v 1.1.1.1 2011/08/24 16:23:39 mmarionn Exp $
//
//



// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//data formats
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <Math/VectorUtil.h>
#include <TMath.h>
#include <TVector2.h>

#include <TH1F.h>
#include <TFile.h>

class CleanMETProducer : public edm::EDProducer {
public:
  explicit CleanMETProducer(const edm::ParameterSet&);
  ~CleanMETProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  

////
//// Main functions ############################################
////
  
private:
  const reco::Vertex defineHardScatterVertex();
  
  //  reco::Vertex* vertexMatching(reco::CandidateRef cand);
  template<class T> bool isHardScatterVertexBC(T cand, const reco::Vertex hsVtx);
  bool isHardScatterVertexJet(const reco::Jet* cand, const reco::Vertex hsVtx);
  
  reco::MET computeCleanMET(const reco::Candidate* refObject);
  reco::MET computeMinMET(reco::MET cleanMET);
  
  TVector2 jetPart(const reco::Candidate* refObject, const reco::Vertex hsVtx);
  TVector2 pfcPart(const reco::Candidate* refObject, const reco::Vertex hsVtx);


////
//// Utility functions ############################################
////
  
private:
  template<class T> const reco::Vertex* vertexMatchingBC(T cand);
  const reco::Vertex* vertexMatchingJet(const reco::Jet* cand);
  float dPhi( float phi1, float phi2 );
  float phi( float x, float y );
  template<class T, class C> float dR( T cand1, C cand2);
  bool isUnClus( reco::Candidate::LorentzVector p4 );
  

////
//// Loading functions ############################################
////

private:
  void loadObjects(  const edm::Event& iEvent );
  void loadPfT1MET(  const edm::Event& iEvent );
  void loadPfJets(   const edm::Event& iEvent );
  void loadPfCands(  const edm::Event& iEvent );
  void loadVertices( const edm::Event& iEvent );
  void loadTracks( const edm::Event& iEvent );

  void loadUserDefVtx( const edm::Event& iEvent);
  void loadUserDefObject(const edm::Event& iEvent);
  void refresh();

////
//// InputTags and Handles
////

private:
  edm::InputTag pfMETCollection_;
  edm::InputTag pfJetCollection_;
  edm::InputTag pfCandidateCollection_;
  edm::InputTag vertexCollection_;
  edm::InputTag trackCollection_;
  edm::InputTag vtxObjectCollection_;
  edm::InputTag refObjectCollection_;

  edm::Handle< reco::METCollection >         pfT1MET_h_;
  edm::Handle< pat::JetCollection >          pfJets_h_;
  edm::Handle< reco::PFCandidateCollection > pfCandidates_h_;
  edm::Handle< reco::VertexCollection >      vertices_h_;
  edm::Handle< reco::TrackCollection >       tracks_h_;

  // UserDef
  edm::Handle< reco::CandidateCollection >   vtxObject_h_;
  edm::Handle< reco::CandidateCollection >   refObject_h_;

  std::string cleanMETName_;
  std::string minMETName_;


////
//// Objects
////

//  reco::Vertex* hsVtx;

  std::vector<reco::PFCandidatePtr> clusPFCs;

  int udvN;
  bool userDefVtx;
  int nro;
  bool debug_;

////
////Thresholds and protectiosn definitions
////

  double jPtHSThr;
  double jBalHSThr;
  double jDirHSThr;
  double jMatchDR;
  double jMatchBal;

};
