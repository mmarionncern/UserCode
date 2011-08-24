// -*- C++ -*-
//
// Package:    CleanMETProducer
// Class:      CleanMETProducer
// 
/**\class CleanMETProducer CleanMETProducer.cc MMarionneau/CleanMETProducer/src/CleanMETProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthieu Pierre Marionneau,27 2-005,+41227673174,
//         Created:  Tue Aug 23 09:36:19 CEST 2011
// $Id: CleanMETProducer.cc,v 1.1.1.1 2011/08/24 16:23:39 mmarionn Exp $
//
//


#include "MMarionneau/CleanMETProducer/interface/CleanMETProducer.h"

#include <iostream>
#include <stdlib.h>

using namespace edm;
using namespace std;

const double myPi       = acos(-1.);


CleanMETProducer::CleanMETProducer(const edm::ParameterSet& iConfig)
{
 //Input
  pfMETCollection_          = iConfig.getParameter<edm::InputTag>("pfMETInput");
  pfJetCollection_          = iConfig.getParameter<edm::InputTag>("pfJETInput");
  pfCandidateCollection_    = iConfig.getParameter<edm::InputTag>("pfCandidateInput");
  vertexCollection_         = iConfig.getParameter<edm::InputTag>("vertexInput");
  trackCollection_          = iConfig.getParameter<edm::InputTag>("trackInput");
  vtxObjectCollection_      = iConfig.getParameter<edm::InputTag>("vtxObjectInput");
  refObjectCollection_      = iConfig.getParameter<edm::InputTag>("refObjectInput");
 
  udvN = iConfig.getUntrackedParameter<int>("nVtxObjectInCol" ,-1);
  userDefVtx = iConfig.getUntrackedParameter<bool>("vtxUserDef" , false);

  nro = iConfig.getUntrackedParameter<int>("nRefObjectInCol" ,-1);

  //Thresholds and protections
  jPtHSThr  = iConfig.getUntrackedParameter<double>("jetPtHSThr");
  jBalHSThr = iConfig.getUntrackedParameter<double>("jetBalHSThr");
  jDirHSThr = iConfig.getUntrackedParameter<double>("jetDirHSThr");
  jMatchDR  = iConfig.getUntrackedParameter<double>("jetMatchDR");
  jMatchBal = iConfig.getUntrackedParameter<double>("jetMatchBal");

  //Output
  cleanMETName_ = iConfig.getParameter<std::string>("nameCleanMET");
  minMETName_ = iConfig.getParameter<std::string>("nameMinMET");

  produces<reco::METCollection>( cleanMETName_ ); 
  produces<reco::METCollection>( minMETName_ ); 

  debug_ = false;
  
}


CleanMETProducer::~CleanMETProducer()
{
 
}


void
CleanMETProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //loading objects
  loadObjects( iEvent );


  //User needs to define precisely the class of the reference object
  // UserDef
  const reco::Candidate* refObject = &(*refObject_h_)[nro];
  

  //MET computation
  reco::MET cleanMET = computeCleanMET( dynamic_cast<const reco::Candidate*>(refObject) );
  reco::MET minMET = computeMinMET( cleanMET );
  
  
  //MET writing in the event
  std::auto_ptr<reco::METCollection> cleanMETcoll;
  cleanMETcoll.reset(new reco::METCollection);
  cleanMETcoll->push_back( cleanMET ) ;
  iEvent.put( cleanMETcoll, cleanMETName_ );  
 
  std::auto_ptr<reco::METCollection> minMETcoll;
  minMETcoll.reset(new reco::METCollection);
  minMETcoll->push_back( minMET ) ;
  iEvent.put( minMETcoll, minMETName_ );  
 
  refresh();
}


void 
CleanMETProducer::beginJob()
{
}


void 
CleanMETProducer::endJob() {

}


void 
CleanMETProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}


void 
CleanMETProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}


void 
CleanMETProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}


void 
CleanMETProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CleanMETProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}






////
//// Main functions ############################################
////


const reco::Vertex
CleanMETProducer::defineHardScatterVertex() {

  if(userDefVtx) { 
    // The physicist controls the vertex determination
    // and have to change the class of the object used
    // as reference for vertex determination
   
    if(udvN != -1 && (size_t)udvN <vtxObject_h_->size() ) {
   
      // UserDef
      reco::CandidateRef cand(vtxObject_h_, udvN );
     // UserDef
      const reco::Vertex* vtmp = vertexMatchingBC<reco::CandidateRef>( cand );
      const reco::Vertex hsVtx = (*vtmp);
      return hsVtx;
    }
    else {

      // UserDef
      reco::CandidateRef cand(vtxObject_h_, 0 );
      // UserDef
      const reco::Vertex* vtmp = vertexMatchingBC<reco::CandidateRef>( cand );
      const reco::Vertex hsVtx = (*vtmp);
      return hsVtx;
    }
    
  }
 else {
   // By default, highest sumPt vertex
   const reco::Vertex hsVtx = ((*vertices_h_)[0]);
   return hsVtx;
 }  

}

template <class T>
bool CleanMETProducer::isHardScatterVertexBC(T cand, const reco::Vertex hsVtx) {

  //This function test the origin of a 
  // basic candidate (all candidates but not a jet)

  bool isMatch=false;

  const reco::Vertex* candVtx = vertexMatchingBC<T>( cand );

  if( candVtx!=NULL ) {

    if(debug_) {cout<<" PFC Vertex ==> "<<candVtx->position()<<"  <->  "<<hsVtx.position()<<"    ";}
    
    if( candVtx->position() == hsVtx.position() ) {
      isMatch=true; 
    }
  }  
  
  return isMatch;
}

bool CleanMETProducer::isHardScatterVertexJet(const reco::Jet* cand, const reco::Vertex hsVtx) {
 
  //This function test the origin of a
  // jet (jet but not another candidate)

  bool isMatch=false;

  const reco::Vertex* candVtx = vertexMatchingJet( cand );
 
  if( candVtx!=NULL ) {

    if(debug_) {cout<<cand->px()<<"  "<<cand->py()
		    <<" Jet Vertex ==> "<<candVtx->position()
		    <<"  <->  "<<hsVtx.position()<<" ->  ";}

    if( candVtx->position() == hsVtx.position() ) {
      if(debug_) cout<< " yes "<<endl;
      isMatch=true; 
    } 
    else{ if(debug_) cout<< " no "<<endl;}
  }
 
  return isMatch;
}



reco::MET
CleanMETProducer::computeCleanMET(const reco::Candidate* refObject) {

  //Compute the pileup cleaned MET
  // Needs to define first the hard scatter vertex
  // and after associate (and reduce) contributions
  // to pileup/hard scatter vertices

  TVector2 p2(0,0);
  
  reco::Vertex hsVtx = defineHardScatterVertex();
  if(debug_) cout<<" ===> Hard Scatter Vertex position "
		 <<cout<<" ==============> "<<hsVtx.position()<<endl;


  // First association jets (needed for indetification
  // of unclustered pf candidates
  p2 += jetPart( refObject, hsVtx );

  // and now pf candidates
  p2 += pfcPart( refObject, hsVtx );

  
  //Create the reco::MET candidate

  math::XYZTLorentzVector p4(0,0,0,0);
  
  p4.SetPxPyPzE(p2.X(), p2.Y(), 0, p2.Mod() );
  
  reco::MET cleanMET(0., p4, hsVtx.position() );

  return cleanMET;
}



reco::MET
CleanMETProducer::computeMinMET( reco::MET cleanMET) {

  // Get the minimum of the pfMET and the cleanMET


  if( !pfT1MET_h_.isValid() ) 
    { cout<<" Error, no jets loaded "<<endl; abort();}
 
  reco::MET pfT1MET = (*pfT1MET_h_)[0];

  if(debug_)  cout<<" ??? clean "<<cleanMET.pt()<<" < "
		  <<pfT1MET.pt()<<" pf -->";
 

 if(cleanMET.pt() < pfT1MET.pt() )
    { 
      if(debug_)  cout<<" yes "<<endl;
      return cleanMET;}
 else 
   {
     if(debug_)  cout<<" no "<<endl;
     return pfT1MET;
   }
 
}

TVector2
CleanMETProducer::jetPart(const reco::Candidate* refObject,  const reco::Vertex hsVtx) {

  //Compute the jet contribution to the cleaned MET
  // First loop on jet : make a list of clustered pf candidates
  // and look for a jet balanced to the reference object and
  // associated to hard scatter vertex
  // Second loop : associate each jet to vertices with the application
  // of protection described in 
  // https://indico.cern.ch/contributionDisplay.py?contribId=5&confId=145789
  // slide 7. 
  // Thresholds and protection parameters can be tuned in the config file

  TVector2 jpart(0,0);
  TVector2 tmp(0,0);

  if( !pfJets_h_.isValid() ) 
    { cout<<" Error, no jets loaded "<<endl; abort();}
  if( !vertices_h_.isValid() ) 
    { cout<<" Error, no vertex loaded "<<endl; abort();}
  
  vector<bool> vtxAssoc(pfJets_h_->size(), false);
  
  bool findjet=false;
  for(int unsigned ij=0; ij<pfJets_h_->size(); ij++ ) {

   pat::JetRef jet(pfJets_h_,ij);
 
   bool isHSVtx = isHardScatterVertexJet( dynamic_cast<const reco::Jet*>(jet.get() ), hsVtx );

   vtxAssoc[ij] = isHSVtx;

   //Match pfCandidates contained in one jet
   vector< reco::PFCandidatePtr > assocPfCs = jet->getPFConstituents();
   for(int unsigned ic=0;ic<assocPfCs.size();ic++) {
     clusPFCs.push_back( assocPfCs[ic] );
   }
   
   
   //Now check jet presence back to back wrt the main candidate
   if( isHSVtx && refObject->pt()> jPtHSThr &&
       jet->pt()/refObject->pt() > jBalHSThr &&
       fabs( dPhi( jet->phi(), refObject->phi() ) ) > jDirHSThr )
     {findjet=true;}
  }

  if(debug_) cout<<" Control clustered particles "<<clusPFCs.size()
		 <<"   -> jet found back-to-back with the main cand ? "
		 <<findjet<<endl;
  
 bool assocJet=false;
 for(int unsigned ij=0;ij<pfJets_h_->size(); ij++) {

   tmp.Set(0.,0.);

   pat::JetRef jet(pfJets_h_,ij);
   assocJet=false;
 
   tmp.Set(jet->px(), jet->py() );
 
   bool jetbal = !findjet && !vtxAssoc[ij]
                 && refObject->pt()> jPtHSThr 
                 && jet->pt()/refObject->pt() > jBalHSThr
                 && fabs( dPhi( jet->phi(), refObject->phi() ) ) > jDirHSThr;

    if( ( vtxAssoc[ij] || jet->pt()>jPtHSThr || jetbal ) )  
      { jpart -= tmp; assocJet=true;}
  
    if(refObject->numberOfDaughters() != 0 && !assocJet) {
      for(int unsigned i=0;i<refObject->numberOfDaughters(); i++) {
	if( ( dR( jet, refObject->daughter(i) )<jMatchDR && jet->pt()/refObject->daughter(i)->pt() > jMatchBal ) )
	  { jpart -= tmp; assocJet=true;}
      }
    }
  
    if(refObject->numberOfDaughters() == 0 && !assocJet) {
      if( ( dR( jet, refObject) <jMatchDR && jet->pt()/refObject->pt() >jMatchBal )  ) 
	{jpart -=  tmp; assocJet=true;}
    }
  
    if( !assocJet)
      {jpart -= tmp/(double)(vertices_h_->size());}
    
    if(debug_) {
      cout<<"\t jet :"<<ij<<"  "<<jet->pt()
	  <<"   "<<jet->eta()<<"   "<<jet->phi()
	  <<" -> "<<vtxAssoc[ij]<<"  "<<jetbal<<"  "
	  <<assocJet<<" & "<<tmp.X()<<"  "<<tmp.Y()<<" --> "
	  <<jpart.X()<<"  "<<jpart.Y()<<" (ref dau) "
	  <<refObject->numberOfDaughters()<<endl;
    }

  }
 
  return jpart;
}




TVector2
CleanMETProducer::pfcPart(const reco::Candidate* refObject, const reco::Vertex hsVtx) {

   //Compute the unclustered pf candidate contribution to the cleaned MET
  // Associate each pfcandidate to vertices with the application
  // of protection described in 
  // https://indico.cern.ch/contributionDisplay.py?contribId=5&confId=145789
  // slide 7. 
  // Thresholds and protection parameters can be tuned in the config file


  TVector2 pfpart(0,0);
  TVector2 tmp(0,0);

 if( !pfCandidates_h_.isValid() ) { cout<<" Error, no pfCandidates loaded "<<endl; abort();}
 if( !pfJets_h_.isValid() ) { cout<<" Error, no jets loaded "<<endl; abort();}
 if( !vertices_h_.isValid() ) { cout<<" Error, no vertex loaded "<<endl; abort();}

 bool assocpf=false; int cnt=0;
 for(int unsigned ip=0; ip<pfCandidates_h_->size(); ip++ ) {
 
   tmp.Set(0.,0.);

   reco::PFCandidateRef pfCRef( pfCandidates_h_, ip );
   assocpf=false;
   tmp.Set( pfCRef->px(), pfCRef->py() );

   //Is Contained in one jet ?
   if( !isUnClus( pfCRef->p4() ) ) { continue; }
   
 
   bool isHSVtx = isHardScatterVertexBC<reco::PFCandidateRef>( pfCRef, hsVtx );
 
   //Now MET
   if( ( isHSVtx && pfCRef->charge()!=0 ) )
     { pfpart -= tmp;assocpf=true;}
   if( refObject->numberOfDaughters() == 0 && !assocpf) {
     if( ( dR(pfCRef, refObject) <jMatchDR && pfCRef->pt()/refObject->pt() >jMatchBal )  ) 
       {pfpart -= tmp; assocpf=true;}
   }
   
   if(!assocpf)
     { pfpart -= tmp/(double)(vertices_h_->size());}
 
   if(debug_) { cout<<" \t "<<ip<<"  "<<assocpf<<"  "<<pfCRef->pdgId()<<"  "
		  <<tmp.X()<<"  "<<tmp.Y()<<"  ---> "
		  <<pfpart.X()<<"   "<<pfpart.Y()<<endl; 
   cnt++;
   }
   
 }
 
 if(debug_) cout<<" Number of pf candidates prcessed "<<cnt<<" / "<<pfCandidates_h_->size()<<endl;
 
    return pfpart;

}


////
//// Utility functions ############################################
////


template <class T>
const reco::Vertex*
CleanMETProducer::vertexMatchingBC( T cand) {

  // Get the 3D closest vertex from one candidate (not a jet)
  
  if( !vertices_h_.isValid() ) 
    {cout<<" Error, no vertex loaded "<<endl; abort();}


  const reco::Vertex* mVertex(0);
  
  float dX;
    float dY;
    float dZ;
    
    float dRtmp=10000;
    
    for(int unsigned iv=0;iv<vertices_h_->size();iv++) {
      
      reco::VertexRef vtx(vertices_h_, iv );
      
      dX =  cand->vx() - vtx->x();
      dY =  cand->vy() - vtx->y();
      dZ =  cand->vz() - vtx->z();
      
      float dR = sqrt(dX*dX + dY*dY + dZ*dZ);
      
      if(dR<dRtmp) {
	mVertex = &((*vertices_h_)[iv]);
	dRtmp = dR;
      }
    }
    return mVertex;
}


const reco::Vertex*
CleanMETProducer::vertexMatchingJet(const reco::Jet* cand) {

  // Get the best vertex for one jet using highest sumPT
  // track-vertex association
  
  const reco::Vertex* mVertex(0);

   
  if(fabs( cand->eta() )>3) {
    if(debug_) cout<<" High eta jet "<<endl;
    return mVertex;
  }    

    map< const reco::Vertex* , float> sumPtAssoc;
    map< const reco::Vertex* , float>::const_iterator itmap;
    int cnt=0;
    float dRtmp = 0.5; float mindr=100;
    //   reco::TrackRefVector trkList = cand->getTrackRefs();
    
    
    for( size_t ii=0; ii<tracks_h_->size(); ii++ )
      {
	reco::TrackRef trk(tracks_h_,ii);
	
	//matching
	float dr = dR( trk, cand);
	if(dr<mindr) mindr =dr;
	if(debug_) cout<<" Track "<<ii<<" -----> "<<dr<<endl;
	if(dr < dRtmp ) {
	  const reco::Vertex* trkVtx = vertexMatchingBC<reco::TrackRef>( trk );
	  if(debug_) cout<<" Track "<<ii<<" -----> Vtx "<<trkVtx<<"   "<<dr<<endl;
	  cnt++; 
	  itmap = sumPtAssoc.find( trkVtx );
	  if( itmap != sumPtAssoc.end() ) {
	    sumPtAssoc[ trkVtx ] += trk->pt();
	    
	  }
	  else {
	    sumPtAssoc[ trkVtx ] = trk->pt();
	  }
	}
      }
  
    float sPttmp=0;
    for(itmap = sumPtAssoc.begin(); itmap != sumPtAssoc.end(); itmap++) {
      if(debug_) cout<<" VtxTrack "<<(*itmap).second<<" -----> "<<(*itmap).first<<endl;
	
      if(sPttmp < (*itmap).second ) {
	sPttmp = (*itmap).second;
	mVertex =  (*itmap).first;
      }
    }
    if(debug_) cout<<" Selected Vertex "<<mVertex<<"   "<<mindr<<endl;
    return mVertex;
    
}
    
 
  



float
CleanMETProducer::dPhi( float phi1, float phi2 ) {

  float phi1_= phi( cos(phi1), sin(phi1) );
  float phi2_= phi( cos(phi2), sin(phi2) );
  float dphi_= phi1_-phi2_;
  if( dphi_ > myPi  ) dphi_-=2.*myPi;
  if( dphi_ < -myPi ) dphi_+=2.*myPi;
  return dphi_;

}

float 
CleanMETProducer::phi( float x, float y )
{
  float phi_ =atan2(y, x);
  return (phi_>=0) ?  phi_ : phi_ + 2*myPi;
}



template<class T, class C>
float CleanMETProducer::dR( T cand1, C cand2) {
  
  float dphi = dPhi( cand1->phi(), cand2->phi() );
  float deta = cand1->eta() - cand2->eta();

  return sqrt(dphi*dphi + deta*deta);
  
}

bool
CleanMETProducer::isUnClus( reco::Candidate::LorentzVector p4 ) {

  //Test the clusterization of a pf candidate

  for(int unsigned ip4=0; ip4<clusPFCs.size(); ip4++) {
    if( (clusPFCs[ip4]->p4()) == p4 ) {
      return false;
    }	
  }
  return true;
}




////
//// Loading functions ############################################
////


void
CleanMETProducer::loadObjects(const edm::Event& iEvent) {
  
  
  loadPfT1MET( iEvent );
  loadPfJets( iEvent );
  loadPfCands( iEvent );
  loadVertices( iEvent );
  loadTracks( iEvent );
 
  loadUserDefVtx( iEvent );
  loadUserDefObject( iEvent );
  
}


void
CleanMETProducer::loadPfT1MET(const edm::Event& iEvent) {
 
  if( pfT1MET_h_.isValid() ) return;
  iEvent.getByLabel( pfMETCollection_, pfT1MET_h_ );
 
}


void
CleanMETProducer::loadPfJets(const edm::Event& iEvent) {
  
   if( pfJets_h_.isValid() ) return;
  iEvent.getByLabel( pfJetCollection_, pfJets_h_ );

}

void
CleanMETProducer::loadPfCands(const edm::Event& iEvent) {

  if( pfCandidates_h_.isValid() ) return;
  iEvent.getByLabel( pfCandidateCollection_, pfCandidates_h_ );
  
}

void
CleanMETProducer::loadVertices(const edm::Event& iEvent) {

  if( vertices_h_.isValid() ) return;
  iEvent.getByLabel( vertexCollection_, vertices_h_ );
  
}

void
CleanMETProducer::loadTracks( const edm::Event& iEvent )
{
  if( tracks_h_.isValid() ) return;
  iEvent.getByLabel( trackCollection_, tracks_h_ );
}


void
CleanMETProducer::loadUserDefVtx(const edm::Event& iEvent) {
  
  if( vtxObject_h_.isValid() ) return;
  iEvent.getByLabel( vtxObjectCollection_, vtxObject_h_ );

}

void
CleanMETProducer::loadUserDefObject(const edm::Event& iEvent) {
  
  if( refObject_h_.isValid() ) return;
  iEvent.getByLabel( refObjectCollection_, refObject_h_ );

}



void
CleanMETProducer::refresh() {

  pfT1MET_h_.clear();
  pfJets_h_.clear();
  pfCandidates_h_.clear();
  vertices_h_.clear();
  tracks_h_.clear();
  vtxObject_h_.clear();
  refObject_h_.clear();
 
  clusPFCs.clear();
 
  
}





//define this as a plug-in
DEFINE_FWK_MODULE(CleanMETProducer);
