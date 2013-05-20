// -*- C++ -*-
//
// Package:    AdvLeptonFilter
// Class:      AdvLeptonFilter
// 
/**\class AdvLeptonFilter AdvLeptonFilter.cc MMarionneau/AdvLeptonFilter/src/AdvLeptonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthieu Pierre Marionneau,8 R-019,+41227675765,
//         Created:  Tue Aug  7 11:56:39 CEST 2012
// $Id: AdvLeptonFilter.cc,v 1.6 2013/05/10 15:09:03 mmarionn Exp $
//
//
#include "MMarionneau/AdvLeptonFilter/interface/AdvLeptonFilter.h"


using namespace edm;
using namespace std;
using namespace cms;



AdvLeptonFilter::AdvLeptonFilter(const edm::ParameterSet& iConfig):
  isolator(),combU()
{

  LeptonCollection_[0]      = iConfig.getParameter<edm::InputTag>("ElectronInput");
  LeptonCollection_[1]      = iConfig.getParameter<edm::InputTag>("MuonInput");
  LeptonCollection_[2]      = iConfig.getParameter<edm::InputTag>("TauInput");
  PhotonCollection_      = iConfig.getParameter<edm::InputTag>("PhotonInput");
//   ZllCollection_[0]      = iConfig.getParameter<edm::InputTag>("ZeeInput");
//   ZllCollection_[1]      = iConfig.getParameter<edm::InputTag>("ZmumuInput");
  PFJetCollection_       = iConfig.getParameter<edm::InputTag>("PFJetInput");
  MetCollection_         = iConfig.getParameter<edm::InputTag>("MetInput");

  sels= iConfig.getUntrackedParameter<vector<string> >("Selections");

  //photon isoaltor initialisation
  isolator.initializePhotonIsolation(kTRUE);
  isolator. setConeSize(0.3);
  phEACh.resize(7);
  phEANe.resize(7);
  phEAPh.resize(7);
  phEACh[0] = 0.012;  phEANe[0] = 0.030; phEAPh[0] = 0.148;
  phEACh[1] = 0.010;  phEANe[1] = 0.057; phEAPh[1] = 0.130;
  phEACh[2] = 0.014;  phEANe[2] = 0.039; phEAPh[2] = 0.112;
  phEACh[3] = 0.012;  phEANe[3] = 0.015; phEAPh[3] = 0.216;
  phEACh[4] = 0.016;  phEANe[4] = 0.024; phEAPh[4] = 0.262;
  phEACh[5] = 0.020;  phEANe[5] = 0.039; phEAPh[5] = 0.260;
  phEACh[6] = 0.012;  phEANe[6] = 0.072; phEAPh[6] = 0.266;

  
  //Parser
  
  parser(sels);
  _selMaps = _selections;

  //Print Selections
  cout<<" =============================== print selections"<<endl;
  for(size_t i=0;i<_selections.size();i++) {
    cout<<" selection : ";

    string selName;
    for(idMIter itS=_selections[i].begin();itS!=_selections[i].end();itS++) {
      cout<<"   "<<(*itS).first<<" x "<<(*itS).second<<"     ";
      
      ostringstream os;
      os << (*itS).second;
      selName += os.str()+(*itS).first+"_";
    }

    counter.push_back( pair<string, int>( selName, 0) );
    cout<<endl;
  }

  //initialize counters
  

}


AdvLeptonFilter::~AdvLeptonFilter()
{ 
}

bool
AdvLeptonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  Nevt_++;
   using namespace edm;
   
 
   iEvent.getByLabel( LeptonCollection_[0], e_h );
   iEvent.getByLabel( LeptonCollection_[1], m_h );
   iEvent.getByLabel( LeptonCollection_[2], t_h );
   iEvent.getByLabel( PhotonCollection_, p_h );
   iEvent.getByLabel( PFJetCollection_, jets_h );
   iEvent.getByLabel( MetCollection_, met_h );

//    edm::Handle< reco::CompositeCandidateCollection >  zee_h;
//    iEvent.getByLabel( ZllCollection_[0], zee_h );
   
//    edm::Handle< reco::CompositeCandidateCollection >  zmm_h;
//    iEvent.getByLabel( ZllCollection_[1], zmm_h );
   
  
   
   //For identification ================================
   
   //Beam spot position
   iEvent.getByLabel( edm::InputTag("offlineBeamSpot"), beamSpot_h_ );
   
   //conversion
   iEvent.getByLabel( edm::InputTag("conversions"), conversions_h_ ); 
   
   //  Vertices
   iEvent.getByLabel( edm::InputTag("offlinePrimaryVertices"), vertices_h_ );

   //rho
   iEvent.getByLabel( edm::InputTag("kt6PFJets:rho"),rho_h_);

   //
   iEvent.getByLabel( edm::InputTag("particleFlow"),pfcs_h_);

   //===================================================


   //variables for DR protection with taus, jets
   vector<TVector3> els;
   vector<TVector3> mus;
   vector<TVector3> taus;
   vector<TVector3> phs;
   TVector3 tmp(0,0,0);


   for( std::vector<pat::Electron>::const_iterator it = e_h->begin(); it != e_h->end(); ++it) {

     //pat::ElectronRef Electron(e_h,ie);
     
     string id = elID(it);
     
     SetCateg("e",(int)it->pt(), id);

     if( it->pt() > 10 ) {
       tmp.SetPtEtaPhi( it->pt(), it->eta(), it->phi() ); 
       els.push_back(tmp);
     }
   }
     
   for(std::vector<pat::Muon>::const_iterator itm=m_h->begin();
       itm !=m_h->end(); ++itm)
     {
       string id= muID( itm );
       
       SetCateg("m",(int)itm->pt(), id);

       if( itm->pt() > 10 ) {
	 tmp.SetPtEtaPhi( itm->pt(), itm->eta(), itm->phi() ); 
	 mus.push_back(tmp);
       }

     }

   //Tau with DR protection against muon and electrons =========================
   TVector3 tauV3(0,0,0);
   for( size_t ii=0; ii<t_h->size(); ii++)
     {
       pat::TauRef Tau( t_h, ii ); 

       tauV3.SetPtEtaPhi(Tau->pt(), Tau->eta(), Tau->phi() );
       bool match = false;

       //if(Tau->pt()>10)
       // cout<<" ====> tau? "<<Tau->pt()<<endl;

       for(size_t ie=0;ie<els.size();ie++) {
	 if(els[ie].DrEtaPhi( tauV3) < 0.3 )
	   {match =true; break;}
       }
       if(match) continue;
       for(size_t im=0;im<mus.size();im++) {
	 if(mus[im].DrEtaPhi( tauV3) < 0.3 )
	   {match =true; break;}
       }
       if(match) continue;
       
     

       string id= tauID( Tau );
     
       //       if(Tau->pt()>10)
       //	 cout<<" ----> tau! "<<Tau->pt()<<"  "<<id<<endl;
  
       SetCateg("t",(int)Tau->pt(), id);

       if( Tau->pt() > 10 ) {
	 tmp.SetPtEtaPhi( Tau->pt(), Tau->eta(), Tau->phi() ); 
	 taus.push_back(tmp);
       }


     }

   
   //Photon==========
   TVector3 phV3(0,0,0);
   for(std::vector<pat::Photon>::const_iterator itp=p_h->begin();
       itp !=p_h->end(); ++itp)
     {
       //const reco::Photon& Photon = *itp; 

       if(itp->pt()< 20 ) continue;

       phV3.SetPtEtaPhi(itp->pt(), itp->eta(), itp->phi() );
       bool match = false;

       for(size_t ie=0;ie<els.size();ie++) {
	 if(els[ie].DrEtaPhi( phV3) < 0.3 )
	   {match =true; break;}
       }
       if(match) continue;

       string id = photonID( itp );

       SetCateg("p",(int)itp->pt(), id);

       tmp.SetPtEtaPhi( itp->pt(), itp->eta(), itp->phi() ); 
       phs.push_back(tmp);

     }

   //Jets==========
   TVector3 jetV3(0,0,0);
   for(size_t ijet=0; ijet<jets_h->size();ijet++)
     {
       pat::JetRef jet(jets_h, ijet); 

       if(jet->pt()< 20 ) continue;

       jetV3.SetPtEtaPhi(jet->pt(), jet->eta(), jet->phi() );
       bool match = false;
       for(size_t ie=0;ie<els.size();ie++) {
	 if(els[ie].DrEtaPhi( jetV3) < 0.3 )
	   {match =true; break;}
       }
       if(match) continue;
       for(size_t im=0;im<mus.size();im++) {
	 if(mus[im].DrEtaPhi( jetV3) < 0.3 )
	   {match =true; break;}
       }
       if(match) continue;
       for(size_t it=0;it<taus.size();it++) {
	 if(taus[it].DrEtaPhi( jetV3) < 0.3 )
	   {match =true; break;}
       }
       if(match) continue; 
       for(size_t ip=0;ip<phs.size();ip++) {
	 if(phs[ip].DrEtaPhi( jetV3) < 0.3 )
	   {match =true; break;}
       }
       if(match) continue;

       //FIXME
       //if(jet->pt()>50 ) cout<<" SUPER JET ============= > "<<jet->pt()<<endl;

       string id = jetID( jet );

       SetCateg("j",(int)jet->pt(), id);

     }

   //MET===========================
   double pfMet = met_h->begin()->et();
   
   SetCateg("h",(int)pfMet, "T");

  
   //===========================================

   //Decision
   if( acceptEvent() )
     return true;
   else
     {
       NFailed_++;
       return false;
     }
}


void 
AdvLeptonFilter::beginJob()
{
 NFailed_=0;Nevt_=0;

 NSingleEl=0;
 NSingleMu=0;
 NSingleTau=0;
 NDoubleEl=0;
 NDoubleMu=0;
 NDoubleTau=0;
 NXElMu=0;
 NXElTau=0;
 NXMuTau=0;
 
}


void 
AdvLeptonFilter::endJob() {
  cout<<" AdvLeptonFilter : number of event failed -> "<<NFailed_<<" on local total number ->"<<Nevt_<<endl;

  cout<<" Detail of selection over "<<Nevt_-NFailed_<<" events (double counting allowed) :"<<endl;
  for(size_t ic=0;ic<counter.size();ic++) {
    cout<<counter[ic].first<<" \t ---> "<<counter[ic].second<<endl;
  }


 //  cout<<" Detail : "<<endl;
//   cout<<"\t --> SingleEl  : "<<NSingleEl<<endl;
//   cout<<"\t --> SingleMu  : "<<NSingleMu<<endl;
//   cout<<"\t --> SingleTau : "<<NSingleTau<<endl;
//   cout<<"\t --> DoubleEl  : "<<NDoubleEl<<endl;
//   cout<<"\t --> DoubleMu  : "<<NDoubleMu<<endl;
//   cout<<"\t --> DoubleTau : "<<NDoubleTau<<endl;
//   cout<<"\t --> XElMu     : "<<NXElMu<<endl;
//   cout<<"\t --> XElTau    : "<<NXElTau<<endl;
//   cout<<"\t --> XMuTau    : "<<NXMuTau<<endl;

}



//IDs ==========================

string 
AdvLeptonFilter::muID(std::vector<pat::Muon>::const_iterator mu) {

  string id="N";

  if( !mu->isGlobalMuon() ) return id;
  //if( !mu->isTrackerMuon() ) return id;
  if( !mu->isPFMuon() ) return id;

  id="V";
  
  reco::TrackRef track_ = mu->globalTrack();
  float dxy = (float)( track_->dxy(beamSpot_h_->position()));
  int trackerHits = track_->hitPattern().trackerLayersWithMeasurement();
  int pixelHits = track_->hitPattern().numberOfValidPixelHits();
  float chi2=mu->globalTrack()->normalizedChi2();
  int muonHits = mu->globalTrack()->hitPattern().numberOfValidMuonHits(); 
  int nStationMatches =mu->numberOfMatchedStations();
  
  if( chi2 > 10 ) return id;
  if( trackerHits < 5) return id;

  id="L";


  if( muonHits<1 ) return id;
  if( nStationMatches<2 ) return id;
  if( pixelHits<1) return id;
  
  id="M";
    
  if( dxy> 1.) return id;
  
  id="T";

  return id;
   
}

string 
AdvLeptonFilter::elID(std::vector<pat::Electron>::const_iterator el) {

  string id="N";

  float chiso = el->chargedHadronIso();
  float nhiso = el->neutralHadronIso();
  float phiso = el->photonIso();


  if(!EgammaCutBasedEleId::PassWP( EgammaCutBasedEleId::VETO, *el , conversions_h_,
				   *(beamSpot_h_.product()), vertices_h_,
				   chiso,phiso,nhiso, *(rho_h_.product()) ) ) {
    return id;
  }
  
  id="V";
  if(!EgammaCutBasedEleId::PassWP( EgammaCutBasedEleId::LOOSE, *el , conversions_h_,
				   *(beamSpot_h_.product()), vertices_h_,
				   chiso,phiso,nhiso, *(rho_h_.product()) ) ) {
    return id;
  }
  id="L";
  if(!EgammaCutBasedEleId::PassWP( EgammaCutBasedEleId::MEDIUM, *el , conversions_h_,
				   *(beamSpot_h_.product()), vertices_h_,
				  chiso,phiso,nhiso, *(rho_h_.product()) ) ) {
    return id;
  }
  id="M";
  if(!EgammaCutBasedEleId::PassWP( EgammaCutBasedEleId::TIGHT, *el , conversions_h_,
				   *(beamSpot_h_.product()), vertices_h_,
				   chiso,phiso,nhiso, *(rho_h_.product()) ) ) {
    return id;
  }
  
  id="T";
  
  return id;
    
}

string 
AdvLeptonFilter::tauID(const pat::TauRef tau) {

  string id="N";

  if(tau->tauID("byVLooseCombinedIsolationDeltaBetaCorr") ) {
    id="V";
  }
  if(tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")) {
    id="L";
  }
  if(tau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")) {
    id="M";
  }
  if(tau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") && 
     tau->tauID("againstMuonLoose") &&
     tau->tauID("againstElectronLoose") ) {
    id="T";
  }
  
  return id;
  
}


string 
AdvLeptonFilter::photonID(std::vector<pat::Photon>::const_iterator ph) {

  string id="N";
  
  float sieie = ph->sigmaIetaIeta();
  float hoe = ph->hcalTowerSumEtConeDR04() + (ph->hadronicOverEm() - ph->hadTowOverEm())*ph->superCluster()->energy()/cosh(ph->superCluster()->eta());
  
  reco::VertexRef vtxRef( vertices_h_, 0 );
  isolator.fGetIsolation( &(*ph),pfcs_h_.product(), vtxRef, vertices_h_);

  int etaB=findPhEtaBin( ph->eta() );
  float eaCh = phEACh[ etaB ];
  float eaPh = phEAPh[ etaB ];
  float eaNe = phEANe[ etaB ];
  float rho = *(rho_h_.product());

  float chiso = max(isolator.getIsolationCharged()-eaCh*rho,(float)0.);
  float phiso = max(isolator.getIsolationPhoton()-eaPh*rho,(float)0.);
  float neiso = max(isolator.getIsolationNeutral()-eaNe*rho,(float)0.);
  float conv = true;//!ConversionTools::hasMatchedPromptElectron(ph->superCluster(), e_h, conversions_h_, (beamSpot_h_.product())->position());

  if(fabs( ph->eta() ) <1.5) {

    if(hoe < 0.05 && sieie<0.012 && conv ) id = "V";
    if(hoe < 0.05 && sieie<0.012 && chiso<2.6 && conv && 
       neiso<3.5+0.04*ph->pt() && phiso<1.3+0.005*ph->pt()) id = "L";
    if(hoe < 0.05 && sieie<0.011 && chiso<1.5 && conv && 
       neiso<1.0+0.04*ph->pt() && phiso<0.7+0.005*ph->pt()) id = "M";
    if(hoe < 0.05 && sieie<0.012 && chiso<0.7 && conv && 
       neiso<0.4+0.04*ph->pt() && phiso<0.5+0.005*ph->pt()) id = "T";
  }
  else {
    if(hoe < 0.05  && sieie<0.034 && conv ) id = "V";
    if(hoe < 0.05 && sieie<0.034 && chiso<2.3 && conv && 
       neiso<2.9+0.04*ph->pt() ) id = "L";
    if(hoe < 0.05 && sieie<0.033 && chiso<1.2 && conv && 
       neiso<1.5+0.04*ph->pt() && phiso<1.0+0.005*ph->pt()) id = "M";
    if(hoe < 0.05 && sieie<0.031 && chiso<0.5 && conv && 
       neiso<1.5+0.04*ph->pt() && phiso<1.0+0.005*ph->pt()) id = "T";
  }


  return id;
  
}

int AdvLeptonFilter::findPhEtaBin( float eta ) {

  if( fabs(eta)<1.0 ) return 0; 
  if( 1.0<fabs(eta) && fabs(eta)<1.479 ) return 1;
  if( 1.479<fabs(eta) && fabs(eta)<2.0 ) return 2;
  if( 2.0<fabs(eta) && fabs(eta)<2.2 ) return 3;
  if( 2.2<fabs(eta) && fabs(eta)<2.3 ) return 4;
  if( 2.3<fabs(eta) && fabs(eta)<2.4 ) return 5;
  if( fabs(eta)>2.4 ) return 6;

  return 0;
}

  
string 
AdvLeptonFilter::jetID(const pat::JetRef jet) {

  string id="V"; //no very loose requirement for jets, enabled by default

  float chFrac = jet->chargedHadronEnergyFraction();
  float nhFrac = jet->neutralHadronEnergyFraction();
  float nhEMFrac = jet->neutralEmEnergyFraction();
  float chEMFrac = jet->chargedEmEnergyFraction();
  int nConst = jet->getPFConstituents().size();
  float chMult = jet->chargedMultiplicity();

  
  if( nhFrac <0.99 && nhEMFrac < 0.99 && nConst > 1 ) id = "L" ;
  if( nhFrac <0.95 && nhEMFrac < 0.95 && nConst > 1 ) id = "M" ;
  if( nhFrac <0.90 && nhEMFrac < 0.90 && nConst > 1 ) id = "T" ;

  if(fabs(jet->eta()) > 2.4 && !(chFrac > 0 && chMult > 0 && chEMFrac<0.99) ) {
    //not identified as valid jet for endcaps
    id = "V";
  }

  return id;
}



//// parsering and selection======================
void  
AdvLeptonFilter::parser(vector<string> sels) {

  //double thresholds[10]={3,5,7,10,15,17,20,25,30,50};

  for(unsigned int is=0;is<sels.size();is++) {
    
    vector< pair<string, int> > blocks;
    
    string sel = sels[is];
    
    string block;
    bool end=false,f=true;
    size_t pos=0,pos1=0;
    int nb=0,nB=0;
    string key;

    vector<string> id;
    vector<vector<idMap> > ids;
    

    //cout<<"============================= line : "<<sel<<endl;

    bool invC=false;

    while(!end) {
    
      if(nB==nb) {
	f=true;
	nb=0;
	id.clear();
      }
      
      pos = sel.find( "_", pos1+1);

      if(pos1==0) {
	block = sel.substr(0,pos);
      }
      else {
	if(pos==(size_t)-1)
	  block = sel.substr(pos1+1 );
	else
	  block = sel.substr(pos1+1, pos-pos1-1 );
      }
      pos1 = pos;
      
      if(f) {
	if(block.substr(0,1)=="!")
	  { nB = 1; invC=true;}
	else
	  nB = atoi(block.substr(0,1).c_str());
	key = block.substr(1,1);
	f=false;
      }
      else {

	if(block.size()==2) {
	  
	  if(invC) {
	    id.push_back( block+"!" );
	    invC=false;
	  }
	  else
	    id.push_back( block+"N" );

	  bool reg=false;
	  for(size_t it=0;it<thrs.size();it++)
	    if(thrs[it] == atoi(block.c_str() ) )
	      reg=true;
	  
	  if(!reg)
	    thrs.push_back( atoi(block.c_str() ) );
	}
	if(block.size()==3) {
	  id.push_back( block );

	  bool reg=false;
	  for(size_t it=0;it<thrs.size();it++)
	    if(thrs[it] == atoi(block.substr(0,2).c_str() ) )
	      reg=true;
	  
	  if(!reg)
	    thrs.push_back( atoi(block.substr(0,2).c_str() ) );
	}
	
	if(nb<nB) {
	  nb++;
	}

	if(nb==nB || pos==(size_t)-1) {//time to build the selection
	  //cout<<nB<<"   "<<key<<"  "<<id.size()<<endl;
	//   for(size_t i=0;i<id.size();i++)
// 	    cout<<" --> "<<id[i]<<endl;
	  
	  ids.push_back( combU.buildSelections(nB, key, id ) );

	  
	}
	
	if(pos==(size_t)-1) break;
      }
      
    }// end sel

    //combine the selections
    MergeSelections( ids );
    
  }//selections


  sort( thrs.begin(), thrs.end() );
  
}


void 
AdvLeptonFilter::SetCateg(string T, int pT, string Id) {

  vector<string> pdg;

  if(T.find("e")!=(size_t)-1) {
    pdg.push_back("e"); //pdg.push_back("l");pdg.push_back("L"); 
  }
  if(T.find("m")!=(size_t)-1) {
    pdg.push_back("m"); //pdg.push_back("l");pdg.push_back("L"); 
  }
  if(T.find("t")!=(size_t)-1) {
    pdg.push_back("t"); //pdg.push_back("l"); 
  }
  if(T.find("h")!=(size_t)-1) {
    pdg.push_back("h");
  }
  if(T.find("g")!=(size_t)-1) {
    pdg.push_back("g");
  }
  if(T.find("j")!=(size_t)-1) {
    pdg.push_back("j");
  }

  vector<string> pTs;
  for(unsigned int i=0;i<thrs.size();i++) {
    if(thrs[i]<pT) {
     
      ostringstream os;
      os << thrs[i];
      pTs.push_back( os.str() );
    }
    else break;
  }

  vector<string> ids;
  if(Id=="T") {
    ids.push_back("N"); ids.push_back("V"); ids.push_back("L");
    ids.push_back("M");
    ids.push_back("T");
  }
  else if(Id=="M") {
     ids.push_back("N"); ids.push_back("V");
     ids.push_back("L");
    ids.push_back("M"); 
  }
  else if(Id=="L") {
     ids.push_back("N"); ids.push_back("V");
    ids.push_back("L");
  }
  else if(Id=="V") {
     ids.push_back("N");
    ids.push_back("V");
  }
  else
    ids.push_back("N");

  ids.push_back("!");
  

  string key;
  for(size_t ip=0;ip<pdg.size();ip++) {
    for(size_t im=0;im<pTs.size();im++) {
      for(size_t ii=0;ii<ids.size();ii++) {
	key = pdg[ip]+pTs[im]+ids[ii];
	
	if(ids[ii]=="!") {
	  for(size_t is=0;is<_selMaps.size();is++) {
	    iter = _selMaps[ is ].find(key);
	    
	    if(iter != _selMaps[ is ].end() ) {
	      
	      _selMaps[ is ][ (*iter).first ] += 1;
	    }
	  }
	}
	else {
	  for(size_t is=0;is<_selMaps.size();is++) {
	    iter = _selMaps[ is ].find(key);
	    
	    if(iter != _selMaps[ is ].end() ) {
	      
	      _selMaps[ is ][ (*iter).first ] -= 1;
	    }
	  }
	  
	}
      }
    }
  }


  
}

bool 
AdvLeptonFilter::acceptEvent() {

  bool acc=false;

  bool debug=false;

  for(size_t i=0;i<_selMaps.size();i++) {
    
    if(debug)
      cout<<"====== new map ===== "<<endl;

    string selName;
    bool accept=true;
    for(iter=_selMaps[i].begin();
	iter!=_selMaps[i].end();iter++) {
      
      if(debug)
	cout<<"\t --> "<<iter->first<<" = "<<iter->second<<endl;

      selName += (*iter).first+"_";
      if( (*iter).second > 0 ) {accept=false; break;}

    }

    if(debug)
      cout<<" map "<<selName<<" accepted ? "<<accept<<endl;

    if(accept) { acc=true; counter[ i ].second += 1; }
  }
 
  //reinit
  _selMaps.clear();
  _selMaps = _selections;

  return acc;

}


void AdvLeptonFilter::MergeSelections(vector<vector<idMap> > inSels) {

  assert(inSels.size()<=5);

  vector<idMap> allSels;
  idMap mMap;
  
  vector<idMap> sel1,sel2,sel3,sel4,sel5;
  
  
  sel1 = inSels[0];
  if(inSels.size()>1)
    sel2 = inSels[1];
  if(inSels.size()>2)
    sel3 = inSels[2];
  if(inSels.size()>3)
    sel4 = inSels[3];
  if(inSels.size()>4) 
    sel5 = inSels[4];
  

  //dummy loop over maximum 5 kind of selections...
  //could be rewritten in a smarter way...
  for(size_t i1=0;i1<sel1.size();i1++) {
    
    if(sel2.size()!=0) {
      for(size_t i2=0;i2<sel2.size();i2++) {
	
	if(sel3.size()!=0) {
	  for(size_t i3=0;i3<sel3.size();i3++) {	
	    
	    if(sel4.size()!=0) {
	      for(size_t i4=0;i4<sel4.size();i4++) {	
		
		if(sel5.size()!=0) {
		  for(size_t i5=0;i5<sel5.size();i5++) {
		    mMap = MergeSels( sel1[i1], sel2[i2] );
		    mMap = MergeSels( mMap, sel3[i3] );
		    mMap = MergeSels( mMap, sel4[i4] );
		    mMap = MergeSels( mMap, sel5[i5] );
		    allSels.push_back( mMap );
		  }
		}
		else { //no sel 5
		  mMap = MergeSels( sel1[i1], sel2[i2] );
		  mMap = MergeSels( mMap, sel3[i3] );
		  mMap = MergeSels( mMap, sel4[i4] );
		  allSels.push_back( mMap );
		}
	      }
	    }
	    else {//no sel 4
	      mMap = MergeSels( sel1[i1], sel2[i2] );
	      mMap = MergeSels( mMap, sel3[i3] );
	      allSels.push_back( mMap );
	    }
	  }
	}
	else {//no sel 3
	  mMap = MergeSels( sel1[i1], sel2[i2] );
	  allSels.push_back( mMap );
	}
      }
    }
    else { //no sel 2
      allSels.push_back( sel1[i1] );
    }
  }

  for(int unsigned is=0;is<allSels.size();is++) {
    bool exists=false;
    for(int unsigned js=0;js<_selections.size();js++) {
      if(CombUtils::exist(_selections[js], allSels[is]) )
      	{ exists=true; break;}
    }
    if(!exists)
      _selections.push_back(allSels[is]);
  }

}


idMap AdvLeptonFilter::MergeSels( idMap s1, idMap s2) {

  idMap o = s1;
  idMIter iter,itS;
  
  for(iter = s2.begin();iter!=s2.end();iter++) {

    itS = s2.find( iter->first );

    if(iter->first.find("!")!=(size_t)-1) {
      o[ iter->first ]=0;
    }
    else {
      if(itS == s2.end() )  
	o[ iter->first ]=iter->second;
      else
	o[ iter->first ]+=iter->second;
    }
  }
  
  return o;
}






//define this as a plug-in
DEFINE_FWK_MODULE(AdvLeptonFilter);
