// -*- C++ -*-
//
// Package:    JetViaTrigPrim
// Class:      JetViaTrigPrim
// 
/**\class JetViaTrigPrim JetViaTrigPrim.cc MMarionneau/JetViaTrigPrim/src/JetViaTrigPrim.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Matthieu Pierre Marionneau
//         Created:  Fri Jul  3 09:12:46 CEST 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "MMarionneau/JetViaTrigPrim/interface/JetViaTrigPrim.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Geometry
#if (CMSSW_COMPAT_VERSION>=210)
#   include "Geometry/Records/interface/CaloGeometryRecord.h"
/**/typedef CaloGeometryRecord MyCaloGeometryRecord;
#else
#   include "Geometry/Records/interface/IdealGeometryRecord.h"
/**/typedef IdealGeometryRecord MyCaloGeometryRecord;
#endif


using namespace cms;
using namespace edm;
using namespace std;


JetViaTrigPrim::JetViaTrigPrim(const edm::ParameterSet& ps):
  collNotFoundWarn_(ps.getUntrackedParameter<bool>("warnIfCollectionNotFound", true)),
  ebDigis_(ps.getParameter<edm::InputTag>("EBDigiCollection"), false,
             collNotFoundWarn_),
  ebSrFlags_(ps.getParameter<edm::InputTag>("EbSrFlagCollection"), false,
             collNotFoundWarn_),
  eeSrFlags_(ps.getParameter<edm::InputTag>("EeSrFlagCollection"), false,
             collNotFoundWarn_),
  triggerTowerMap_(0)
{

  trigPrimProducer_ = ps.getParameter<string>("trigPrimProducer");
  trigPrimCollection_ = ps.getParameter<string>("trigPrimCollection");
  outputFile_ = ps.getUntrackedParameter<string>("outputFile", "");
  
}


JetViaTrigPrim::~JetViaTrigPrim()
{
 
  
}


void
JetViaTrigPrim::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
 
  readAllCollections(event); //Reading collections
  const EcalTrigPrimDigiCollection* TPColl_;
  TPColl_ = (getTrigPrims(event)); //Read TP collection


  for(EcalTrigPrimDigiCollection::const_iterator it = TPColl_->begin(); //Loop over TPs
      it != TPColl_->end(); ++it){
    const EcalTriggerPrimitiveDigi& trigPrim = *it;

    EBSrFlagCollection::const_iterator srf  //Find associated SR flag
      = ebSrFlags_->find(trigPrim.id() );

  
    int flag = srf->value() & ~EcalSrFlag::SRF_FORCED_MASK; //Value of SR flag, removing forced bit

    if(flag==EcalSrFlag::SRF_ZS1) { // If TT zero suppressed
      
      vector<DetId> tpCrystals =  RetrieveCrystals( trigPrim.id() );
      
      for(int unsigned it=0;it<tpCrystals.size();it++)  { //Matching "geometrical" crystals with crystals not zero suppressed 
	
	const EBDigiCollection::const_iterator digi = ebDigis_->find( static_cast<const EBDetId&>(tpCrystals[it]) );
	
	if(digi != ebDigis_->end()) //Finding crystals not suppressed
	  {
	    //... Do what you want with crystals
	  }

	//If you don't want to match crystals like this you can retrieve crystals appearing in the jet you are studying
	//and matching them with the same algo. 


      }//end crystal loop

    }//end if ZS
    else if(flag==EcalSrFlag::SRF_FULL) { // If TT fully read out
      
       vector<DetId> tpCrystals =  RetrieveCrystals( trigPrim.id() );
    
       //... Do what you want with crystals and TPs.
    }//end if FRO

  }//end TP loop
}// end function





void 
JetViaTrigPrim::beginJob(const edm::EventSetup& setup)
{
  
  edm::ESHandle<EcalTrigTowerConstituentsMap> hTriggerTowerMap;
  setup.get<IdealGeometryRecord>().get(hTriggerTowerMap);
  triggerTowerMap_ = hTriggerTowerMap.product();
  
}


void 
JetViaTrigPrim::endJob() {
}


void JetViaTrigPrim::readAllCollections(const edm::Event& event) {
  ebSrFlags_.read(event); //Read SRFlags collections
  eeSrFlags_.read(event);
  ebDigis_.read(event);
}

//Filling TP collection
const EcalTrigPrimDigiCollection*
JetViaTrigPrim::getTrigPrims(const edm::Event& event)
{
  edm::Handle<EcalTrigPrimDigiCollection> hTPDigis;
  event.getByLabel(trigPrimProducer_, trigPrimCollection_, hTPDigis);
  return hTPDigis.product();
}


//Return crystals Contained in a trigger tower.
const vector<DetId> JetViaTrigPrim::RetrieveCrystals(const EcalTrigTowerDetId& ttid ) {

  return triggerTowerMap_->constituentsOf( ttid ); 
}


//define this as a plug-in
DEFINE_FWK_MODULE(JetViaTrigPrim);
