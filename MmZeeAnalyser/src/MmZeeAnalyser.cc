// -*- C++ -*-
//
// Package:    MmZeeAnalyser
// Class:      MmZeeAnalyser
// 
/**\class MmZeeAnalyser MmZeeAnalyser.cc mmarionn/MmZeeAnalyser/src/MmZeeAnalyser.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Matthieu Marionneau
//         Created:  Mon Nov 10 14:59:45 CET 2008
// $Id: MmZeeAnalyser.cc,v 1.2 2008/11/21 17:13:04 mmarionn Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/FWLite/interface/Event.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "../interface/MmZeeAnalyser.h"

using namespace edm;
using namespace std;
using namespace cms;


// Constructors

MmZeeAnalyser::MmZeeAnalyser(const edm::ParameterSet& iConfig):

  //trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")),
  electronCollection_(iConfig.getUntrackedParameter<edm::InputTag>("electronCollection_"))
{
   //now do what ever initialization is needed

  /*Electron_charge = book1D("ChargeElectron",
			   "Charge Distribution",
			   5,-1.5,1.5);
  */

}

// Destructors

MmZeeAnalyser::~MmZeeAnalyser()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MmZeeAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //Initialisation
   for(int i=0;i<5;i++)
     Candidate[i] = 0;

   
   //Electron Recognition
   
   cout << "analyzing new event " <<endl;
   
    using reco::GsfElectronCollection;
   
   edm::Handle<GsfElectronCollection> gsfElectrons;
   iEvent.getByLabel(electronCollection_, gsfElectrons);
   edm::LogInfo("")<<"\n\n =================> Treating event "<<iEvent.id()<<" Number of electrons "<<gsfElectrons.product()->size();
   
   //Evaluation of the NumberElectron_ value    
   NumberElectrons(gsfElectrons);
  
   for(GsfElectronCollection::const_iterator itElec = gsfElectrons->begin();
       itElec != gsfElectrons->end();
       itElec++) 
     {
       cout<<" test "<<itElec->isElectron()<<"  "<<itElec->gsfTrack().index() 
	   <<"   "<<itElec->charge()<<"  "<<itElec->energy()<<"  "<<itElec->pt()  <<endl;
       
       
       // fill(Electron_charge,itElec->charge());
     }
   
   cout<<" Number of electrons  "<<NumberElectron_<<" ; Event  "<<iEvent.id()<<endl;

   ZeeCandidate(gsfElectrons);

  /*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
  */


}


// ------------ method called once each job just before starting event loop  ------------
void 
MmZeeAnalyser::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MmZeeAnalyser::endJob() {
}


//--------------------- Finding best "Z to ee" candidate in the event
void 
MmZeeAnalyser::ZeeCandidate(const  edm::Handle<reco::GsfElectronCollection>& gsfElectrons)
{
  using reco::GsfElectronCollection;

  int numCandidate=0;

 

  if(NumberElectron_<2)
    {
      cout<<" No Z reconstruction possibility in this event "<<endl;
      return;
    }

  else
    {
      for(GsfElectronCollection::const_iterator itElec = gsfElectrons->begin();
	  itElec != gsfElectrons->end(); itElec++) 
	{
	  if(itElec->charge() == -1){
	    
	  

	    for(GsfElectronCollection::const_iterator itPosi = gsfElectrons->begin();
		itPosi != gsfElectrons->end(); itPosi++) 
	    {
	      if(itPosi->charge() == 1)
		{
		  Candidate[numCandidate]= itElec->gsfTrack().index() 
		                         + 1000*(itPosi->gsfTrack().index());
		  numCandidate++;
		
		}
	    } 
	  }
     	}
     

      for(int i=0;i<numCandidate;i++)
	{
	
	  cout<<" Candidat Zee : Positron "<<Candidate[i]%1000
	      <<"  ; Electron "<<Candidate[i]-Candidate[i]*1000
	      <<endl;
	}
    }
     
}


//----------------------- Method called to looking for the number of electrons
void
MmZeeAnalyser::NumberElectrons(const  edm::Handle<reco::GsfElectronCollection>& gsfElectrons) 
{
  using reco::GsfElectronCollection;
  NumberElectron_ = 0;
  int tmp_number=0;

  for(GsfElectronCollection::const_iterator itElec = gsfElectrons->begin();
      itElec != gsfElectrons->end();  itElec++) 
    tmp_number++;

  NumberElectron_ = tmp_number;
}








//define this as a plug-in
DEFINE_FWK_MODULE(MmZeeAnalyser);
