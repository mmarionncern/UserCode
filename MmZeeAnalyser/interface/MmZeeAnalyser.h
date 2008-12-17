#ifndef GsfElectronDataAnalyzer_h
#define GsfElectronDataAnalyzer_h



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

#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

class DQMStore;
class MonitorElement;


using namespace edm;
using namespace std;
using namespace cms;
//using reco::GsfElectronCollection;

 //Z mass (GeV)
static const double mZ0= 91.19;


class MmZeeAnalyser : public edm::EDAnalyzer
{
 public: //  -----------  functions ---------------
  MmZeeAnalyser(const edm::ParameterSet&);
  ~MmZeeAnalyser();


 public: //  -----------  member data --------------





 private: //  -----------  member data --------------
 
  edm::InputTag electronCollection_;

 

  //Number of electron in the event
  int NumberElectron_;

  //Zee Candidates
  int Candidate[5];

  //Histogrammes
  MonitorElement* Electron_charge;


 private: //  -----------  functions ---------------
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // General functions
  void NumberElectrons(const  edm::Handle<reco::GsfElectronCollection>& gsfElectrons);

  // Zee analyse functions
  void ZeeCandidate(const  edm::Handle<reco::GsfElectronCollection>& gsfElectrons);

 // Wrapper to fill methods of DQM monitor elements.
   
  void fill(MonitorElement* me, float x){
    if(me) me->Fill(x);
  }
  void fill(MonitorElement* me,float x, float yw){
    if(me) me->Fill(x, yw);
  }
  void fill(MonitorElement* me,float x, float y, float zw){
    if(me) me->Fill(x, y, zw);
  }
  void fill(MonitorElement* me,float x, float y, float z, float w){
    if(me) me->Fill(x, y, z, w);
  }

 
  // Wrappers to the book methods of the DQMStore DQM
  //  histogramming interface.
  
  MonitorElement* book1D(const std::string& name,
			 const std::string& title, int nbins,
			 double xmin, double xmax);
 
  MonitorElement* book2D(const std::string& name,
			 const std::string& title,
			 int nxbins, double xmin, double xmax,
			 int nybins, double ymin, double ymax);

  MonitorElement* bookProfile(const std::string& name,
			      const std::string& title,
			      int nbins, double xmin, double xmax);
  
  MonitorElement* bookProfile2D(const std::string& name,
				const std::string& title,
				int nbinx, double xmin, double xmax,
				int nbiny, double ymin, double ymax,
				const char* option = "");


  
 
};



#endif
