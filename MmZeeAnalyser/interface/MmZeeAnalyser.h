#ifndef GsfElectronDataAnalyzer_h
#define GsfElectronDataAnalyzer_h



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <utility>

class DQMStore;
class MonitorElement;


using namespace edm;
using namespace std;
using namespace cms;
//using reco::GsfElectronCollection;

 //Z mass (GeV)
//static const double mZ0= 91.19;
//static const double myPi=acos(-1);


class MmZeeAnalyser : public edm::EDAnalyzer
{
 public: //  -----------  functions ---------------
  MmZeeAnalyser(const edm::ParameterSet&);
  ~MmZeeAnalyser();


 public: //  -----------  member data --------------





 private: //  -----------  member data --------------
 
 
  InputTag electronCollectionTag_;
  InputTag ZCandCollectionTag_;
  InputTag electronGSFCollectionTag_;
  InputTag genParticlesTag_;
  InputTag trackCollectionTag_;
 

  //Number of electron in the event
  int NumberElectron_;

  int EventNumber_;
  int GoodCandidate_;
  int NoZCandidate_;

  //Zee Candidates
  // int Candidate[5];

  //Histogrammes
  MonitorElement* MCTruth_Mass;
  MonitorElement* Mass_Z_HCL;
  MonitorElement* Mass_Z_LCL;

  MonitorElement* Zmass_bySC;
  MonitorElement* Mass_no_MCT;
  MonitorElement* Mass_no_MCT2;

  MonitorElement* Pt_Z_HCL;
  MonitorElement* Eta_Z_HCL;
  MonitorElement* Phi_Z_HCL;
  
  MonitorElement* LevelConf_Z;
  MonitorElement* MCFalse_LevelConf_Z;
  MonitorElement* MCFalse_CL_ColS_Z;

  MonitorElement* PtvsdPhi_Z;

  MonitorElement* Pt_dau_0;
  MonitorElement* Eta_dau_0;
  MonitorElement* Phi_dau_0;
  MonitorElement* Class_dau_0 ;
  MonitorElement* deltaPtMC_dau_0;
  MonitorElement* Pt_dau_1;
  MonitorElement* Eta_dau_1;
  MonitorElement* Phi_dau_1;
  MonitorElement* Class_dau_1;
  MonitorElement* deltaPtMC_dau_1;
 
  MonitorElement* dPhi_dau;
  MonitorElement* EOP_electron;

  MonitorElement* deltaMassMC_Z;

  MonitorElement* Eta_MCT;
  MonitorElement* Pt_MCT;
  MonitorElement* Eta_vs_E_MCT;
  MonitorElement* Nelectron;
  MonitorElement* No_Zreco_code;
  MonitorElement* MCFalse_Z_code;

  MonitorElement* SuperCluster_phiWidth;
  MonitorElement* SuperCluster_etaWidth;
 
 private: //  -----------  functions ---------------
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // General functions

  bool MCTruth(const reco::GenParticleCollection& Genparticles,
	       const reco::CompositeCandidate& Zee,
	       double deltaPt[2], double& deltaZMass, int& code);

  void ElectronTruth(const edm::View<pat::Electron>& ElectronColl,
		     const reco::GsfTrackCollection& TrackColl,
		     const reco::GenParticleCollection& Genparticles,
		     int& code);

  // Zee analyse functions
  vector<reco::CompositeCandidate> ClassZeeCandidate(const  reco::CompositeCandidateCollection& ZeeCandidates);

  //Zee Confidence level calculation with electron conditions
  int ZeeConfLvlFromElec(const reco::CompositeCandidate& Zee);

  //Electron analyze functions
  void ElectronFromZAnalysis(const reco::CompositeCandidate& Zee);

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
  
 /** List of enabled histograms. Special name "all" is used to indicate
   * all available histograms.
   */
  std::set<std::string> histList_;

  /** When true, every histogram is enabled.
   */
  bool allHists_;

  /** Histogram directory PATH in DQM or within the output ROOT file
   */
  std::string histDir_;
  
  /** List of available histograms. Filled by the booking methods.
   * key: name, value: title.
   */
  std::map<std::string, std::string> availableHistList_;

  /** Register a histogram in the available histogram list and check if
   * the histogram is enabled. Called by the histogram booking methods.
   * @return true if the histogram is enable, false otherwise
   */
  bool registerHist(const std::string& name, const std::string& title);

  /** Prints the list of available histograms
   * (registered by the registerHist method), including disabled one.
   */
  void printAvailableHists();

  ///Histogramming interface
  DQMStore* dbe_;
  
  ///Output file for histograms
  std::string outputFile_;
  
  
 
};



#endif
