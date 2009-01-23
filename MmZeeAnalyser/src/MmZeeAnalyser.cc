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
// $Id: MmZeeAnalyser.cc,v 1.3 2008/12/17 09:13:21 mmarionn Exp $
//
//
#include "MMarionneau/MmZeeAnalyser/interface/MmZeeAnalyser.h"

// system include files
#include <memory>
#include <string.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/FWLite/interface/Event.h"

#include "PhysicsTools/CandUtils/interface/CandMatcher.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

//#include"Utilitaires.cc"

using namespace edm;
using namespace std;
using namespace cms;



// Constructors

MmZeeAnalyser::MmZeeAnalyser(const edm::ParameterSet& iConfig):

 
  // electronCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronCollection")),
  ZCandCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("ZCandidateCollection"))
  //ZRefCandCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("ZCandidateCollection")),
  // ElectronMCTruthTag_(iConfig.getUntrackedParameter<edm::InputTag>("ElectronMCTruthMap")),
  //ZeeMCTruthTag_(iConfig.getUntrackedParameter<edm::InputTag>("ZeeMCTruthCollection"))
{
  // DQM ROOT output
  outputFile_ = iConfig.getUntrackedParameter<string>("outputFile", ""); 
 if(outputFile_.size() != 0){
    LogInfo("OutputInfo") << " Zee Analysis Task histograms will be saved to '"
			  << outputFile_.c_str() << "'";
  } else{
    LogInfo("OutputInfo") << " Zee Analysis Task histograms will NOT be saved";
  }

  // get hold of back-end interface
 dbe_ = edm::Service<DQMStore>().operator->();
  
  /* if(verbose_){
    dbe_->setVerbose(1);
  } else{
    dbe_->setVerbose(0);
    }*/
  
  // if(verbose_) dbe_->showDirStructure();
  
 dbe_->setCurrentFolder(histDir_);

  vector<string>
    hists(iConfig.getUntrackedParameter<vector<string> >("histograms",
						    vector<string>(1, "all")));
  
  for(vector<string>::iterator it = hists.begin();
      it!=hists.end(); ++it) histList_.insert(*it);
  if(histList_.find("all") != histList_.end()) allHists_ = true;

  



  Mass_Z_HCL = book1D("MassofZforHCL",
		     "mass of Z reconstructed with CL of 4, 3 or 2;"
		     "Mass;Nevts",
		     100, 0., 200.);

  Mass_Z_LCL = book1D("MassofZforLCL",
		     "mass of Z reconstructed with CL of 1;"
		     "Mass;Nevts",
		     100, 0., 200.);

  Pt_Z_HCL = book1D("PtofZforHCL",
		    "Pt of Z reconstructed with CL of 4, 3 or 2;"
		    "Pt;Nevts",
		    100, 0., 200.);
  Eta_Z_HCL = book1D("EtaofZforHCL",
		     "Eta of Z reconstructed with CL of 4, 3 or 2;"
		     "Eta;Nevts",
		     80, -4., 4.);
  Phi_Z_HCL = book1D("PhiofZforHCL",
		     "Phi of Z reconstructed with CL of 4, 3 or 2;"
		     "Phi;Nevts",
		     64, -3.2, 3.2);
  
  LevelConf_Z = book1D("LevelConfforZ",
			"Confidence Level of Z;"
		       "CL;Nevts",
		       5,-0.5,4.5);

  MCFalse_LevelConf_Z = book1D("MCFalseLevelConfforZ",
			"Confidence Level of Z if MC->False;"
		       "CL;Nevts",
		       5,-0.5,4.5);
  
  PtvsdPhi_Z = book2D("PtVsdPhiforZ",
		      "Pt vs dphi of Z;"
		      "Pt;dPhi",
		      100,0.,200.,90,-1.,181.);
  
  Pt_dau_0 = book1D("Ptofdaughter_0",
		     "Pt of the daughter 1;"
		     "Pt;Nevts",
		     100, 0., 200.);
 
  Eta_dau_0 = book1D( "Etaofdaughter_0",
		       "Eta of the daughter 1;"
		       "Eta;Nevts",
			  80, -4., 4.);
 
  Phi_dau_0 = book1D("Phiofdaughter_0",
		     "Phi of the daughter 1;"
		      "Phi;Nevts",
		      180, -181., 181.);
 
  Class_dau_0 = book1D("Classificationofdaughter_0",
		       "Classification of the daughter 1;"
		       "Classification;Nevts",
		       5, -0.5, 4.5);
  
  deltaPtMC_dau_0 = book1D("deltaPtofdaughter_0",
			   "#deltaPt mc/pat of the daughter 1;"
			   "Pt;Nevts",
			   100, 0., 20.);
  
  Pt_dau_1 = book1D("Ptofdaughter_1",
		    "Pt of the daughter 2;"
		    "Pt;Nevts",
		    100, 0., 200.);
  
  Eta_dau_1 = book1D( "Etaofdaughter_1",
		      "Eta of the daughter 2;"
		      "Eta;Nevts",
		      80, -4., 4.);
  
  Phi_dau_1 = book1D("Phiofdaughter_1",
		     "Phi of the daughter 2;"
		     "Phi;Nevts",
		      180, -181., 181.);
  
  Class_dau_1 = book1D("Classificationofdaughter_1",
		       "Classification of the daughter 2;"
		       "Classification;Nevts",
		       5, -0.5, 4.5);

  deltaPtMC_dau_1 = book1D("deltaPtofdaughter_1",
			   "deltaPt mc/pat of the daughter 2;"
			   "#deltaPt;Nevts",
			   100, 0., 20.);
  
  dPhi_dau =  book1D("dPhidaughters",
		     "#Delta#Phi between daughters ;"
		     "dPhi;Nevts",
		     90, -1., 181);

  deltaMassMC_Z = book1D("deltaMassMC_Z",
			 "#DeltaMass mc/pat of the Z candidate;"
			 "#DeltaM;NEvts",
			 80,0.,4.);

}

// Destructors

MmZeeAnalyser::~MmZeeAnalyser()
{
 if(outputFile_.size()!=0) dbe_->save(outputFile_);
}

// ------------ method called to for each event  ------------
void
MmZeeAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  EventNumber_++;
  // using namespace edm;
  using namespace reco;
  using namespace pat;
   
  //Collection reading    
  /* edm::Handle<edm::View<pat::Electron> > ElectronsHandle;
   iEvent.getByLabel(electronCollectionTag_, ElectronsHandle);
   edm::View<pat::Electron> Electrons= *ElectronsHandle;*/

   edm::Handle<reco::CompositeCandidateCollection> ZeeCandidatesH;
   iEvent.getByLabel(ZCandCollectionTag_, ZeeCandidatesH);
   reco::CompositeCandidateCollection ZeeCandidates = *ZeeCandidatesH;
   
   
   vector<reco::CompositeCandidate> ClassifiedZColl;
   ClassifiedZColl = ClassZeeCandidate(ZeeCandidates);
   
   for(int unsigned iZee=0;iZee<ClassifiedZColl.size();iZee++)
     {
       fill(LevelConf_Z, ZeeConfLvlFromElec(ClassifiedZColl[iZee]));
       /* cout<<"Classification result  "
	   <<ZeeConfLvlFromElec(ClassifiedZColl[iZee])<<"   "
	   <<ClassifiedZColl[iZee].mass();*/
       
       double deltaPt[2] = {0,0};
       double deltaZMass = 0;
       
       bool mctruth =  MCTruth(ClassifiedZColl[iZee], deltaPt, deltaZMass);
       
       if(mctruth)
	 {
	   fill(deltaMassMC_Z,deltaZMass);
	   fill(deltaPtMC_dau_0,deltaPt[0]);
	   fill(deltaPtMC_dau_1,deltaPt[2]);
	 }
       else 
	 {
	   fill(MCFalse_LevelConf_Z, ZeeConfLvlFromElec(ClassifiedZColl[iZee]));
	   cout<<"Warning, Zee is not validated by MCTruth : "
	      <<" #Evt = "<<EventNumber_<<" ; Numero of Z  "<<iZee+1<<endl;
	 cout<<"Classification of Z "<<ZeeConfLvlFromElec(ClassifiedZColl[iZee])
	     <<" mass "<<ClassifiedZColl[iZee].mass()<<endl;
	 }
       if(ZeeConfLvlFromElec(ClassifiedZColl[0])>1)
	 {

	   ElectronFromZAnalysis(ClassifiedZColl[0]);

	   fill(Mass_Z_HCL, ClassifiedZColl[0].mass());
	   fill(Pt_Z_HCL, ClassifiedZColl[0].pt()); 
	   fill(Eta_Z_HCL, ClassifiedZColl[0].eta()); 
	   fill(Phi_Z_HCL, ClassifiedZColl[0].phi()); 
	 }
       else
	 fill(Mass_Z_LCL,ClassifiedZColl[0].mass());
     
     }
}


void 
MmZeeAnalyser::beginJob(const edm::EventSetup&)
{
  EventNumber_=0;

}


void 
MmZeeAnalyser::endJob() {
}


//--------------------- classing "Z to ee" candidate in the event
 vector<reco::CompositeCandidate>
MmZeeAnalyser::ClassZeeCandidate(const  reco::CompositeCandidateCollection& ZeeCandidates)
{
  vector<reco::CompositeCandidate> ClassifiedZColl ;

  if(ZeeCandidates.size()==0)
    {
      cout<<" Warning, no reconstructed Z in this event"<<endl;
      return ClassifiedZColl;
    }
  else if(ZeeCandidates.size()==1)
    {
      if(ZeeConfLvlFromElec(ZeeCandidates[0]) < 1 )
	{
	  cout<<"Warning, Confidence level on electrons is really bad"<<endl;
	  cout<<"but alone candidate"<<endl;
	  ClassifiedZColl.push_back(ZeeCandidates[0]);
	}
      else
	{
	  ClassifiedZColl.push_back(ZeeCandidates[0]);
	}
    }
  else
    {
      int iZee=0; int goodZCand=0;
      for(reco::CompositeCandidateCollection::const_iterator itZee = ZeeCandidates.begin(); itZee != ZeeCandidates.end(); itZee++, iZee++)
	{
	  const reco::CompositeCandidate& Zee = *itZee;
	
	  if(ZeeConfLvlFromElec(Zee) < 1 ){
	    cout<<"Warning, Confidence level on electrons is really bad"<<endl;
	    cout<<" candidate not classified "<<endl;
	    continue;
	  }
	  else
	    {
	      ClassifiedZColl.push_back(Zee);
	      goodZCand++;
	    }
	 
	}

      for(int i=0;i<goodZCand;i++){
	for(int j=0;j<goodZCand;j++)
	  {
	    if(i != j)
	      if( ZeeConfLvlFromElec(ClassifiedZColl[i])
		  > ZeeConfLvlFromElec(ClassifiedZColl[j]) )
		{
		  const reco::CompositeCandidate Zee_tmp = ClassifiedZColl[i];
		  ClassifiedZColl[i] = ClassifiedZColl[j];
		  ClassifiedZColl[j] = Zee_tmp;
		}
	      else if( ZeeConfLvlFromElec(ClassifiedZColl[i])
		       == ZeeConfLvlFromElec(ClassifiedZColl[j]) )
		{
		  double dm1 = fabs(ClassifiedZColl[i].mass()-mZ0);
		  double dm2 = fabs(ClassifiedZColl[j].mass()-mZ0);
		  
		  if(dm1 < dm2 )
		    {
		   const reco::CompositeCandidate Zee_tmp = ClassifiedZColl[i];
		   ClassifiedZColl[i] = ClassifiedZColl[j];
		   ClassifiedZColl[j] = Zee_tmp;
		    }
		}
	  }
      }//fin algo de tri des Z
    }
  return ClassifiedZColl;
}

 


int 
MmZeeAnalyser::ZeeConfLvlFromElec(const reco::CompositeCandidate& Zee){

  float dau_Class[2] = {-1,-1};
 
  if(Zee.numberOfDaughters()!=2)
    {
      cout<<" Bad reconstruction : number of daughter not egal to two "<<endl;
      return -1;
    }

 const pat::Electron *originalElectron1 = dynamic_cast<const pat::Electron *>(Zee.daughter(0)->masterClone().get());
 const pat::Electron *originalElectron2 = dynamic_cast<const pat::Electron *>(Zee.daughter(1)->masterClone().get());

 dau_Class[0] = originalElectron1->leptonID("robust") + 
                2*originalElectron1->leptonID("tight") ;

 dau_Class[1] = originalElectron2->leptonID("robust") +
                2*originalElectron2->leptonID("tight") ;

 //Filling classification histos
 fill(Class_dau_0,dau_Class[0]);
 fill(Class_dau_1,dau_Class[1]);


  if(dau_Class[0]==3 && dau_Class[1]==3) //Two electrons are tights. ***
    return 4;
  else if(dau_Class[0]==1 && dau_Class[1]==3) //One electron is tight
    return 3;
  else if(dau_Class[1]==1 && dau_Class[0]==3) // and one is robust. ***
    return 3;
  else if(dau_Class[1]==1 && dau_Class[0]==1) //Two electrons are robust. ***
    return 2;
  else if(dau_Class[1]==0 && ( dau_Class[0]==3 || dau_Class[0]==1)) //One elec
    return 1;
  else if(dau_Class[0]==0 && ( dau_Class[1]==3 || dau_Class[1]==1))// is bad. *
    return 1;
  else if(dau_Class[0]==0 && dau_Class[1]==0) //The two are bad. ***
    return 0;
  else
    {
      cout<<endl;
      cout<<" Bad recognition, there is no two electron in this Z candidate"<<endl;
      return -1;
    }

}

//---------------------- Histos Method --------------------

MonitorElement* MmZeeAnalyser::book1D(const std::string& name, const std::string& title, int nbins, double xmin, double xmax){
  if(!registerHist(name, title)) return 0; //this histo is disabled
  MonitorElement* result = dbe_->book1D(name, title, nbins, xmin, xmax);
  if(result==0){
    throw cms::Exception("Histo")
      << "Failed to book histogram " << name;
  }
  return result;
}

MonitorElement* MmZeeAnalyser::book2D(const std::string& name, const std::string& title, int nxbins, double xmin, double xmax, int nybins, double ymin, double ymax){
  if(!registerHist(name, title)) return 0; //this histo is disabled
  MonitorElement* result = dbe_->book2D(name, title, nxbins, xmin, xmax,
					nybins, ymin, ymax);
  if(result==0){
    throw cms::Exception("Histo")
      << "Failed to book histogram " << name;
  }
  return result;
}

MonitorElement* MmZeeAnalyser::bookProfile(const std::string& name, const std::string& title, int nbins, double xmin, double xmax){
  if(!registerHist(name, title)) return 0; //this histo is disabled
  MonitorElement* result = dbe_->bookProfile(name, title, nbins, xmin, xmax,
					     0, 0, 0);
  if(result==0){
    throw cms::Exception("Histo")
      << "Failed to book histogram " << name;
  }
  return result;
}

MonitorElement* MmZeeAnalyser::bookProfile2D(const std::string& name, const std::string& title, int nbinx, double xmin, double xmax, int nbiny, double ymin, double ymax, const char* option){
  if(!registerHist(name, title)) return 0; //this histo is disabled
  MonitorElement* result
    = dbe_->bookProfile2D(name,
			  title,
			  nbinx, xmin, xmax,
			  nbiny, ymin, ymax,
			  0, 0, 0,
			  option);
  if(result==0){
    throw cms::Exception("Histo")
      << "Failed to book histogram " << name;
  }
  return result;
}

bool MmZeeAnalyser::registerHist(const std::string& name,
				 const std::string& title){
  availableHistList_.insert(pair<string, string>(name, title));
  return allHists_ || histList_.find(name)!=histList_.end();
}


//MC validation

bool MmZeeAnalyser::MCTruth(const reco::CompositeCandidate& Zee,
			    double deltaPt[2],
			    double deltaZMass){
  using namespace reco;
  
  bool mctruth=false;

  const pat::Electron *originalElectron1 = dynamic_cast<const pat::Electron *>(Zee.daughter(0)->masterClone().get());
  const pat::Electron *originalElectron2 = dynamic_cast<const pat::Electron *>(Zee.daughter(1)->masterClone().get());
  
  
   reco::GenParticleRef match1 = originalElectron1->genParticleRef();
   reco::GenParticleRef match2 = originalElectron2->genParticleRef();
   
   if( (match1->pdgId()==originalElectron1->pdgId())
       && (match2->pdgId()==originalElectron2->pdgId()) )
     {
       if( (match1->numberOfMothers()>0) && (match2->numberOfMothers()>0))
	 {
	   const Candidate * mother1 = match1->mother();
	   const Candidate * mother2 = match2->mother();
	   if( (mother1->numberOfMothers()>0) && (mother2->numberOfMothers()>0))
	     {
	       const Candidate * mother_lvl2_1 = mother1->mother();
	       const Candidate * mother_lvl2_2 = mother2->mother();
	       
	       if((mother_lvl2_1->pdgId()==23) && (mother_lvl2_2->pdgId()==23))
		 {
		   deltaPt[0] = abs(originalElectron1->pt() - match1->pt());
		   deltaPt[1] = abs(originalElectron2->pt() - match2->pt());
		   deltaZMass = abs(Zee.mass()-mother_lvl2_1->mass());
		   mctruth=true;
		 }
	       else
		 mctruth=false;
	     }
	   else mctruth=false;
	 }
       else mctruth=false;
     }
   else mctruth=false;
   
   if(!mctruth)
     {
       deltaPt[0]=-1;
       deltaPt[1]=-1;
       deltaZMass=-1;  
     }
   return mctruth;
}

void MmZeeAnalyser::ElectronFromZAnalysis(const reco::CompositeCandidate& Zee){

  using namespace reco;

  const pat::Electron *originalElectron1 = dynamic_cast<const pat::Electron *>(Zee.daughter(0)->masterClone().get());
  const pat::Electron *originalElectron2 = dynamic_cast<const pat::Electron *>(Zee.daughter(1)->masterClone().get());
  
  fill( Pt_dau_0,originalElectron1->pt());
  fill( Eta_dau_0,originalElectron1->eta());
  fill( Phi_dau_0,originalElectron1->phi()*180/myPi);

  fill( Pt_dau_1,originalElectron2->pt());
  fill( Eta_dau_1,originalElectron2->eta());
  fill( Phi_dau_1,originalElectron2->phi()*180/myPi);

  return;
}
  


//define this as a plug-in
DEFINE_FWK_MODULE(MmZeeAnalyser);
