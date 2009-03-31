#ifndef GsfElectronDataAnalyzer_cc
#define GsfElectronDataAnalyzer_cc

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
// $Id: MmZeeAnalyser.cc,v 1.7 2009/01/28 13:12:38 mmarionn Exp $
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
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/FWLite/interface/Event.h"

#include "PhysicsTools/CandUtils/interface/CandMatcher.h"
#include "RecoEgamma/EgammaTools/interface/HoECalculator.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "MMarionneau/MmZeeAnalyser/Utilitaires/Utilitaires_Geom.cc"


using namespace edm;
using namespace std;
using namespace cms;



// Constructors

MmZeeAnalyser::MmZeeAnalyser(const edm::ParameterSet& iConfig):

  electronCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronCollection")),
  ZCandCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("ZCandidateCollection")),
  genParticlesTag_(iConfig.getUntrackedParameter<edm::InputTag>("genParticleCollection")),
  trackCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("trackCollection"))
{

  //Analyse type Definition
  analyse_type_ = iConfig.getUntrackedParameter<string>("analyseType");

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
  
 dbe_->setCurrentFolder(histDir_);

  vector<string>
    hists(iConfig.getUntrackedParameter<vector<string> >("histograms",
						    vector<string>(1, "all")));
  
  for(vector<string>::iterator it = hists.begin();
      it!=hists.end(); ++it) histList_.insert(*it);
  if(histList_.find("all") != histList_.end()) allHists_ = true;

  
  MCTruth_Mass =  book1D("MassofZMC",
			 "mass of Z from MCTruth;"
			 "Mass;Nevts",
			 200, 0., 200.);


  Mass_Z_HCL = book1D("MassofZforHCL",
		     "mass of Z reconstructed with CL of 4, 3 or 2;"
		     "Mass;Nevts",
		     200, 0., 200.);

  Mass_Z_LCL = book1D("MassofZforLCL",
		     "mass of Z reconstructed with CL of 1;"
		     "Mass;Nevts",
		     200, 0., 200.);

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

  MCFalse_LevelConf_Z = book2D("MCFalseLevelConfforZ_EColS",
			"Confidence Level of Z if MC->False;"
		       "CL;NbElectron",
		       5,-0.5,4.5,11,-0.5,10.5);

  MCFalse_CL_ColS_Z = book2D("MCFalseLevelConfforZ_ColS",
			"Confidence Level of Z if MC->False;"
		       "CL;Place dans la liste",
		       5,-0.5,4.5,11,-0.5,10.5);
  
  PtvsdPhi_Z = book2D("PtVsdPhiforZ",
		      "Pt vs dphi of Z;"
		      "Pt;dPhi",
		      100,0.,200.,92,-3.,183.);
  
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
		      182, -183., 183.);
 
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
		      182, -183., 183.);
  
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
		     92, -3., 183);

  deltaMassMC_Z = book1D("deltaMassMC_Z",
			 "#DeltaMass mc/pat of the Z candidate;"
			 "#DeltaM;NEvts",
			 160,-40.,40.);
  
  EOP_electron = book1D("EOP_electron",
			"Energy on momentum for electrons;"
			"EOP (GeV);Nevts",
			100,0.5,4.5);
  

  Eta_MCT = book2D( "Eta_MCT",
		    "Eta of electron;"
		    "Eta;Nevts",
		    120, -6., 6.,120, -6., 6.);
  
  Pt_MCT = book2D("Pt_MCT",
		  "Pt of electron;"
		  "Pt1;Pt2",
		  100, 0., 50.,100, 0., 50.);

  Eta_vs_E_MCT = book2D("Eta_cs_E",
			"",
			120,-6., 6.,100,0.,200.);

  Nelectron = book1D("Nelec",
		     "",
		     120,0.,6);

  No_Zreco_code = book1D("No_Zreco_code",
			 "Code corresponding to the non-tagged Z;"
			 "Code;Nevts",
			 11,-0.5,10.5);

  MCFalse_Z_code = book1D("MCFalse_Z_code",
			 "Code corresponding to the Z without MCT;"
			 "Code;Nevts",
			 11,-1.5,9.5);

  SuperCluster_etaWidth = book1D( "etaWidth",
				  "EtaWidth of supercluster;"
				  "Eta;Nevts",
				  100, 0., 0.2);

  SuperCluster_phiWidth = book1D( "phiWidth",
				  "PhiWidth of supercluster;"
				  "Phi;Nevts",
				  100, 0., 0.2);

  Zmass_bySC = book1D("MassofZbySC",
		     "mass of Z reconstructed with the SC data;"
		     "Mass;Nevts",
		     200, 0., 200.);

  Mass_no_MCT = book1D("MassofZnoMCT_saved",
		     "mass of Z reconstructed for no MCT;"
		     "Mass;Nevts",
		     200, 0., 200.);

 Mass_no_MCT2 = book1D("MassofZnoMCT_unsaved",
		     "mass of Z reconstructed for no MCT;"
		     "Mass;Nevts",
		     200, 0., 200.);

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
   edm::Handle<edm::View<pat::Electron> > ElectronsHandle;
   iEvent.getByLabel(electronCollectionTag_, ElectronsHandle);
   edm::View<pat::Electron> Electrons= *ElectronsHandle;

   edm::Handle<reco::CompositeCandidateCollection> ZeeCandidatesH;
   iEvent.getByLabel(ZCandCollectionTag_, ZeeCandidatesH);
   reco::CompositeCandidateCollection ZeeCandidates = *ZeeCandidatesH;

   edm::Handle<reco::GsfTrackCollection> GsfTrackH;
   iEvent.getByLabel(trackCollectionTag_, GsfTrackH);
   reco::GsfTrackCollection Tracks = *GsfTrackH;

    edm::Handle<reco::GenParticleCollection> GenParticlesH;
   iEvent.getByLabel(genParticlesTag_, GenParticlesH);
   reco::GenParticleCollection GenParticles = *GenParticlesH;
   

   //liste classifiee des candidates  
   vector<reco::CompositeCandidate> ClassifiedZColl;
   //liste des nuiveaux de confiance de ces candidats
   vector< int > ClassifiedZCL;

   ClassifiedZColl = ClassZeeCandidate(ZeeCandidates,Tracks);

   if(analyse_type_ =="TDR")
     for(int unsigned iZee=0; iZee<ClassifiedZColl.size();iZee++)
       ClassifiedZCL.push_back( ZeeConfLvlFromTDRElec(ClassifiedZColl[iZee]) );
   else if(analyse_type_ =="EWK")
     for(int unsigned iZee=0; iZee<ClassifiedZColl.size();iZee++)
       ClassifiedZCL.push_back( ZeeConfLvlFromEWKElec(ClassifiedZColl[iZee],Tracks) );

     //Matching MCTruth
   if(ClassifiedZColl.size()==0)
     {
       int code=-1;
       
       ElectronTruth(Electrons,Tracks,GenParticles,code);
       
       fill(No_Zreco_code,code);
       NoZCandidate_++;
     }
   else{
     double deltaPt[2] = {0,0};
     double deltaZMass = 0;
     
     int code=-10;
     
     bool mctruth =  MCTruth(GenParticles,
			     ClassifiedZColl[0],
			     deltaPt, deltaZMass,code);
     if(mctruth)
       {
	 GoodCandidate_++;
	 fill(deltaMassMC_Z,deltaZMass);
	 fill(deltaPtMC_dau_0,deltaPt[0]);
	 fill(deltaPtMC_dau_1,deltaPt[1]);
	 
	 fill(LevelConf_Z, ClassifiedZCL[0]);
	 
	 //cout<< ClassifiedZCL[0]<<endl;
       }
     else 
       {
	 // if(iZee==0)
	 fill(MCFalse_LevelConf_Z, ClassifiedZCL[0] ,Electrons.size());
	 fill(MCFalse_CL_ColS_Z,
	      ClassifiedZCL[0],
	      1);
	 fill(MCFalse_Z_code,code);
       }
     
     if( ClassifiedZCL[0] >1)
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
  GoodCandidate_=0;
  NoZCandidate_=0;
}


void 
MmZeeAnalyser::endJob() {

  cout<<endl;
  cout<<"Number of event : "<<EventNumber_<<endl;
  cout<<"Number of selected Candidates "<<EventNumber_-NoZCandidate_<<endl;
  cout<<"Number of good candidates (MCT) : "<<GoodCandidate_<<endl;
  cout<<" total efficiency "<<(float)GoodCandidate_/EventNumber_
      <<" partial efficiency "<<(float)GoodCandidate_/(EventNumber_-NoZCandidate_)<<endl;
  cout<<"Reconstruction loose "<<NoZCandidate_<<" ->% "
      <<(float)NoZCandidate_/EventNumber_<<endl;
  

}


//--------------------- classing "Z to ee" candidate in the event
 vector<reco::CompositeCandidate>
MmZeeAnalyser::ClassZeeCandidate(const  reco::CompositeCandidateCollection& ZeeCandidates,
				 const reco::GsfTrackCollection Tracks)
{
  vector<reco::CompositeCandidate> ClassifiedZColl ;
  if(analyse_type_ == "TDR"){
    if(ZeeCandidates.size()==0) //pas de candidat Z
      {
	// Warning, no reconstructed Z in this event
	return ClassifiedZColl;
      }
    else if(ZeeCandidates.size()==1) // un seul candidat Z
      {
	if(ZeeConfLvlFromTDRElec(ZeeCandidates[0]) < 1 )
	  {
	    //Warning, Confidence level on electrons is really bad
	    //but alone candidate
	    ClassifiedZColl.push_back(ZeeCandidates[0]);
	  }
	else
	  {
	    ClassifiedZColl.push_back(ZeeCandidates[0]);
	  }
      }
    else     //plusieurs candidats Z, algo de tri necessaire
      {      // base sur le niveau de confiance des electrons, puis sur la masse
	int iZee=0; int goodZCand=0;
	for(reco::CompositeCandidateCollection::const_iterator itZee = ZeeCandidates.begin(); itZee != ZeeCandidates.end(); itZee++, iZee++)
	  {
	    const reco::CompositeCandidate& Zee = *itZee;
	    
	    if(ZeeConfLvlFromTDRElec(Zee) < 1 )
	      {//Warning, Confidence level on electrons is really bad
	      // candidate not classified
		continue;}
	    else
	      {ClassifiedZColl.push_back(Zee);
	      goodZCand++;}
	  }
	for(int i=0;i<goodZCand;i++){
	  for(int j=0;j<goodZCand;j++)
	    {
	      if(i != j)
		if( ZeeConfLvlFromTDRElec(ClassifiedZColl[i])
		    > ZeeConfLvlFromTDRElec(ClassifiedZColl[j]) )
		  {const reco::CompositeCandidate Zee_tmp = ClassifiedZColl[i];
		  ClassifiedZColl[i] = ClassifiedZColl[j];
		  ClassifiedZColl[j] = Zee_tmp;
		  }
		else if( ZeeConfLvlFromTDRElec(ClassifiedZColl[i])
			 == ZeeConfLvlFromTDRElec(ClassifiedZColl[j]) )
		  {
		    double dm1 = fabs(ClassifiedZColl[i].mass()-mZ0);
		    double dm2 = fabs(ClassifiedZColl[j].mass()-mZ0);
		    
		    if(dm1 < dm2 )
		      {const reco::CompositeCandidate Zee_tmp=ClassifiedZColl[i];
		      ClassifiedZColl[i] = ClassifiedZColl[j];
			ClassifiedZColl[j] = Zee_tmp;
		      }
		  }
	    }
	}//fin algo de tri des Z
      }
  }
  else if (analyse_type_ == "EWK")
    {
      if(ZeeCandidates.size()==0) //pas de candidat Z
	{ //Warning, no reconstructed Z in this event
	  return ClassifiedZColl;cout<<" youpiyoupa_1_2_1!  "<<endl;}
      else if(ZeeCandidates.size()==1) // un seul candidat Z
	{
	  if(ZeeConfLvlFromEWKElec(ZeeCandidates[0],Tracks) < 1 )
	    { //Warning, Confidence level on electrons is really bad
	      //but alone candidate
	      ClassifiedZColl.push_back(ZeeCandidates[0]);}
	  else
	    {ClassifiedZColl.push_back(ZeeCandidates[0]);}
	}
      else     //plusieurs candidats Z, algo de tri necessaire
	{      // base sur le niveau de confiance des electrons, puis sur la masse
	  int iZee=0; int goodZCand=0;
	  for(reco::CompositeCandidateCollection::const_iterator itZee = ZeeCandidates.begin(); itZee != ZeeCandidates.end(); itZee++, iZee++)
	    {
	      const reco::CompositeCandidate& Zee = *itZee;
	      if(ZeeConfLvlFromEWKElec(Zee,Tracks) < 1 )
		{//Warning, Confidence level on electrons is really bad
		  //candidate not classified
		  continue;}
	      else
		{ ClassifiedZColl.push_back(Zee);
		goodZCand++; }
	    }
	  for(int i=0;i<goodZCand;i++){
	    for(int j=0;j<goodZCand;j++){
	      if(i != j)
		if( ZeeConfLvlFromEWKElec(ClassifiedZColl[i], Tracks)
		    > ZeeConfLvlFromEWKElec(ClassifiedZColl[j], Tracks) )
		  {
		    const reco::CompositeCandidate Zee_tmp = ClassifiedZColl[i];
		    ClassifiedZColl[i] = ClassifiedZColl[j];
		    ClassifiedZColl[j] = Zee_tmp;
		  }
		else if( ZeeConfLvlFromEWKElec(ClassifiedZColl[i], Tracks)
			 == ZeeConfLvlFromEWKElec(ClassifiedZColl[j], Tracks) )
		  {
		    double dm1 = fabs(ClassifiedZColl[i].mass()-mZ0);
		    double dm2 = fabs(ClassifiedZColl[j].mass()-mZ0);
		    
		    if(dm1 < dm2 )
		      {const reco::CompositeCandidate Zee_tmp=ClassifiedZColl[i];
		      ClassifiedZColl[i] = ClassifiedZColl[j];
		      ClassifiedZColl[j] = Zee_tmp;
		      }
		  }
	    }
	  }
	}//fin algo de tri des Z
    }
   return ClassifiedZColl;
}

//Electron Classification using TDR criteria
int 
MmZeeAnalyser::ZeeConfLvlFromTDRElec(const reco::CompositeCandidate& Zee){

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

//Electron Classification using EWK note criteria       *********************** doesn't work yet ************************************
int MmZeeAnalyser::ZeeConfLvlFromEWKElec(const reco::CompositeCandidate& Zee,
					 const reco::GsfTrackCollection& TrackColl){
  float EOP1=0,EOP2=0;
  float IsoTrack1=0.,IsoTrack2=0.;
  int NbTrack1=0,NbTrack2=0;
  bool EOP[2] = {false,false};
  bool HOE[2] = {false,false};
  bool IsoTrack[2] = {false,false};
  bool NTrack[2] = {false,false};
  
  bool Electron[2]={false,false};
  
  if(Zee.numberOfDaughters()!=2)
    {
      cout<<" Bad reconstruction : number of daughter not egal to two "<<endl;
      return -1;
    }
  const pat::Electron *originalElectron1 = dynamic_cast<const pat::Electron *>(Zee.daughter(0)->masterClone().get());
  const pat::Electron *originalElectron2 = dynamic_cast<const pat::Electron *>(Zee.daughter(1)->masterClone().get());
  
  reco::SuperClusterRef SC1 = originalElectron1->superCluster();
  reco::SuperClusterRef SC2 = originalElectron2->superCluster();
  reco::GsfTrackRef Trk1 = originalElectron1->gsfTrack();
  reco::GsfTrackRef Trk2 = originalElectron2->gsfTrack();
  
  EOP1 = SC1->rawEnergy()/Trk1->outerP();
  EOP2 = SC2->rawEnergy()/Trk2->outerP();
 
  if(EOP1<3.0 && EOP1>0.75) {EOP[0]=true;} //Condition sur EOP
  if(EOP2<3.0 && EOP2>0.75) {EOP[1]=true;}
   
  for(reco::GsfTrackCollection::const_iterator itTrack = TrackColl.begin();
      itTrack !=TrackColl.end(); itTrack++)
    {
      if(itTrack->pt()>2 && itTrack->found()>5)
	{ 
	  float Dr1 = mmDeltaR(itTrack->eta(),itTrack->phi(),
			 originalElectron1->eta(),originalElectron1->phi(),
			 itTrack->pt(),originalElectron1->pt(), false);
	  float Dr2 = mmDeltaR(itTrack->eta(),itTrack->phi(),
			  originalElectron2->eta(),originalElectron2->phi(),
			  itTrack->pt(),originalElectron2->pt(), false);
	  if(Dr1<0.15) {NbTrack1++;}
	  if(Dr2<0.15) {NbTrack2++;}
	  
	  if(Dr1<0.2)
	    {IsoTrack1 += itTrack->pt() -
	       ConversionEpt(SC1->rawEnergy(),SC1->eta());}
	  if(Dr2<0.2)
	    {IsoTrack2 += itTrack->pt() -
	       ConversionEpt(SC2->rawEnergy(),SC1->eta());}
	}
    }
 
  IsoTrack1 /= ConversionEpt(SC1->rawEnergy(),SC1->eta());
  IsoTrack2 /= ConversionEpt(SC2->rawEnergy(),SC1->eta());
  
  if(NbTrack1<3) {NTrack[0]=true;} //Condition sur le nombre de traces
  if(NbTrack2<3) {NTrack[1]=true;}
  if(IsoTrack1<0.34) {IsoTrack[0]=true;} //Condition sur l'isolation des
  if(IsoTrack1<0.34) {IsoTrack[1]=true;} //electrons
   
  /* double EEnergy=0,HEnergy=0;
  cout<<" SC1 "<<SC1->caloID().detectors()<<endl;
  for(reco::basicCluster_iterator itCluster=SC1->clustersBegin();
      itCluster!=SC1->clustersEnd(); itCluster++)
    { cout<<" IDCLuster1 "<<(*itCluster)->caloID().detectors()
	  <<" energy "<<(*itCluster)->energy()<<endl;
      if((*itCluster)->caloID().detector(reco::CaloID::DET_ECAL_BARREL) || 
	 (*itCluster)->caloID().detector(reco::CaloID::DET_ECAL_ENDCAP)  )
	{ EEnergy += (*itCluster)->energy(); cout<<" eenergy "<<EEnergy<<endl;}
      if((*itCluster)->caloID().detector(reco::CaloID::DET_HCAL_BARREL) || 
	 (*itCluster)->caloID().detector(reco::CaloID::DET_HCAL_ENDCAP)  )
	{ HEnergy += (*itCluster)->energy();cout<<" henergy "<<HEnergy<<endl;}
    }
  if( (HEnergy/EEnergy)<0.08) {HOE[0]=true;} //Condition sur HOE
  cout<<" HOE1 "<<HEnergy/EEnergy<<" ";
  EEnergy=0;HEnergy=0;
  for(reco::basicCluster_iterator itCluster=SC2->clustersBegin();
      itCluster!=SC2->clustersEnd(); itCluster++)
    {  cout<<" IDCLuster2 "<<(*itCluster)->caloID().detectors()
	   <<" energy "<<(*itCluster)->energy()<<endl;
      if((*itCluster)->caloID().detector(reco::CaloID::DET_ECAL_BARREL) || 
	 (*itCluster)->caloID().detector(reco::CaloID::DET_ECAL_ENDCAP)  )
	{ EEnergy += (*itCluster)->energy();cout<<" eenergy2 "<<EEnergy<<endl;}
      if((*itCluster)->caloID().detector(reco::CaloID::DET_HCAL_BARREL) || 
	 (*itCluster)->caloID().detector(reco::CaloID::DET_HCAL_ENDCAP)  )
	{ HEnergy += (*itCluster)->energy();cout<<" henergy2 "<<HEnergy<<endl;}
    }
  cout<<" HOE2 "<<HEnergy/EEnergy<<endl;
  if( (HEnergy/EEnergy)<0.08) {HOE[1]=true;} //Condition sur HOE*/

 if( originalElectron1->hadronicOverEm()<0.08) {HOE[0]=true;}
 if( originalElectron2->hadronicOverEm()<0.08) {HOE[1]=true;}

 for(int i=0;i<2;i++)  //Validation des electrons
    {if(HOE[i] && EOP[i] && NTrack[i] && IsoTrack[i]) Electron[i]=true;}
  if(Electron[0] && Electron[1]) // electrons  valides -> 4
    return 4;
  else if( (!Electron[0] && Electron[1]) || //Un electron valide
	   (Electron[0] && !Electron[1])  )  //un autre non reconnu -> 1
    return 1;
  else return 0; //Aucun electron reconnu
    
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

// Zee MC validation

bool MmZeeAnalyser::MCTruth(const reco::GenParticleCollection& MCElectronColl,
			    const reco::CompositeCandidate& Zee,
			    double deltaPt[2],
			    double& deltaZMass,int& code){
  using namespace reco;
  
  bool mctruth=false;

  code=-2;

  const pat::Electron *originalElectron1 = dynamic_cast<const pat::Electron *>(Zee.daughter(0)->masterClone().get());
  const pat::Electron *originalElectron2 = dynamic_cast<const pat::Electron *>(Zee.daughter(1)->masterClone().get());
  
  reco::GenParticleRef match1;
  reco::GenParticleRef match2;
  
    if(originalElectron1->genParticlesSize()!=0)
    match1 = originalElectron1->genParticleRef();
    else return false;
 
    if(originalElectron2->genParticlesSize()!=0)
      match2 = originalElectron2->genParticleRef();
    else return false;
    if(!match1.isNull() && !match2.isNull()){
    if( (match1->pdgId()==originalElectron1->pdgId())
	&& (match2->pdgId()==originalElectron2->pdgId()) )
      {
	if( (match1->numberOfMothers()!=0) && (match2->numberOfMothers()!=0))
	  {
	    const Candidate * mother1 = match1->mother();
	    const Candidate * mother2 = match2->mother();
	    if( (mother1->numberOfMothers()!=0)&&
		(mother2->numberOfMothers()!=0) )
	      {
		const Candidate * mother_lvl2_1 = mother1->mother();
		const Candidate * mother_lvl2_2 = mother2->mother();
		if((mother_lvl2_1->pdgId()==23)&&(mother_lvl2_2->pdgId()==23))
		  {
		    deltaPt[0] = abs(originalElectron1->pt() - match1->pt());
		    deltaPt[1] = abs(originalElectron2->pt() - match2->pt());
		    deltaZMass = Zee.mass()-mother_lvl2_1->mass();
		    mctruth=true; code=0;
		    fill(MCTruth_Mass, mother_lvl2_1->mass());
		  }
		else
		  mctruth=false;
	      }
	    else {code=3;mctruth=false;}
	  }
	else mctruth=false;
      }
    else {code=2;mctruth=false;}
    }
    else {
      code=-1; mctruth=false;
    

      bool cand1=false,cand2=false;
      double ZMCmass=-1;

      for(reco::GenParticleCollection::const_iterator itMC = MCElectronColl.begin(); itMC != MCElectronColl.end(); itMC++)
	{
	  const Candidate * mother = itMC->mother();
	 	    
	  if(abs(itMC->pdgId())==11 && mother->pdgId()==23)
	    {if(match1.isNull()){
		double Dr=mmDeltaR(originalElectron1->eta(),
				   originalElectron1->phi(),
				   itMC->eta(),itMC->phi(),
				   originalElectron1->pt(),itMC->pt(),true);
		if(Dr<0.4 && originalElectron1->pdgId()==itMC->pdgId() )
		  { cand1 = true;}
	       }
	    else cand1=true;
	    
	    if(match2.isNull()){
	      double Dr=mmDeltaR(originalElectron2->eta(),
				 originalElectron2->phi(),
				 itMC->eta(),itMC->phi(),
				 originalElectron1->pt(),itMC->pt(),true);
	      if(Dr<0.4 && originalElectron2->pdgId()==itMC->pdgId())
		{cand2 = true;}
	      }
	      else cand2 =true;

	    ZMCmass=mother->mass();
	    }
	}

      if(cand1 && cand2)  {
	mctruth=true;code=8;  fill(Mass_no_MCT,Zee.mass());
	deltaZMass = Zee.mass()-ZMCmass;
      }
      else   fill(Mass_no_MCT2,Zee.mass());
    }
    if(!mctruth)
      {
	deltaPt[0]=-1;
	deltaPt[1]=-1;
	deltaZMass=-1;  
      }
    return mctruth;
}


//Electron MC validation
void MmZeeAnalyser::ElectronTruth(const edm::View<pat::Electron>& ElectronColl,
				  const reco::GsfTrackCollection& TrackColl,
				  const reco::GenParticleCollection& MCElectronColl, int& code){
  
  using namespace reco;
  int Nelec=-10;
  double eta1=0;
  double eta2=0;
  double pt1=0,pt2=0;
  double mass=0;

  code=10;

  for(reco::GenParticleCollection::const_iterator itMC = MCElectronColl.begin(); itMC != MCElectronColl.end(); itMC++)
    {
      if(abs(itMC->pdgId())==11)
	{Nelec++;
	if(itMC->numberOfMothers()!=0){
	  const Candidate * mother = itMC->mother();
	  if(mother->pdgId()==23)
	    {	    
	    mass=mother->mass();
	    if(itMC->pdgId()==-11)
	      {eta1=itMC->eta();
	      pt1=itMC->pt();
	      }
	    if(itMC->pdgId()==11)
	      {eta2=itMC->eta();
	      pt2=itMC->pt();
	      }
	    }
	}
	else cout<<"pas de mère"<<endl;
	}
    }
  fill(Pt_MCT,pt1,pt2);
  fill(Eta_MCT,eta1,eta2);
  if(abs(eta1)>2.45 || abs(eta2)>2.45 || 
     (eta1<1.55 && eta1>1.45) || (eta1>-1.55 && eta1<-1.45)  || 
     (eta2<1.55 && eta2>1.45) || (eta2>-1.55 && eta2<-1.45)  )
    {code=1;}
  else if(pt1<10 || pt2<10) {code=2; }
  else if(mass<30) {code=3; }
  else if(ElectronColl.size()>=2){
  for(int unsigned it =0; it < ElectronColl.size(); it++)
    {const pat::Electron itE = ElectronColl[it];
    for(reco::GenParticleCollection::const_iterator itMC = MCElectronColl.begin(); itMC != MCElectronColl.end(); itMC++)
      {
	if(abs(itMC->pdgId())==11)
	  {double Dr=mmDeltaR(itE.eta(),itE.phi(),itMC->eta(),itMC->phi(),itE.pt(),itMC->pt(),false);
	  if(Dr<0.5 && itE.pdgId()!=itMC->pdgId())
	    {code=4;}
	  }
      }
    }
  }
  else if(code != 4){
    for(reco::GenParticleCollection::const_iterator itMC = MCElectronColl.begin(); itMC != MCElectronColl.end(); itMC++)
      { int NTrack=0;
      const Candidate * mother = itMC->mother();
      if(abs(itMC->pdgId())==11 && mother->pdgId()==23){
	for(reco::GsfTrackCollection::const_iterator itTrack = TrackColl.begin(); itTrack != TrackColl.end();itTrack++)
	  {double Dr=mmDeltaR(itTrack->outerEta(),itTrack->outerPhi(),itMC->eta(),itMC->phi(),itTrack->outerPt(),itMC->pt(),false);
	  if(Dr<0.5)
	    {NTrack++;
	    if(itTrack->outerP()/itMC->p()<0.8) {code=5;fill(Nelectron,abs(itMC->eta()));} //pas sûr de moi sur ce coup là
	    }
	  }
	if(NTrack>5) {code=6;}     
      }
      }
  }
  return;
}


//Filling Electrons Histos
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

  double Dphi = abs((originalElectron1->phi()*180/myPi)
		    -(originalElectron2->phi()*180/myPi));

  fill(dPhi_dau,Dphi);
  fill( PtvsdPhi_Z,Zee.pt(),Dphi);

  double EOP[2] ={-10,-10};
  
  reco::SuperClusterRef SC1 = originalElectron1->superCluster();
  reco::SuperClusterRef SC2 = originalElectron2->superCluster();

  reco::GsfTrackRef Trk1 = originalElectron1->gsfTrack();
  reco::GsfTrackRef Trk2 = originalElectron2->gsfTrack();

  fill(SuperCluster_etaWidth,SC1->etaWidth());
  fill(SuperCluster_phiWidth,SC1->phiWidth());

  if(SC1.isNull()) cout<<" SC1 null ";
  if(SC2.isNull()) cout<<" SC2 null ";
  if(Trk1.isNull()) cout<<" Trk1 null ";
  if(Trk2.isNull()) cout<<" Trk2 null ";
 
  EOP[0] = SC1->rawEnergy()/Trk1->outerP();
  EOP[1] = SC2->rawEnergy()/Trk2->outerP();

  double RecoEnergySC[2]={-1,-1};
  RecoEnergySC[0]=(sin(ConversionEtaTheta(SC1->eta())) +
		   cos(ConversionEtaTheta(SC1->eta()))/cos(SC1->phiWidth()) )*SC1->rawEnergy();
  RecoEnergySC[1]=(sin(ConversionEtaTheta(SC2->eta())) +
		 cos(ConversionEtaTheta(SC2->eta()))/cos(SC2->phiWidth()) )*SC2->rawEnergy();
  
  float Cluster1[3] = {SC1->x(),SC1->y(),SC1->z()};
  float Cluster2[3] = {SC2->x(),SC2->y(),SC2->z()};

  double ZmassSC = -1;
  ZmassSC = sqrt(2*RecoEnergySC[0]*RecoEnergySC[1]*
		 (1-cos(Angle_2vecteurs(Cluster1,Cluster2,true))) );

  fill(EOP_electron,EOP[0]);fill(EOP_electron,EOP[1]);
  fill(Zmass_bySC,ZmassSC);
  return;
}
  


//define this as a plug-in
DEFINE_FWK_MODULE(MmZeeAnalyser);


#endif
