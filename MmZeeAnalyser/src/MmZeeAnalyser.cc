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
// $Id: MmZeeAnalyser.cc,v 1.8 2009/03/31 15:01:29 mmarionn Exp $
//
//
#include "MMarionneau/MmZeeAnalyser/interface/MmZeeAnalyser.h"
//#include "MMarionneau/MmZeeAnalyser/interface/MmZeeAnalysis.h"
//#include "MMarionneau/MmZeeAnalyser/interface/MmZmumuAnalysis.h"

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

#include "DataFormats/PatCandidates/interface/Lepton.h"
#include <DataFormats/PatCandidates/interface/Particle.h>

#include "PhysicsTools/CandUtils/interface/CandMatcher.h"
#include "RecoEgamma/EgammaTools/interface/HoECalculator.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

//#include "MMarionneau/MmZeeAnalyser/Utilitaires/Utilitaires_Geom.cc"


using namespace edm;
using namespace std;
using namespace cms;



// Constructors

MmZeeAnalyser::MmZeeAnalyser(const edm::ParameterSet& iConfig):

  electronCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronCollection")),
  muonCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("muonCollection")),
  ZeeCandCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("ZeeCandidateCollection")),
  ZmumuCandCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("ZmumuCandidateCollection")),
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



  DecMod_[0]="Zee";
  DecMod_[1]="Zmumu";

  for(int imode=0;imode<2;imode++)
    {
      string DM = DecMod_[imode];
  
    
  MCTruth_Mass[DM] =  book1D("MassofZMC_"+DM,
			 "mass of Z from MCTruth;"
			 "Mass;Nevts",
			 200, 0., 200.);


  Mass_Z_HCL[DM] = book1D("MassofZforHCL_"+DM,
		     "mass of Z reconstructed with CL of 4, 3 or 2;"
		     "Mass;Nevts",
		     200, 0., 200.);

  Mass_Z_LCL[DM] = book1D("MassofZforLCL_"+DM,
		     "mass of Z reconstructed with CL of 1;"
		     "Mass;Nevts",
		     200, 0., 200.);

  Pt_Z_HCL[DM] = book1D("PtofZforHCL_"+DM,
		    "Pt of Z reconstructed with CL of 4, 3 or 2;"
		    "Pt;Nevts",
		    100, 0., 200.);
  Eta_Z_HCL[DM] = book1D("EtaofZforHCL_"+DM,
		     "Eta of Z reconstructed with CL of 4, 3 or 2;"
		     "Eta;Nevts",
		     80, -4., 4.);
  Phi_Z_HCL[DM] = book1D("PhiofZforHCL_"+DM,
		     "Phi of Z reconstructed with CL of 4, 3 or 2;"
		     "Phi;Nevts",
		     64, -3.2, 3.2);
  
  LevelConf_Z[DM] = book1D("LevelConfforZ_"+DM,
		       "Confidence Level of Z;"
		       "CL;Nevts",
		       5,-0.5,4.5);

  MCFalse_LevelConf_Z[DM] = book2D("MCFalseLevelConfforZ_EColS_"+DM,
			"Confidence Level of Z if MC->False;"
		       "CL;NbLepton",
		       5,-0.5,4.5,11,-0.5,10.5);

  MCFalse_CL_ColS_Z[DM] = book2D("MCFalseLevelConfforZ_ColS_"+DM,
			"Confidence Level of Z if MC->False;"
		       "CL;Place dans la liste",
		       5,-0.5,4.5,11,-0.5,10.5);
  
  PtvsdPhi_Z[DM] = book2D("PtVsdPhiforZ_"+DM,
		      "Pt vs dphi of Z;"
		      "Pt;dPhi",
		      100,0.,200.,92,-3.,183.);
  
  Pt_dau_0[DM] = book1D("Ptofdaughter_0_"+DM,
		     "Pt of the daughter 1;"
		     "Pt;Nevts",
		     100, 0., 200.);
 
  Eta_dau_0[DM] = book1D( "Etaofdaughter_0_"+DM,
		       "Eta of the daughter 1;"
		       "Eta;Nevts",
			  80, -4., 4.);
 
  Phi_dau_0[DM] = book1D("Phiofdaughter_0_"+DM,
		     "Phi of the daughter 1;"
		      "Phi;Nevts",
		      182, -183., 183.);
 
  Class_dau_0[DM] = book1D("Classificationofdaughter_0_"+DM,
		       "Classification of the daughter 1;"
		       "Classification;Nevts",
		       5, -0.5, 4.5);
  
  deltaPtMC_dau_0[DM] = book1D("deltaPtofdaughter_0_"+DM,
			   "#deltaPt mc/pat of the daughter 1;"
			   "Pt;Nevts",
			   100, 0., 20.);
  
  Pt_dau_1[DM] = book1D("Ptofdaughter_1_"+DM,
		    "Pt of the daughter 2;"
		    "Pt;Nevts",
		    100, 0., 200.);
  
  Eta_dau_1[DM] = book1D( "Etaofdaughter_1_"+DM,
		      "Eta of the daughter 2;"
		      "Eta;Nevts",
		      80, -4., 4.);
  
  Phi_dau_1[DM] = book1D("Phiofdaughter_1_"+DM,
		     "Phi of the daughter 2;"
		     "Phi;Nevts",
		      182, -183., 183.);
  
  Class_dau_1[DM] = book1D("Classificationofdaughter_1_"+DM,
		       "Classification of the daughter 2;"
		       "Classification;Nevts",
		       5, -0.5, 4.5);

  deltaPtMC_dau_1[DM] = book1D("deltaPtofdaughter_1_"+DM,
			   "deltaPt mc/pat of the daughter 2;"
			   "#deltaPt;Nevts",
			   100, 0., 20.);
  
  dPhi_dau[DM] =  book1D("dPhidaughters_"+DM,
		     "#Delta#Phi between daughters ;"
		     "dPhi;Nevts",
		     92, -3., 183);

  deltaMassMC_Z[DM] = book1D("deltaMassMC_Z_"+DM,
			 "#DeltaMass mc/pat of the Z candidate;"
			 "#DeltaM;NEvts",
			 160,-40.,40.);
  

  

  Eta_MCT[DM] = book2D( "Eta_MCT_"+DM,
		    "Eta of lepton;"
		    "Eta;Nevts",
		    120, -6., 6.,120, -6., 6.);
  
  Pt_MCT[DM] = book2D("Pt_MCT_"+DM,
		  "Pt of lepton;"
		  "Pt1;Pt2",
		  100, 0., 50.,100, 0., 50.);

  Eta_vs_E_MCT[DM] = book2D("Eta_cs_E_"+DM,
			"",
			120,-6., 6.,100,0.,200.);

  Nlepton[DM] = book1D("Nlepton_"+DM,
		     "",
		     120,0.,6);

  No_Zreco_code[DM] = book1D("No_Zreco_code_"+DM,
			 "Code corresponding to the non-tagged Z;"
			 "Code;Nevts",
			 11,-0.5,10.5);

  MCFalse_Z_code[DM] = book1D("MCFalse_Z_code_"+DM,
			 "Code corresponding to the Z without MCT;"
			 "Code;Nevts",
			 11,-1.5,9.5);

  Mass_no_MCT[DM] = book1D("MassofZnoMCT_saved_"+DM,
		     "mass of Z reconstructed for no MCT;"
		     "Mass;Nevts",
		     200, 0., 200.);

 Mass_no_MCT2[DM] = book1D("MassofZnoMCT_unsaved_"+DM,
		     "mass of Z reconstructed for no MCT;"
		     "Mass;Nevts",
		     200, 0., 200.);


}


 SuperCluster_etaWidth = book1D( "etaWidth",
				  "EtaWidth of supercluster;"
				  "Eta;Nevts",
				  100, 0., 0.2);

  SuperCluster_phiWidth = book1D( "phiWidth",
				  "PhiWidth of supercluster;"
				  "Phi;Nevts",
				  100, 0., 0.2);

  EOP_electron = book1D("EOP_electron_Zee",
			"Energy on momentum for leptons;"
			"EOP (GeV);Nevts",
			100,0.5,4.5);
  
  Zmass_bySC = book1D("MassofZbySC",
		      "mass of Z reconstructed with the SC data;"
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

  cout<<" new event "<<endl;

 vector< reco::CompositeCandidateCollection > ZCandidates;
 // vector< edm::View< pat::Particle> > Leptons;

  //Collection reading    
   edm::Handle<edm::View<pat::Electron> > ElectronsHandle;
   iEvent.getByLabel(electronCollectionTag_, ElectronsHandle);
   edm::View<pat::Electron> Electrons =  *ElectronsHandle ;
   // Leptons.push_back(  *ElectronsHandle );

   edm::Handle<edm::View<pat::Muon> > MuonsHandle;
   iEvent.getByLabel(muonCollectionTag_, MuonsHandle);
   edm::View<pat::Muon> Muons =  *MuonsHandle ;
   //Leptons.push_back(  *MuonsHandle );

   edm::Handle<reco::CompositeCandidateCollection> ZeeCandidatesH;
   iEvent.getByLabel(ZeeCandCollectionTag_, ZeeCandidatesH);
   ZCandidates.push_back( *ZeeCandidatesH );

   edm::Handle<reco::CompositeCandidateCollection> ZmumuCandidatesH;
   iEvent.getByLabel(ZmumuCandCollectionTag_, ZmumuCandidatesH);
   ZCandidates.push_back( *ZmumuCandidatesH);

   edm::Handle<reco::GsfTrackCollection> GsfTrackH;
   iEvent.getByLabel(trackCollectionTag_, GsfTrackH);
   reco::GsfTrackCollection Tracks = *GsfTrackH;

    edm::Handle<reco::GenParticleCollection> GenParticlesH;
   iEvent.getByLabel(genParticlesTag_, GenParticlesH);
   reco::GenParticleCollection GenParticles = *GenParticlesH;
   
  

   //liste classifiee des candidats  
   //  vector<reco::CompositeCandidate> ClassifiedZColl;
   map <string, vector<reco::CompositeCandidate> > ClassifiedZColl;
   //liste des niveaux de confiance de ces candidats
   map <string, vector< int > > ClassifiedZCL;
   vector< int > ClassifiedZCL_tmp;
   //cout<<"test_1"<<endl;
   for(int imode=0;imode<2;imode++)
     {
   

       string DM =DecMod_[imode];
      
       ClassifiedZColl[ DM ] = ClassZCandidate(ZCandidates[imode],Tracks, DM);
      
       if(analyse_type_ =="TDR")
	     for(int unsigned iZ=0; iZ<ClassifiedZColl[ DM ].size();iZ++)
	       ClassifiedZCL_tmp.push_back( ZConfLvlFromTDR(ClassifiedZColl[DM][iZ], DM) );

       ClassifiedZCL[ DM ] = ClassifiedZCL_tmp;

       //Matching MCTruth
       if(ClassifiedZColl[DM].size()==0)
	 {
	   int code=-1;
	   
	   LeptonTruth( Electrons,Muons,Tracks,GenParticles,code, DM);
	   
	   fill(No_Zreco_code[DM],code);
	   NoZCandidate_[imode]++;
	 }
       else{
	 double deltaPt[2] = {0,0};
	 double deltaZMass = 0;
	 
	 int code=-10;
	  
	 bool mctruth =  MCTruth(GenParticles,
				 ClassifiedZColl[DM][0],
				 deltaPt, deltaZMass,code, DM);

	 LeptonFromZAnalysis(ClassifiedZColl[DM][0],DM);

	 if(mctruth)
	   {
	     GoodCandidate_[imode]++;
	     fill(deltaMassMC_Z[DM],deltaZMass);
	     fill(deltaPtMC_dau_0[DM],deltaPt[0]);
	     fill(deltaPtMC_dau_1[DM],deltaPt[1]);
	     
	     fill(LevelConf_Z[DM], ClassifiedZCL[DM][0]);
	    
	   }
	 else 
	   {
	   if(DM==DecMod_[0])
	     fill(MCFalse_LevelConf_Z[DM], ClassifiedZCL[DM][0] ,Electrons.size());
	   if(DM==DecMod_[1])
	     fill(MCFalse_LevelConf_Z[DM], ClassifiedZCL[DM][1] ,Muons.size());
	   fill(MCFalse_CL_ColS_Z[DM],  ClassifiedZCL[DM][0], 1);
	   fill(MCFalse_Z_code[DM],code);
	  
	   }
	 
	 if( ClassifiedZCL[DM][0] >1)
	   {
	    
	     
	     fill(Mass_Z_HCL[DM], ClassifiedZColl[DM][0].mass());
	     fill(Pt_Z_HCL[DM], ClassifiedZColl[DM][0].pt()); 
	     fill(Eta_Z_HCL[DM], ClassifiedZColl[DM][0].eta()); 
	     fill(Phi_Z_HCL[DM], ClassifiedZColl[DM][0].phi()); 
	    
	   }
	 else
	   {
	     fill(Mass_Z_LCL[DM],ClassifiedZColl[DM][0].mass());
	   }
       }
     }
}

void 
MmZeeAnalyser::beginJob(const edm::EventSetup&)
{
  EventNumber_=0;
  GoodCandidate_[0]=0;
  GoodCandidate_[1]=0;
  NoZCandidate_[0]=0;
  NoZCandidate_[1]=0;
}


void 
MmZeeAnalyser::endJob() {

  

  cout<<endl;
  cout<<"Number of event : "<<EventNumber_<<endl;
  cout<<"Number of selected Candidates "<<EventNumber_-NoZCandidate_[0]<<endl;
  cout<<"Number of good candidates (MCT) : "<<GoodCandidate_[0]<<endl;
  cout<<" total efficiency "<<(float)GoodCandidate_[0]/EventNumber_
      <<" partial efficiency "<<(float)GoodCandidate_[0]/(EventNumber_-NoZCandidate_[0])<<endl;
  cout<<"Reconstruction loose "<<NoZCandidate_[0]<<" ->% "
      <<(float)NoZCandidate_[0]/EventNumber_<<endl;
  

}


//--------------------- classing Z candidate in the event -----------------------
 vector<reco::CompositeCandidate>
MmZeeAnalyser::ClassZCandidate(const  reco::CompositeCandidateCollection& ZCandidates,
			     const reco::GsfTrackCollection Tracks, string channel)
{
  vector<reco::CompositeCandidate> ClassifiedZColl ;
  if(analyse_type_ == "TDR"){
    if(ZCandidates.size()==0) //pas de candidat Z
      {
	// Warning, no reconstructed Z in this event
	return ClassifiedZColl;
      }
    else if(ZCandidates.size()==1) // un seul candidat Z
      {
	if(ZConfLvlFromTDR(ZCandidates[0],channel) < 1 )
	  {
	    //Warning, Confidence level on electrons is really bad
	    //but alone candidate
	    ClassifiedZColl.push_back(ZCandidates[0]);
	  }
	else
	  {
	    ClassifiedZColl.push_back(ZCandidates[0]);
	  }
      }
    else     //plusieurs candidats Z, algo de tri necessaire
      {      // base sur le niveau de confiance des electrons, puis sur la masse
	int iZ=0; int goodZCand=0;
	for(reco::CompositeCandidateCollection::const_iterator itZ = ZCandidates.begin(); itZ != ZCandidates.end(); itZ++, iZ++)
	  {
	    const reco::CompositeCandidate& Z = *itZ;
	    
	    if(ZConfLvlFromTDR(Z,channel) < 1 )
	      {//Warning, Confidence level on electrons is really bad
	      // candidate not classified
		continue;}
	    else
	      {ClassifiedZColl.push_back(Z);
	      goodZCand++;}
	  }
	for(int i=0;i<goodZCand;i++){
	  for(int j=0;j<goodZCand;j++)
	    {
	      if(i != j)
		if( ZConfLvlFromTDR(ClassifiedZColl[i],channel)
		    > ZConfLvlFromTDR(ClassifiedZColl[j],channel) )
		  {const reco::CompositeCandidate Z_tmp = ClassifiedZColl[i];
		  ClassifiedZColl[i] = ClassifiedZColl[j];
		  ClassifiedZColl[j] = Z_tmp;
		  }
		else if( ZConfLvlFromTDR(ClassifiedZColl[i],channel)
			 == ZConfLvlFromTDR(ClassifiedZColl[j], channel) )
		  {
		    double dm1 = fabs(ClassifiedZColl[i].mass()-mZ0);
		    double dm2 = fabs(ClassifiedZColl[j].mass()-mZ0);
		    
		    if(dm1 < dm2 )
		      {const reco::CompositeCandidate Z_tmp=ClassifiedZColl[i];
		      ClassifiedZColl[i] = ClassifiedZColl[j];
			ClassifiedZColl[j] = Z_tmp;
		      }
		  }
	    }
	}//fin algo de tri des Z
      }
  }
 
   return ClassifiedZColl;
}

// ---------- Lepton Classification using TDR criteria -------------
int MmZeeAnalyser::ZConfLvlFromTDR(const reco::CompositeCandidate& Z, string channel){

  // MmZeeAnalysis  ZeeAnalysis;
  // MmZmumuAnalysis  ZmumuAnalysis;

  if(channel=="Zmumu")
    return ZmumuConfLvlFromTDR(Z);
  else if(channel=="Zee")
    return ZeeConfLvlFromTDR(Z);
  else return 0;
}


//--------------- Z MC validation ----------------------
bool MmZeeAnalyser::MCTruth(const reco::GenParticleCollection& MCLeptonColl,
			    const reco::CompositeCandidate& Z,
			    double deltaPt[2],
			    double& deltaZMass,int& code,string channel){
  using namespace reco;
  
  bool mctruth=false;

  code=-2;
   
  reco::GenParticleRef match1;
  reco::GenParticleRef match2;


  const reco::Candidate &cand1 = *Z.daughter(0)->masterClone();
  const reco::Candidate &cand2 = *Z.daughter(1)->masterClone();
  
  if (cand1.isMuon() && cand2.isMuon()) {
    const pat::Muon & Mu1 = dynamic_cast<const pat::Muon &>(cand1);
    const pat::Muon & Mu2 = dynamic_cast<const pat::Muon &>(cand2);
    
    if(Mu1.genParticlesSize()!=0)
      {
      match1 = Mu1.genParticleRef();}
    else return false;
    if(Mu2.genParticlesSize()!=0)
      match2 = Mu2.genParticleRef();
    else return false;
    
    
  } else if (cand1.isElectron() && cand2.isElectron()) {
    
    const pat::Electron & Elec1 = dynamic_cast<const pat::Electron &>(cand1);
    const pat::Electron & Elec2 = dynamic_cast<const pat::Electron &>(cand2);
    
    if(Elec1.genParticlesSize()!=0)
      {
      match1 = Elec1.genParticleRef();}
    else return false;
    if(Elec2.genParticlesSize()!=0)
      match2 = Elec2.genParticleRef();
    else return false;
  }
  else{ cout<<"Error, Z shaped with l and l' "<<endl;    return false;}

    if(!match1.isNull() && !match2.isNull()){
    if( (match1->pdgId()==cand1.pdgId())
	&& (match2->pdgId()==cand2.pdgId()) )
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
		    deltaPt[0] = abs(cand1.pt() - match1->pt());
		    deltaPt[1] = abs(cand2.pt() - match2->pt());
		    deltaZMass = Z.mass()-mother_lvl2_1->mass();
		    mctruth=true; code=0;
		    fill(MCTruth_Mass[channel], mother_lvl2_1->mass());
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
    

      bool candMC1=false,candMC2=false;
      double ZMCmass=-1;

      for(reco::GenParticleCollection::const_iterator itMC = MCLeptonColl.begin(); itMC != MCLeptonColl.end(); itMC++)
	{
	  const Candidate * mother = itMC->mother();
	 	    
	  if(abs(itMC->pdgId())==11 && mother->pdgId()==23)
	    {if(match1.isNull()){
		double Dr=mmDeltaR(cand1.eta(),
				   cand1.phi(),
				   itMC->eta(),itMC->phi(),
				   cand1.pt(),itMC->pt(),true);
		if(Dr<0.4 && cand1.pdgId()==itMC->pdgId() )
		  { candMC1 = true;}
	       }
	    else candMC1=true;
	    
	    if(match2.isNull()){
	      double Dr=mmDeltaR(cand2.eta(),
				 cand2.phi(),
				 itMC->eta(),itMC->phi(),
				 cand1.pt(),itMC->pt(),true);
	      if(Dr<0.4 && cand2.pdgId()==itMC->pdgId())
		{candMC2 = true;}
	      }
	      else candMC2 =true;

	    ZMCmass=mother->mass();
	    }
	}

      if(candMC1 && candMC2)  {
	mctruth=true;code=8;  fill(Mass_no_MCT[channel],Z.mass());
	deltaZMass = Z.mass()-ZMCmass;
      }
      else   fill(Mass_no_MCT2[channel],Z.mass());
    }
    if(!mctruth)
      {
	deltaPt[0]=-1;
	deltaPt[1]=-1;
	deltaZMass=-1;  
      }
    return mctruth;
}


//------------------ Lepton MC validation ----------------------------------
void MmZeeAnalyser::LeptonTruth(const edm::View<pat::Electron>& ElectronColl,
				const edm::View<pat::Muon>& MuonColl,
				const reco::GsfTrackCollection& TrackColl,
				const reco::GenParticleCollection& MCParticleColl,
				int& code, string channel){
  
  int pdgID=0;
  if(channel==DecMod_[0]) pdgID=11;
  if(channel==DecMod_[1]) pdgID=13;

  using namespace reco;
  int Nlep=-10;
  double eta1=0;
  double eta2=0;
  double pt1=0,pt2=0;
  double mass=0;

  code=10;

  for(reco::GenParticleCollection::const_iterator itMC = MCParticleColl.begin(); itMC != MCParticleColl.end(); itMC++)
    {
      if(abs(itMC->pdgId())==pdgID)
	{Nlep++;
	if(itMC->numberOfMothers()!=0){
	  const Candidate * mother = itMC->mother();
	  if(mother->pdgId()==23)
	    {	    
	    mass=mother->mass();
	    if(itMC->pdgId()==-pdgID)
	      {eta1=itMC->eta();
	      pt1=itMC->pt();
	      }
	    if(itMC->pdgId()==pdgID)
	      {eta2=itMC->eta();
	      pt2=itMC->pt();
	      }
	    }
	}
	else cout<<"pas de mère"<<endl;
	}
    }
  fill(Pt_MCT[channel],pt1,pt2);
  fill(Eta_MCT[channel],eta1,eta2);
  if(abs(eta1)>2.45 || abs(eta2)>2.45 || 
     (eta1<1.55 && eta1>1.45) || (eta1>-1.55 && eta1<-1.45)  || 
     (eta2<1.55 && eta2>1.45) || (eta2>-1.55 && eta2<-1.45)  )
    {code=1;}
  else if(pt1<10 || pt2<10) {code=2; }
  else if(mass<30) {code=3; }
  else if(channel==DecMod_[0] && ElectronColl.size()>=2){   //Electrons
  for(int unsigned it =0; it < ElectronColl.size(); it++)
    {const pat::Electron itE = ElectronColl[it];
    for(reco::GenParticleCollection::const_iterator itMC = MCParticleColl.begin(); itMC != MCParticleColl.end(); itMC++)
      {
	if(abs(itMC->pdgId())==pdgID)
	  {double Dr=mmDeltaR(itE.eta(),itE.phi(),itMC->eta(),itMC->phi(),itE.pt(),itMC->pt(),false);
	  if(Dr<0.5 && itE.pdgId()!=itMC->pdgId())
	    {code=4;}
	  }
      }
    }
  }
  else if(channel==DecMod_[1] && MuonColl.size()>=2){    //Muons
  for(int unsigned it =0; it < MuonColl.size(); it++)
    {const pat::Muon itE = MuonColl[it];
    for(reco::GenParticleCollection::const_iterator itMC = MCParticleColl.begin(); itMC != MCParticleColl.end(); itMC++)
      {
	if(abs(itMC->pdgId())==pdgID)
	  {double Dr=mmDeltaR(itE.eta(),itE.phi(),itMC->eta(),itMC->phi(),itE.pt(),itMC->pt(),false);
	  if(Dr<0.5 && itE.pdgId()!=itMC->pdgId())
	    {code=4;}
	  }
      }
    }
  }
  else if(code != 4){
    for(reco::GenParticleCollection::const_iterator itMC = MCParticleColl.begin(); itMC != MCParticleColl.end(); itMC++)
      { int NTrack=0;
      const Candidate * mother = itMC->mother();
      if(abs(itMC->pdgId())==pdgID && mother->pdgId()==23){
	for(reco::GsfTrackCollection::const_iterator itTrack = TrackColl.begin(); itTrack != TrackColl.end();itTrack++)
	  {double Dr=mmDeltaR(itTrack->outerEta(),itTrack->outerPhi(),itMC->eta(),itMC->phi(),itTrack->outerPt(),itMC->pt(),false);
	  if(Dr<0.5)
	    {NTrack++;
	    if(itTrack->outerP()/itMC->p()<0.8) {code=5;fill(Nlepton[channel],abs(itMC->eta()));} //pas sûr de moi sur ce coup là
	    }
	  }
	if(NTrack>5) {code=6;}     
      }
      }
  }
  return;
}

// ---------------------- Lepton Variables ------------------------------
void MmZeeAnalyser::LeptonFromZAnalysis(const reco::CompositeCandidate& Z,
					  string channel){
  using namespace reco;

 const reco::Candidate &cand1 = *Z.daughter(0)->masterClone();
 const reco::Candidate &cand2 = *Z.daughter(1)->masterClone();
  
  fill( Pt_dau_0[channel],cand1.pt());
  fill( Eta_dau_0[channel],cand1.eta());
  fill( Phi_dau_0[channel],cand1.phi()*180/myPi);

  fill( Pt_dau_1[channel],cand2.pt());
  fill( Eta_dau_1[channel],cand2.eta());
  fill( Phi_dau_1[channel],cand2.phi()*180/myPi);

  double Dphi = abs((cand1.phi()*180/myPi)
		    -(cand2.phi()*180/myPi));

  fill(dPhi_dau[channel],Dphi);
  fill( PtvsdPhi_Z[channel],Z.pt(),Dphi);

  if(channel==DecMod_[0])
    {
      const pat::Electron & Elec1 = dynamic_cast<const pat::Electron &>(cand1);
      const pat::Electron & Elec2 = dynamic_cast<const pat::Electron &>(cand2);

      double EOP[2] ={-10,-10};
      
      reco::SuperClusterRef SC1 = Elec1.superCluster();
      reco::SuperClusterRef SC2 = Elec2.superCluster();
      
      reco::GsfTrackRef Trk1 = Elec1.gsfTrack();
      reco::GsfTrackRef Trk2 = Elec2.gsfTrack();
      
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
    }

  return;
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


//---------------- Particular Functions (no general)

//Electron Classification using TDR criteria
int 
MmZeeAnalyser::ZeeConfLvlFromTDR(const reco::CompositeCandidate& Zee){

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
 //fill(Class_dau_0,dau_Class[0]);
 //fill(Class_dau_1,dau_Class[1]);

 if(dau_Class[0]==3 && dau_Class[1]==3) //Two electrons are tights. ***
   return 4;
 else if(dau_Class[0]==1 && dau_Class[1]==3) //One electron is tight
   return 3;
 else if(dau_Class[1]==1 && dau_Class[0]==3) // and one is robust. ***
   return 3;
 else if(dau_Class[1]==1 && dau_Class[0]==1) //Two electrons are robust. ***
   return 2;
 else if(dau_Class[1]==0 && ( dau_Class[0]==3 || dau_Class[0]==1)) //One electron
   return 1;
 else if(dau_Class[0]==0 && ( dau_Class[1]==3 || dau_Class[1]==1))// is bad. ***
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

int 
MmZeeAnalyser::ZmumuConfLvlFromTDR(const reco::CompositeCandidate& Zmumu){

  float dau_Class[2] = {-1,-1};
 
  if(Zmumu.numberOfDaughters()!=2)
    {
      cout<<" Bad reconstruction : number of daughter not egal to two "<<endl;
      return -1;
    }
  const pat::Muon *originalMuon1 = dynamic_cast<const pat::Muon *>(Zmumu.daughter(0)->masterClone().get());
  const pat::Muon *originalMuon2 = dynamic_cast<const pat::Muon *>(Zmumu.daughter(1)->masterClone().get());
 

  /* dau_Class[0] = originalMuon1->leptonID("robust") + 
                2*originalMuon1->leptonID("tight") ;

 dau_Class[1] = originalMuon2->leptonID("robust") +
 2*originalMuon2->leptonID("tight") ;*/

  float global[2]={0,0};
  
  global[0] = (int)originalMuon1->isGlobalMuon() + IsoMuon(originalMuon1);
  global[1] = (int)originalMuon2->isGlobalMuon() + IsoMuon(originalMuon2);

 //Filling classification histos
 // fill(Class_dau_0,dau_Class[0]);
 //fill(Class_dau_1,dau_Class[1]);

 if(dau_Class[0]==3 && dau_Class[1]==3) //Two global isolated muons. ***
   return 4;
 else if(dau_Class[0]==1 && dau_Class[1]==3) //One global isolated Muon
   return 3;
 else if(dau_Class[1]==1 && dau_Class[0]==3) // and one global muon. ***
   return 3;
 else if(dau_Class[1]==1 && dau_Class[0]==1) //Two global muons. ***
   return 2;
 else if(dau_Class[1]==0 && ( dau_Class[0]==3 || dau_Class[0]==1)) //One muon is
   return 1;
 else if(dau_Class[0]==0 && ( dau_Class[1]==3 || dau_Class[1]==1))// not global and not isolated. ***
   return 1;
  else if(dau_Class[0]==0 && dau_Class[1]==0) //The two are bad. ***
    return 0;
 else
   {
     cout<<endl;
     cout<<" Bad recognition, there is no two Muons in this Z candidate"<<endl;
     return -1;
   }

}



int 
MmZeeAnalyser::IsoMuon(const pat::Muon*& muon) {

  reco::TrackRef Trk = muon->globalTrack();
  reco::MuonIsolation Iso = muon->isolationR03();

  float TrackIso = Iso.sumPt; 
  float CaloIso = Iso.hoEt + Iso.hadEt ;
  
  if(TrackIso < 0.92 && Trk->normalizedChi2() < 10 && CaloIso < 5)
    return 2;
  else return 0;

}


//define this as a plug-in
DEFINE_FWK_MODULE(MmZeeAnalyser);


#endif



