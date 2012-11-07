#include "LQUseTree.hh"

#include <iomanip>
#include <TLorentzVector.h>
#include <boost/timer.hpp>
#include <boost/progress.hpp>

#include <TH1I.h>
#include <TObjArray.h>

using namespace std;

ClassImp(LQUseTree)


LQUseTree::LQUseTree():
UseTree()
{
}

void
LQUseTree::AddVariables(bool skim) {
  skimming = skim;
}

void
LQUseTree::FillMETTree() {

  //Add Variables to be computed

  PrepareHistograms();
  
  //===========================================================================================================
  //End Declaration, now preparation
  histoManager.PrepareHistos(name);
  histoManager.Prepare2DHistos(name);
  histoManager.PrepareProfiles(name);

  vector<vector<vector<float> > > Wghts;
  // if(ShapeWeight)
  //  Wghts = GetFitWeight();
  
  vector<int> NumberEntries(nt+1,0);

  //sélection 
  for(int i=0;i<nt+1;i++) {
    
    _isData = (i==nt);

    //Déclaration des variables

    vc.ResetVariables();
    histoManager.InitStatus();
    
    //Skimming;
    // TFile* oFile=NULL;
//     TTree* skimtree=NULL;

    //Global variables ====== (saving time)
    
    //====================================

    float Weight=1;
    
    cout<<" beginning tree  " <<name[i]<<"   "<<tChains[i]->GetEntries()<<endl;
    //Boucle et sélection, remplissage histos
   
    int ent = tChains[i]->GetEntries();

    if(skimming && ent==0/*(NoData && _isData)*/ ) {
      _storeTuple=true;
      FinalizeSkimming();
    }

    //protection against empty trees
    if(ent==0) continue;

    int runM=0;
    
    boost::progress_display show_progress( ent );
    for(int ie=-1;ie<ent;ie++) { //-1 for the initialization
      ++show_progress;

    
	 
      if(ie==0) { //first scan done, now intitialize stuff
	vc.FinalizeInit();
	vc.BuildTree( tChains[i], skimming );
	histoManager.StartFilling();
      }
      
      if(skimming) {
 	FinalizeSkimming();
    	InitSkimming(i, ie, ent);
      }

      tChains[i]->GetEntry(ie); 
      
      if(i!=nt ) {
	Weight = GetWeight(i, ie);
	if(ie==-1) Weight =0;
      }
      else {
	Weight = 1;
	if(ie==-1) Weight =0;
      }

      if(vc.getUI("run") > runM )
	runM = vc.getUI("run");

      if(i==nt)
      	{
      
      	  bool doubleCount=false;
	  std::pair<int,int> tmp(vc.getUI("run"),vc.getUI("event"));
	  EventIter = Events.find( tmp );
	  if(EventIter != Events.end() ) {
	    doubleCount=true;
	    //abort(); ?? FIXME -> no abort by default
	  }
      	  if(doubleCount || (EventFilter && vc.getUI("run")>EventNum ) )
      	    {  continue; }
	  string t1(""),t2("");
	  std::pair<string,string> tmp2( t1, t2 );
	  
	  Events[ tmp ] = tmp2;
	  EvtsInFile.push_back(vc.getUI("event") );
   
	  //}
      	}

      
      //Skimming procedure
      //=====================================================================================
      
      if(_isData) {
	bool filters=passJetMETFilters();
	if(!MakeCut<bool>(filters, true, "=", i , "filters", Weight)) continue;
      }
      else {
	MakeCut<bool>(true, true, "=", i , "filters", Weight);
      }

      bool hlt = passHLT();
      if(!MakeCut<bool>(hlt, true, "=", i , "HLT", Weight)) continue;
   
      int NTaus = vc.getSize("HPSTauEta");
      int NEls  = vc.getSize("ElectronPt");

      int nLooseEl=0;
      int nGTau=0;
      int nGEl=0;
      int vtxId=-1;

      vc.tau.SetPtEtaPhiM(0.0001,0,0,1.776);
      vc.el.SetPtEtaPhiM(0.0001,0,0,0.0005);

      for(int ie=0;ie<NEls;ie++) {
	
	if(SimpleCut<bool>(passLooseElectronID(ie), false, "=") ) 
	  nLooseEl++;
	
	if(SimpleCut<bool>(passElectronID(ie), false, "=") ) continue;

      	if(vc.el.Pt() < vc.getD("ElectronPt", ie) ) {
	  vc.el.SetPtEtaPhiM( vc.getD("ElectronPt",ie), vc.getD("ElectronEta",ie), vc.getD("ElectronPhi",ie), 0.0005 );
	  nGEl++;
	  vtxId = vc.getI("ElectronVtxIndex",ie);
	}
      }
      
      for(int it=0;it<NTaus;it++) {

	if(SimpleCut<bool>(passTauID(it), false, "=") ) continue;
	//if(SimpleCut<bool>(passTauVtx(it), false, "=") ) continue;
	
	if(SimpleCut<float>(dR( vc.getD("HPSTauEta") ,vc.el.Eta(), vc.getD("HPSTauPhi"), vc.el.Phi()), 0.5 ,">")) continue; 
	if(vc.tau.Pt() < vc.getD("HPSTauPt", it) ) {
	  vc.tau.SetPtEtaPhiM( vc.getD("HPSTauPt"), vc.getD("HPSTauEta"), vc.getD("HPSTauPhi"), 1.776 );
	  nGTau++;
	}
      }

      //cout<<nGTau<<"  "<<nGEl<<endl;

      if(!MakeCut<bool>( nGEl>0 , true, "=", i , "preselelectron", Weight)) continue;
      if(!MakeCut<bool>( nGTau>0 , true, "=", i , "preseltau", Weight)) continue;

      // end filling ==========================================================
      
      if(skimming && ie!=-1) {
	//cout<<skimTree<<endl;
	skimTree->Fill();
      }    

      if(ie!=-1)
	NumberEntries[i]++;
     
    }//End events
   
    //skimming for last dataset
    if(skimming && _isData ) {
      FinalizeSkimming();
    }

    //Pointer deletion =====
  //   if(skimming) {
//       cout<<" End Writing, check pointers "<<endl;
//       cout<<skimtree<<"     "<<oFile<<endl;
//       //delete skimtree;
//       //delete oFile;
//       cout<<" End deletion... "<<endl;
//     }
    //=====================

  } //End datasets
 
  //cout<<" End Filling , S="<<S<<" ;  B= "<<B<<endl;
  cout<<" detail "<<endl;
  for(int i=0;i<nt+1;i++) {
    cout<<" --->  "<<name[i]<<"  "<< NumberEntries[i]<<endl;
  }
  
  Histos = histoManager.GetHistos();
  Histos2D = histoManager.GetHistos2D();
  Profiles = histoManager.GetProfiles();
  
  for(size_t i=0;i<Histos.size();i++)
    Weighted.push_back(false);
  
  cout<<"Histos "<<Histos.size()<<"   "<<Histos2D.size()<<endl;
  
}



void
LQUseTree::PrepareDatasets() {
  
  colors.push_back(kBlack); //End Line
 
  for(int unsigned i=0;i<datasets.size();i++) {
    if(useXSect) {
      GetNProcEvent(datasets[i]);
    }
    FillAddWeight(datasets[i]);
    //Type.push_back(datasets[i]);
  }
  name.push_back("data");

  for(int i=0;i<nt;i++) {
    cout<<name[nt-1-i]<<endl; 
    Sorder.push_back(name[nt-1-i]);
  }
  
  cout<<" End dataset preparation : nt="<<nt<<endl;

}


//========================================================================================

//Old function for PU contamination correction

TH1F* 
LQUseTree::GetVtxCorrection(string obs,int nds, int bin) {

  //Protection
  if(NVert!=1) {
    cout<<" Warning, no requirement of only 1 Vtx event !!!!!"<<endl;
    abort();
  }
  
  float contamination = 0.069;
  float Systcontamination = 0.02;

  int Obs = histoManager.FindNVar(obs);
  int Ndata= histoManager.access(Obs,nds);
  TH1F* dataTMP = (TH1F*)Histos[ Ndata ]->Clone();

  float Int = dataTMP->Integral(0, dataTMP->GetNbinsX()+1);
  
  float Nevts = Int*contamination;
  Obs = histoManager.FindNVar(obs+"Control");
  Ndata= histoManager.access(Obs,nds);
  TH1F* dataTMP2 = (TH1F*)Histos[ Ndata ]->Clone();
 
  TH1F* Syst = (TH1F*)Histos[ Ndata ]->Clone();
  float NevtsSyst = Int*Systcontamination;

  Int = dataTMP2->Integral(0, dataTMP2->GetNbinsX()+1);

  dataTMP2->Scale( Nevts/Int );
  dataTMP2->Rebin(bin);
  dataTMP->Rebin(bin);

  Syst->Scale( NevtsSyst/Int );
  Syst->Rebin(bin);
 
  VtxCorSyst.clear();
  if(VtxCorSyst.size()==0)
    for(int i=1;i< Syst->GetNbinsX();i++)
      VtxCorSyst.push_back(Syst->GetBinContent(i) );
  
  cout<<" Contamination  "<<contamination<<" ;  "<<Nevts<<"  ;  "<<Int<<"   ;  "<<Nevts/Int<<"   "<<endl;
  if(Draw3on1)
    getRatio =true;
  if(!getRatio) {
    string nameO =  (histoManager.FindLeg( observable )).first;
    histoManager.RatioDistributionsNorm( (TH1F*)dataTMP->Clone(),(TH1F*)dataTMP2->Clone(), nameO, RangeX);
    c2->cd();
    getRatio=true;
  }
  return dataTMP2;

}


//Old function for on the fly reweighting
vector<vector<vector<float> > >  
LQUseTree::GetFitWeight() {

  vector<int> NumberEntries(nt+1,0);
  
  cout<<" ==================================================== "<<endl;
  cout<<" ================= En cours de weight =============== "<<endl;

  //Needed by the NN
  vector<double> inputVec(8,0);

  //sÃ©lection
  for(int i=0;i<nt+1;i++) {

    //Déclaration des variables
    
    int nVertex;
    tChains[i]->SetBranchAddress("nVertex", &nVertex);

    cout<<" beginning tree " <<name[i]<<"   "<<tChains[i]->GetEntries()<<endl;
    //Boucle et sélection, remplissage histos
    int ent = tChains[i]->GetEntries();
    
    float Weight;

    for(int ie=0;ie<ent;ie++) {
      
      // if(ie%1000==0)
      // 	fprintf(stdout, "Done : %d / %d\r", ie, ent );

      if(i!=nt)
	Weight = GetWeight(i, ie);
      else
	Weight = 1;
      
      tChains[i]->GetEntry(ie);      

      histoManager.fill("NVertexControl",i,nVertex, Weight );

      
    }//End events
  }// End datasets
 
    
  Histos = histoManager.GetHistos();
  Histos2D = histoManager.GetHistos2D();
    
  vector<vector<vector<float> > > Wghts(5,vector<vector<float> >(2,vector<float>(2,1)));
 
  return Wghts;
    
}


void
LQUseTree::AddTreeLoading( string type, vector<string> AddTNames ) {

  vector<TChain*> tmpch;
  
 //Preparation des TChains
   for(int i=0;i<nt+1;i++) {
     TChain* ctmp = new TChain(AddTNames[0].c_str());
     tmpch.push_back(ctmp);
   }
   
   cout<<" debut lecture "<<reposi<<endl;
   TFile* datafile;
   if(!NoData) {
     for(size_t i=0;i<data.size();i++) {
       if(data[i]!="") {
	 string NameF = "data/"+reposi+"/"+data[i]+".root"; 
	 if(reposi.find(":")!=(size_t)-1) NameF=reposi+"/"+datasets[i]+".root";

	 datafile = TFile::Open(NameF.c_str()); //new TFile(NameF.c_str(), "READ");
	 if(datafile==NULL) { cout<<" No such file "<<histoManager.Name<<endl; return;}
	 if(AddTNames.size()!=1) {
	   for(size_t it=0; it<AddTNames.size();it++) {
	     TTree* tmptree = (TTree*)datafile->Get( (AddTNames[it]).c_str() );
	     if(tmptree != NULL ) {
	       tmpch[nt]->Add( (NameF+"/"+AddTNames[it]).c_str()); }
	     delete tmptree;
	   }
	 }
	 else{
	   tmpch[nt]->Add(NameF.c_str());
	 }
	 
	 datafile->Close();
       }
     }
   }
   
   for(size_t i=0;i<datasets.size();i++) {
    string NameF = "data/"+reposi+"/"+datasets[i]+".root";
    if(reposi.find(":")!=(size_t)-1) NameF=reposi+"/"+datasets[i]+".root";	     

    datafile =  TFile::Open(NameF.c_str());//new TFile(NameF.c_str(), "READ");
    if(AddTNames.size()!=1) {
      for(size_t it=0; it<AddTNames.size();it++) {
	TTree* tmptree = (TTree*)datafile->Get( (AddTNames[it]).c_str() );
	if(tmptree != NULL ) {
	  tmpch[ dt[i] ]->Add( (NameF+"/"+AddTNames[it]).c_str()); }
	delete tmptree;
      }
    }
    else{
      tmpch[ dt[i] ]->Add(NameF.c_str());
    }
    datafile->Close();
  }

  AddChains[ type ] = tmpch;
  cout<<" Additionnal tree loaded : "<<type<<endl;

}





void 
LQUseTree::SetAdditionnalTrees(string s, vector<string> vs ) {
  AddTreeLoading( s, vs );
}


TTree*
LQUseTree::GetAdditionnalTree(string type, int ds ) {

  map<string, vector<TChain*> >::const_iterator iter;

  iter = AddChains.find(type);

  return (TTree*)(((*iter).second)[ds]);

}


//=================================================================================
// User functions




void 
LQUseTree::GetNProcEvent(string dataset) {

  string NameF = "data/"+reposi+"/"+dataset+".root";
  if(reposi.find(":")!=(size_t)-1) NameF=reposi+"/"+dataset+".root";	 
  TFile* file = TFile::Open( NameF.c_str() );
  TH1I* htmp = (TH1I*)file->Get("rootTupleTree/nEventProc");
  internalNEvt[ dataset ] = htmp->Integral(0,1001);
  
  delete htmp;
  file->Close();
}


void
LQUseTree::PrepareHistograms() {
  
  //
  histoManager.AddVariable("qT",400,0,400,"p_{T}(#mu#tau) [GeV]","pT");
  histoManager.AddVariable("mass",1600,0,400,"M_{#mu#tau} [GeV]","mass");
  histoManager.AddVariable("Tmass",1600,0,400,"M_{T}(#mu,#slash{E}_{T}) [GeV]","Tmass");

  histoManager.AddVariable("ptl1",400,0,200,"p_{T}(#mu_{1}) [GeV]","ptl1");
  histoManager.AddVariable("ptl2",400,0,200,"p_{T}(#mu_{2}) [GeV]","ptl2");
  histoManager.AddVariable("etal1",120,-3,3,"#eta(#mu_{1}) ","etal1");
  histoManager.AddVariable("etal2",120,-3,3,"#eta(#mu_{2}) ","etal2");
  histoManager.AddVariable("phil1",128,0,TMath::Pi()*2,"#phi(#mu_{1}) [rad]","phil1");
  histoManager.AddVariable("phil2",128,0,TMath::Pi()*2,"#phi(#mu_{2}) [rad]","phil2");

  histoManager.AddVariable("njet",21,-0.5,20.5,"Jet multiplicity","nJet");

  //Vertex
  histoManager.AddVariable("NVertex",50,0,50,"number of vertices","Vertex");
  
}



bool
LQUseTree::passLooseElectronID(size_t ie) {

  float eta = vc.getD("ElectronSCEta",ie);
  float pt = vc.getD("ElectronPt",ie);
  
  if(SimpleCut<float>(fabs(eta), 2.1, ">") ||
     SimpleCut<float>(fabs(eta), 1.4444,"[]",1.56) ) return false;
  if(SimpleCut<float>(pt, 25, "<")) return false;

  //Identification
  bool isEB= fabs(eta)<1.5;
  float hoe = vc.getD("ElectronHoE",ie);
  float sieie = vc.getD("ElectronSigmaIEtaIEta",ie);
  float deta = vc.getD("ElectronDeltaEtaTrkSC",ie);
  float dphi = vc.getD("ElectronDeltaPhiTrkSC",ie);
  
  if(SimpleCut<float>(hoe,isEB?0.12:0.10, ">")) return false;
  if(SimpleCut<float>(sieie,isEB?0.01:0.03, ">")) return false;
  if(SimpleCut<float>(fabs(deta),isEB?0.007:0.009, ">")) return false;
  if(SimpleCut<float>(fabs(dphi),isEB?0.15:0.10, ">")) return false;
 
 //COnversion
  float convP=0.;// vc.getD("ElectronConvFitProb",ie);
  int nMsHit = vc.getI("ElectronMissingHitsEG",ie);

  if(SimpleCut<float>( convP, 1e-6, ">")) return false;
  if(SimpleCut<int>( nMsHit, 1, ">")) return false;


  //Isolation
//   float rho = vc.getD("rhoJets");
//   float chIso = vc.getD("ElectronPFChargedHadronIso04");
//   float nhIso = vc.getD("ElectronPFNeutralHadronIso04");
//   float phIso = vc.getD("ElectronPFPhotonIso04");

//   float isocorr = chIso + max(nhIso+phIso - rho * Aeff(eta), (float)0.);

//   if(isocorr/pt > 0.15) return false;

  return true;


}

bool
LQUseTree::passElectronID(size_t ie) {

  float eta = vc.getD("ElectronSCEta",ie);
  float pt = vc.getD("ElectronPt",ie);
  
  if(SimpleCut<float>(fabs(eta), 2.1, ">") ||
     SimpleCut<float>(fabs(eta), 1.4444,"[]",1.56) ) return false;
 
  if(SimpleCut<float>(pt, 25, "<")) return false;
 
  //Identification
  bool isEB= fabs(eta)<1.5;
  float hoe = vc.getD("ElectronHoE",ie);
  float sieie = vc.getD("ElectronSigmaIEtaIEta",ie);
  float deta = vc.getD("ElectronDeltaEtaTrkSC",ie);
  float dphi = vc.getD("ElectronDeltaPhiTrkSC",ie);
  
  if(SimpleCut<float>(hoe,isEB?0.12:0.10, ">")) return false;
  if(SimpleCut<float>(sieie,isEB?0.01:0.03, ">")) return false;
  if(SimpleCut<float>(fabs(deta),isEB?0.004:0.007, ">")) return false;
  if(SimpleCut<float>(fabs(dphi),isEB?0.06:0.03, ">")) return false;
 
  //COnversion
  float convP=0.;// vc.getD("ElectronConvFitProb",ie);
  int nMsHit = vc.getI("ElectronMissingHitsEG",ie);

  if(SimpleCut<float>( convP, 1e-6, ">")) return false;
  if(SimpleCut<int>( nMsHit, 1, ">")) return false;


  //Isolation
//   float rho = vc.getD("rhoJets");
//   float chIso = vc.getD("ElectronPFChargedHadronIso04");
//   float nhIso = vc.getD("ElectronPFNeutralHadronIso04");
//   float phIso = vc.getD("ElectronPFPhotonIso04");

//   float isocorr = chIso + max(nhIso+phIso - rho * Aeff(eta), (float)0.);

//   if(isocorr/pt > 0.15) return false;

  return true;

}


float 
LQUseTree::Aeff( float eta) {

  if(_isData) {
    if(fabs(eta)<1.0) return 0.19;
    if(fabs(eta)>1.0 && fabs(eta)<1.479) return 0.25;
    if(fabs(eta)>1.479 && fabs(eta)<2.0) return 0.12;
    if(fabs(eta)>2.0 && fabs(eta)<2.2) 	return 0.21;
    if(fabs(eta)>2.2 && fabs(eta)<2.3) 	return 0.27;
    if(fabs(eta)>2.3 && fabs(eta)<2.4) 	return 0.44;
    if(fabs(eta)>2.4) 	return 0.52; 
  }
  else { //To be fixed //FIXME
    if(fabs(eta)<1.0) return 0.19;
    if(fabs(eta)>1.0 && fabs(eta)<1.479) return 0.25;
    if(fabs(eta)>1.479 && fabs(eta)<2.0) return 0.12;
    if(fabs(eta)>2.0 && fabs(eta)<2.2) 	return 0.21;
    if(fabs(eta)>2.2 && fabs(eta)<2.3) 	return 0.27;
    if(fabs(eta)>2.3 && fabs(eta)<2.4) 	return 0.44;
    if(fabs(eta)>2.4) 	return 0.52; 

  }

  return 0.;

}

bool
LQUseTree::passTauID(size_t it) {

  if(SimpleCut<double>(vc.getD("HPSTauEta", it), 2.1, ">")) return false;
  if(SimpleCut<double>(vc.getD("HPSTauPt", it), 25, "<")) return false;
  
  // cout<<" héhé "<<endl;

  if(SimpleCut<double>(vc.getD("HPSTauAgainstElectronMVADiscr", it),0.,"=")) return false;
  if(SimpleCut<double>(vc.getD("HPSTauAgainstMuonLooseDiscr", it),0.,"=")) return false;
  if(SimpleCut<double>(vc.getD("HPSTauLooseCombinedIsolationDeltaBetaCorrDiscr", it),0.,"=")) return false;
  

  return true;
  
}

bool 
LQUseTree::passHLT() {


  bool pass=false;
  int run = vc.getUI("run");
  if(!_isData) run=-1;

  for(size_t ih=0;ih<vc.getSize("HLTInsideDatasetTriggerNames");ih++) {
    
    if( run >=190456 && run <=193661) {
      pass = triggerDecision(ih, "HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20");//||
	    //  triggerDecision(ih, "HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20L1Jet") ||
// 	     triggerDecision(ih, "HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau22L1Jet") ||
// 	     triggerDecision(ih, "HLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20") ||
    // triggerDecision(ih, "HLT_Ele22_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20") ;
  
    }
    else if( run>=193834 && run<=203002 ) {
      pass = //triggerDecision(ih, "HLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20") ||
	triggerDecision(ih, "HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20");
      
    }
    else{
      pass = triggerDecision(ih, "HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20");
    }
    //cout<<vc.getS("HLTInsideDatasetTriggerNames",ih)<<"   "<<pass<<endl;

    if(pass) return true;

  }

  return false;
}

bool
LQUseTree::triggerDecision(size_t ihlt, string ref) {

  if(SimpleCut<size_t>(vc.getS("HLTInsideDatasetTriggerNames",ihlt).find(ref),-1,"=") ) return false;
  //cout<<" line found "<<vc.getS("HLTInsideDatasetTriggerNames",ihlt)<<endl;
  if(SimpleCut<int>(vc.getI("HLTInsideDatasetTriggerPrescales", ihlt),1,"!=")) return false;
  //cout<<"\t ---------> "<<vc.getI("HLTInsideDatasetTriggerPrescales", ihlt)<<endl;
  if(SimpleCut<bool>(vc.getB("HLTInsideDatasetTriggerDecisions",ihlt),true,"!=")) return false;
  //cout<<"\t ---------> "<<vc.getB("HLTInsideDatasetTriggerDecisions",ihlt)<<endl;
  return true;
}

bool
LQUseTree::passJetMETFilters() {

  bool hbhe = vc.getB("passHBHENoiseFilter");
  bool trkFail = !vc.getB("isTrackingFailure");
  bool bmScrap = !vc.getB("isBeamScraping");
  bool isPV = vc.getB("isPrimaryVertex");
  bool bmHalo = vc.getB("passBeamHaloFilterLoose");
  bool hcalLas = !vc.getB("passHcalLaserEventFilter");
  bool ecalLas = vc.getB("passEcalLaserCorrFilter");
  bool ecalDC = !vc.getB("passEcalDeadCellBoundaryEnergyFilter");
  bool ecalTP = vc.getB("passEcalDeadCellTriggerPrimitiveFilter");
  bool ecalMk = vc.getB("passEcalMaskedCellDRFilter");
  bool ecalSC = !vc.getB("passBadEESupercrystalFilter");
  bool caloBnd = vc.getB("passCaloBoundaryDRFilter");

 //   cout<<hbhe <<"  "<< trkFail <<"  "<< bmScrap <<"  "<< isPV <<"  "<< bmHalo <<"  "<<
//      hcalLas <<"  "<< ecalLas <<"  "<< ecalDC /*<<"  "<< ecalMk*/ <<"  "<< 
//      ecalSC <<"  "<< caloBnd<<endl;

  return (hbhe && trkFail && bmScrap && isPV && bmHalo &&
	  hcalLas && ecalLas && ecalDC /*&& ecalMk*/ && 
	  ecalSC && caloBnd);

}
