#include "LeptoQuarkUseTree.hh"

#include <iomanip>
#include <TLorentzVector.h>
// #include <boost/timer.h>
// #include <boost/progress.hpp>

#include <TH1I.h>
#include <TObjArray.h>

//#include "/home/mmarionn/Documents/CMS/LeptoQuark/KKCode/ForReweighting/GetPUWeight.C"

using namespace std;

ClassImp(LeptoQuarkUseTree)


LeptoQuarkUseTree::LeptoQuarkUseTree():
UseTree()
{
  
  //skimming
  skimming=false; 

}

void
LeptoQuarkUseTree::AddVariables( bool skim ) {
  
  skimming=skim;

}

void
LeptoQuarkUseTree::FillLQTree() {

  //Add Variables to be computed

  PrepareHistograms();

  // string flagId="";
  // if(UsePdgId && pdgId==13) flagId="Mu";
  // if(UsePdgId && pdgId==11) flagId="Elec";

  //===========================================================================================================
  //End Declaration, now preparation
  histoManager.PrepareHistos(name);
  histoManager.Prepare2DHistos(name);
 

  vector<vector<vector<float> > > Wghts;
  // if(ShapeWeight)
  //  Wghts = GetFitWeight();
  
  string Epart="EE";
  if(EcalP==0)
    Epart = "EB";
  if(EcalP==2)
    Epart="Combined";
  if(EcalP==3)
    Epart="EBEE for Z studies";

  cout<<" Ecal Part "<<Epart<<endl;

  cout<<" check cuts EB : id "<<_IDcuts[0][0]<<"  "
      <<_IDcuts[0][1]<<"  "<<_IDcuts[0][2]<<"  "<<_IDcuts[0][3]<<endl;
  cout<<" check cuts EB : iso "<<_Isocuts[0][0]<<"  "
      <<_Isocuts[0][1]<<"  "<<_Isocuts[0][2]<<endl;

  cout<<" check cuts EE : id "<<_IDcuts[1][0]<<"  "
      <<_IDcuts[1][1]<<"  "<<_IDcuts[1][2]<<"  "<<_IDcuts[1][3]<<endl;
  cout<<" check cuts EE : iso "<<_Isocuts[1][0]<<"  "
      <<_Isocuts[1][1]<<"  "<<_Isocuts[1][2]<<endl;

  int S=0,B=0;
 
  vector<int> NumberEntries(nt+1,0);


  //sélection 
  for(int i=0;i<nt+1;i++) {
    
    //Déclaration des variables
    VarClassLQ vc;
    histoManager.InitStatus();
    
    //Skimming;
    TFile* oFile;
    TTree* skimtree;

    //Global variables ====== (saving time)
   
    //====================================

    vector<string> bTagEvent;
    


    float Weight=1;
    
    //    float Weight2=1;

    cout<<" beginning tree  " <<name[i]<<"   "<<tChains[i]->GetEntries()<<endl;
    //Boucle et sélection, remplissage histos
    //int EP=0;
    int ent = tChains[i]->GetEntries();

    //protection against empty trees
    if(ent==0) continue;

    int runM=0;

    //    boost::progress_display show_progress( ent );
    for(int ie=-1;ie<ent;ie++) { //-1 for the initialization
      // ++show_progress;
	 
      if(ie==0) {
	//	ie += 7080;
	vc.FinalizeInit();
	vc.BuildTree( tChains[i], 0 );
	if(skimming) {
	 oFile = new TFile( ("Skimming/"+name[i]+"_skim.root").c_str(),"RECREATE");
	  tChains[i]->LoadTree(0);
	  skimtree = tChains[i]->CloneTree(0);
	}
	histoManager.StartFilling();
      }
      
      if(ie%1000==0)
	fprintf(stdout, "Done : %d / %d\r", ie, ent );
      

      tChains[i]->GetEntry(ie); 
      
    
      if(i!=nt  ) {
	if(name[i].substr(0,2)!="LQ")
	  Weight = Lumi*vc.getD("weight");//*pileupWeight_Fall11_InTime( vc.getI("nPU",1) ); //cout<<Weight<<endl;
	else
	  {
	    Weight = Lumi*0.5/36000 ;//* pileupWeight_Fall11_InTime( vc.getI("nPU",1) );
	  }
	//cout<<"  "<<vc.getD("weight")<<"   "<<vc.getI("nPU",1)<<"   "<<pileupWeight_Fall11_InTime( vc.getI("nPU",1) )<<endl;
	// if(ie<10)
	//   cout<<ie<<"   "<<Weight<<endl;
	if(ie==-1) Weight =0;
      }
      else
	Weight = 1;
      
      // isEl=false;
      // isMu=false;
      
      // cout<<Weight<<endl;
      //   bTagEvent.clear();
    
      // if((ie+1)%1000==0)
      //cout<<" event "<<ie+1<<endl;
      
      if(vc.getUI("run") > (size_t)runM )
	runM = vc.getUI("run");

 

      if(i==nt)
      	{
      
      	  bool doubleCount=false;
	  std::pair<int,int> tmp(vc.getUI("run"),vc.getUI("event"));
	  EventIter = Events.find( tmp );
	  //	  cout<<ie<<"   "<<vc.run<<"   "<<vc.event<<endl;
	  if(EventIter != Events.end() ) {
	    doubleCount=true;//cout<<ie<<" ============>>> found  double event  "<<vc.getUI("run")<<"   "<<vc.getUI("event")<<"   "
	    //	 <<vc.getD("PFTauPt")<<" first tau pt? "<<endl;
	    abort();
	  }
      	  if(doubleCount || (EventFilter && vc.getUI("run")>(size_t)EventNum ) )
      	    {  continue; }
	  string t1(""),t2("");
	  std::pair<string,string> tmp2( t1, t2 );
	  
	  Events[ tmp ] = tmp2;
	  EvtsInFile.push_back(vc.getUI("event") );
   
	  //}
      	}

 
      //All events
      //=====================================================================================
      
      
      MakeCut<bool>(true,true,"=", i, "Beginning", Weight);
  
    
      //HLT
      bool passHLT = HLTRequirement( &vc, i );
      //cout<<" pass HLT "<<passHLT<<"   "<<ie<<endl;
      if(!MakeCut<bool>(passHLT, true, "=",i,"HLT", Weight )) continue;
  
      //Electron
      bool hasEle = ElectronID( &vc );

      if(!MakeCut<bool>(hasEle, true, "=",i,"electron requirement", Weight )) continue;
    
      //Vertex
      bool passVtx = VtxID( &vc );
      
      if(!MakeCut<bool>(passVtx, true, "=",i,"vtx requirement", Weight )) continue;
      //  cout<<" event "<<ie+1<<endl;
      //tau
      bool hasTau = TauID( &vc );

      if(!MakeCut<bool>(hasTau, true, "=",i,"tau requirement", Weight )) continue;
    
       
      //Jets
      vc.noBTag =false; //Def, event is btag
      bool hasJets = JetID( &vc );

       if(!MakeCut<bool>(hasJets, true, "=",i,"jet requirement", Weight )) continue;
       
   

       //Var filling 
	vc.dijet = vc.jet1 + vc.jet2;
      	vc.tjet1 = vc.tau + vc.jet1;
      	vc.ljet1 = vc.lep + vc.jet1;
      	vc.tjet2 = vc.tau + vc.jet2;
      	vc.ljet2 = vc.lep + vc.jet2;

      	vc.tljet1 = vc.tau + vc.lep + vc.jet1;
      	vc.tljet2 = vc.tau + vc.lep + vc.jet2;
      	vc.tjj = vc.tau + vc.jet1 + vc.jet2;
      	vc.ljj = vc.lep + vc.jet1 + vc.jet2;
      	vc.tljj = vc.tau + vc.lep + vc.jet1 + vc.jet2;


	FillPlots( (&vc), "Presel", i, Weight );
	
	//skimming
	if(skimming) {
	  skimtree->Fill();
	}

	//Advanced selection ===============================

	bool advKine = vc.lep.Pt() > ptLep && vc.tau.Pt() > ptTau;
	
	if( !MakeCut<bool>(advKine, true, "=", i,"Adv Kine Sel",Weight) ) continue;
	
	FillPlots( (&vc), "AdvKine", i, Weight );
	
	//Charge selection
	int lepChr = vc.getI("ElectronCharge",vc.lepL);
	int tauChr = vc.getI("PFTauCharge",vc.tauL);
	bool charge = chrReq!=0?(lepChr*tauChr==chrReq):true;
	//	cout<<ie<<"  "<<chrReq<<"   "<<lepChr<<"   "<<ch<<endl;//<<tauChr<<"   "<<endl;//<<charge<<endl;
	//	cout<<ie<<"   "<<chrReq<<"   "<<vc.getI("ElectronCharge",vc.lepL)*vc.getI("PFTauCharge",vc.tauL)<<"   "<<charge<<"        "<<endl;
	if( !MakeCut<bool>(charge, true, "=", i,"Charge",Weight) ) continue;
	//	cout<<" pass "<<endl;
	FillPlots( (&vc), "Charge", i, Weight );	

	//Advanced Vtx selection

	bool advVtxMatchZ = abs(vc.getD("ElectronVtxDistZ",vc.lepL)) <0.2;
	bool advVtxMatchR = sqrt( pow(vc.getD("ElectronVtxDist3D",vc.lepL),2) - pow(vc.getD("ElectronVtxDistZ",vc.lepL), 2)) <0.045;
	
	//	cout<<vc.getD("ElectronVtxDistZ",vc.lepL)<<"   "<<sqrt( pow(vc.getD("ElectronVtxDist3D",vc.lepL),2) - pow(vc.getD("ElectronVtxDistZ",vc.lepL), 2))<<endl;

	if( !MakeCut<bool>( (advVtxMatchZ && advVtxMatchR), true, "=", i,"Adv Vtx Match",Weight) ) continue;

	FillPlots( (&vc), "VtxMatch", i, Weight );
	
	
	//btagging
	if (!MakeCut<bool>( vc.noBTag, false, "=",i, "b tagging", Weight ) ) continue;

	FillPlots( (&vc), "BTag", i, Weight );

	//Nominal isolation
	
	if( !MakeCut<int>( vc.getI("PFTauLooseCombIsolation", vc.tauL), 1, "=" , i, "Nom Iso", Weight ) ) continue;
	
	FillPlots( (&vc), "TauIso", i, Weight );

	//===================================================================== end selection
	    
	    
      if(name[i]=="ZZ #rightarrow 2l2#nu")
	{ S++; 
	}
      else if(name[i]!="data")
	{B++; }
     
     
      //  cout<<" end event "<<endl;
      NumberEntries[i]++;
   
  }//End events
  //cout<<" Run Max "<<runM<<endl;
  //Pointer deletion =====
  
      //skimming
      if(skimming) {
	skimtree->Write();
	oFile->Close();
      }

    //=====================


  } //End datasets
 
  cout<<" End Filling , S="<<S<<" ;  B= "<<B<<endl;
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
LeptoQuarkUseTree::PrepareDatasets() {
  
  colors.push_back(kBlack); //End Line
 
  for(int unsigned i=0;i<datasets.size();i++) {
    if(useXSect) {
      GetNProcEvent(datasets[i]);
    }
    FillAddWeight(datasets[i]);
    Type.push_back(datasets[i]);
  }
  name.push_back("data");

  for(int i=0;i<nt;i++) {
    cout<<name[nt-1-i]<<endl; 
    Sorder.push_back(name[nt-1-i]);
  }
  
  cout<<" End dataset preparation : nt="<<nt<<endl;

}


TH1F* 
LeptoQuarkUseTree::GetVtxCorrection(string obs,int nds, int bin) {

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



vector<vector<vector<float> > >  
LeptoQuarkUseTree::GetFitWeight() {

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
      
      if(ie%1000==0)
	fprintf(stdout, "Done : %d / %d\r", ie, ent );

      if(i!=nt)
	Weight = GetWeight(i, ie);// special for LQ ...
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
 
float LeptoQuarkUseTree::Rapidity(float p, float mass, float pt, float eta) {

  //compute pz
  float pz = sqrt(p*p-pt*pt);

  //compute E
  float E = sqrt(p*p + mass*mass);

  //and now rapidity
  float y = 0.5 * log( (E+pz) / (E-pz) );

  return y;

}


int
LeptoQuarkUseTree::MuonClassification( double pt, double eta, int isID) {
  int muID = 0;
  
  // if( !( isGlob==1 && isTrk==1 && fabs(dxy)<0.2 && mStations>=2 && mPixHit>=1 && mTrkHit>=11) )
  //   { return 0; }
  
  if( fabs(eta)<2.4) {
    muID += 1;}

  if( pt > 20 ) {
    muID += 2;
  }
  
  if( isID==1 )
    muID += 4;

  return muID;
}

int
LeptoQuarkUseTree::ElectronClassification( double pt, double eta, double sieie, double deta, double dphi, double hoe ) {

  int elID=0;

  if(pt>20) elID+=2;

  if( fabs(eta) < 1.45 || 
      ( fabs(eta)>1.54 && fabs(eta)<2.5 ) )
    elID+=1;

    if( fabs(eta) <1.5) {

    if( sieie >= 0.01 ) return elID;
    if( fabs(deta) >= 0.005 ) return elID;
    if( fabs(dphi) >= 0.027 ) return elID;
    //   if( hoe >= 0.025 ) return elID;

  }
  else {
    
    if( sieie >= 0.031 ) return elID;
    if( fabs(deta) >= 0.006 ) return elID;
    if( fabs(dphi) >= 0.021 ) return elID;
    //    if( hoe >= 0.025 ) return elID;
    
  }
    
    elID +=4;
    return elID;
}

int
LeptoQuarkUseTree::TauClassification(double pt, double eta) {
  
  int tauID = 0;
  
  if( fabs(eta) < 5 ) {
    tauID += 1;
  }

  if( pt > 20 )
    tauID += 2;

  return tauID;
}



int
LeptoQuarkUseTree::JetClassification( double pt, double eta, double phi,
				      double l1Eta, double l1Phi,
				      double l2Eta, double l2Phi ) {

  int jetId=0;
  if( pt > 30 ) {
    jetId += 1; }
  if( fabs( eta ) < 4.5 ) {
    jetId += 2;
  }

  float dR1 = UseTree::dR( eta, l1Eta, phi, l1Phi );
  float dR2 = UseTree::dR( eta, l2Eta, phi, l2Phi );
  //cout<<eta<<"   "<<phi<<"   "<<l1Eta<<"    "<<l1Phi<<"    "<<fabs(eta-l1Eta)<<"   "<<dPhi(phi,l1Phi)<<"    "<<dR1<<endl;
  if(dR1 > 0.5 && dR2 > 0.5) {
    jetId += 4; }

  return jetId;

}



void
LeptoQuarkUseTree::AddTreeLoading( string type, vector<string> AddTNames ) {

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
	 datafile = new TFile(NameF.c_str(), "READ");
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
   
   //   int nTmp=0; 
   //int ndt=0;
  for(size_t i=0;i<datasets.size();i++) {
    string NameF = "data/"+reposi+"/"+datasets[i]+".root";

    // if(ndt!=dt[i])
    //   nTmp=0;
    //ndt=dt[i];
    
    datafile = new TFile(NameF.c_str(), "READ");
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
LeptoQuarkUseTree::SetAdditionnalTrees(string s, vector<string> vs ) {
  AddTreeLoading( s, vs );
}


TTree*
LeptoQuarkUseTree::GetAdditionnalTree(string type, int ds ) {

  map<string, vector<TChain*> >::const_iterator iter;

  iter = AddChains.find(type);

  return (TTree*)(((*iter).second)[ds]);

}



bool 
LeptoQuarkUseTree::isMatch(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) {

  
  if( dR(eta1,eta2,phi1,phi2) < 0.2 && 
      ( (fabs(pt1)-fabs(pt2))/fabs(pt2) >0 ) && (fabs(pt1)-fabs(pt2))/fabs(pt2) <0.6 )
    return true;
  else
    return false;

}


void 
LeptoQuarkUseTree::GetNProcEvent(string dataset) {

  string NameF = "data/"+reposi+"/"+dataset+".root";
  TFile* file = TFile::Open( NameF.c_str() );
  
  TH1I* htmp = (TH1I*)file->Get("EventCounter");
  
  internalNEvt[ dataset ] = htmp->GetBinContent(1);
  
  delete htmp;
  file->Close();
}


void
LeptoQuarkUseTree::FillPlots( VarClassLQ* vc, string anlvl, int i, float Weight ) {

  anlvl = "_"+anlvl;

  vc->dilep = vc->tau + vc->lep;

  float taujet1dR = dR(vc->tau.Eta(), vc->jet1.Eta(), vc->tau.Phi(), vc->jet1.Phi() ); 
  float taujet2dR = dR(vc->tau.Eta(), vc->jet2.Eta(), vc->tau.Phi(), vc->jet2.Phi() ); 
  float lepjet1dR = dR(vc->lep.Eta(), vc->jet1.Eta(), vc->lep.Phi(), vc->jet1.Phi() ); 
  float lepjet2dR = dR(vc->lep.Eta(), vc->jet2.Eta(), vc->lep.Phi(), vc->jet2.Phi() ); 

  float taujet1dPhi = dPhi(vc->tau.Phi(), vc->jet1.Phi() )*180/TMath::Pi() ;
  float taujet2dPhi = dPhi(vc->tau.Phi(), vc->jet2.Phi() )*180/TMath::Pi() ;
  float lepjet1dPhi = dPhi(vc->lep.Phi(), vc->jet1.Phi() )*180/TMath::Pi() ;
  float lepjet2dPhi = dPhi(vc->lep.Phi(), vc->jet2.Phi() )*180/TMath::Pi() ;

  float taujet1dPhilepjet2 = dPhi( (vc->tjet1).Phi(), (vc->ljet2).Phi() )*180/TMath::Pi();
  float taujet2dPhilepjet1 = dPhi( (vc->tjet2).Phi(), (vc->ljet1).Phi() )*180/TMath::Pi();

  float taujet1dRlepjet2 = dR( (vc->tjet1).Eta(), (vc->ljet2).Eta(), (vc->tjet1).Phi(), (vc->ljet2).Phi() );
  float taujet2dRlepjet1 = dR( (vc->tjet2).Eta(), (vc->ljet1).Eta(), (vc->tjet2).Phi(), (vc->ljet1).Phi() );

  float taulepDPhi = UseTree::dPhi(vc->tau.Phi(), vc->lep.Phi() )*180/TMath::Pi() ;
  float taulepDEta =  vc->tau.Eta() - vc->lep.Eta();

  float jjDPhi = UseTree::dPhi(vc->jet1.Phi(), vc->jet2.Phi() )*180/TMath::Pi() ;
  float jjDEta =  vc->jet1.Eta() - vc->jet2.Eta();

  float taulepjetjetDPhi = UseTree::dPhi( vc->dilep.Phi(), vc->dijet.Phi() )*180/TMath::Pi();
  float taulepjetjetDEta = ( vc->dilep.Rapidity() - vc->dijet.Rapidity() );

  histoManager.fill("PtTau"+anlvl,i, vc->tau.Pt() ,Weight);
  histoManager.fill("PtLep"+anlvl,i, vc->lep.Pt(),Weight);
  histoManager.fill("EtaTau"+anlvl,i, vc->tau.Eta(),Weight);
  histoManager.fill("EtaLep"+anlvl,i, vc->lep.Eta(),Weight);
  histoManager.fill("PhiTau"+anlvl,i, vc->tau.Phi()*180/TMath::Pi(),Weight);
  histoManager.fill("PhiLep"+anlvl,i, vc->lep.Phi()*180/TMath::Pi(),Weight);
  


  histoManager.fill("TLMass"+anlvl,i,vc->dilep.M(),Weight);
  histoManager.fill("TLPt"+anlvl,i,vc->dilep.Pt(),Weight);
  histoManager.fill("TLY"+anlvl,i,vc->dilep.Rapidity(),Weight);
  histoManager.fill("TLPhi"+anlvl,i,vc->dilep.Phi()*180/TMath::Pi(),Weight);

  histoManager.fill("TLDPhi"+anlvl,i,taulepDPhi,Weight);
  histoManager.fill("TLDEta"+anlvl,i,taulepDEta,Weight);


  // histoManager.fill("TrkIsoT"+anlvl,i, (*MuonTrkIso)[muLead] ,Weight);
  // histoManager.fill("TrkIsoM2"+anlvl,i, (*MuonTrkIso)[muTrail] ,Weight);
  // histoManager.fill("EcalIsoM1"+anlvl,i, (*MuonEcalIso)[muLead] ,Weight);
  // histoManager.fill("EcalIsoM2"+anlvl,i, (*MuonEcalIso)[muTrail] ,Weight);
  // histoManager.fill("HcalIsoM1"+anlvl,i, (*MuonHcalIso)[muLead] ,Weight);
  // histoManager.fill("HcalIsoM2"+anlvl,i, (*MuonHcalIso)[muTrail] ,Weight);
  // histoManager.fill("IsoM1"+anlvl,i, ( (*MuonTrkIso)[muLead] + (*MuonEcalIso)[muLead] + (*MuonHcalIso)[muLead]) / (*MuonPt)[muLead],Weight);
  // histoManager.fill("IsoM2"+anlvl,i, ( (*MuonTrkIso)[muTrail] + (*MuonEcalIso)[muTrail] + (*MuonHcalIso)[muTrail]) / (*MuonPt)[muTrail],Weight);
 
  histoManager.fill("PtJ1"+anlvl,i, vc->jet1.Pt(),Weight);
  histoManager.fill("PtJ2"+anlvl,i, vc->jet2.Pt(),Weight);
  histoManager.fill("EtaJ1"+anlvl,i, vc->jet1.Eta(),Weight);
  histoManager.fill("EtaJ2"+anlvl,i, vc->jet2.Eta(),Weight);
  histoManager.fill("PhiJ1"+anlvl,i, vc->jet1.Phi()*180/TMath::Pi(),Weight);
  histoManager.fill("PhiJ2"+anlvl,i, vc->jet2.Phi()*180/TMath::Pi(),Weight);

  histoManager.fill("JJMass"+anlvl,i,vc->dijet.M(),Weight);
  histoManager.fill("JJPt"+anlvl,i,vc->dijet.Pt(),Weight);
  histoManager.fill("JJY"+anlvl,i,vc->dijet.Rapidity(),Weight);
  histoManager.fill("JJPhi"+anlvl,i,vc->dijet.Phi()*180/TMath::Pi(),Weight);

  histoManager.fill("JJDPhi"+anlvl,i,jjDPhi,Weight);
  histoManager.fill("JJDEta"+anlvl,i,jjDEta,Weight);
 
   histoManager.fill("J1BTag"+anlvl,i,vc->getD("PFJetTrackCountingHighEffBTag",vc->jetLead), Weight);
  histoManager.fill("J2BTag"+anlvl,i,vc->getD("PFJetTrackCountingHighEffBTag",vc->jetTrail), Weight);

  histoManager.fill("TLJJDPhi"+anlvl,i,taulepjetjetDPhi,Weight);
  histoManager.fill("TLJJDY"+anlvl,i,taulepjetjetDEta,Weight);
  histoManager.fill("TLJJMass"+anlvl,i,vc->tljj.M(),Weight);
  histoManager.fill("TLJJPt"+anlvl,i,vc->tljj.Pt(),Weight);
  histoManager.fill("TLJJPhi"+anlvl,i,vc->tljj.Phi()*180/TMath::Pi(),Weight);
  histoManager.fill("TLJJY"+anlvl,i,vc->tljj.Rapidity(),Weight);

  histoManager.fill("ST"+anlvl,i,vc->jet1.Pt() + vc->jet2.Pt() + vc->tau.Pt() + vc->lep.Pt(),Weight);

  histoManager.fill("TauJJMass"+anlvl,i,vc->tjj.M(),Weight);
  histoManager.fill("LepJJMass"+anlvl,i,vc->ljj.M(),Weight);
  histoManager.fill("TLJ1Mass"+anlvl,i,vc->tljet1.M(),Weight);
  histoManager.fill("TLJ2Mass"+anlvl,i,vc->tljet2.M(),Weight);
      
  histoManager.fill("TauJ1Mass"+anlvl,i,vc->tjet1.M(),Weight);
  histoManager.fill("LepJ1Mass"+anlvl,i,vc->ljet1.M(),Weight);
  histoManager.fill("TauJ2Mass"+anlvl,i,vc->tjet2.M(),Weight);
  histoManager.fill("LepJ2Mass"+anlvl,i,vc->ljet2.M(),Weight);

  histoManager.fill("TauJ1DR"+anlvl,i,taujet1dR,Weight);
  histoManager.fill("TauJ2DR"+anlvl,i,taujet2dR,Weight);
  histoManager.fill("LepJ1DR"+anlvl,i,lepjet1dR,Weight);
  histoManager.fill("LepJ2DR"+anlvl,i,lepjet2dR,Weight);
    
  histoManager.fill("TauJ1DPhi"+anlvl,i,taujet1dPhi,Weight);
  histoManager.fill("TauJ2DPhi"+anlvl,i,taujet2dPhi,Weight);
  histoManager.fill("LepJ1DPhi"+anlvl,i,lepjet1dPhi,Weight);
  histoManager.fill("LepJ2DPhi"+anlvl,i,lepjet2dPhi,Weight);

  histoManager.fill("TauJ1DPhiLepJ2"+anlvl,i,taujet1dPhilepjet2,Weight);
  histoManager.fill("TauJ2DPhiLepJ1"+anlvl,i,taujet2dPhilepjet1,Weight);
    
  histoManager.fill("TauJ1DRLepJ2"+anlvl,i,taujet1dRlepjet2,Weight);
  histoManager.fill("TauJ2DRLepJ1"+anlvl,i,taujet2dRlepjet1,Weight);

  float metLTDPhi = UseTree::dPhi(vc->getD("PFMETPhi"), vc->dilep.Phi() )*180/TMath::Pi() ;
  float metJJDPhi = UseTree::dPhi(vc->getD("PFMETPhi"), vc->dijet.Phi() )*180/TMath::Pi() ;
  float mettauDPhi = UseTree::dPhi(vc->getD("PFMETPhi"), vc->tau.Phi() )*180/TMath::Pi() ;
  float metlepDPhi = UseTree::dPhi(vc->getD("PFMETPhi"), vc->lep.Phi() )*180/TMath::Pi() ;
  float metjet1DPhi = UseTree::dPhi(vc->getD("PFMETPhi"), vc->jet1.Phi() )*180/TMath::Pi() ;
  float metjet2DPhi = UseTree::dPhi(vc->getD("PFMETPhi"), vc->jet2.Phi() )*180/TMath::Pi() ;
  
  histoManager.fill("MET"+anlvl,i, vc->getD("PFMET"), Weight);
  histoManager.fill("METTLDPhi"+anlvl,i,metLTDPhi, Weight);
  histoManager.fill("METJJDPhi"+anlvl,i,metJJDPhi, Weight);
  histoManager.fill("METTauDPhi"+anlvl,i,mettauDPhi, Weight);
  histoManager.fill("METLepDPhi"+anlvl,i,metlepDPhi, Weight);
  histoManager.fill("METJ1DPhi"+anlvl,i,metjet1DPhi, Weight);
  histoManager.fill("METJ2DPhi"+anlvl,i,metjet2DPhi, Weight);




}

void
LeptoQuarkUseTree::PrepareHistograms() {

  vector<string> flavlvls(1,"");
  flavlvls[0]="";
  // flavlvls[1]="Mu";
  // flavlvls[2]="Elec";

  vector<string> anlvls;

  for(size_t is=0;is<flavlvls.size();is++) {

  anlvls.push_back("Presel"+flavlvls[is]);
  anlvls.push_back("AdvKine"+flavlvls[is]);
  anlvls.push_back("Charge"+flavlvls[is]);
  anlvls.push_back("VtxMatch"+flavlvls[is]);
  anlvls.push_back("BTag"+flavlvls[is]);
  anlvls.push_back("TauIso"+flavlvls[is]);
  }
  
  for(unsigned int ia=0;ia<anlvls.size();ia++) {

    string anlvl = "_"+anlvls[ia];


    histoManager.AddVariable("PtTau"+anlvl,200,0,400,"p_{T} #tau_{had} [GeV]","PtTau"+anlvl);
    histoManager.AddVariable("PtLep"+anlvl,200,0,400,"p_{T} l_{#tau} [GeV]","PtLep"+anlvl);
    histoManager.AddVariable("EtaTau"+anlvl,100,-3,3,"#eta #tau_{had} ","EtaTau"+anlvl);
    histoManager.AddVariable("EtaLep"+anlvl,100,-3,3,"#eta l_{#tau} ","EtaLep"+anlvl);
    histoManager.AddVariable("PhiTau"+anlvl,128,-180,180,"#phi #tau_{had} [#circ]","PhiTau"+anlvl);
    histoManager.AddVariable("PhiLep"+anlvl,128,-180,180,"#phi l_{#tau} [#circ]","PhiLep"+anlvl);

    histoManager.AddVariable("TLMass"+anlvl,500,0,500,"M_{#tau_{had} l_{#tau}} [GeV]","TLMass"+anlvl);
    histoManager.AddVariable("TLPt"+anlvl,400,0,400,"q_{T} [GeV]","TLPt"+anlvl);
    histoManager.AddVariable("TLY"+anlvl,100,-5,5,"Y_{#tau_{had}l_{#tau}}","TLY"+anlvl);
    histoManager.AddVariable("TLPhi"+anlvl,100,-3,3,"#phi_{#tau_{had}l_{#tau}}","TLPhi"+anlvl);

    histoManager.AddVariable("TLDPhi"+anlvl,128,-180,180,"#Delta#phi (#tau_{had}, l_{#tau}) [#circ]","TLDPhi"+anlvl);
    histoManager.AddVariable("TLDEta"+anlvl,100,-5,5,"#Delta#eta ({#tau}_{had}, l_{#tau}) ","TLDEta"+anlvl);


    // histoManager.AddVariable("TrkIsoTau"+anlvl,200,0,5,"TrkIso/p_{T} #mu_{1}","TrkIsoTau"+anlvl);
    // histoManager.AddVariable("TrkIsoLep"+anlvl,200,0,5,"TrkIso/p_{T} #mu_{2}","TrkIsoTau"+anlvl);
    // histoManager.AddVariable("EcalIsoTau"+anlvl,200,0,5,"EcalIso/p_{T} #mu_{1}","EcalIsoTau"+anlvl);
    // histoManager.AddVariable("EcalIsoLep"+anlvl,200,0,5,"EcalIso/p_{T} #mu_{2}","EcalIsoLep"+anlvl);
    // histoManager.AddVariable("HcalIsoTau"+anlvl,200,0,5,"HcalIso/p_{T} #mu_{1}","HcalIsoTau"+anlvl);
    // histoManager.AddVariable("HcalIsoLep"+anlvl,200,0,5,"HcalIso/p_{T} #mu_{2}","HcalIsoLep"+anlvl);
    // histoManager.AddVariable("IsoTau"+anlvl,200,0,5,"#Sigma iso/p_{T} #mu_{1}","IsoTau"+anlvl);
    // histoManager.AddVariable("IsoLep"+anlvl,200,0,5,"#Sigma iso/p_{T} #mu_{2}","IsoLep"+anlvl);

    //Jets
    histoManager.AddVariable("PtJ1"+anlvl,400,0,400,"p_{T} j_{1} [GeV]","PtJ1"+anlvl);
    histoManager.AddVariable("PtJ2"+anlvl,400,0,400,"p_{T} j_{2} [GeV]","PtJ2"+anlvl);
    histoManager.AddVariable("EtaJ1"+anlvl,100,-5,5,"#eta j_{1} ","EtaJ1"+anlvl);
    histoManager.AddVariable("EtaJ2"+anlvl,100,-5,5,"#eta j_{2} ","EtaJ2"+anlvl);
    histoManager.AddVariable("PhiJ1"+anlvl,128,-180,180,"#phi j_{1} [#circ]","PhiJ1"+anlvl);
    histoManager.AddVariable("PhiJ2"+anlvl,128,-180,180,"#phi j_{2} [#circ]","PhiJ2"+anlvl);
  
    histoManager.AddVariable("JJMass"+anlvl,500,0,500,"M_{jj} [GeV]","JJMass"+anlvl);
    histoManager.AddVariable("JJPt"+anlvl,400,0,400,"p_{T}(jj) [GeV]","JJPt"+anlvl);
    histoManager.AddVariable("JJY"+anlvl,100,-5,5,"Y_{jj}","JJY"+anlvl);
    histoManager.AddVariable("JJPhi"+anlvl,100,-3,3,"Y_{jj}","JJPhi"+anlvl);
  
    histoManager.AddVariable("JJDPhi"+anlvl,128,-180,180,"#Delta#phi_{jj} [#circ]","JJDPhi"+anlvl);
    histoManager.AddVariable("JJDEta"+anlvl,100,-5,5,"#Delta#eta_{jj} ","JJDEta"+anlvl);
  
    histoManager.AddVariable("J1BTag"+anlvl,440,-10,100,"TCHE j1","J1BTag"+anlvl);
    histoManager.AddVariable("J2BTag"+anlvl,440,-10,100,"TCHE j2","J2BTag"+anlvl);

    //Jets + leptons
    histoManager.AddVariable("TLJJDPhi"+anlvl,128,-180,180,"#Delta#phi (jj , #tau_{had}l_{#tau}) [#circ]","TLJJDPhi"+anlvl);
    histoManager.AddVariable("TLJJDY"+anlvl,100,-5,5,"#Delta Y (jj , #tau_{had}l_{#tau} )","TLJJDY"+anlvl);
    histoManager.AddVariable("TLJJMass"+anlvl,1000,0,1000,"M_{jj#tau_{had}l_{#tau}} [GeV]","TLJJMass"+anlvl);
    histoManager.AddVariable("TLJJPt"+anlvl,400,0,400,"p_{T} (jj#tau_{had}l_{#tau}) [GeV]","TLJJPt"+anlvl);
    histoManager.AddVariable("TLJJPhi"+anlvl,128,-180,180,"#phi (jj#tau_{had}l_{#tau}) [#circ]","TLJJPhi"+anlvl);
    histoManager.AddVariable("TLJJY"+anlvl,100,-5,5," Y (jj#tau_{had}l_{#tau})","TLJJY"+anlvl);

    histoManager.AddVariable("ST"+anlvl,100,0,1000,"S_{T} [GeV]","ST"+anlvl);

    histoManager.AddVariable("TauJJMass"+anlvl,1000,0,1000,"M_{jj#tau_{had}} [GeV]","TauJJMass"+anlvl);
    histoManager.AddVariable("LepJJMass"+anlvl,1000,0,1000,"M_{jjl_{#tau}} [GeV]","LepJJMass"+anlvl);
    histoManager.AddVariable("TLJ1Mass"+anlvl,1000,0,1000,"M_{j1#tau_{had}l_{#tau}} [GeV]","TLJ1Mass"+anlvl);
    histoManager.AddVariable("TLJ2Mass"+anlvl,1000,0,1000,"M_{j2#tau_{had}l_{#tau}} [GeV]","TLJ2Mass"+anlvl);

    histoManager.AddVariable("TauJ1Mass"+anlvl,1000,0,1000,"M_{j1#tau_{had}} [GeV]","TauJ1Mass"+anlvl);
    histoManager.AddVariable("LepJ1Mass"+anlvl,1000,0,1000,"M_{j1l_{#tau}} [GeV]","LepJ1Mass"+anlvl);
    histoManager.AddVariable("TauJ2Mass"+anlvl,1000,0,1000,"M_{j2#tau_{had}} [GeV]","TauJ2Mass"+anlvl);
    histoManager.AddVariable("LepJ2Mass"+anlvl,1000,0,1000,"M_{j2l_{#tau}} [GeV]","LepJ2Mass"+anlvl);

    histoManager.AddVariable("TauJ1DR"+anlvl,100,0,6,"#Delta R(#tau_{had}j1)","TauJ1DR"+anlvl);
    histoManager.AddVariable("TauJ2DR"+anlvl,100,0,6,"#Delta R(#tau_{had}j2)","TauJ2DR"+anlvl);
    histoManager.AddVariable("LepJ1DR"+anlvl,100,0,6,"#Delta R(l_{#tau}j1)","LepJ1DR"+anlvl);
    histoManager.AddVariable("LepJ2DR"+anlvl,100,0,6,"#Delta R(l_{#tau}j2)","LepJ2DR"+anlvl);

    histoManager.AddVariable("TauJ1DPhi"+anlvl,128,-180,180,"#Delta#phi(#tau_{had},j1) [#circ]","TauJ1DPhi"+anlvl);
    histoManager.AddVariable("TauJ2DPhi"+anlvl,128,-180,180,"#Delta#phi(#tau_{had},j2) [#circ]","TauJ2DPhi"+anlvl);
    histoManager.AddVariable("LepJ1DPhi"+anlvl,128,-180,180,"#Delta#phi(l_{#tau},j1) [#circ]","LepJ1DPhi"+anlvl);
    histoManager.AddVariable("LepJ2DPhi"+anlvl,128,-180,180,"#Delta#phi(l_{#tau},j2) [#circ]","LepJ2DPhi"+anlvl);

    histoManager.AddVariable("TauJ1DPhiLepJ2"+anlvl,128,-180,180,"#Delta#phi(#tau_{had}j1, l_{#tau}j2) [#circ]","TauJ1DPhiLepJ2"+anlvl);
    histoManager.AddVariable("TauJ2DPhiLepJ1"+anlvl,128,-180,180,"#Delta#phi(#tau_{had}j2, l_{#tau}j1) [#circ]","TauJ2DPhiLepJ1"+anlvl);
  
    histoManager.AddVariable("TauJ1DRLepJ2"+anlvl,100,0,6,"#Delta R(#tau_{had}j1, l_{#tau}j2)","TauJ1DRLepJ2"+anlvl);
    histoManager.AddVariable("TauJ2DRLepJ1"+anlvl,100,0,6,"#Delta R(#tau_{had}j2, l_{#tau}j1)","TauJ2DRLepJ1"+anlvl);

   

    //MET plots
    histoManager.AddVariable("MET"+anlvl,400,0,400,"#slash{E}_{T} [GeV]","MET"+anlvl);
    histoManager.AddVariable("METTLDPhi"+anlvl,128,-180,180,"#Delta#phi (#tau_{had}l_{#tau}, #slash{E}_{T}) [#circ]","METMMDPhi"+anlvl);
    histoManager.AddVariable("METJJDPhi"+anlvl,128,-180,180,"#Delta#phi (jj, #slash{E}_{T}) [#circ]","METJJDPhi"+anlvl);
    histoManager.AddVariable("METTauDPhi"+anlvl,128,-180,180,"#Delta#phi (#tau_{had, #slash{E}_{T}) [#circ]","METTauDPhi"+anlvl);
    histoManager.AddVariable("METLepDPhi"+anlvl,128,-180,180,"#Delta#phi (l_{#tau}, #slash{E}_{T}) [#circ]","METLepDPhi"+anlvl);
    histoManager.AddVariable("METJ1DPhi"+anlvl,128,-180,180,"#Delta#phi (j1, #slash{E}_{T}) [#circ]","METJ1DPhi"+anlvl);
    histoManager.AddVariable("METJ2DPhi"+anlvl,128,-180,180,"#Delta#phi (j2, #slash{E}_{T}) [#circ]","METJ2DPhi"+anlvl);


    //Lepton plots
    histoManager.AddVariable("TLCharge"+anlvl,3,-1.5,1.5,"C_{#tau}#timesC_{l}","TLCharge");


  }

  //===============================================================================================================

 //Misc
  histoManager.AddVariable("TauID",11,-0.5,10.5,"TauID","TauID");
  histoManager.AddVariable("ElecID",11,-0.5,10.5,"ElecID","ElecID");
  histoManager.AddVariable("MuID",11,-0.5,10.5,"MuID","MuID");
  histoManager.AddVariable("JetID",11,-0.5,10.5,"JetID","JetID");
  histoManager.AddVariable("JetMult",11,-0.5,10.5,"JetMult","JetMult");

 

  histoManager.AddVariable("HLT24Prescale",51,-0.5,50.5,"HLT24Prescale","HLT24Prescale");
  histoManager.AddVariable("HLT15Prescale",51,-0.5,50.5,"HLT24Prescale","HLT15Prescale");

 //Vertex
  histoManager.AddVariable("NVertex",20,0,20,"number of vertices","Vertex");
  histoManager.AddVariable("NVertexControl",20,0,20,"number of vertices","Vertex");
  
  





}





void 
LeptoQuarkUseTree::DisableBranchs(TTree* tree) {
  
  TObjArray* branchs =  tree->GetListOfBranches();
    string name;
    //cout<<" size "<<branchs->GetEntries()<<endl;
    for(int ib=0;ib<branchs->GetEntries();ib++) {
    
      name = (string)( ((*branchs)[ib])->GetName());
      cout<<name<<"   "<<((TBranch*)((*branchs)[ib]))->GetClassName()<<endl;
      tree->SetBranchStatus( name.c_str() , 0);
    }
    //  delete branchs;
    
}





bool
LeptoQuarkUseTree::HLTRequirement(VarClassLQ* vc, int nds) {

  if( vc->getSize("HLTResults") == 0) return false;
 
  if(nds==nt) { //data
    if ( vc->getI("HLTResults",39)==1
	 ||vc->getI("HLTResults",42)==1
	 ||vc->getI("HLTResults",43)==1
	 ||vc->getI("HLTResults",44)==1
	 ||vc->getI("HLTResults",45)==1
	 ||vc->getI("HLTResults",48)==1
	 ||vc->getI("HLTResults",51)==1
	 ||vc->getI("HLTResults",52)==1)
      return true;
  }
  else { //MC
    if ( vc->getI("HLTResults",38)==1
	 ||vc->getI("HLTResults",41)==1
	 ||vc->getI("HLTResults",42)==1
	 ||vc->getI("HLTResults",43)==1
	 ||vc->getI("HLTResults",44)==1
	 ||vc->getI("HLTResults",47)==1
	 ||vc->getI("HLTResults",50)==1
	 ||vc->getI("HLTResults",51)==1)
      return true;

  }

  return false;

}


bool
LeptoQuarkUseTree::ElectronID(VarClassLQ* vc) {

  vc->v_electrons_vbtf95.clear();
  vc->v_electrons.clear();

  size_t nEl = vc->getSize( "ElectronPt" );

  for(size_t ie=0;ie<nEl;ie++) {


    if( SimpleCut<double>(vc->getD("ElectronPt", ie ),20, "<" ) ) continue;
    
    bool ooEcal=false;
    if( fabs( vc->getD("ElectronEta", ie ) ) > 2.1 || 
	( fabs( vc->getD("ElectronEta", ie ) )< 1.566 
	  && fabs( vc->getD("ElectronEta", ie ) )> 1.4442 ) )
      ooEcal =true;
 
    if( SimpleCut<bool>( ooEcal, true, "=" ) )
      continue;


    float see = vc->getD("ElectronSigmaEtaEta",ie);
    float deta = vc->getD("ElectronDeltaEtaTrkSC",ie);
    float dphi = vc->getD("ElectronDeltaPhiTrkSC",ie);
    float hoe = vc->getD("ElectronHoE",ie);
    float iso = vc->getD("ElectronPFRelIso",ie);
    float edist = vc->getD("ElectronDist",ie);
    float dcot = vc->getD("ElectronDCotTheta",ie);

    bool quality_EB_loose = fabs(vc->getD("ElectronEta",ie)) <= 1.4442
      && see < 0.01
      && fabs(deta) < 0.007
      && fabs(dphi) < 0.8
      && hoe < 0.15
      && iso < 0.3;

    bool quality_EE_loose = fabs(vc->getD("ElectronEta",ie)) >= 1.566
      && see < 0.03
      && fabs(deta) < 0.01
      && fabs(dphi) < 0.7
      && hoe < 0.15
      && iso <0.3;
  
    // cout<<" el "<<ie<<"     "<<quality_EB_loose<<"    "<<quality_EE_loose<<"    "<<vc->getD("ElectronEta",ie)<<"    "<<see<<"   "<<deta<<"   "<<dphi<<"    "<<hoe<<"   "<<iso<<endl;

    if (quality_EB_loose || quality_EE_loose) 
      vc->v_electrons_vbtf95.push_back(ie);

    // General quality for conversion
    if ( SimpleCut<int>(vc->getI("ElectronMissingHits",ie), 0, ">" ) ) continue;
    bool quality_conv = fabs( edist ) < 0.02 && fabs( dcot ) < 0.02;
 
    //If the quality is satisfied, then electron is from conversion and needs to be rejected 
    if (SimpleCut<bool> ( quality_conv, true, "=" ) ) continue;
	
    // bool quality_EB_lowpt = fabs(vc->getD("ElectronEta",ie)) <= 1.4442
    //   && vc->getD("ElectronSigmaEtaEta",ie) < 0.01
    //   && vc->getD("ElectronDeltaEtaTrkSC",ie) < 0.004
    //   && vc->getD("ElectronDeltaPhiTrkSC",ie) < 0.03
    //   && (vc->getD("ElectronFBrem",ie) > 0.15 
    //       || 
    //       (fabs(vc->getD("ElectronEta",ie)) < 1. && vc->getD("ElectronEOverPin",ie) >0.95))
    //   && vc->getD("ElectronHoE",ie) < 0.025
    //   && vc->getD("ElectronPFRelIso",ie) < 0.1;
	
 
    bool quality_EB_highpt = fabs(vc->getD("ElectronEta",ie)) <= 1.4442
      && see < 0.01
      && fabs(deta) < 0.004
      && fabs(dphi) < 0.06
      && hoe < 0.04
      && iso < 0.1;
	
    // bool quality_EE_lowpt = fabs(vc->getD("ElectronEta",ie)) >= 1.566
    // && vc->getD("ElectronSigmaEtaEta",ie) < 0.03
    // && vc->getD("ElectronDeltaEtaTrkSC",ie) < 0.005
    // && vc->getD("ElectronDeltaPhiTrkSC",ie) < 0.02
    // && vc->getD("ElectronHoE",ie) < 0.1
    // && vc->getD("ElectronPFRelIso",ie) < 0.1;
	
    bool quality_EE_highpt = fabs(vc->getD("ElectronEta",ie)) >= 1.566
      && see < 0.03
      && fabs(deta) < 0.007
      && fabs(dphi) < 0.03
      && hoe < 0.1
      && iso < 0.1;
	
  
    if ( SimpleCut<bool> ( (quality_EB_highpt||quality_EE_highpt), true, "!=") )
      continue;
  
      vc->v_electrons.push_back(ie);
  }// End of loop over electrons

  // Move to next event if zero electron is found that passes the selection
  if ( SimpleCut<int>(vc->v_electrons.size(), 0, "=" ) ) return false;
 
  // Reject the event if there is an electron with OS and WP95 ID
  bool secondElOS = false;
  for (int unsigned i=0; i>vc->v_electrons_vbtf95.size(); ++i) {
    if ((vc->getI("ElectronCharge",vc->v_electrons[0])*vc->getI("ElectronCharge",vc->v_electrons_vbtf95[i]))<0)
      secondElOS = true;
  }

  if ( SimpleCut<bool>(secondElOS,true,"=") ) return false;

  vc->lep.SetPtEtaPhiM( vc->getD("ElectronPt",vc->v_electrons[0]),
			vc->getD("ElectronEta",vc->v_electrons[0]),
			vc->getD("ElectronPhi",vc->v_electrons[0]), 0. );

  vc->lepL = vc->v_electrons[0];

  return true;

}



bool
LeptoQuarkUseTree::VtxID(VarClassLQ* vc) {

  int nGoodVertices = 0;
  vc->zpv = -100;
  for(unsigned int ivtx = 0; ivtx != vc->getSize("VertexZ"); ++ivtx) {
    //double chi2 = vc->getD("VertexChi2",ivtx);
    double ndf  = vc->getD("VertexNDF",ivtx);
    double z    = vc->getD("VertexZ",ivtx);
    double rho  = vc->getD("VertexRho",ivtx);
    
    bool goodVtx = !vc->getB("VertexIsFake",ivtx) 
      && ndf > 4 
      && fabs(z) < 24 
      && fabs(rho)< 2;

    if (  SimpleCut<bool>(goodVtx,false,"=") ) continue;
    {
      ++nGoodVertices;
      if ((int)ivtx==vc->getI("ElectronVtxIndex", vc->lepL) ) 
	vc->zpv = z;
    }
  }
  if (nGoodVertices==0) return false;
  
  return true;

}


bool
LeptoQuarkUseTree::TauID(VarClassLQ* vc) {


  vector<int> v_taus;
  for (int unsigned itau=0; itau< vc->getSize("PFTauPt"); ++itau) {
	
    if( SimpleCut<double>(vc->getD("PFTauPt",itau), 20, "<" ) ) continue;
    if( SimpleCut<double>(fabs( vc->getD("PFTauEta",itau) ), 2.3, ">=") ) continue;
    //EB-EE cracks are removde by anti-e discriminator
 	
    int pfdecay = vc->getI("PFTauDecayModeFinding",itau);
    int iso = vc->getI("PFTauVeryLooseCombIsolation",itau);
    int eldis = vc->getI("PFTauAgainstElectronTight",itau);
    int mudis = vc->getI("PFTauAgainstMuonLoose",itau);
    int elmva = vc->getI("PFTauAgainstElectronMVA",itau);
    float pftauvtx = vc->getD("PFTauVertexZ",itau);

    bool quality = pfdecay == 1
      && iso == 1
      && eldis == 1
      && mudis == 1
      && elmva == 1
    && fabs( pftauvtx - vc->zpv) < 0.2;
    if ( SimpleCut<bool>(quality,true,"!=") ) continue;	
	
   
    //Check that the tau candidate does is away from any electron from looseelectrons
    //by at least dR=0.5
    bool tauCloseToEl = false;
    // cout<<" EV "<<"   "<<vc->v_electrons_vbtf95.size()<<endl;// if(!histoManager.GetInitStatus() )abort();
    for (int unsigned i=0; i<vc->v_electrons_vbtf95.size(); ++i){	  
      int i_el = vc->v_electrons_vbtf95[i];
      //   cout<<" ---> "<<vc->getD("PFTauEta",itau)<<"   "<<vc->getD("ElectronEta",i_el)<<"    "<<dR( vc->getD("PFTauEta",itau), vc->getD("ElectronEta",i_el) , vc->getD("PFTauPhi",itau), vc->getD("ElectronPhi",i_el) )<<endl;
      if ( dR( vc->getD("PFTauEta",itau), vc->getD("ElectronEta",i_el) ,
	       vc->getD("PFTauPhi",itau), vc->getD("ElectronPhi",i_el) ) < 0.5)
	tauCloseToEl = true;
    }
	
    // if ( SimpleCut<bool>(tauCloseToEl,true,"=") ) continue;             
	
    v_taus.push_back(itau);
  }//End of loop over taus

  if (v_taus.size()==0) return false;

  vc->tauL = v_taus[0];

  vc->tau.SetPtEtaPhiM( vc->getD("PFTauPt",v_taus[0]),
			vc->getD("PFTauEta",v_taus[0]),
			vc->getD("PFTauPhi",v_taus[0]), 0. );
  
  return true;

}



bool
LeptoQuarkUseTree::JetID(VarClassLQ* vc) {

  bool bt1=false, bt2=false; 

  vector<int> v_pfjets; vector<int> v_pfjetsBT;
  for (int unsigned ijet=0; ijet<vc->getSize("PFJetPt"); ++ijet) {
	
    if( SimpleCut<double>(vc->getD("PFJetPt", ijet ),30, "<") ) continue;
    if( SimpleCut<double>(fabs(vc->getD("PFJetEta", ijet ) ) ,2.4, ">" ) ) continue;
    
    if ( SimpleCut<int>(vc->getI("PFJetPassLooseID",ijet), 1, "!=" ) ) continue;
	
    
    if ( SimpleCut<float>( dR( vc->getD("PFJetEta", ijet ), vc->tau.Eta(), vc->getD("PFJetPhi", ijet ), vc->tau.Phi() ),
			   0.5, "<" ) )  continue;
	 
    if ( SimpleCut<float>( dR( vc->getD("PFJetEta", ijet ), vc->lep.Eta(), vc->getD("PFJetPhi", ijet ), vc->lep.Phi() ),
			   0.5, "<") )  continue;

    double tchel = vc->getD("PFJetTrackCountingHighEffBTag", ijet);
    
    if(!bt1 && tchel>1.7 ) bt1=true;
    if(bt1 && !bt2 && tchel>1.7 ) bt2=true;

    v_pfjets.push_back(ijet);
    if(tchel > 1.7 )
      v_pfjetsBT.push_back(ijet);
  }//End of loop over PF jets

  if( SimpleCut<unsigned int>(v_pfjets.size(),2, "<") ) return false;

  int jetLead = v_pfjets[0];
  int jetTrail = v_pfjets[1];
  vc->jetLead = jetLead;
  vc->jetTrail = jetTrail;
  
  if(applyBtag1) {
    if( v_pfjetsBT.size() == 0 ) vc->noBTag=true;
    else {
      vc->bJetLead = v_pfjetsBT[0];
      jetLead = jetLead>v_pfjetsBT[0]?jetLead:v_pfjetsBT[0];
      jetTrail = jetLead>v_pfjetsBT[0]?v_pfjetsBT[0]:jetTrail;
    }
  }
  if(applyBtag2) {
    if( v_pfjetsBT.size() <= 1 ) vc->noBTag=true;
    else {
      vc->bJetTrail = v_pfjetsBT[1];
      jetLead = v_pfjetsBT[0];
      jetTrail = v_pfjetsBT[1];
    }
  }
  
  
  vc->jet1.SetPtEtaPhiE( vc->getD("PFJetPt",jetLead),vc->getD("PFJetEta",jetLead),vc->getD("PFJetPhi",jetLead),vc->getD("PFJetEnergy",jetLead) );
  vc->jet2.SetPtEtaPhiE( vc->getD("PFJetPt",jetTrail),vc->getD("PFJetEta",jetTrail),vc->getD("PFJetPhi",jetTrail),vc->getD("PFJetEnergy",jetTrail) );

  return true;

}
