#include "METUseTree.hh"

#include <iomanip>
#include <TLorentzVector.h>
#include <boost/timer.hpp>
#include <boost/progress.hpp>

#include <TH1I.h>
#include <TObjArray.h>

using namespace std;

ClassImp(METUseTree)


METUseTree::METUseTree():
UseTree()
{
  
  LoadPUWeights();
  LoadVtxReweight();
  
  // addMETname("mva___mvaIdLoose_","mvaIdLoose");
  // addMETname("mva___mvaIdMedium_","mvaIdMedium");
  // addMETname("mva___mvaIdTight_","mvaIdTight");
  // addMETname("atlasAssoc___","atlasAssoc");
  // addMETname("atlas_10_","atlas10");
  // addMETname("atlas_15_","atlas15");
  // addMETname("atlas_20_","atlas20");
  // addMETname("atlas_25_","atlas25");
  // addMETname("atlas_30_","atlas30");
  // addMETname("atlas_35_","atlas35");
  // addMETname("basic___","basic");
  addMETname("pf___","pfT1");
  addMETname("pf___","pfT1Phi"); //phi
  addMETname("pf_pfMet__","pf");
  addMETname("pf_pfMet__","pfPhi");
  addMETname("pf_pfType1p0CorrectedMet__","pfT01");
  addMETname("pf_pfType1p0CorrectedMet__","pfT01Phi"); //Phi
  // addMETname("pf_pfType1p0p2CorrectedMet__","pfT012");
  // addMETname("pf_pfType1p0p2CorrectedMet__","pfT012Phi"); //Ph
 
  //FIXME
  addMETname("pf__Smear_","pfSmearT1"); //pfType1CorrectedMet
  addMETname("pf__Smear_","pfSmearT1Phi"); //Phi
  addMETname("pf_pfType1p0CorrectedMet_Smear_","pfSmearT01");
  addMETname("pf_pfType1p0CorrectedMet_Smear_","pfSmearT01Phi"); //phi
  addMETname("pat_patType1CorrectedPFMet__","patT1");
  addMETname("pat_patType1CorrectedPFMet__","patT1Phi");
 // addMETname("pat_patType1p2CorrectedPFMet__","patT12");
// addMETname("pat_patType1CorrectedPFMetElectronEnDown__","patT1ElED");
// addMETname("pat_patType1CorrectedPFMetElectronEnUp__","patT1ElEU");
// addMETname("pat_patType1CorrectedPFMetJetEnDown__","patT1JED");
// addMETname("pat_patType1CorrectedPFMetJetEnUp__","patT1JEU");
// addMETname("pat_patType1CorrectedPFMetMuonEnDown__","patT1MED");
// addMETname("pat_patType1CorrectedPFMetTauEnDown__","patT1JEU");
// addMETname("pat_patType1CorrectedPFMetTauEnUp__","patT1TED");
// addMETname("pat_patType1CorrectedPFMetUnclusteredEnDown__","patT1UED");
// addMETname("pat_patType1CorrectedPFMetUnclusteredEnUp__","patT1UEU");
  
//uncertainties
  addUnc("pat_patType1CorrectedPFMetJetEnDown", "JES", "Do");
  addUnc("pat_patType1CorrectedPFMetJetEnUp", "JES", "Up");
  // addUnc("pat_patType1CorrectedPFMetJetResDown", "JER", "Do");
  // addUnc("pat_patType1CorrectedPFMetJetResUp", "JER", "Up");
  //addUnc("pat_patType1CorrectedPFMetMuonEnDown", "MES", ""); //symetrized by default...
  //  addUnc("pat_patType1CorrectedPFMetMuonEnUp", "MES", "Up");
  addUnc("pat_patType1CorrectedPFMetUnclusteredEnDown", "UES", "Do");
  addUnc("pat_patType1CorrectedPFMetUnclusteredEnUp", "UES", "Up");


}

void  
METUseTree::LoadPUWeights() {
  
  TFile * file = new TFile("/home/mmarionn/Documents/CMS/METStudies/PUReduction/puWeightsS10.root","READ"); //puWeightsS10
  puweights = (TH1F*)file->Get("pileup");

}

float 
METUseTree::SearchWeight(float trueNint) {
  return puweights->GetBinContent( puweights->GetXaxis()->FindBin( trueNint)  );
}

void
METUseTree::AddVariables() {

}

void
METUseTree::FillMETTree() {

  //Add Variables to be computed

  PrepareHistograms();

  // string flagId="";
  // if(UsePdgId && pdgId==13) flagId="Mu";
  // if(UsePdgId && pdgId==11) flagId="Elec";

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
    
    //Déclaration des variables
    VarClassLQ vc;
    histoManager.InitStatus();
    
    //Skimming;
    //TFile* oFile;
    //TTree* skimtree;

    //Global variables ====== (saving time)
    TVector2 qT;
    TVector2 met;

    float sumEt;

    float smearB = 0.011; //0.011
    float smearE = 0.011; //0.015
    //====================================

    float Weight=1;
    
    cout<<" beginning tree  " <<name[i]<<"   "<<tChains[i]->GetEntries()<<endl;
    //Boucle et sélection, remplissage histos
    //int EP=0;
    int ent = tChains[i]->GetEntries();

    //protection against empty trees
    if(ent==0) continue;

    int runM=0;
    
    boost::progress_display show_progress( ent );
    for(int ie=-1;ie<ent;ie++) { //-1 for the initialization
      ++show_progress;

      //if(ie>10000) continue;
	 
      if(ie==0) {
	//	ie += 7080;
	vc.FinalizeInit();
	vc.BuildTree( tChains[i], 0 );
	// if(skimming) {
	//  oFile = new TFile( ("Skimming/"+name[i]+"_skim.root").c_str(),"RECREATE");
	//   tChains[i]->LoadTree(0);
	//   skimtree = tChains[i]->CloneTree(0);
	// }
	histoManager.StartFilling();
      }
      
      // if(ie%1000==0)
      // 	fprintf(stdout, "Done : %d / %d\r", ie, ent );
      

      tChains[i]->GetEntry(ie); 
      
      if(i!=nt ) {
	Weight = SearchWeight( vc.getF("trueNint") )*GetWeight(i, ie)*GetVtwW(vc.getI("nVtx") );  //FIXME
	//if(GetVtwW(vc.getI("nVtx") ) > 10000 || GetVtwW(vc.getI("nVtx") ) ==0 ) cout<<vc.getI("nVtx")<<"   "<<GetVtwW(vc.getI("nVtx"))<<endl;
	if(ie==-1) Weight =0;
      }
      else {
	Weight = 1;
	if(ie==-1) Weight =0;
      }

      if(vc.getI("run") > (size_t)runM )
	runM = vc.getI("run");

      if(i==nt)
      	{
      
      	  bool doubleCount=false;
	  std::pair<int,int> tmp(vc.getI("run"),vc.getI("event"));
	  EventIter = Events.find( tmp );
	  if(EventIter != Events.end() ) {
	    doubleCount=true;
	    //abort(); ?? FIXME
	  }
      	  if(doubleCount || (EventFilter && vc.getI("run")>(size_t)EventNum ) )
      	    {  continue; }
	  string t1(""),t2("");
	  std::pair<string,string> tmp2( t1, t2 );
	  
	  Events[ tmp ] = tmp2;
	  EvtsInFile.push_back(vc.getI("event") );
   
	  //}
      	}

      
 
      //All events
      //=====================================================================================
 
      if(!MakeCut<float>( vc.getF("massZ"), 70., "[]", i , "mass Cut", Weight, 120. )) continue;
      
      //Filling the METs
      
      bool isD= (i==nt);
      
      //FIXME
      //   if(!MakeCut<float>( vc.getF("ptZ"), 10., "[]", i, "qt cut",Weight,20.)) continue;
      // if(!MakeCut<float>( vc.getF("eta_m1"), -1.3, "]![", i, "eta1 cut",Weight,1.3)) continue;
      // if(!MakeCut<float>( vc.getF("eta_m2"), -1.3, "]![", i, "eta2 cut",Weight,1.3)) continue;

      
      if(i==nt) {
	l1.SetPtEtaPhiM(vc.getF("pt_m1"), vc.getF("eta_m1"), vc.getF("phi_m1"), 0.105);
	l2.SetPtEtaPhiM(vc.getF("pt_m2"), vc.getF("eta_m2"), vc.getF("phi_m2"), 0.105);
	
	l1 *= (fabs(vc.getF("eta_m1"))>1.3)?1.000824:1.001305;
	l2 *= (fabs(vc.getF("eta_m2"))>1.3)?1.000824:1.001305;
      }
      else {
	
	l1.SetPtEtaPhiM(vc.getF("pt_m1"), vc.getF("eta_m1"), vc.getF("phi_m1"), 0.105);
	l2.SetPtEtaPhiM(vc.getF("pt_m2"), vc.getF("eta_m2"), vc.getF("phi_m2"), 0.105);

	l1 *= rnd.Gaus(1,(fabs(vc.getF("eta_m1"))<1.3)?smearB:smearE); 
	l2 *= rnd.Gaus(1,(fabs(vc.getF("eta_m2"))<1.3)?smearE:smearB); 

	//uncertainty on mass and qT from muon pt measurement : conservative 2%
	l1U = l1;
	l2U = l2;
	l1U *= rnd.Gaus(1,0.02);
	l2U *= rnd.Gaus(1,0.02);
	histoManager.fillUnc("qT", "MES",i, (l1U+l2U).Pt(), Weight, "" );
	histoManager.fillUnc("mass", "MES",i, (l1U+l2U).M(), Weight, "" );

      }

      Zv4 = l1+l2;
      
      histoManager.fill("qT", i, Zv4.Pt() /*vc.getF("ptZ")*/ , Weight);
      histoManager.fill("mass", i, Zv4.M() /*vc.getF("massZ")*/ , Weight);

      //transverse mass to check the W+jet contribution
      float MT=sqrt( 2* vc.getF("pt_m1")*vc.getF("pf___pt")*(1-cos(dPhi( vc.getF("phi_m1"),vc.getF("pf___phi")))));

      histoManager.fill("Tmass", i, MT , Weight);

      qT.SetMagPhi( vc.getF("ptZ"), vc.getF("phiZ") ); //old Z
      
      histoManager.fill("ptl1",i,l1.Pt()/*vc.getF("pt_m1")*/,Weight);
      histoManager.fill("ptl2",i,l2.Pt()/*vc.getF("pt_m2")*/,Weight);
      histoManager.fill("etal1",i,l1.Eta()/*vc.getF("eta_m1")*/,Weight);
      histoManager.fill("etal2",i,l2.Eta()/*vc.getF("eta_m2")*/,Weight);
      histoManager.fill("phil1",i,l1.Phi()/*vc.getF("phi_m1")*/,Weight);
      histoManager.fill("phil2",i,l2.Phi()/*vc.getF("phi_m2")*/,Weight);

      histoManager.fill("njet",i,vc.getI("njet"),Weight);
      
      // MET ==========================================================
      //

      for(size_t im=0;im<mets.size();im++) {
	
      	//exceptions...
      	if(SimpleCut<bool>(isD,true,"=")) {
	  bool isExp=false;
      	  bool test = mets[im].second== "pf_pfType1CorrectedMet_Smear_";
      	  if(SimpleCut<bool>(test,true,"=") )
	    { met.SetMagPhi( vc.getF("pf___pt"), vc.getF("pf___phi") );isExp=true;
	      sumEt=vc.getF("pf___sumEt");}
	  
	  test = mets[im].second == "pf_pfType1PhiCorrectedMet_Smear_";
	  if(SimpleCut<bool>(test,true, "=") )
	    { met.SetMagPhi( vc.getF("pf_pfType1PhiCorrectedMet__pt"), vc.getF("pf_pfType1PhiCorrectedMet__phi") );isExp=true;
	      sumEt = vc.getF("pf_pfType1PhiCorrectedMet__sumEt");}
	  
	  test = mets[im].second == "pf_pfType1p0CorrectedMet_Smear_";
	  if(SimpleCut<bool>(test,true, "=") )
	    { met.SetMagPhi( vc.getF("pf_pfType1p0CorrectedMet__pt"), vc.getF("pf_pfType1p0CorrectedMet__phi") );isExp=true;
	      sumEt = vc.getF("pf_pfType1p0CorrectedMet__sumEt");}
	  
	  test = mets[im].second == "pf_pfType1p0PhiCorrectedMet_Smear_";
	  if(SimpleCut<bool>(test,true, "=") )
	    { met.SetMagPhi( vc.getF("pf_pfType1p0PhiCorrectedMet__pt"), vc.getF("pf_pfType1p0PhiCorrectedMet__phi") );isExp=true;
	      sumEt = vc.getF("pf_pfType1p0PhiCorrectedMet__sumEt");}

	  test = mets[im].second == "pf__Smear_";
      	  if(SimpleCut<bool>(test,true, "=") )
	    { met.SetMagPhi( vc.getF("pf___pt"), vc.getF("pf___phi") ); isExp=true;sumEt = vc.getF("pf___sumEt");}

	  if(!isExp) {
	    met.SetMagPhi( vc.getF(mets[im].second+"pt"), vc.getF(mets[im].second+"phi") );
	    sumEt = vc.getF(mets[im].second+"sumEt");
	  }
      	}
	else { //MC
	  met.SetMagPhi( vc.getF(mets[im].second+"pt"), vc.getF(mets[im].second+"phi") );
	  sumEt = vc.getF(mets[im].second+"sumEt");
	}

	met += Zv4.Vect().XYvector() - qT; 

      	if( (mets[im].first.find("Phi")!=(size_t)-1) ) {
      	  met = phiCorrection(met,vc.getI("nVtx"),isD,vc.getI("run"));
	  //sumEt = vc.getF(mets[im].second+"sumEt"); 
      	}

      	fillMETHistos(mets[im].first, met, qT, vc.getI("nVtx"), sumEt/*vc.getF(mets[im].second+"sumEt")*/ , i, Weight );
	fillMETUnc(mets[im].first, met, qT, i, Weight, &vc);

	// //Jet mult splitting
	{
	  ostringstream os;
	  os << min( (int)(qT.Mod()/50),4);//min(vc.getI("njet"),4);
	  fillMETHistos(mets[im].first+"NJ"+os.str(), met, Zv4.Vect().XYvector(), vc.getI("nVtx"), sumEt /*vc.getF(mets[im].second+"sumEt")*/ , i, Weight );
	}

	// //nvtx splitting
	// if(vc.getI("njet")==2 ){
	//   ostringstream os;
	//   os << min(vc.getI("nVtx"),19);
	//   fillMETHistos(mets[im].first+"NV"+os.str(), met, Zv4.Vect().XYvector(), vc.getI("nVtx"), sumEt , i, Weight );
	// }
	
      }

      //
      // MET ==========================================================
      

      histoManager.fill("NVertex", i, vc.getI("nVtx") , Weight);

      NumberEntries[i]++;
   
    }//End events
  //cout<<" Run Max "<<runM<<endl;
  //Pointer deletion =====
  
      //skimming
      // if(skimming) {
      // 	skimtree->Write();
      // 	oFile->Close();
      // }

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
METUseTree::PrepareDatasets() {
  
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


TH1F* 
METUseTree::GetVtxCorrection(string obs,int nds, int bin) {

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
METUseTree::GetFitWeight() {

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
METUseTree::AddTreeLoading( string type, vector<string> AddTNames ) {

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
METUseTree::SetAdditionnalTrees(string s, vector<string> vs ) {
  AddTreeLoading( s, vs );
}


TTree*
METUseTree::GetAdditionnalTree(string type, int ds ) {

  map<string, vector<TChain*> >::const_iterator iter;

  iter = AddChains.find(type);

  return (TTree*)(((*iter).second)[ds]);

}


void 
METUseTree::GetNProcEvent(string dataset) {

  string NameF = "data/"+reposi+"/"+dataset+".root";
  TFile* file = TFile::Open( NameF.c_str() );
  
  TH1I* htmp = (TH1I*)file->Get("all/all/nEvtProc__run");
  
  internalNEvt[ dataset ] = htmp->Integral(0,1001);
  
  delete htmp;
  file->Close();
}


void
METUseTree::PrepareHistograms() {

  for(size_t im=0;im<mets.size();im++) 
    prepareMETHistos(mets[im].first);
  
  for(size_t im=0;im<mets.size();im++) {
    for(int ij=0;ij<5;ij++) {
      ostringstream os;
      os<<ij;
      prepareMETHistos(mets[im].first+"NJ"+os.str());
    }
  }

  for(size_t im=0;im<mets.size();im++) {
    for(int in=0;in<20;in++) {
      ostringstream os;
      os<<in;
      prepareMETHistos(mets[im].first+"NV"+os.str());
    }
  }


  //Z
  histoManager.AddVariable("qT",400,0,400,"q_{T} [GeV]","qT");
  histoManager.AddVariable("mass",1600,0,400,"M_{#mu#mu} [GeV]","mass");
  histoManager.AddVariable("Tmass",1600,0,400,"M_{T}(#mu,#slash{E}_{T}) [GeV]","mass");

  histoManager.AddVariable("ptl1",400,0,200,"p_{T}(#mu_{1}) [GeV]","ptl1");
  histoManager.AddVariable("ptl2",400,0,200,"p_{T}(#mu_{2}) [GeV]","ptl2");
  histoManager.AddVariable("etal1",120,-3,3,"#eta(#mu_{1}) ","etal1");
  histoManager.AddVariable("etal2",120,-3,3,"#eta(#mu_{2}) ","etal2");
  histoManager.AddVariable("phil1",128,0,TMath::Pi()*2,"#phi(#mu_{1}) [rad]","phil1");
  histoManager.AddVariable("phil2",128,0,TMath::Pi()*2,"#phi(#mu_{2}) [rad]","phil2");

  histoManager.AddVariable("njet",21,-0.5,20.5,"Jet multiplicity","nJet");

  //Vertex
  histoManager.AddVariable("NVertex",50,0,50,"number of vertices","Vertex");
  histoManager.AddVariable("NVertexControl",50,0,50,"number of vertices","Vertex");

}


void
METUseTree::prepareMETHistos(string name) {

  histoManager.AddVariable( name+"MET",400,0,400,"#slash{E}_{T} [GeV]", name+"MET" );
  histoManager.AddVariable( name+"Phi",128,0,TMath::Pi()*2,"#phi(#slash{E}_{T}) [GeV]", name+"Phi" );
  histoManager.AddVariable( name+"X",400,-200,200,"#slash{E}_{T,X} [GeV]", name+"X" );
  histoManager.AddVariable( name+"Y",400,-200,200,"#slash{E}_{T,Y} [GeV]", name+"Y" );
  
  histoManager.AddVariable( name+"SumEt",1000,0,2000,"#Sum {E}_{T} [GeV]", name+"SumEt" );
 

  histoManager.AddVariable( name+"Upara",800,-600,200,"u_{||} [GeV]", name+"Upara" );
  histoManager.AddVariable( name+"Uperp",400,-200,200,"u_{#perp}   [GeV]", name+"Uperp" );
  histoManager.AddVariable( name+"redUpara",400,-200,200,"u_{||}+q_{T} [GeV]", name+"redUpara" );
  
  histoManager.AddProfVariable( name+"RespvsNvtx",40,0,40,"number of vertices","<u_{||}/q_{T}> [GeV]");

  histoManager.AddProfVariable( name+"XvsNvtx",40,0,40,"number of vertices","<#slash{E}_{T,X}> [GeV]");
  histoManager.AddProfVariable( name+"YvsNvtx",40,0,40,"number of vertices","<#slash{E}_{T,Y}> [GeV]");

  histoManager.AddProfVariable( name+"XvsSumEt",40,0,40,"#sum E_{T} [GeV] ","<#slash{E}_{T,X}> [GeV]");
  histoManager.AddProfVariable( name+"YvsSumEt",40,0,40,"#sum E_{T} [GeV] ","<#slash{E}_{T,Y}> [GeV]");
}

void
METUseTree::addMETname(string tn, string name) {
  pair<string,string> p(name, tn);
  mets.push_back( p );
}

void
METUseTree::addUnc(string tn, string name, string dir) {
  pair<string,string> p(name, dir);
  uncMap[ tn ]=p;
}


void 
METUseTree::fillMETHistos(string name, TVector2 met, TVector2 qT, int nvtx, float sumEt, int i, float Weight) {

  TVector2 u(0,0);
  u -= (met +qT);

  histoManager.fill( name+"MET",i, met.Mod() , Weight );  
  histoManager.fill( name+"Phi",i, met.Phi() , Weight );  

  histoManager.fill( name+"X",i, met.X() , Weight );  
  histoManager.fill( name+"Y",i, met.Y() , Weight ); 
  histoManager.fill( name+"SumEt",i, sumEt , Weight ); 

  float upara = (float)(( qT*u)/qT.Mod());
  u = u.Rotate( TMath::Pi()/2);
  float uperp = (float)(( qT*u)/qT.Mod());
  histoManager.fill( name+"Upara",i, upara , Weight );
  histoManager.fill( name+"Uperp",i, uperp , Weight );
  histoManager.fill( name+"redUpara",i, upara+qT.Mod() , Weight ); 

  if(qT.Mod()>10)
    histoManager.fillProf( name+"RespvsNvtx", i, nvtx, -upara/qT.Mod(), Weight );

  histoManager.fillProf( name+"XvsNvtx", i, nvtx, met.X(), Weight );
  histoManager.fillProf( name+"YvsNvtx", i, nvtx, met.Y(), Weight );

  histoManager.fillProf( name+"XvsSumEt", i, sumEt, met.X(), Weight );
  histoManager.fillProf( name+"YvsSumEt", i, sumEt, met.Y(), Weight );

}

void
METUseTree::fillMETUnc(string name, TVector2 met, TVector2 qT, int i,float Weight, VarClassLQ* vc) {

  if( (name!="pfT01Phi" && name!="patT1Phi"&& name!="pfSmearT01Phi") || nt==i ) return;

  mettmp.SetMagPhi(vc->getF("pat_patType1CorrectedPFMet__pt"),vc->getF("pat_patType1CorrectedPFMet__phi"));
  //cout	<<mettmp.Mod()<<"   "<<mettmp.Phi()<<endl;
  if( (name.find("Phi")!=(size_t)-1) ) {
    mettmp = phiCorrection(mettmp, vc->getI("nVtx"),false, vc->getI("run"));
    
  }
  
  for(itUM=uncMap.begin();itUM!=uncMap.end();itUM++) {
    mettmpUnc.SetMagPhi(vc->getF(itUM->first+"__pt"), vc->getF(itUM->first+"__phi"));
    if( (name.find("Phi")!=(size_t)-1) ) {
      mettmpUnc = phiCorrection(mettmpUnc, vc->getI("nVtx"),false, vc->getI("run"));
    }
    metUnc = mettmpUnc.Mod()/mettmp.Mod()*met.Mod();
    angleUnc = mettmpUnc.Phi()-mettmp.Phi() + met.Phi(); //works only on 0 to 2 pi

    // cout<<metUnc<<"  "<<angleUnc<<" / "
    // 	<<met.Mod()<<"   "<<met.Phi()<<" /  "
    // 	<<mettmp.Mod()<<"   "<<mettmp.Phi()<<" /  "
    // 	<<mettmpUnc.Mod()<<"   "<<mettmpUnc.Phi()<<endl;

    mettmpUnc.SetMagPhi(metUnc,angleUnc);
    
  
    TVector2 u(0,0);
    u -= (mettmpUnc + qT ); 
    
    histoManager.fillUnc(name+"MET", itUM->second.first,i, mettmpUnc.Mod(), Weight, itUM->second.second );
    histoManager.fillUnc(name+"Phi", itUM->second.first,i, mettmpUnc.Phi(), Weight, itUM->second.second );

    histoManager.fillUnc(name+"X", itUM->second.first,i, mettmpUnc.X(), Weight, itUM->second.second );
    histoManager.fillUnc(name+"Y", itUM->second.first,i, mettmpUnc.X(), Weight, itUM->second.second ); 

    float upara = (float)(( qT*u)/qT.Mod());
    u = u.Rotate( TMath::Pi()/2);
    float uperp = (float)(( qT*u)/qT.Mod());
    histoManager.fillUnc( name+"Upara", itUM->second.first,i, upara , Weight, itUM->second.second );
    histoManager.fillUnc( name+"Uperp", itUM->second.first,i, uperp , Weight, itUM->second.second );
    histoManager.fillUnc( name+"redUpara", itUM->second.first,i, upara+qT.Mod() , Weight, itUM->second.second ); 
  }
  
}


TVector2 
METUseTree::phiCorrection(TVector2 met, int Nvtx, bool isD, int run) {

  float corX,corY;
  //  cout<<" tapoué "<<endl;
  //Old corrections 2012A
  if(isD) {
    // //2012A
    // corX = +3.54233e-01 + 2.65299e-01*Nvtx;
    // corY = +1.88923e-01 - 1.66425e-01*Nvtx;

    //pat 52x
    // corX = -0.087244 +0.305351*Nvtx ;
    // corY = -0.0156588 -0.20106*Nvtx ;

    //2012AB
    // corX = +1.68804e-01 + 3.37139e-01*Nvtx;
    // corY = -1.72555e-01 - 1.79594e-01*Nvtx;

    //2012ABC
    // corX = 0.122578 + 0.332148*Nvtx  -0.07908;
    // corY = 0.0103295 - 0.227043*Nvtx -0.1135;

    //CHristian corrections 2012ABC
    corX = +2.87340e-01 + 3.29813e-01*Nvtx; 
    corY = -2.27938e-01 - 1.71272e-01*Nvtx;

  }
  else {
    //2012A
    // corX = -2.99576e-02 - 6.61932e-02*Nvtx;
    // corY = +3.70819e-01 - 1.48617e-01*Nvtx;

    //2012AB
    // corX = +2.22335e-02 - 6.59183e-02*Nvtx;
    // corY = +1.52720e-01 - 1.28052e-01*Nvtx;

    //2012AB MM
    // corX = 0.160696 - 0.0778255*Nvtx -0.0629;
    // corY = 0.680142  -0.154409*Nvtx -0.17836;

    //CHristian corrections 2012ABC
    corX = +8.72683e-02 - 1.66671e-02*Nvtx;
    corY = +1.86650e-01 - 1.21946e-01*Nvtx;

  }

  // if(isD) {
  //   if(run<193621) { //2012A
  //   corX = +0.07439 +0.3374 *Nvtx;
  //   corY = +0.4649 -0.2518 *Nvtx;
  //   }
  //   else if(run< 197770) { // 2012B
  //     corX = -0.4144 +0.3911 *Nvtx;
  //     corY = -0.0322 -0.2661 *Nvtx;
  //   }
  //   else {//2012 C
  //     corX = -0.2959 +0.3806 *Nvtx;
  //     corY = 0.08905 -2.596 *Nvtx;
  //   }
  // }
  // else {
  //   corX = +0.09972 - 0.003793*Nvtx;
  //   corY = +0.4434 - 0.2117*Nvtx;
  // }


  


  TVector2 corMET(0,0);
  corMET.Set( met.X() - corX, met.Y() -corY );

  // cout<<met.X()<<"  "<<corX<<"  "<<corMET.X()<<" :::: "
  //     <<met.Y()<<"  "<<corY<<"  "<<corMET.Y()<<endl;

  return corMET;

}




// dummy vertex reweighting
void METUseTree::LoadVtxReweight() {

  TH1F *ratioNVertex__1 = new TH1F("ratioNVertex__1","92",50,0,50);
  ratioNVertex__1->SetBinContent(1,4.352382);
   ratioNVertex__1->SetBinContent(2,0.912782);
   ratioNVertex__1->SetBinContent(3,0.9450765);
   ratioNVertex__1->SetBinContent(4,0.9245538);
   ratioNVertex__1->SetBinContent(5,0.9409372);
   ratioNVertex__1->SetBinContent(6,0.9554647);
   ratioNVertex__1->SetBinContent(7,0.9564364);
   ratioNVertex__1->SetBinContent(8,0.962751);
   ratioNVertex__1->SetBinContent(9,0.9730481);
   ratioNVertex__1->SetBinContent(10,0.9775884);
   ratioNVertex__1->SetBinContent(11,0.987715);
   ratioNVertex__1->SetBinContent(12,0.9920352);
   ratioNVertex__1->SetBinContent(13,1.002308);
   ratioNVertex__1->SetBinContent(14,1.001181);
   ratioNVertex__1->SetBinContent(15,1.006853);
   ratioNVertex__1->SetBinContent(16,1.009103);
   ratioNVertex__1->SetBinContent(17,1.014458);
   ratioNVertex__1->SetBinContent(18,1.01718);
   ratioNVertex__1->SetBinContent(19,1.020113);
   ratioNVertex__1->SetBinContent(20,1.02148);
   ratioNVertex__1->SetBinContent(21,1.029908);
   ratioNVertex__1->SetBinContent(22,1.027123);
   ratioNVertex__1->SetBinContent(23,1.028488);
   ratioNVertex__1->SetBinContent(24,1.038353);
   ratioNVertex__1->SetBinContent(25,1.028244);
   ratioNVertex__1->SetBinContent(26,1.032906);
   ratioNVertex__1->SetBinContent(27,1.041499);
   ratioNVertex__1->SetBinContent(28,1.061431);
   ratioNVertex__1->SetBinContent(29,1.049645);
   ratioNVertex__1->SetBinContent(30,1.027846);
   ratioNVertex__1->SetBinContent(31,1.042164);
   ratioNVertex__1->SetBinContent(32,1.079537);
   ratioNVertex__1->SetBinContent(33,1.151077);
   ratioNVertex__1->SetBinContent(34,1.081402);
   ratioNVertex__1->SetBinContent(35,1.116077);
   ratioNVertex__1->SetBinContent(36,1.025585);
   ratioNVertex__1->SetBinContent(37,1.196205);
   ratioNVertex__1->SetBinContent(38,1.240716);
   ratioNVertex__1->SetBinContent(39,0.8237543);
   ratioNVertex__1->SetBinContent(40,0.9180249);
   ratioNVertex__1->SetBinContent(41,0.7386603);
   ratioNVertex__1->SetBinContent(42,0.4086069);
   ratioNVertex__1->SetBinContent(43,0.3347449);
   ratioNVertex__1->SetBinContent(44,0.3824325);
   ratioNVertex__1->SetBinContent(45,0.2143044);
   ratioNVertex__1->SetBinContent(48,13.12211);


   vtxW = (TH1F*)ratioNVertex__1->Clone();
   vtxW->Reset("ICEM");
   for(int i=0;i<ratioNVertex__1->GetNbinsX()+2;i++) {
     vtxW->SetBinContent(i, ratioNVertex__1->GetBinContent(i));
    //  if(ratioNVertex__1->GetBinContent(i)!=0)
   //     vtxW->SetBinContent(i, 1./ratioNVertex__1->GetBinContent(i));
   //   else
   //     vtxW->SetBinContent(i, 0.);
   }
}

float METUseTree::GetVtwW(int nvtx) {
  //if(nvtx==44)cout<<vtxW->GetBinContent( vtxW->GetXaxis()->FindBin( (float)nvtx)+1 )<<endl;
  //cout<<vtxW->GetBinContent( vtxW->GetXaxis()->FindBin( (float)nvtx) )<<endl;
  if(nvtx>= vtxW->GetNbinsX()) return 0;//nvtx=vtxW->GetNbinsX();
  return vtxW->GetBinContent( vtxW->GetXaxis()->FindBin( (float)nvtx)+1 );
}
