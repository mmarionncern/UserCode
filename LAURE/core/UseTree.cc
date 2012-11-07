#include "UseTree.hh"

#include <iomanip>

using namespace std;


ClassImp(UseTree)

UseTree::UseTree() :
savedMC(0)
{
  
  METType=0;
  MEtType="";
  METcut =0;
  PTcut =0;
  invCut=false;

  Bin=2;

  RangeY[0]=0.02;  RangeY[1]= 100;
  RangeX[0]=-100;  RangeX[1]= 100; 

  unbinned=true;
  response=false;

  Suffix.push_back("Pf");
  Suffix.push_back("Tc");
  Suffix.push_back("Calo");
  Suffix.push_back("CaloTypeI");
  Suffix.push_back("CaloTypeII");
  Suffix.push_back("PfTypeI");

  AddSyst = false;

  IDs.push_back("sigieie");
  IDs.push_back("Deta");
  IDs.push_back("Dphi");
  IDs.push_back("HoE");
  
  Isos.push_back("TrackIso");
  Isos.push_back("EcalIso");
  Isos.push_back("HcalIso");

  TGaxis::SetMaxDigits(3);

  nt=-1;
  
  treenames.push_back("METCommVariables");
  
  saveFile = false;

  ShowDMCRatio= false;
  //  savedMC = new TH1F(NULL);

  inAcc=true;
  useAccForEff =false;


  //Plotstyle by default not defstyle
  defStyle =false;


  _isData=false;

  //Skimming

  wDs="";
  _storeTuple=false;
  oFile = NULL;
  skimTree = NULL;

}




void UseTree::cmsPrel(double intLumi) {
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);

  //  latex.SetTextAlign(31); // align right
  // latex.DrawLatex(0.90,0.96,"#sqrt{s} = 7 TeV");

  if (intLumi > 0.) {
    latex.SetTextAlign(31); // align right

    if(intLumi > 1000 )
      latex.DrawLatex(0.886,0.862,Form(" %.1f fb^{-1}  at   #sqrt{s} = 8 TeV",intLumi/1000.));
    else
      latex.DrawLatex(0.886,0.862,Form(" %.0f pb^{-1}  at   #sqrt{s} = 8 TeV",intLumi));
  }
  latex.SetTextAlign(11); // align left
  latex.DrawLatex(0.54,0.96,"CMS preliminary 2012");
}


void UseTree::Delete() {

  // delete c2;
  // for(int unsigned i=0;i< files.size();i++) {
  //   delete files[i];
  //   delete trees[i];
  // }
  
  for(int unsigned i=0;i<tChains.size();i++)
    delete tChains[i];
  
  /* for(int unsigned i=0;i<Histos.size();i++)
    delete Histos[i];
  for(int unsigned i=0;i<Histos2D.size();i++)
    delete Histos2D[i];
  
  

  delete leg;
  delete leg2D;
  */

}


bool UseTree::Cbool(bool skip, bool bname) {

  if(!skip)
    return true;
  else{
    if(bname)
      return false;
    else return true;
  }

}

void UseTree::ConfigureIDIso(int iso, int id) {
  srand(time(NULL)); 
  int IDtmp= 0;
  int ISOtmp= 0;
  
  bool noIso=false;
  bool noID=false;

  if(iso ==70)
    ISOtmp =0;
  else if(iso ==80)
    ISOtmp =1;
  else if(iso ==90)
    ISOtmp =2;
  else if(iso ==95)
    ISOtmp =3;
  else { cout<<"Warning, ISO not coded, no ISO cuts "<<endl; noIso=true;}

  if(id ==70)
    IDtmp =0;
  else if(id ==80)
    IDtmp =1;
  else if(id ==90)
    IDtmp =2;
  else if(id ==95)
    IDtmp =3;
  else { cout<<"Warning, ID not coded, no ID cuts "<<endl; noID=true;}
   
  if(!noIso)
    for(int j=0;j<2;j++)
      for(int i=0;i<3;i++) {
	_Isocuts[j][i] = _IsoCuts[j][ISOtmp][i];
      }
  if(!noID)
    for(int j=0;j<2;j++)
      for(int i=0;i<4;i++) {
	_IDcuts[j][i]= _IDCuts[j][IDtmp][i];
      }
  cout<<" Using Predef cuts : iso "<<iso<<"  ; id "<<id<<endl;
  ISO = iso;
  ID =id;
}

int 
UseTree::GetNbyName(string sname) {
  int N=nt-1;

  for(int unsigned i=0;i<name.size();i++) {
    if(sname==name[i]) {
      N=i;
      break;
    }
  }
  return N;

}

void UseTree::ConfigurePlots(int Binning,int AddBinBkg, bool overflow, bool underflow, double rangeY[2], double rangeX[2],string ytitle, int XDiv[3], int YDiv[3], float markerS, float lineW,vector<int> lines,vector<string> sevObs,  bool Draw3o1,bool sucolor,bool switchrms, string rmsopt,bool FillProfile,bool basicProf, bool ShowRatio, bool logY, bool addSyst ) {

  Bin = Binning;
  BinBkg = AddBinBkg;
  OverFlowBin = overflow;
  UnderFlowBin = underflow;
  ShowDMCRatio = ShowRatio;
  if(UnderFlowBin)
    OverFlowBin=true;
  RangeY[0] = rangeY[0];
  RangeY[1] = rangeY[1];
  RangeX[0] = rangeX[0];
  RangeX[1] = rangeX[1];
 
  YTitle = ytitle;

  logYScale=logY;

  AddSyst = addSyst;

  for(int i=0;i<3;i++) {
    Xdiv[i] = XDiv[i];
    Ydiv[i] = YDiv[i];
  }

  Draw3on1=Draw3o1;
  
  Lines = lines;

  SuperColor =sucolor;

  SevObs = sevObs;
  /* for(int i=0;i<3;i++)
     cout<<" test "<<SevObs[i]<<endl;*/

  switchRMS = switchrms;
  FillProf = FillProfile;
  BasicProf = basicProf;
  RMSOpt = rmsopt;

  MarkerSize= markerS;
  LineWidth = lineW;

}

void UseTree::ConfigureEndAnalysis(bool ForceCut,int iso, int id,float IdCuts[2][4],
				   float IsoCuts[2][3],float PtCut, float MetCut, 
				   string MetType,float MTcut,int nvert,string Normalisation) {

  ConfigureIDIso(iso, id);

  if(ForceCut) {
    for(int j=0;j<2;j++) {
      for(int i=0;i<4;i++){
	if(i<3){
	  if(IsoCuts[j][i]!=-1) {
	    cout<<" Cut Forced "<<Isos[j][i]<<"  "<<_IsoCuts[j][i]<<endl;
	    _Isocuts[j][i]= IsoCuts[j][i];
	  }
	}
	if(IdCuts[j][i]!=-1) {
	  cout<<" Cut Forced "<<IDs[j][i]<<"  "<<_IDCuts[j][i]<<endl;
	  _IDcuts[j][i] = IdCuts[j][i];
	}
      }
    }
  }
  
  NVert=nvert;
  MTCut=MTcut;
  METcut = MetCut;
  PTcut = PtCut;
  MEtType = MetType;
  if(MetType == "pfMET")
    METType = 0;
  else if(MetType == "caloMET")
    METType= 1;
  else if(MetType == "tcMET")
    METType= 2;
 else if(MetType == "rcMET")
    METType= 6;
 else if(MetType == "cT2MET")
    METType= 8;
  else
    {cout<<"Warning, no such MET "
	 <<MetType
	 <<" ; using pfMET by default"<<endl; }
  
  ParserNorm(Normalisation);


}

void UseTree::ConfigureData(bool eventfilter, int eventnum,bool MCOnly ) {
  EventFilter=eventfilter;
  EventNum=eventnum;

  NoData = MCOnly;

}



void UseTree::DrawLine(string variable, int WP, int EP,int rangeY) {

  int nWP =3;

  bool lineOK=true;

  if(WP ==70)
    nWP =0;
  else if(WP ==80)
    nWP=1;
  else if(WP ==90)
    nWP=2;
  else if(WP ==95)
    nWP=3;
  else { cout<<"No such WP "<<endl;lineOK=false; return;}
  
  float XLine=0;

  if(variable == "sigieie") 
    XLine = _IDCuts[ EP ][ nWP ][0];
  else if(variable == "Deta") 
    XLine = _IDCuts[ EP ][ nWP ][1];
  else if(variable == "Dphi") 
    XLine = _IDCuts[ EP ][ nWP ][2];
  else if(variable == "HoE") 
    XLine = _IDCuts[ EP ][ nWP ][3];
  
  else if(variable == "TrackIso") 
    XLine = _IsoCuts[ EP ][ nWP ][0];
  else if(variable == "EcalIso") 
    XLine = _IsoCuts[ EP ][ nWP ][1];
  else if(variable == "HcalIso") 
    XLine = _IsoCuts[ EP ][ nWP ][2];
  else { cout<<" No such observable, check name; No limit line drawn"<<endl;lineOK=false; return;}

  float YLinem = rangeY*0.75;
 
  float YLineM = rangeY*0.90;

  cout<<" Line "<<XLine<<"  "<<YLinem<<"  "<<YLineM<<endl;
  

  TArrow* line = new TArrow(XLine,YLinem,XLine,YLineM);
  line->SetLineColor(2);
  line->SetLineWidth(2);
  line->SetArrowSize(0.05);
  line->SetAngle(60);
  line->SetOption("<");
  if(lineOK)
    line->Draw("same");

}

void UseTree::SetTreeName(string tname) {
  
  if(treenames[0] == "METCommVariables")
    { treenames.clear(); }
  
  treenames.push_back( tname );

}

void UseTree::EndAnalysis( vector<string> Data, string Observable, 
			   int ecalP, bool ConvRej, bool VetoE, bool Invc,
			   bool SkipCut, string SkipVar, bool NoDeltaCor,
			   string QCDtype,string Wtype, string rep) {

  observable=Observable;
  EcalP = ecalP;
  
  vetoE = VetoE;
  invCut = Invc;
  skipCut = SkipCut;
  Nm1Var = SkipVar;
  noDeltaCor = NoDeltaCor;
  convRej = ConvRej;

   //Put Datasets
  QCDType = QCDtype;
  WType = Wtype;
  data =Data;
  reposi = rep;

  nt += 1; //to get exactly the good number
 
  PrepareDatasets();
 

  //Config analysis , number of datatype (data not with it)
  histoManager.ConfigAnalysis(nt);


  //Preparation des TChains
   for(int i=0;i<nt+1;i++) {
     TChain* ctmp = new TChain(treenames[0].c_str());
     tChains.push_back(ctmp);
   }
   
   cout<<" debut lecture "<<reposi<<endl;
   string path(getenv("UseTree"));
   TFile* datafile;
   if(!NoData) {
     for(size_t i=0;i<data.size();i++) {
       if(data[i]!="") {
	 string NameF = path+"/data/"+reposi+"/"+data[i]+".root"; 
	 if(reposi.find(":")!=(size_t)-1) NameF=reposi+"/"+data[i]+".root";	 

	 cout<< NameF <<endl;
	 datafile = TFile::Open(NameF.c_str());// new TFile(NameF.c_str(), "READ");
	 if(datafile==NULL) { cout<<" No such file "<<histoManager.Name<<endl; return;}
	 if(treenames.size()!=1) {
	   for(size_t it=0; it<treenames.size();it++) {
	     TTree* tmptree = (TTree*)datafile->Get( (treenames[it]).c_str() );
	     if(tmptree != NULL ) {
	       cout <<"Data Tree n"<<it<<"  --> "<<NameF+"/"+treenames[it]
		    <<"   "<<tmptree->GetEntries()<<"   "<<tChains[ nt ]->GetEntries()<<endl;
	       tChains[nt]->Add( (NameF+"/"+treenames[it]).c_str()); }
	     delete tmptree;
	   }
	 }
	 else{
	   tChains[nt]->Add(NameF.c_str());
	 }
	 cout<<NameF<<" entries -->  "<<tChains[nt]->GetEntries()<<endl;
	 if(tChains[nt]->GetEntries()==0)
	   { cout<<"No data in data file"<< path+"/data/"+reposi+"/"+data[i]+".root"<<endl;}
       
	 datafile->Close();
       }
     }
   }
   
  for(int i=0;i<nt;i++) {
    vector<std::pair<string,int> > tmp;
    ChanNum.push_back(tmp);
  }

  int nTmp=0; int ndt=0;
  for(size_t i=0;i<datasets.size();i++) {
    string NameF = path+"/data/"+reposi+"/"+datasets[i]+".root";
    if(reposi.find(":")!=(size_t)-1) NameF=reposi+"/"+datasets[i]+".root";

    if(ndt!=dt[i])
      nTmp=0;
    ndt=dt[i];
    
    datafile = TFile::Open(NameF.c_str()); //new TFile(NameF.c_str(), "READ");
    if(treenames.size()!=1) {
      for(size_t it=0; it<treenames.size();it++) {
	TTree* tmptree = (TTree*)datafile->Get( (treenames[it]).c_str() );
	if(tmptree != NULL ) {
	  tChains[ dt[i] ]->Add( (NameF+"/"+treenames[it]).c_str()); }
	delete tmptree;
	cout <<" Tree n"<<it<<"  --> "<<NameF+"/"+treenames[it]<<"   "<<tChains[ dt[i] ]->GetEntries()<<endl;
      }
    }
    else{
      tChains[ dt[i] ]->Add(NameF.c_str());
    }
    datafile->Close();

    cout<<" Adding "<<reposi+"/"+datasets[i]+".root"<<"  to "<<name[ dt[i] ]
    	<<"   :  nEvt "<<tChains[ dt[i] ]->GetEntries()<<endl;
    
    FillEventMap(datasets[i],nTmp, dt[i] );
    
     nTmp = tChains[ dt[i] ]->GetEntries();
  }
  for(int i=0;i<nt;i++) {
    FillEventMap("End",tChains[i]->GetEntries(), i );
  }

  cout<<" data ---> "<<tChains[ nt ]->GetEntries()<<endl;

  for(int i=0;i<nt+1;i++) {
    NEntries.push_back(tChains[i]->GetEntries());
  }

  //Trees charg√©s et combin√©s

  //Prepare renorm shapes
  for(int i=0;i<nt;i++) {
    shape.push_back(NULL);
  }
}


void
UseTree::SavePlot(bool sfile) {
  saveFile = sfile;
}


TTree* UseTree::LoadChain(int dnum) {

  TChain* ctmp = new TChain(treenames[0].c_str()); //
  TFile* tmp;

  cout<<" Loading chain, take care !!!!!!!!!!!!!!!!!!!!!!! "<<endl;
    for(int i=0;i<nt;i++) {
      vector<std::pair<string,int> > tmp;
      ChanNum.push_back(tmp);
    }

    int nTmp=0; int ndt=0;
    for(size_t i=0;i<datasets.size();i++) {
      string NameF = "data/"+reposi+"/"+datasets[i]+".root";
      if(reposi.find(":")!=(size_t)-1) NameF=reposi+"/"+datasets[i]+".root";

      if(dt[i]!=dnum) continue; //Select only good DS

      if(ndt!=dt[i])
	nTmp=0;
      ndt=dt[i];
      
      tmp=TFile::Open(NameF.c_str());
      //  ctmp->Add(NameF.c_str());
      if(treenames.size() ==1 ) {
	ctmp = (TChain*)tmp->Get(treenames[0].c_str());
      }
      else {
	ctmp->Clear();
	for(size_t it=0;it<treenames.size();it++) {
	 
	  tmp=TFile::Open(NameF.c_str());
	  ctmp->Add( (NameF+"/"+treenames[it]).c_str());
	  cout<< "several contributions to "<<datasets[i]<<" cumulated-> "<<ctmp->GetEntries()<<endl;
	}
      }

      cout<<" Adding "<<reposi+"/"+datasets[i]+".root"<<"  to "<<name[ dnum ]
	  <<"   :  nEvt "<<ctmp->GetEntries()<<endl;
    
      FillEventMap(datasets[i],nTmp, dnum );
     
      nTmp += ctmp->GetEntries();
    }
    for(int i=0;i<nt;i++) {
      FillEventMap("End",ctmp->GetEntries(), i );
    }
    //  }

    delete ctmp;
    tmp->Close();
    ctmp =NULL;
  //NEntries.push_back(ctmp->GetEntries());

  //delete ctmp;
  
  return ctmp;
}

void UseTree::SaveAllPlots(string path, bool erase) {

  Draw3on1=false;
  defStyle=true;
  SuperColor=false;
  ShowDMCRatio =true;
 
  if(erase) {
    TString command_ = TString("rm -r ") + path+"/"; 
    system( command_.Data());
  }

  vector<string> vars = histoManager.GetVariables(); 
  for(unsigned int iv=0;iv<vars.size();iv++) {

    string obs = vars[iv];
    
    c2->cd();
    
    vector<TPad*> pad_;
    vector<TPad*> padlow_;
    if(!ShowDMCRatio) {
      pad_  = PreparePads(1);
    }
    else {
      vector<vector<TPad*> > allpads = PreparePadsWithRatio(1);
      pad_  = allpads[0];
      padlow_ = allpads[1];
    }


    pad_[0]->Draw();
    pad_[0]->cd();
    if(logYScale)
      pad_[0]->SetLogy(1);

    leg=new TLegend(0.59,0.54,0.88,0.796);
    Plot1DHisto( obs, true, false, 0, 0 );
    
 
    leg->Draw("same");
    if(!NoData)
      cmsPrel(Lumi);

    pad_[0]->RedrawAxis();
    c2->cd();
    
    if(ShowDMCRatio) {
      padlow_[0]->Draw();
      padlow_[0]->cd();
      PlotDataMCRatio(obs, false);
      padlow_[0]->RedrawAxis();
      c2->cd();
    }

    c2->RedrawAxis();
    
    histoManager.SavePlots( path, c2);
    
    delete pad_[0];
    delete padlow_[0];
    //  delete c2;
   
  }
  
  defStyle=false;
  
}


void UseTree::PlotDistribution(string obs, bool Fill2DHistos,bool FillProfile ) {

  cout<<" Now plot!!! ratio ? "<<ShowDMCRatio<<endl;
  //Ok histo filled, find observable
  if(!Fill2DHistos && !FillProfile ) {
    
    if(!Draw3on1) {

      leg=new TLegend(0.59,0.54,0.88,0.796);
      Wline = 2;

      vector<TPad*> pad_;
      vector<TPad*> padlow_;
      if(!ShowDMCRatio) {
	pad_  = PreparePads(1);
      }
      else {
	vector<vector<TPad*> > allpads = PreparePadsWithRatio(1);
	pad_  = allpads[0];
	padlow_ = allpads[1];
      }
    
      int color=0,color2=0;
      if(SuperColor) {
      
	if(obs.substr(0,2)=="ca") {
	  color=400;   //kYellow;
	  color2=800;  //kOrange+7;
	}
	else if(obs.substr(0,2)=="pf") {
	  color=38;
	  color2=854;  //kAzure-6;
	}
	else if(obs.substr(0,2)=="tc") {
	  color=896;
	  color2=904;    //kPink+4;	
	}

      }
      pad_[0]->Draw();
      pad_[0]->cd();
      Plot1DHisto(obs, true, false,color,color2);
   
      leg->Draw("same");
      if(!NoData)
	cmsPrel(Lumi);

      if(logYScale)
	pad_[0]->SetLogy(1);

      pad_[0]->RedrawAxis();
      c2->cd();

      if(ShowDMCRatio) {
	padlow_[0]->Draw();
	padlow_[0]->cd();
	PlotDataMCRatio(obs, false);
	padlow_[0]->RedrawAxis();
	c2->cd();
      }

      c2->RedrawAxis();
    
    }
    else { //combine  Canvas
   
      leg=new TLegend(0.62,0.53,0.94,0.9);
      Wline = 1;

      vector<TPad*> pad_;
      vector<TPad*> padlow_;
      if(!ShowDMCRatio) {
	pad_  = PreparePads(3);
      }
      else {
	vector<vector<TPad*> > allpads = PreparePadsWithRatio(3);
	pad_  = allpads[0];
	padlow_ = allpads[1];
      }
    
    
      //Canvas 1************************************************************
      {
	pad_[2]->Draw();
	pad_[2]->cd();
	leg=new TLegend(0.666,0.579,0.939,0.912);


	int color=0,color2=0;
	if(SuperColor) {
	  color=38;
	  color2=kAzure-6;
	  //	leg=new TLegend(0.52,0.64,0.87,0.85);
	  //	leg=new TLegend(0.58,0.70,0.93,0.91);
	  leg=new TLegend(0.576,0.596,0.927,0.807);
	}

	Plot1DHisto(SevObs[2], true, true,color,color2);
	//Latest Pads********************************
	TPaveText* LumiPad;//=new TPaveText(0.56,0.75,0.92,0.89,"NDC");
	if( SevObs[2]=="pfRecoilPara" ){ // bottom-left
	  LumiPad = new TPaveText(0.417,0.868,0.949,0.919,"NDC");
	}
	else if( SevObs[2]=="pfRecoilPerp" ||
		 SevObs[2]=="pfdPhiRecoil" ){ // top-left
	  LumiPad = new TPaveText(0.2548779, 0.722797, 0.582312,
				  0.87,"NDC");//0.927762
	}
	else{ // bottom-right
	  LumiPad = new TPaveText(0.417,0.868,0.949,0.919,"NDC");
	  //	LumiPad = new TPaveText(0.56,0.50,0.92,0.64,"NDC");
	}
	LumiPad->SetFillColor(0);
	LumiPad->SetFillStyle(0);
	LumiPad->SetLineColor(0);
	LumiPad->SetTextSize(0.05);
	//      LumiPad->AddText("#sqrt{s} = 7 TeV      ");
	LumiPad->AddText("");
	LumiPad->AddText("");
	if(Lumi > 1000 )
	  LumiPad->AddText(Form("%.1f fb^{-1} at  #sqrt{s} = 8 TeV",Lumi/1000.));
	else
	  LumiPad->AddText(Form("%.0f pb^{-1} at  #sqrt{s} = 8 TeV",Lumi));
	if(!NoData)
	  LumiPad->Draw("same");
	if(SuperColor)
	  leg->Draw("same");

	if(logYScale)
	  pad_[2]->SetLogy(1);
	
	pad_[2]->RedrawAxis();
      
	if(ShowDMCRatio) {
	  c2->cd();
	  padlow_[2]->Draw();
	  padlow_[2]->cd();
	  PlotDataMCRatio(SevObs[2], true);
	  padlow_[2]->RedrawAxis();
	  c2->cd();
	}

	c2->cd(3)->Update();
	c2->cd(6)->Update();

      }

      //Canvas 2 **********************************************************************

      {
	leg=new TLegend(0.666,0.579,0.939,0.912);
	c2->cd();
	pad_[1]->Draw();
	pad_[1]->cd();

	int color=0,color2=0;
	if(SuperColor) {
	  color=896;
	  color2=kPink+4;	
	  //leg=new TLegend(0.52,0.64,0.87,0.85);
	  //      	leg=new TLegend(0.58,0.70,0.93,0.91);
	  leg=new TLegend(0.576,0.596,0.927,0.807);
	}
	Plot1DHisto(SevObs[1], true, true,color,color2);
	// if(SuperColor)
	leg->Draw("same");

	if(logYScale)
	  pad_[1]->SetLogy(1);


	pad_[1]->RedrawAxis();  

	if(ShowDMCRatio) {
	  c2->cd();
	  padlow_[1]->Draw();
	  padlow_[1]->cd();
	  PlotDataMCRatio(SevObs[1], true);
	  padlow_[1]->RedrawAxis();
	  c2->cd();
	}
	c2->cd(3)->Update();
	c2->cd(5)->Update();

      }
      //Canvas 3 *************************************************************

      {
	c2->cd();
	pad_[0]->Draw();
	pad_[0]->cd();
	leg=new TLegend(0.666,0.579,0.939,0.912);

	int color=0,color2=0;
	if(SuperColor) {
	  color=kYellow;
	  color2=kOrange-3;
	  //leg=new TLegend(0.52,0.64,0.87,0.85);
	  leg=new TLegend(0.576,0.596,0.927,0.807);
	  //      	leg=new TLegend(0.58,0.70,0.93,0.91);
	}

	Plot1DHisto(SevObs[0], true, false,color,color2);
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.055);
	latex.SetTextAlign(11); // align left
	if(!NoData)
	  latex.DrawLatex(0.437,0.96,"CMS preliminary 2012");
	if(SuperColor)
	  leg->Draw("same");

	if(logYScale)
	  pad_[0]->SetLogy(1);
	
	pad_[0]->RedrawAxis();

	if(ShowDMCRatio) {
	  c2->cd();
	  padlow_[0]->Draw();
	  padlow_[0]->cd();
	  PlotDataMCRatio(SevObs[0], false );
	  padlow_[0]->RedrawAxis();
	  c2->cd();
	}

	c2->cd(1)->Update();
	c2->cd(4)->Update();
      }

    }
  }  //2D Variable *****************************
  else if(Fill2DHistos) {
    ShowDMCRatio =false; //FIXME
    cout<<" Switching to 2D plots ******************************* "<<endl;
  
    if(!Draw3on1) {
      
      int color=38;
      if(obs.substr(0,2)=="ca")
	color = kOrange+7;
      else if(obs.substr(0,2)=="tc")
	color = 896;
      
      Wline = 2;
      leg2D=new TLegend(0.70,0.62,0.9,0.82);
      
      vector<TPad*> pad_ = PreparePads(1);
      
      pad_[0]->Draw();
      pad_[0]->cd();
      
      if(!switchRMS)
	Plot2DHisto(obs,true, false);
      else
	PlotRMS(obs,true, false,color);
      cmsPrel(Lumi);
      leg2D->Draw("same");
      c2->cd();
      c2->RedrawAxis();
    }

    //********************* 3 Canvas ************************************
    else {
      Wline = 1;
      leg2D=new TLegend(0.49,0.53,0.81,0.9);

      vector<TPad*> pad_ = PreparePads(3);
      
      //Canvas 1 *****************************************************************************
      { 

	int color=38; //38

	pad_[2]->Draw();
	pad_[2]->cd();
	if(!switchRMS)
	  Plot2DHisto(SevObs[2],false, true);
	else {
	  PlotRMS(SevObs[2],true, true,color);
	  leg2D->Draw("same");
	}
	//Latest Pads********************************
	TPaveText* LumiPad=new TPaveText(0.55,0.75,0.91,0.89,"NDC");
      
	LumiPad->SetFillColor(0);
	LumiPad->SetFillStyle(0);
	LumiPad->SetLineColor(0);
	LumiPad->SetTextSize(0.05);
	/*	LumiPad->AddText("#sqrt{s} = 7 TeV      ");
	LumiPad->AddText("");
	LumiPad->AddText("");*/
	if(Lumi > 1000 )
	  LumiPad->AddText(Form("%.1f fb^{-1} at  #sqrt{s} = 8 TeV",Lumi/1000.));
	else
	  LumiPad->AddText(Form("%.0f pb^{-1} at  #sqrt{s} = 8 TeV",Lumi));
	LumiPad->Draw("same");
	//	pad_[2]->RedrawAxis();

      }
    
      //Canvas 2 ****************************************************************
	
      {
	int color=896; //896

	c2->cd();
	pad_[1]->Draw();
	pad_[1]->cd();
 	if(!switchRMS)
	  Plot2DHisto(SevObs[1],true, true);
	else {
	  leg2D=new TLegend(0.49,0.53,0.81,0.9);
	  PlotRMS(SevObs[1],true, true,color);
	}
	leg2D->Draw("same");
	//	pad_[1]->RedrawAxis();
      }
     
      //Canvas 3 *******************************************************
      {
	int color=kOrange+7; //kOrange+7

	c2->cd();
	pad_[0]->Draw();
	pad_[0]->cd();
	if(!switchRMS)
	  Plot2DHisto(SevObs[0],false, false);
	else {
	  leg2D=new TLegend(0.49,0.53,0.81,0.9);
	  PlotRMS(SevObs[0],true, false,color);
	  leg2D->Draw("same");
	}
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.055);
	latex.SetTextAlign(11); // align left
	latex.DrawLatex(0.20,0.96,"CMS preliminary 2012");
	//	pad_[0]->RedrawAxis();

      }
      
    }//End multicanvas
  } //End 2D histos
  else if(FillProfile) {
    ShowDMCRatio =false; //FIXME
    if(!Draw3on1) {
      Wline = 2;
      leg2D=new TLegend(0.60,0.63,0.9,0.80);
      
      vector<TPad*> pad_ = PreparePads(1);

      int color=38;
      if(obs.substr(0,2)=="ca")
	color = kOrange+7;
      else if(obs.substr(0,2)=="tc")
	color = 896;

   
     pad_[0]->Draw();
     pad_[0]->cd();
     
      PlotProfile(obs,true, false,color);
      cmsPrel(Lumi);
      leg2D->Draw("same");

      c2->cd();
      c2->RedrawAxis();

    }
 //********************* 3 Canvas ************************************
    else {
      Wline = 1;
      leg2D=new TLegend(0.49,0.53,0.81,0.9);
         
      vector<TPad*> pad_ = PreparePads(3);
    
      //Canvas 1 *****************************************************************************
      { 
	c2->cd();
	pad_[2]->Draw();
	pad_[2]->cd();
	int color=38; //38
	leg2D=new TLegend(0.58,0.20,0.91,0.34);
	PlotProfile(SevObs[2],true, true,color);
	leg2D->Draw("same");
	//Latest Pads********************************
	TPaveText* LumiPad=new TPaveText(0.56,0.38,0.92,0.52,"NDC");

	LumiPad->SetFillColor(0);
	LumiPad->SetFillStyle(0);
	LumiPad->SetLineColor(0);
	LumiPad->SetTextSize(0.05);
	/*	LumiPad->AddText("#sqrt{s} = 7 TeV      ");
	LumiPad->AddText("");
	LumiPad->AddText("");*/
	if(Lumi > 1000 )
	  LumiPad->AddText(Form("%.1f fb^{-1} at  #sqrt{s} = 8 TeV",Lumi/1000.));
	else
	  LumiPad->AddText(Form("%.0f pb^{-1} at  #sqrt{s} = 8 TeV",Lumi));
	//	LumiPad->AddText(Form(" %.0f pb^{-1}  at #sqrt{s} = 8 TeV",Lumi));
	LumiPad->Draw("same");
	pad_[2]->RedrawAxis();

      }
    
      //Canvas 2 ****************************************************************
	
      {
	c2->cd();
	pad_[1]->Draw();
	pad_[1]->cd();
	int color=896; //896
	leg2D=new TLegend(0.58,0.20,0.91,0.34);
	PlotProfile(SevObs[1],true, true,color);
	leg2D->Draw("same");
	pad_[1]->RedrawAxis();
      }
     
      //Canvas 3 *******************************************************
      {
	c2->cd();
	pad_[0]->Draw();
	pad_[0]->cd();
	int color=kOrange+7; //kOrange+7
	leg2D=new TLegend(0.58,0.20,0.91,0.34);
	PlotProfile(SevObs[0],true, false,color);
	leg2D->Draw("same");
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.055);
	latex.SetTextAlign(11); // align left
	latex.DrawLatex(0.20,0.96,"CMS preliminary 2012");
	pad_[0]->RedrawAxis();

      }
      
    }//End multicanvas
  } //End Profile

} //End macro





void UseTree::Plot1DHisto(string obs, bool FillLeg, bool removelabel, int color, int color2) {

  vector<int> tmpcol(colors.size(), 0);

  if(defStyle) {
    Bin=5;
    BinBkg = 2;
    OverFlowBin = false;
    UnderFlowBin = false;
  
    YTitle = "number of events ";
    
  }


  if(SuperColor) {
    for(int i=0;i<nt-1;i++)
      tmpcol[i] = color2;
    tmpcol[nt-1] = color;
    tmpcol[nt] = colors[nt];
  }
  else
    for(int i=0;i<nt+1;i++)
      tmpcol[i]=colors[i];

  int Obs = histoManager.FindNVar(obs);
  
  if(Obs==-1) {
    cout<<" Be careful, no such observable "<<obs<<" , MET by default "<<endl;
    return;
  }
  else
    cout<<" Observable "<<obs<<" found "<<endl;
  vector<float> templates = histoManager.GetTemplate(Obs);
  cout<<" ajajaj "<<endl;
  XTitle = (histoManager.FindLeg( obs )).first;

  //  THStack* stackedObs = new THStack("stackedM",obs.c_str());
  TH1F* totalObs =new TH1F("totalObs","Observable",(int)templates[0],templates[1],templates[2]);
  TH1F* totalObsStep2 =new TH1F("totalObsS2","Observable",(int)templates[0],templates[1],templates[2]);
  
  TH1F* emptyHisto= new TH1F("bidon","bidon",(int)templates[0],templates[1],templates[2]);

  vector<TH1F*> Clones;
  vector<TH1F*> ClonesStack;

  leg->SetFillColor(0);
  leg->SetLineColor(1);
  
  map<string,float> weights;
  weights = Reweight(Obs, true);
  
  for(int i=0;i<nt;i++) {
  
    int Nhisto = histoManager.access(Obs, i);
    Clones.push_back( (TH1F*)Histos[Nhisto]->Clone() );
    Clones[i]->Rebin(Bin);

    cout<<name[i]<<" weight "<<weights[ name[i] ]<<" Nevt  "<<Clones[i]->Integral(0,Clones[i]->GetNbinsX()+1)<<endl;

    if(i==0) {
      totalObsStep2 = (TH1F*)Clones[i]->Clone();
    }
    else {
      totalObsStep2->Add(Clones[i]);
    }
    // stackedObs->Add( Clones[i],"ah" );
    
  }
 
  RooHist* roohist;
  TGraphAsymmErrors* DataGraph(0);
  TH1F* CloneData(0);
  
  int Ndata= histoManager.access(Obs,nt);
  cout<<" Observable number -> "<<Ndata<<endl;
  if(!NoData) {
    CloneData = (TH1F*)Histos[ Ndata ]->Clone();
    CloneData->Rebin(Bin);
  }
  cout<<" QCD ratio weights ---> "<<"   "<<QCDweight<<" / "<<normQCDWeight<<"  =  "<<QCDweight/normQCDWeight<<endl;
    
  float UncerWeight=1;

  if(Remove2V && !NoData) {
    TH1F* V2Cor = GetVtxCorrection(obs, nt,Bin);
    // cout<<V2Cor->Integral()<<endl;
    CloneData->Add(V2Cor,-1);
    cout<<" Integral after cleaning "<<CloneData->Integral(0, CloneData->GetNbinsX()+1)<<endl;
  }
  if(Remove2V) {
    int nn=0;
    if( reposi != "Wen")
      nn= nt;
    for(int i=nn;i<nt;i++) {
      TH1F* V2Cor = GetVtxCorrection(obs, i,Bin); //nt-i
      Clones[i]->Add(V2Cor,-1);
      for(int j=0;j<Clones[i]->GetNbinsX()+2;j++) {
	if(Clones[i]->GetBinContent(j)<0)
	  Clones[i]->SetBinContent(j,0);
      }
      cout<<" Integral after cleaning "<< Clones[i]->Integral(0, Clones[i]->GetNbinsX()+1)<<endl;
    }
  }
  float NEWK =0;
  float NQCD =0;
  float NData =0;
  if(reposi=="Wen" && !NoData) {
    NEWK = Clones[0]->Integral(0, Clones[4]->GetNbinsX()+1) + Clones[4]->Integral(0, Clones[4]->GetNbinsX()+1) + Clones[3]->Integral(0, Clones[3]->GetNbinsX()+1);
   NQCD = Clones[1]->Integral(0, Clones[1]->GetNbinsX()+1) + Clones[2]->Integral(0, Clones[2]->GetNbinsX()+1);
   NData = CloneData->Integral(0, CloneData->GetNbinsX()+1);
  }
  else {
    if(Remove2V || ShapeWeight) {
      FitR=false;

      if(vetoE)
	{
	  NData = CloneData->Integral(0, CloneData->GetNbinsX()+1); 
	  NEWK = Clones[nt-1]->Integral(0, Clones[nt-1]->GetNbinsX()+1);cout<<" NData****** "<<NData<<"  "<<NEWK<<endl;
	  Clones[nt-1]->Scale(NData/NEWK);
	}
      else
	{
	  /*NEWK = Clones[0]->Integral(0, Clones[0]->GetNbinsX()+1) + Clones[4]->Integral(0, Clones[4]->GetNbinsX()+1) + Clones[3]->Integral(0, Clones[3]->GetNbinsX()+1);
	    NQCD = Clones[2]->Integral(0, Clones[2]->GetNbinsX()+1) + Clones[1]->Integral(0, Clones[1]->GetNbinsX()+1);
	    NData = CloneData->Integral(0, CloneData->GetNbinsX()+1);
	    float R = NData/(NEWK+NQCD);
	  */
	  if(Norm) { //FIXME
	    NData = CloneData->Integral(0, CloneData->GetNbinsX()+1);
	    
	    float totint=0;
	    for(int i=0;i<nt;i++) {
	      totint +=Clones[i]->Integral(0, Clones[i]->GetNbinsX()+1);
	      cout<<i<<"   "<<totint<<endl;
	    }
	    
	    for(int i=0;i<nt;i++) {
	      float integ = Clones[i]->Integral(0, Clones[i]->GetNbinsX()+1);
	      cout<<" integ "<<integ<<"   "<<totint<<"   "<<NData<<"  "<<NData/totint *integ<<endl;
	      if(integ!= 0)
		Clones[i]->Scale(NData/totint);
	    }
	    UncerWeight = NData/totint;
	  }
	  /*
	    Clones[nt-1]->Scale(R*NEWK/Clones[4]->Integral(0, Clones[4]->GetNbinsX()+1)); cout<<R*NEWK<<endl;
	  Clones[nt-2]->Scale(R*NEWK/Clones[3]->Integral(0, Clones[3]->GetNbinsX()+1));
	  Clones[0]->Scale(R*NEWK/Clones[0]->Integral(0, Clones[0]->GetNbinsX()+1));
	  Clones[1]->Scale(R*NQCD/Clones[1]->Integral(0, Clones[1]->GetNbinsX()+1));
	  Clones[2]->Scale(R*NQCD/Clones[2]->Integral(0, Clones[2]->GetNbinsX()+1));*/
	}

    }
  }


  if(Norm && !ShapeWeight) {
    
    if(vetoE)
      {
	if(!NoData) {
	  NData = CloneData->Integral(0, CloneData->GetNbinsX()+1); 
	  NEWK = Clones[nt-1]->Integral(0, Clones[nt-1]->GetNbinsX()+1);cout<<" NData****** "<<NData<<"  "<<NEWK<<endl;
	  Clones[nt-1]->Scale(NData/NEWK);
	}
      }
    else {
      
      if(NoData)
	NData = 1; //to be checked
      else
	NData = CloneData->Integral(0, CloneData->GetNbinsX()+1);
      
      float totint=0;
      for(int i=0;i<nt;i++) {
	totint +=Clones[i]->Integral(0, Clones[i]->GetNbinsX()+1);
	cout<<i<<"   "<<totint<<endl;
      }

      for(int i=0;i<nt;i++) {
	float integ = Clones[i]->Integral(0, Clones[i]->GetNbinsX()+1);
	cout<<" integ "<<integ<<"   "<<totint<<"   "<<NData<<"  "<<NData/totint *integ<<endl;
	if(integ!= 0)
	  Clones[i]->Scale(NData/totint);
      }
    }

    UncerWeight = NData/Clones[nt-1]->Integral(0, Clones[nt-1]->GetNbinsX()+1);

  }
 
  if(FitR && !MultiWeight && !NoData) {
    
    //Get MT distribution and determine a************************
    TString obsMT;
    if(obs.substr(0,4)=="pfT1")
      obsMT="pfT1";
    else if(obs.substr(0,2)=="tc")
      obsMT="tc";
    else if(obs.substr(0,6)=="caloT1")
      obsMT="caloT1";
    else if(obs.substr(0,6)=="caloT2")
      obsMT="caloT2";
    else if(obs.substr(0,2)=="ca")
      obsMT="calo";
    else
      obsMT="pf";
    vector<float> a = GetAlpha(obsMT+"MET");
    //***********************************************************

    
    float NormEWK2 = a[0];
    float NormQCD2 = a[1];


    float R = NData/(NormEWK2*NEWK+NormQCD2*NQCD);
    
    for(int i=0;i<nt;i++) {
      if(Clones[i]->Integral()!=0) {

      if(i>2 || i==0)
	Clones[i]->Scale( NormEWK2*R );
      else 
	Clones[i]->Scale( NormQCD2*R );
      UncerWeight =  NormEWK2*R ;
      float integ = Clones[i]->Integral(0, Clones[i]->GetNbinsX()+1 ); 
      cout<<" ---> After reweight "<<name[i] << " -> "<<integ<<endl;
      }
    }
  }
  else if(MultiWeight &&  !NoData) {
    float R = NData/(NEWK+NQCD);
    for(int i=0;i<nt;i++) {
      if(Clones[i]->Integral()!=0) {
	if(i>1)
	  Clones[i]->Scale( R );
	else 
	  Clones[i]->Scale( R );

	UncerWeight = R ;
	float integ = Clones[i]->Integral(0, Clones[i]->GetNbinsX()+1 ); 
	cout<<" ---> After reweight "<<name[i] << " -> "<<integ<<endl;
      }
    }
  }
    

  float tot2=0;
 
  for(int i=0;i<nt;i++) { 
    
    ClonesStack.push_back( (TH1F*)Clones[i]->Clone() );
    
    Clones[i]->SetFillColor(tmpcol[i]);
    ClonesStack[i]->SetFillColor(tmpcol[i]);
   
    if(i==nt-1) {
      if(!SuperColor) {
	Clones[i]->SetLineColor(tmpcol[i]);
	Clones[i]->SetLineWidth(LineWidth);
	Clones[i]->SetMarkerStyle(0);
	ClonesStack[i]->SetLineColor(tmpcol[i]);
	ClonesStack[i]->SetLineWidth(LineWidth);
	ClonesStack[i]->SetMarkerStyle(0);
      }
      else  {
	Clones[i]->SetLineColor(color);	
	ClonesStack[i]->SetLineColor(color);	
      }
    }
    else
      if(!SuperColor) {
	Clones[i]->SetLineColor(tmpcol[i]);
	ClonesStack[i]->SetLineColor(tmpcol[i]);
      }
      else {
	Clones[i]->SetLineColor(color2);
    	ClonesStack[i]->SetLineColor(color2);
      }
 
    Clones[i]->SetLineWidth(LineWidth);
    Clones[i]->SetMarkerColor(tmpcol[i]);
    Clones[i]->GetYaxis()->SetLabelSize(0);
    Clones[i]->GetYaxis()->SetTitleSize(0);
    Clones[i]->GetYaxis()->SetLabelOffset(1000);
    Clones[i]->GetYaxis()->SetTitleOffset(1000);

    if(!defStyle) {

      Clones[i]->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
      if(!OverFlowBin && !UnderFlowBin)
	Clones[i]->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1] );
      else {
	if(OverFlowBin)
	  Clones[i]->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1] + Clones[i]->GetBinWidth(1)  -0.001);
	if(UnderFlowBin)
	  Clones[i]->GetXaxis()->SetRangeUser(RangeX[0] - Clones[i]->GetBinWidth(1) ,RangeX[1] + Clones[i]->GetBinWidth(1) -0.001 );
      }

    }
    else {
      if(!logYScale)
	Clones[i]->GetYaxis()->SetRangeUser(0.1, Clones[i]->GetBinContent(Clones[i]->GetMaximumBin()) );
      else
	Clones[i]->GetYaxis()->SetRangeUser(0.5, Clones[i]->GetBinContent(Clones[i]->GetMaximumBin())*2 );
    }

    Clones[i]->GetYaxis()->SetNdivisions(0,0,0);
    Clones[i]->GetXaxis()->SetNdivisions(0,0,0);
 
    ClonesStack[i]->SetLineWidth(LineWidth);
    ClonesStack[i]->SetMarkerColor(tmpcol[i]);
    ClonesStack[i]->GetYaxis()->SetLabelSize(0);
    ClonesStack[i]->GetYaxis()->SetTitleSize(0);
    ClonesStack[i]->GetYaxis()->SetLabelOffset(1000);
    ClonesStack[i]->GetYaxis()->SetTitleOffset(1000);

    if(!defStyle) {
      ClonesStack[i]->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
      if(!OverFlowBin && !UnderFlowBin)
	ClonesStack[i]->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1] );
      else {
	if(OverFlowBin)
	  ClonesStack[i]->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1] + ClonesStack[i]->GetBinWidth(1) -0.001);
	if(UnderFlowBin)
	  ClonesStack[i]->GetXaxis()->SetRangeUser(RangeX[0] - ClonesStack[i]->GetBinWidth(1) ,RangeX[1] + ClonesStack[i]->GetBinWidth(1) -0.001);
      }
    }
    else {
      if(!logYScale)
	Clones[i]->GetYaxis()->SetRangeUser(0.1, Clones[i]->GetBinContent(Clones[i]->GetMaximumBin()) );
      else
	Clones[i]->GetYaxis()->SetRangeUser(0.5, Clones[i]->GetBinContent(Clones[i]->GetMaximumBin())*2 );
    }

    ClonesStack[i]->GetYaxis()->SetNdivisions(0,0,0);
    ClonesStack[i]->GetXaxis()->SetNdivisions(0,0,0);
  
    //background stacks
    if(i<nt-1)
      for(int j=0;j<i;j++) {
    	ClonesStack[i]->Add(Clones[j],1);
      }
    tot2 += Clones[i]->Integral(0,Clones[i]->GetNbinsX()+1);
  }  

  //Additionnal binning for background
  if(nt!=0) {

    if(BinBkg!=1) {
      for(int i=0;i<nt-1;i++) {
	ClonesStack[i]->Rebin(BinBkg);
	ClonesStack[i]->Scale( 1./BinBkg);
      }
    }

    //Last Clone
    for(int ib=0;ib<ClonesStack[nt-1]->GetNbinsX()+2;ib++) {
      float center = ClonesStack[nt-1]->GetBinCenter(ib);
      float val = ClonesStack[nt-1]->GetBinContent(ib);
      //and the background
      float addval = 0;
      if(nt>1)
	addval = ClonesStack[nt-2]->GetBinContent( ClonesStack[nt-2]->FindBin(center) ) ;
    
    
      ClonesStack[nt-1]->SetBinContent(ib, val + addval );
    }
    cout<<" end stacking "<<endl;
    if( !NoData ) {
      cout<<" Checking the number of events before overflowing : "<<CloneData->Integral(0,1000)<<" <-> "<<ClonesStack[nt-1]->Integral(0,1000)<<endl;
    }   

  }

 

  //Overflow bin
  vector<TH1F*> ovflowHists;
  vector<TH1F*> unflowHists;

  if(OverFlowBin) {
    float back=0;
    for(int i=0;i<nt;i++) {
    
    //prepare ovf histos;
    ovflowHists.push_back( (TH1F*)Clones[i]->Clone() );
    ovflowHists[i]->Reset();

    int ovbin = Clones[i]->FindBin(RangeX[1]) ;
    float ovinteg = Clones[i]->Integral(ovbin,2000);
    if(i==nt-1)
      Clones[i]->SetBinContent(ovbin, ovinteg);
    else
      Clones[i]->SetBinContent(ovbin, ovinteg*BinBkg);
    int ovbinstack = ClonesStack[i]->FindBin(RangeX[1]) ;
    float ovintegstack = ClonesStack[i]->Integral(ovbinstack,2000);
    //  cout<<i<<" before ovbin "<<ovbin<<"   "<<ovinteg<<" stack  "<<ovbinstack<<"  "<<ovintegstack<<"    "<< Clones[i]->GetBinContent(2000)<<endl;
  
    if(i<nt-1 ) {
      if(BinBkg!=1)
	ovintegstack *= (double)BinBkg;
    }
    else
      ovintegstack = back + ovinteg;
    
    ClonesStack[i]->SetBinContent(ovbinstack , ovintegstack);
    if(i==0 && !NoData) //do data
      {
	ovinteg = CloneData->Integral(ovbin,2000);
	CloneData->SetBinContent(ovbin, ovinteg);
      }

     ovflowHists[i]->SetBinContent(ovbin, ovintegstack);

     //     cout<<i<<" ovbin "<<ovbin<<"   "<<ovinteg<<" stack  "<<ovbinstack<<"  "<<ovintegstack<<"    "<< Clones[i]->GetBinContent(2000)<<endl;

    if(i==nt-2)
      back = ovintegstack;

    }  
  }

  if(UnderFlowBin) {
    float back=0;
    for(int i=0;i<nt;i++) {
      
    //prepare ovf histos;
    unflowHists.push_back( (TH1F*)Clones[i]->Clone() );
    unflowHists[i]->Reset();

    int ovbin = Clones[i]->FindBin(RangeX[0]) -1 ;
    float ovinteg = Clones[i]->Integral(-1,ovbin);
    if(i==nt-1)
      Clones[i]->SetBinContent(ovbin, ovinteg);
    else
      Clones[i]->SetBinContent(ovbin, ovinteg*BinBkg);
    int ovbinstack = ClonesStack[i]->FindBin(RangeX[0]) -1 ;
    float ovintegstack = ClonesStack[i]->Integral(-1,ovbinstack);
    cout<<i<<" before unbin "<<ovbin<<"   "<<ovinteg<<" stack  "<<ovbinstack<<"  "<<ovintegstack<<"    "<< Clones[i]->GetBinContent(2000)<<endl;
  
    if(i<nt-1 ) {
      if(BinBkg!=1)
	ovintegstack *= (double)BinBkg;
    }
    else
      ovintegstack = back + ovinteg;
    
    ClonesStack[i]->SetBinContent(ovbinstack , ovintegstack);
    if(i==0 && !NoData) //do data
      {
	ovinteg = CloneData->Integral(-1,ovbin);
	CloneData->SetBinContent(ovbin, ovinteg);
      }

     unflowHists[i]->SetBinContent(ovbin, ovintegstack);

     cout<<i<<" unbin "<<ovbin<<"   "<<ovinteg<<" stack  "<<ovbinstack<<"  "<<ovintegstack<<"    "<< Clones[i]->GetBinContent(2000)<<endl;

    if(i==nt-2)
      back = ovintegstack;

    }  
  }

  //Now total obs, cross check of stacked clones
  if(nt!=0) {
    totalObs->Rebin(Bin);
    for(int ib=0;ib<ClonesStack[nt-1]->GetNbinsX()+2;ib++) {
      float val = ClonesStack[nt-1]->GetBinContent(ib);
      //float center = totalObs->GetBinCenter(ib);
    
      /*    for(int i=0;i<nt-1;i++) {
	    TH1F* tmp = (TH1F*)Clones[i]->Clone();
	    if(BinBkg!=1) {
	    tmp->Rebin(BinBkg);
	    tmp->Scale( 1./BinBkg);
	    }
	    val += tmp->GetBinContent( tmp->FindBin(center) ) ;
	    }*/
        
      totalObs->SetBinContent(ib, val);
    }

    //plotting
    totalObs->SetFillColor(0);
    totalObs->SetLineColor(tmpcol[nt]);
    totalObs->SetLineWidth(LineWidth);
  }
  
  savedMC = (TH1F*)totalObs->Clone();

  if(removelabel) 
    emptyHisto->GetYaxis()->SetLabelOffset(1000);

  if(!NoData) {
 
    if(FillLeg)
      leg->AddEntry(Histos[ Ndata ],"data","pl");
    
    TFile* out(0);
    if(saveFile) {
      out= new TFile((obs+".root").c_str(),"RECREATE");
      CloneData->SetName("data");
      CloneData->Write();
      if(nt!=0) {
	
	string tts = name[nt-1];
	if(tts.size() > 16)
	  tts = name[nt-1].substr(0,16);

	Clones[nt-1]->SetName( tts.c_str() );
	Clones[nt-1]->Write();
	//Clones[1]->SetName("datauncor");
	if(nt>1) {
	  Clones[nt-2]->SetName( "otherproc" );
	  Clones[nt-2]->Write();
	  for(int i=0;i<nt-1;i++) {
	    tts = name[i];
	    if(tts.size() > 16)
	      tts = name[i].substr(0,16);
	    // cout<<" Saving plot "<< tts.c_str()<<endl;
	    Clones[i]->SetName( tts.c_str() );
	    Clones[i]->Write();
	  }
	
	}

	//Clones[1]->Write();
	totalObs->SetName("totalObs");
	totalObs->Write();
      }
    }
    roohist =new RooHist((*CloneData));
    cout<<endl<<" ==================== "<<endl;
    cout<<" RMS of data : "<<CloneData->GetRMS()<<" +- "<<CloneData->GetRMSError()<<endl<<endl;
    cout<<" Mean of data : "<<CloneData->GetMean()<<" +- "<<CloneData->GetMeanError()<<endl;
    cout<<" ==================== "<<endl<<endl;
    if(saveFile) {
      out->Close();
    }
    float DatNorm=1;
    if(NormPlot)
      DatNorm = 1./CloneData->Integral(0,1000);
    
    //    cout<<" DatNorm "<<DatNorm<<"    "<<CloneData->Integral(0,1000)<<endl;
    if(CloneData->GetEntries()!=0) {
      int Nn0=0;
      vector<double> vY;
      vector<double> vX;
      vector<double > veY;
      vector<double > veX;
      vector<double> tmp(0,2);
  
      for(int ip=0;ip<roohist->GetN();ip++) {
	double Y,X;
      
	roohist->GetPoint(ip,X,Y);
    
	if(Y!=0) {
	  Nn0++;
      
	  vY.push_back(Y);
	  vX.push_back(X);
	  veX.push_back( roohist->GetErrorXlow(ip) );
	  veX.push_back( roohist->GetErrorXhigh(ip) );
	  veY.push_back( roohist->GetErrorYlow(ip)*DatNorm );
	  veY.push_back( roohist->GetErrorYhigh(ip)*DatNorm );
	}
      }
      DataGraph=new TGraphAsymmErrors(Nn0);
      for(int ip=0;ip<Nn0;ip++) {
	DataGraph->SetPoint(ip,vX[ip],vY[ip]*DatNorm);
	DataGraph->SetPointError(ip,0,0/*veX[ip*2],veX[ip*2+1]*/,veY[ip*2],veY[ip*2+1]);
	//	cout<<" Errors "<<ip<<"  "<<vX[ip]<<"  --->  "<<veY[ip*2]<<"   "<<veY[ip*2+1]<<endl;
      }
      if(RangeX[1]> 110)
	DataGraph->SetMarkerSize(MarkerSize);//0.9
      else
	DataGraph->SetMarkerSize(MarkerSize); 
      DataGraph->SetMarkerStyle(20);
      DataGraph->SetMarkerColor(1);
      DataGraph->SetLineColor(1);
      DataGraph->SetLineWidth(LineWidth);
  
    }
    if(nt!=0) {
    totalObs->GetYaxis()->SetTitleOffset(1000);
    totalObs->GetYaxis()->SetLabelOffset(1000);
    totalObs->GetYaxis()->SetLabelSize(0);
    totalObs->GetYaxis()->SetTitleSize(0);
    }
    if(CloneData->GetEntries()!=0) {
      DataGraph->GetYaxis()->SetLabelSize(0);
      DataGraph->GetYaxis()->SetTitleSize(0);
    }
    if(removelabel) {
      if(CloneData->GetEntries()!=0)
	DataGraph->GetYaxis()->SetLabelOffset(1000);
      emptyHisto->GetYaxis()->SetLabelOffset(1000);
    }
    
    //Cleaning Removing 2V
    if(Remove2V) {
      double x,y;
      for(int k=0;k<DataGraph->GetN();k++)
	{
	
	  DataGraph->GetPoint(k,x,y);
	  //cout<<" x "<<x<<"  y "<<y<<endl;
	  if(y <= RangeY[0])  {
	    DataGraph->SetPointError(k,0,0,0,0);
	    //   DataGraph->RemovePoint(k);
	  }
	}
    }
    
  }

  if(!defStyle) {
    
    if(RangeX[1]>=1000) {
      cout<<"------------------------->   "<<XTitle.find("GeV")<<endl;
	if(XTitle.find("GeV") != (size_t)-1 ) {
	
	XTitle.replace( XTitle.find("GeV") , 3, "TeV" );
	cout<<XTitle<<endl;
      }
    }
    emptyHisto->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
    if(!OverFlowBin && !UnderFlowBin) {
      emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1] );
    }
    else {
      if(OverFlowBin)
	emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1] + Clones[0]->GetBinWidth(1) -0.001);
      if(UnderFlowBin)
	emptyHisto->GetXaxis()->SetRangeUser(RangeX[0] - Clones[0]->GetBinWidth(1) ,RangeX[1] + Clones[0]->GetBinWidth(1) -0.001);
    
      cout<<" limits "<<RangeX[0] - Clones[0]->GetBinWidth(1) <<"    "<<RangeX[1] + Clones[0]->GetBinWidth(1)<<endl;
    }

  }
  else {

    //MAX value along Y axis
    float YMax = ClonesStack[nt-1]->GetBinContent(ClonesStack[nt-1]->GetMaximumBin()) *1.1;
    if( !NoData )
      if( CloneData->GetBinContent(CloneData->GetMaximumBin()) *1.1 > YMax )
	YMax = CloneData->GetBinContent(CloneData->GetMaximumBin())*1.1;
  
    if(logYScale)
      YMax *=2;

    if(NormPlot) {
      YMax=0;
      float tmp;
      for(int i=0;i<nt;i++) {
	tmp = Clones[i]->GetBinContent(Clones[i]->GetMaximumBin())/Clones[i]->Integral(0,1000) *1.1;
	if(tmp>YMax)
	  YMax = tmp;

	if(logYScale)
	  YMax *=2;
      }
    }
    
    emptyHisto->GetYaxis()->SetRangeUser( logYScale?0.9:0.1, YMax );
    emptyHisto->GetXaxis()->SetRangeUser(templates[1],templates[2] );

    if(templates[2]>=1000) {
      if(XTitle.find("GeV") != (size_t)-1 ) {

	XTitle.replace( XTitle.find("GeV") , 3, "TeV" );
      }
    }

  }
  
  /// configure the binning of the Y title =======
  float bin = Clones[nt-1]->GetXaxis()->GetBinWidth(1);
  ostringstream os;
  os<<bin;
  string addbin=" / "+(os.str()).substr(0,3)+" GeV";
  ///=============================================


  emptyHisto->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
  emptyHisto->GetYaxis()->SetNdivisions(Ydiv[0],Ydiv[1],Ydiv[2]);
  emptyHisto->GetXaxis()->SetTitle(XTitle.c_str());
  emptyHisto->GetYaxis()->SetTitle((YTitle+addbin).c_str());
  
  if(RangeY[1] >= 1000 && !removelabel ) {
    emptyHisto->GetYaxis()->SetTitleOffset(1.5);
    if( !NoData )
      CloneData->GetYaxis()->SetTitleOffset(1.5);
    //   stackedObs->GetYaxis()->SetTitleOffset(1.6);
    if(nt!=0) {
      totalObs->GetYaxis()->SetTitleOffset(1.5);
      for(int i=0;i<nt;i++) {
	Clones[i]->GetYaxis()->SetTitleOffset(1.5);
	ClonesStack[i]->GetYaxis()->SetTitleOffset(1.5);
      }
    }
  }

  string legopt="f";
  if(NormPlot)
    legopt = "l";

  if(FillLeg && nt!=0) {
    for(int i=1;i<nt+1;i++) {
      if(!SuperColor)
	leg->AddEntry(Clones[nt-i],name[nt-i].c_str(),legopt.c_str());
    }
    //    if(SuperColor) {
      // leg->AddEntry(Clones[nt-1],name[nt-1].c_str(),legopt.c_str());
      // leg->AddEntry(Clones[0],"background",legopt.c_str());
      if(AddSyst) {
	TH1F* unctmp=new TH1F("unctmp","",2,0,1);
	unctmp->SetMarkerStyle(1);
	unctmp->SetMarkerColor(1);
	unctmp->SetLineColor(1);
	unctmp->SetFillStyle(3001);
	unctmp->SetFillColor(1);
	leg->AddEntry(unctmp,"uncertainties",legopt.c_str());
      }
      // }
  }

  if(!NoData && CloneData->GetEntries()!=0) {
    double x,y;
   
    for(int k=0;k<DataGraph->GetN();k++)
      {
	DataGraph->GetPoint(k,x,y);
	//	cout<<x<<"   "<<y<<"  "<<RangeY[0]<<endl;
	if(y <= RangeY[0]  ) {
	  DataGraph->SetPointError(k,0,0,0,0);
	  // DataGraph->RemovePoint(knt);
	  //	  cout<<"removing "<<endl;
	}
      }
  }
 
  if(data.size() != 0 && !NoData && CloneData->GetEntries()!=0) {
    emptyHisto->Draw();
    if(!NormPlot) {

      for(int i=nt-1;i>-1;i--) {
	ClonesStack[i]->Draw("same hist");
      }

      // stackedObs->Draw("ahsame hist"); 
   
      //OverFlow bin
      if(OverFlowBin){
	for(int i=nt-1;i>-1;i--)
	  ovflowHists[i]->Draw("same hist");
      }      
      if(UnderFlowBin){
	for(int i=nt-1;i>-1;i--)
	  unflowHists[i]->Draw("same hist");
      }    
      if(nt!=0) {
	totalObs->Draw("same hist");
      }
    if(AddSyst/*  && reposi!="OfficialZee"*/)
      AddSystematics(obs, totalObs,1./*UncerWeight*/);
    }
    else {
      for(int i=0;i<nt;i++) {
	Clones[i]->SetFillColor(0);
	Clones[i]->Scale(1./Clones[i]->Integral(0,1000) );
	Clones[i]->Draw("same hist");
      }
    }
    if(CloneData->GetEntries()!=0 && DataGraph->GetN()!=0) {
      DataGraph->Draw("PZ"); 
      DataGraph->SetName("data");
    }
  }
  else {
    emptyHisto->Draw();
    if(!NormPlot) {
      for(int i=nt-1;i>-1;i--) {
	ClonesStack[i]->Draw("same hist");
      }
      // stackedObs->Draw("AHsame hist");
     
      totalObs->Draw("same hist");
   
      if(AddSyst/*  && reposi!="OfficialZee"*/)
	AddSystematics(obs, totalObs,1./*UncerWeight*/);
   
    }
    else {
      for(int i=0;i<nt;i++) {
	Clones[i]->SetFillColor(0);
	Clones[i]->Scale(1./Clones[i]->Integral(0,1000) );
	Clones[i]->Draw("same hist");
      }
    }
    
   
  }

  
  //Blanco
  //	TPave* pv = new TPave(1000.94,0.25,1000.97,0.31);
  TPave* pv = new TPave(1016.06,-0.8,1070.2,120.7);
  //TPave* pv = new TPave(0.9722,0.403579,1.00676,0.527054,4,"NDC");
  pv->SetFillColor(0);
  pv->SetShadowColor(0);
  pv->SetLineColor(0);
  pv->Draw();


  if(EcalP!=2)
    for(size_t i=0;i<Lines.size();i++)
	DrawLine(obs,Lines[i],EcalP,RangeY[1]);

  if(!NoData && CloneData->GetEntries()!=0 && nt!=0) {
    cout<<" Checking the number of events before reweight : "<<CloneData->Integral(0,1000)<<" <-> "<<totalObsStep2->Integral(0,1000)<<endl;
    cout<<" Checking the number of events : "<<CloneData->Integral(0,1000)<<" <-> "<<totalObs->Integral(0,1000)<<endl;

    if(OverFlowBin)
      cout<<" Checking in reduced range : "<<CloneData->Integral(0, CloneData->FindBin(RangeX[1] ) )
	  <<" <-> "<<totalObs->Integral(0, totalObs->FindBin(RangeX[1] ) )<<endl;

    /*  cout<<"***** Data Statistics numbers "<<endl;
    cout<<" Mean ----> "<<CloneData->GetMean()<<endl;
    cout<<" RMS ----> "<<CloneData->GetRMS()<<endl;
    
    cout<<"***** Kolmogorov Test "<<endl;
    KolmogorovTest(totalObs,CloneData,obs);
    // Chi2Test(totalObs,CloneData,observable);
    
    cout<<" ****** 2e passage" <<endl;
    vector<TH1F*> vt = HistoManager::ReduceHistoToNoEmptyBin(CloneData,totalObs);
    KolmogorovTest(vt[1],vt[0],obs);*/
  }
  
  //Rapport S/B
  float s=0,b=0;
  for(int i=0;i<nt;i++) {

    // cout<<" name "<<name[i]<<"  "<< Clones[i]->Integral(0,1000)<<endl;

    if(i==0)
      s = Clones[i]->Integral(0,1000);
    else
      b += Clones[i]->Integral(0,1000); 
  }
  cout<<" Rapport S/sqrt(B) ----> "<<s/sqrt(b)<<"       ("<<s<<"/"<<b<<")"<<endl;


}


void UseTree::Plot2DHisto(string obs, bool fillleg, bool removelabel) {

  int Obs = histoManager.FindNVar2D(obs);
      
      if(Obs==-1) {
	cout<<" Be careful, no such observable "<<obs<<endl;
	return;
      }

      map<string,float> weights;
      weights = Reweight(Obs, false);
     
      leg2D->SetFillColor(0);
      leg2D->SetLineColor(1);

      int Ndata= histoManager.access2D(Obs,nt);
      cout<<obs<<"    "<<Obs<<"   "<<Ndata<<endl;
      // TGraph* dataPointWhite=  (TGraph*)Histos2D[Ndata]->Clone();
      vector<TH2F*> cloneH;
      TH2F* cloneData;
      //  vector<TGraph*> graphs;
      XTitle = (histoManager.FindLeg2D( obs )).first;
      if(!Draw3on1)
	YTitle = (histoManager.FindLeg2D( obs )).second;

      vector<float> temp = histoManager.GetTemplate2D(Obs);
      TH2F* emptyHisto= new TH2F(("bidon"+obs).c_str(),"bidon",(int)temp[0],temp[1],temp[2],(int)temp[3],temp[4],temp[5]);
      cout<<temp[0]<<"   "<<temp[1]<<"  "<<RangeY[0]<<endl;
      emptyHisto->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
      emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);
      emptyHisto->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
      emptyHisto->GetYaxis()->SetNdivisions(Ydiv[0],Ydiv[1],Ydiv[2]);
      emptyHisto->GetXaxis()->SetTitle(XTitle.c_str());
      emptyHisto->GetYaxis()->SetTitle(YTitle.c_str());
 
      emptyHisto->SetFillColor(0);
      emptyHisto->SetLineColor(0);
      emptyHisto->Fill(0.,0.);
      //Ellipse

      TMarker* mark = new TMarker();
      mark->SetMarkerStyle(4);
      mark->SetMarkerColor(kGreen+2);
      cout<<" NObs "<<Ndata<<"   "<<Histos2D.size()<<endl;
      
      cloneData = (TH2F*)Histos2D[Ndata]->Clone();
      cloneData->SetMarkerColor(kGreen+2);
      cloneData->SetLineColor(kWhite);
      cloneData->SetMarkerStyle(20);
      cloneData->SetMarkerSize(MarkerSize);
      cloneData->SetFillColor(kWhite);
      if(fillleg) {
	leg2D->AddEntry(cloneData,"data","p"); }

      for(int i=0;i<nt;i++) {
	int Nhisto = histoManager.access2D(Obs,nt-1-i);
	TH2F* tmp = (TH2F*)Histos2D[ Nhisto ]->Clone();

	cloneH.push_back(tmp);
	cloneH[i]->SetMarkerColor(colors[nt-i]);
	cloneH[i]->SetMarkerStyle(8);
	cloneH[i]->SetFillColor(colors[nt-i]);
	cloneH[i]->SetLineColor(colors[nt-i]);
	cloneH[i]->SetMarkerSize(MarkerSize);

      }  
   
      for(int i=0;i<nt;i++) {
	//Mise en forme
	cloneH[i]->SetMarkerColor(colors[i]);
	cloneH[i]->SetFillColor(colors[i]);
	cloneH[i]->SetLineColor(colors[i]);
	cloneH[i]->SetMarkerStyle(8);
	cloneH[i]->SetMarkerSize(MarkerSize);
      }
   
      //Mise en forme donnÈes

	if(fillleg) {
	for(int i=1;i<nt+1;i++) {
	leg2D->AddEntry(cloneH[nt-i],name[nt-i].c_str(),"f"); }
      }
  

      /*  dataPointWhite->SetMarkerStyle(20);
	  dataPointWhite->SetMarkerSize(0.4);
	  dataPointWhite->SetMarkerColor(kWhite);*/
    
      
      // vector<TH2F*> Histos2d =  histoManager.ConvertAll(graphs,obs,name,XTitle);
      vector <TH2F*> Ordered;
      for(int i=0;i<nt;i++) {
	for(int j=0;j<nt;j++) {
	  if(Sorder[i]==name[j]) {
	    cout<<Sorder[i]<<"   "<<name[j]<<endl;
	    Ordered.push_back(cloneH[i]);
	  }
	}
      }
      /* float a=0.5;
	 if(FitR) {
	
	 //Get MT distribution and determine a************************
	 vector<TH1F*> ClonesMT;
	 string obsMT;
	 if(observable.substr(0,2)=="pf")
	 obsMT="pfMET";
	 else if(observable.substr(0,2)=="tc")
	 obsMT="tcMET";
	 else if(observable.substr(0,6)=="caloT1")
	 obsMT="caloT1MET";
	 else if(observable.substr(0,6)=="caloT2")
	 obsMT="caloT2MET";
	 else
	 obsMT="caloMET";
	
	 cout<<" test "<<obsMT<<endl;
	int ObsMT = histoManager.FindNVar(obsMT);
	cout<<" ob "<<Obs<<"  "<<ObsMT<<endl;
	//  weights = Reweight(ObsMT, true);
	for(int i=0;i<nt;i++) {
	  int NhistoMT = histoManager.access(ObsMT, i);
	  cout<<"  greoureo "<<NhistoMT<<endl;
	  ClonesMT.push_back( (TH1F*)Histos[NhistoMT]->Clone() );
	  cout<<" coin "<<endl;
	  ClonesMT[i]->Rebin(2);
	}
	int NdataMT= histoManager.access(ObsMT,nt);
	TH1F* CloneDataMT = (TH1F*)Histos[ NdataMT ]->Clone();
	CloneDataMT->Rebin(2);
	a= FitReweight(ClonesMT, CloneDataMT);
	cout<<" Best a "<<a<<" ; 0.585 as reference  -> r="<<a/0.585<<endl;
	//=============================================================
	
	float NormQCD = Ordered[3]->Integral() + Ordered[2]->Integral();
	float NormEWK = Ordered[0]->Integral() + Ordered[1]->Integral();
	float NormData = dataPointWhite->GetN();
	for(int i=0;i<4;i++) {
	  cout<<" test "<<weights[ name[i] ]<<endl;
	  if(i<2) //bck
	    weights[ name[i] ] *= (1-a)*NormData/NormQCD;
	  else //signal
	    weights[ name[i] ] *= a*NormData/NormEWK;
	  cout<<" testend "<<weights[ name[i] ]<<endl;
	}
      }
      */
    
      vector<TH2F*> Stacked = histoManager.Stack2DHistos(Ordered,Sorder,weights);
    
      if(removelabel) {
	cout<<" remove label!! "<<obs<<endl;
	for(int i=0;i<nt;i++) {
	  Stacked[i]->GetYaxis()->SetLabelOffset(1000);
	  Stacked[i]->GetYaxis()->SetTitleOffset(1000);
	}
	emptyHisto->GetYaxis()->SetLabelOffset(1000);
      }
      // else {
   
      for(int i=0;i<nt;i++) {
	Stacked[i]->Rebin2D(Bin,Bin);
	Stacked[i]->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
	Stacked[i]->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);
	Stacked[i]->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
	Stacked[i]->GetYaxis()->SetNdivisions(Ydiv[0],Ydiv[1],Ydiv[2]);
	Stacked[i]->GetXaxis()->SetTitle(XTitle.c_str());
	Stacked[i]->GetYaxis()->SetTitle(YTitle.c_str());
      }
      //  }
    
      emptyHisto->Draw("box");
      emptyHisto->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
      emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);
 
      for(int i=0;i<nt;i++) {

	Stacked[i]->SetFillColor(colors[nt-1-i]);
	if(i==0)
	  Stacked[i]->Draw("box");
	else
	  Stacked[i]->Draw("boxsame");
      }
      Histos2D[Ndata]->Draw("boxsame");
      // dataPointWhite->Draw("Psame");
    
}



void UseTree::PlotProfile(string obs, bool fillleg, bool removelabel,int color) {

  cout<<" Starting profile plotting"<<endl;

  int Obs = histoManager.FindNVarP(obs);
      
      if(Obs==-1) {
	cout<<" Be careful, no such observable "<<obs<<endl;
	return;
      }
      cout<<" begin reweight"<<endl;
      // map<string,float> weights;
      // if(nt!=0) {
      // 	weights = Reweight(Obs, false);
      // }
      cout<<" end reweight"<<endl;
      leg2D->SetFillColor(0);
      leg2D->SetLineColor(1);

      int Ndata= histoManager.accessProf(Obs,nt);
      cout<<" histo val "<<Ndata<<endl;
      //  TGraph* dataPointWhite=  (TGraph*)Histos2D[Ndata]->Clone();
      vector<TProfile*> cloneP;
      TProfile* dataProf;
      //vector<TGraph*> graphs;
      XTitle = (histoManager.FindLegP( obs )).first;
      if(!Draw3on1)
	YTitle = (histoManager.FindLegP( obs )).second;
      cout<<" title "<<XTitle<<endl;
      TH1F* emptyHisto;
      vector<float> temp = histoManager.GetTemplateP(Obs);

      double* varBin = histoManager.GetVarTemplateP(obs);
      if(varBin==NULL ) {
	emptyHisto = new TH1F(("bidon"+obs).c_str(),"bidon",(int)temp[0], temp[1], temp[2] );
      }
      else
	{
	  emptyHisto = new TH1F(("bidon"+obs).c_str(),"bidon", (int)temp[0] , 0,1000 );
	}

      emptyHisto->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
      emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);
      emptyHisto->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
      emptyHisto->GetYaxis()->SetNdivisions(Ydiv[0],Ydiv[1],Ydiv[2]);
      emptyHisto->GetXaxis()->SetTitle(XTitle.c_str());
      emptyHisto->GetYaxis()->SetTitle(YTitle.c_str());
      if(!Draw3on1) {
	emptyHisto->GetYaxis()->SetTitleOffset(1.1);
	if(RangeY[1]>= 1000)
	  emptyHisto->GetYaxis()->SetTitleOffset(1.60); 
      }

      emptyHisto->SetFillColor(0);
      emptyHisto->SetLineColor(0);
      
      //Ellipse
      TMarker* mark = new TMarker();
      mark->SetMarkerStyle(4);
      mark->SetMarkerColor(kGreen+2);
     
     
      //prepare order of histos
      for(int i=0;i<nt;i++) {
	int Nhisto = histoManager.accessProf(Obs,nt-1-i);
	cloneP.push_back((TProfile*)Profiles[Nhisto]->Clone() );
      }
      
      float wtot = 0;
      for(int i=1;i<nt;i++) {
	wtot = GetWeight(nt-1-i,0);
      }

      TProfile* MC(0);
      if(nt!=0) {
	MC =(TProfile*)cloneP[0]->Clone();
      // cout<<" number of events "<<MC->GetEntries()<<"  "<<MC->GetEntries()*weights[ name[nt-1] ]<<endl;
      for(int i=1;i<nt;i++) {
	// cout<<" number of events "<<cloneP[i]->GetEntries()<<"  "<<name[nt-1-i]<<"   "<<weights[ name[nt-1-i] ]<<"  "<<cloneP[i]->GetEntries()*weights[ name[nt-1-i] ]<<endl;
	
	MC->Add((TProfile*)cloneP[i],GetWeight(nt-1-i,0)/wtot );
		//weights[ name[nt-1-i] ]/weights[ name[nt-1]] );
      }
      }      

      dataProf = (TProfile*)Profiles[Ndata]->Clone();
      
      if(removelabel) {
	cout<<" remove label!! "<<obs<<endl;
	if(nt!=0) {
	  MC->GetYaxis()->SetLabelOffset(1000);
	}
	emptyHisto->GetYaxis()->SetLabelOffset(1000);
      }
      emptyHisto->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
      emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);

      if(nt!=0) {
	MC->Rebin(Bin);
        MC->SetLineColor(color);
	MC->SetLineWidth(LineWidth);
	MC->SetMarkerColor(color);
	MC->SetMarkerStyle(25);
      }
	emptyHisto->Draw();
    
      //MC->SetMarkerSize(1.5);
      if(BasicProf) {
	if(nt!=0) {
	  MC->Draw("same");
	}
	dataProf->SetLineWidth(LineWidth);
	dataProf->SetMarkerSize(MarkerSize);
	dataProf->Draw("same");
      }
      
      if(fillleg) {
	leg2D->AddEntry(dataProf,"data","lp");
	leg2D->AddEntry(MC,"simulation","lp");
      }

      TLine* lineR=new TLine(RangeX[0],1,RangeX[1],1);
      lineR->SetLineColor(kGreen-1);
      lineR->SetLineWidth(2);
      lineR->SetLineStyle(2);
      lineR->Draw("same");

      //Save the plots
      TFile out( (obs.substr(0,2)+"Resp.root").c_str(),"RECREATE");
      dataProf->SetName( (obs.substr(0,2)+"data").c_str() );
      dataProf->Write();
      if(nt!=0) {
	MC->SetName( (obs.substr(0,2)+"mc").c_str() );
	MC->Write();
      }
      out.Close();


      /*
      vector<TGraph*> grou;
      if(!BasicProf) {
	cout<<" Best Profile plots *************************** "<<endl;
	int BinVar[12]={0,1,10,20,30,40,50,60,80,100,150,200};
	grou = histoManager.ConvertGraphToErrorProfile(graphs[0],11,BinVar);
	cout<<" end Best Profile plots *************************** "<<endl;
      }
      
          
      
      MC->Rebin(Bin);
      emptyHisto->Draw();
      MC->SetLineColor(color);
      MC->SetLineWidth(2);
      MC->SetMarkerColor(color);
      MC->SetMarkerStyle(25);
      //MC->SetMarkerSize(1.5);
      if(BasicProf)
	MC->Draw("same");
	
      vector<vector<double> > Datavalues;
      vector<double> tmp;
      for(int i=0;i<n;i++)
	Datavalues.push_back(tmp);
      
      //  TProfile* dataProf=new TProfile(("data"+observable).c_str(),"",n,VarBin);
      TProfile* dataProf=new TProfile(("data"+obs).c_str(),"",11,0,11);
      dataProf->SetErrorOption(RMSOpt.c_str() );
      TH2F* Data(0);
      TPolyLine* polline(0);

      double x,y;
      for(int i=0;i<Histos2D[Ndata]->GetN();i++) {
	Histos2D[Ndata]->GetPoint(i,x,y);
	if(x!=0 && y!=0) { {
	  dataProf->Fill(x,y);

	  for(int k=0;k<n;k++)
	    if(VarBin[k]<=x && VarBin[k+1]>x )
	      Datavalues[k].push_back(y);
	  }
	}
      }
    

      if(BasicProf) {
	//	dataProf->Rebin(Bin);
	dataProf->SetLineWidth(2);
	dataProf->SetMarkerSize(1.2);
	dataProf->Draw("same");

	if(response) {
	  double BinVar[27] = {0,10,20,30,40,50,60,80,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000};
	  MC = histoManager.ConvertGraphToResponsePlotVarBin(graphs[0],26,BinVar,name[0], histoManager.access2D(Obs,0) );  //FIXME
	  dataProf = histoManager.ConvertGraphToResponsePlotVarBin(Histos2D[Ndata],26,BinVar,name[nt], Ndata );
	 
	  if(switchRMS) {
	   
	    MC = (TProfile*)histoManager.ConvertGraphToRMSPlotVarBin(graphs[0],26,BinVar,obs, histoManager.access2D(Obs,0), false );
	    dataProf = (TProfile*)histoManager.ConvertGraphToRMSPlotVarBin(Histos2D[Ndata],26,BinVar,obs, Ndata, Remove2V);
	    // cout<<dataProf->GetBinError(0)<<"  "<<dataProf->GetBinError(1)<<endl;
	  }
	  //  c2->cd();
	  emptyHisto->Draw();
	  MC->SetLineColor(color);
	  MC->SetLineWidth(2);
	  MC->SetMarkerColor(color);
	  MC->SetMarkerStyle(25);
	  dataProf->SetLineWidth(2);
	  dataProf->SetMarkerSize(1.2);
	  MC->Draw("same");
	  dataProf->Draw("same");
	}
	
      }
      else {
        double BinX[20]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,125,150,200};
      double BinY[21]={0,0.157,0.314,0.471,0.628,0.785,0.942,
		       1.099,1.256,1.413,1.57,1.727,1.884,2.041,
		       2.198,2.355,2.512,2.669,2.826,2.983,3.14};

      Data = histoManager.ConvertGraphToHisto(dataPointWhite,19,BinX,20,BinY,"greu");

      dataPointWhite->SetMarkerStyle(20);
      dataPointWhite->SetMarkerSize(0.6);
      // dataPointWhite->Draw("P");
      gStyle->SetPalette(156);
     
      Data->SetFillColor(color);
      //dataProf->Draw("same");

      vector<TGraph*> gra;
      gra.push_back(grou[1]);
      gra.push_back(grou[2]);

      polline = PolyLineFromGraph(gra);
     
      polline->SetFillColor(kYellow-9);
      //polline->Draw("Fsame");
    
     
      grou[0]->SetLineColor(kRed+2);
      grou[1]->SetLineColor(kYellow-9);
      grou[2]->SetLineColor(kYellow-9);
      grou[0]->SetLineWidth(2);
      grou[1]->SetLineWidth(2);
      grou[2]->SetLineWidth(2);
      grou[1]->SetLineStyle(2);
      grou[2]->SetLineStyle(2);
      grou[0]->SetMarkerStyle(0);
      grou[1]->SetMarkerStyle(0);
      grou[2]->SetMarkerStyle(0);

      if(!Draw3on1)
	grou[0]->SetLineWidth(4);
      Data->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);
      Data->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
      if(removelabel) {
	Data->GetYaxis()->SetTitleSize(0);
	Data->GetYaxis()->SetLabelSize(0);
      }
      else
	Data->GetYaxis()->SetTitle(YTitle.c_str());
      Data->GetXaxis()->SetTitle(XTitle.c_str());
      Data->Draw("box");
      polline->Draw("Fsame");
      Data->Draw("samebox");
      grou[0]->Draw("Lsame");
    
      //  grou[1]->Draw("L");
      //  grou[2]->Draw("L");
   
      // emptyHisto->Draw("same");
      //c2->RedrawAxis();
      }

      if(fillleg) {
	if(BasicProf)
	  leg2D->AddEntry(dataProf,"data","lp");
	else
	  leg2D->AddEntry(Data,"data","f");
	if(!switchRMS) {
	  if(BasicProf )
	    leg2D->AddEntry(MC,"Z #rightarrow ee MC","lp");
	  else {
	    leg2D->AddEntry(grou[0],"Z #rightarrow ee MC","l");
	    leg2D->AddEntry(polline,"C.L. 68%","f");
	  }
	}
      }
      if(switchRMS)
	leg2D->AddEntry(MC,"Z #rightarrow ee MC","lp");

       if(switchRMS && !response){
	TH1F* MC2=new TH1F("cloneMC","",MC->GetNbinsX(),0,200);

	leg2D->AddEntry(MC2,"Pythia6 MC","pl");
	for(int i=0;i<MC->GetNbinsX();i++) {
	  //float kurt = Kurt(MCvalues[i],MC->GetBinContent(i) , MC->GetBinError(i));
	  MC2->SetBinContent(i, MC->GetBinError(i) );
	  MC2->SetBinError(i, 0 );
	}
	//TH1F* data2=new TH1F("clonedata","",n,VarBin);
	TH1F* data2=new TH1F("clonedata","",11,0,11);
	for(int i=0;i<dataProf->GetNbinsX();i++) {
	  if(dataProf->GetBinEntries(i)!=0) {
	    float kurt = Kurt(Datavalues[i],dataProf->GetBinContent(i) , dataProf->GetBinError(i));
 	    data2->SetBinContent(i,dataProf->GetBinError(i) );
	    data2->SetBinError(i, kurt );
	    cout<<Datavalues[i].size()<<"  "<<i<<"  "<<kurt<<endl;
	  }
	}


	MC2->SetLineColor(color);
	MC2->SetLineWidth(2);
	MC2->SetMarkerColor(color);
	MC2->SetMarkerStyle(29);
	MC2->SetMarkerSize(1.2);
	data2->SetLineWidth(2);
	data2->SetMarkerSize(1.2);
	emptyHisto->Draw();
	MC2->Draw("samePE");
	

	
	//data2->Draw("samePE");
       }
       */
      
}

void UseTree::PlotRMS(string obs, bool fillleg, bool removelabel,int color) {

  int Obs = histoManager.FindNVar2D(obs);
      
  if(Obs==-1) {
    cout<<" Be careful, no such observable "<<obs<<endl;
    return;
  }

      map<string,float> weights;
      weights = Reweight(Obs, false);
     
      leg2D->SetFillColor(0);
      leg2D->SetLineColor(1);

      int Ndata= histoManager.access2D(Obs,nt);
     
      XTitle = (histoManager.FindLeg2D( obs )).first;
      if(!Draw3on1)
	YTitle = (histoManager.FindLeg2D( observable )).second;

      vector<float> temp = histoManager.GetTemplate2D(Obs);
      cout<<" Template RMS "<<temp[0]<<"   "<<temp[1]<<"   "<<temp[2]<<"  mc used  "<<name[nt-1]<<endl;
      TH1F* emptyHisto= new TH1F(("bidon"+obs).c_str(),"bidon",(int)temp[0],temp[1],temp[2]);
      emptyHisto->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
      emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);
      emptyHisto->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
      emptyHisto->GetYaxis()->SetNdivisions(Ydiv[0],Ydiv[1],Ydiv[2]);
      emptyHisto->GetXaxis()->SetTitle(XTitle.c_str());
      emptyHisto->GetYaxis()->SetTitle(YTitle.c_str());
      emptyHisto->SetFillColor(0);
      emptyHisto->SetLineColor(0);
      emptyHisto->Fill(0.,0.);
      
      int Nhisto = histoManager.access2D(Obs,nt-1);
      TH2F* HSignal = (TH2F*)Histos2D[ Nhisto ]->Clone();
      TH2F* HData = (TH2F*)Histos2D[ Ndata ]->Clone();
      
      int Nbin=26;
      //double BinVar[27] = {0,10,20,30,40,50,60,80,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000};
      double BinVar[27] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,350,400,450,500};
      TH1D* RMSMC = histoManager.ConvertHistoToRMSPlotVarBin(HSignal,Nbin,BinVar,obs,Remove2V,false);  c2->cd();
      TH1D* RMSData = histoManager.ConvertHistoToRMSPlotVarBin(HData,Nbin,BinVar,obs,Remove2V,false);  c2->cd();


      if(fillleg) {
	leg2D->AddEntry(RMSData,"data","lp");
	leg2D->AddEntry(RMSMC,"Z #rightarrow ee MC","lp");
      }

      if(removelabel) {
	cout<<" remove label!! "<<obs<<endl;
	RMSMC->GetYaxis()->SetLabelOffset(1000);
	emptyHisto->GetYaxis()->SetLabelOffset(1000);
      }
     
      RMSData->SetLineWidth(LineWidth);
      RMSData->SetMarkerSize(MarkerSize);
      RMSMC->SetLineColor(color);
      RMSMC->SetLineWidth(LineWidth);
      RMSMC->SetMarkerColor(color);
      RMSMC->SetMarkerStyle(25);
      
      emptyHisto->Draw();
      RMSMC->Draw("same");
      RMSData->Draw("same");

      //   TF1* f1fit = new TF1( ("f1"+obs).c_str(),"sqrt( pow( sqrt(x+[1])*[0]+[2]  ,2 ) )",30,100);

      //      RMSData->Fit(("f1"+obs).c_str(),"R");


      //Save the plots
      TFile out( (obs.substr(0,2)+"Reso.root").c_str(),"RECREATE");
      RMSData->SetName( (obs.substr(0,2)+"data").c_str() );
      RMSData->Write();
      RMSMC->SetName( (obs.substr(0,2)+"mc").c_str() );
      RMSMC->Write();
      out.Close();


}


void UseTree::PlotEffRej(string Obs, string ssig, string Obs2, string Obs3 ) {

  bool dv=false;
  if(Obs2!="")
    dv=true;

  TCanvas* c4=new TCanvas("c4","Test",300,300,600,600);
  c4->cd();
  leg = new TLegend(0.70,0.62,0.9,0.82);
  leg->SetLineColor(0);
  leg->SetFillColor(0);

  vector<string> Observables;
  Observables.push_back(Obs);
  Observables.push_back(Obs2);
  if(Obs3!="")
    Observables.push_back(Obs3);
  
  vector<int> Colors;
  Colors.push_back(kTeal-6);
  Colors.push_back(kRed+1);
  Colors.push_back(kSpring+7);
  
  vector<int> Markers;
  Markers.push_back(20);
  Markers.push_back(22);
  Markers.push_back(21);

  int N = GetNbyName(ssig);

  vector<TGraphErrors*> EffRejs;

  for(int unsigned i=0;i<Observables.size();i++) {

    TH1F* signal(0);
    TH1F* totalObs(0);
    
    int nvar = histoManager.FindNVar((Observables[i]));
   
    //  Reweight(nvar, true);

    for(int j=0;j<nt;j++) {
     
      int Nhisto = histoManager.access(nvar, j);
      
      if(j==N) {
	signal = (TH1F*)Histos[Nhisto]->Clone();
      }
      else {
	if(totalObs==NULL) {
	  totalObs = (TH1F*)Histos[Nhisto]->Clone();
	}
	else {
	  totalObs->Add(Histos[Nhisto]);
	}
      }
    }//End dataset loop
    
    EffRejs.push_back( histoManager.EffRej( signal, totalObs) );
    EffRejs[i]->GetXaxis()->SetTitle("signal efficiency" );
    EffRejs[i]->GetYaxis()->SetTitle("background efficiency" );
    EffRejs[i]->GetYaxis()->SetRangeUser(0.00001,1);
    EffRejs[i]->SetMarkerStyle(Markers[i]);
    EffRejs[i]->SetMarkerColor(Colors[i]);
    EffRejs[i]->SetLineColor(Colors[i]);
    leg->AddEntry(EffRejs[i], (Observables[i]).c_str(),"pl" );

    if(i==0)
      EffRejs[i]->Draw("APL");
    else
      EffRejs[i]->Draw("PL");

    delete signal;
    delete totalObs;

    if(!dv){
      break;
    }

  }
  leg->Draw("same");
}

void UseTree::PlotSigRej(string Obs, string ssig, string Obs2, string Obs3) {


  bool dv=false;
  if(Obs2!="")
    dv=true;
  
  TCanvas* c4=new TCanvas("c5","Test",300,300,600,600);
  c4->cd();
  leg = new TLegend(0.70,0.62,0.9,0.82);
  leg->SetLineColor(0);
  leg->SetFillColor(0);

  vector<string> Observables;
  Observables.push_back(Obs);
  Observables.push_back(Obs2);
  if(Obs3!="")
    Observables.push_back(Obs3);
  
  vector<int> Colors;
  Colors.push_back(kTeal-6);
  Colors.push_back(kRed+1);
  Colors.push_back(kSpring+7);
  
  vector<int> Markers;
  Markers.push_back(20);
  Markers.push_back(22);
  Markers.push_back(21);

  int N = GetNbyName(ssig);
  
  vector<TGraphErrors*> SigRejs;

  for(int unsigned i=0;i<Observables.size();i++) {

    TH1F* signal(0);
    TH1F* totalObs(0);
    
    int nvar = histoManager.FindNVar((Observables[i]));
    
    for(int j=0;j<nt;j++) {
     
      int Nhisto = histoManager.access(nvar, j);
      
      if(j==N) {
	signal = (TH1F*)Histos[Nhisto]->Clone();
      }
      else {
	if(totalObs==NULL) {
	  totalObs = (TH1F*)Histos[Nhisto]->Clone();
	}
	else {
	  totalObs->Add(Histos[Nhisto]);
	}
      }
    }//End dataset loop
    
    SigRejs.push_back( histoManager.SigniErr( signal, totalObs) );
    SigRejs[i]->GetXaxis()->SetTitle("Signal efficiency" );
    SigRejs[i]->GetYaxis()->SetTitle("S/B" );
    SigRejs[i]->SetMarkerStyle(Markers[i]);
    SigRejs[i]->SetMarkerColor(Colors[i]);
    SigRejs[i]->SetLineColor(Colors[i]);
    SigRejs[i]->GetYaxis()->SetRangeUser(0,60);
    leg->AddEntry(SigRejs[i], (Observables[i]).c_str(),"pl" );

    if(i==0)
      SigRejs[i]->Draw("APL");
    else
      SigRejs[i]->Draw("PL");

    delete signal;
    delete totalObs;

    if(!dv){
      break;
    }

  }

  leg->Draw("same");
}



map<string,float> UseTree::Reweight(int Obs, bool OneDim) {
 
  map<string,float> weights;
  //QCD Weight *************************************************
  double nEWK= 0;
  double nGJ=0;
  double nData=0;

  if(!MultiWeight) {

  if(OneDim) {

    for(int i=0;i<nt+1;i++) {
      
      int Nhisto = histoManager.access(Obs, i);
      
      if("EWK"==name[i] || "W"==name[i].substr(0,1) )
	nEWK += Histos[Nhisto]->Integral(0,Histos[Nhisto]->GetNbinsX()+1); 
      else if("#gamma+jet"==name[i])
	nGJ = Histos[Nhisto]->Integral(0,Histos[Nhisto]->GetNbinsX()+1); 
      else if(i==nt)
	nData= Histos[Nhisto]->Integral(0,Histos[Nhisto]->GetNbinsX()+1);
    }

  }
  else {

      for(int i=0;i<nt+1;i++) {
	int Nhisto = histoManager.access2D(Obs, i);
	
	if("EWK"==name[i] || "W"==name[i].substr(0,1) )
	  nEWK += Histos2D[Nhisto]->GetEntries(); 
	else if("#gamma+jet"==name[i])
	  nGJ = Histos2D[Nhisto]->GetEntries(); 
	else if(i==nt)
	  nData= Histos2D[Nhisto]->GetEntries();
      }
  }

  QCDweight = (nData - nEWK - nGJ);
  //************************************************************
  
  for(int i=0;i<nt;i++) {
    
    //Binning and coloring :s
    if(OneDim) {
      
      int Nhisto = histoManager.access(Obs, i);

      if(Weighted[Nhisto]) continue; //Skip histos already weighted
      
      float Weight =1;
      if("QCD"==name[i] ) {
	cout<<" NQcd "<<Histos[Nhisto]->Integral(0,Histos[Nhisto]->GetNbinsX()+1)<<endl;
	QCDweight = QCDweight/ Histos[Nhisto]->Integral(0,Histos[Nhisto]->GetNbinsX()+1);
	cout<<QCDweight*Histos[Nhisto]->Integral()<<"  + "<<nEWK+nGJ<<" = "<<QCDweight*Histos[Nhisto]->Integral()+nEWK<<" <> "<<nData<<endl;
	Weight = QCDweight;
	//normQCDWeight = Lumi;
	if(!QCDAutoWeight) {
	  Weight = 1;
	}
      }
      else
	Weight = 1;
      if("QCD"==name[i] ) {
	if(QCDType=="EM")
	  Histos[Nhisto]->Scale(Weight); 
	else
	  Histos[Nhisto]->Scale(Weight); 
      }
      
      Weighted[Nhisto] = true;
      weights[ name[i] ] = Weight;
      
    }
    else { //2D 
      int Nhisto = histoManager.access2D(Obs, i);
      
      float Weight =1;
      if("QCD"==name[i] ) {
	cout<<" NQcd "<<Histos2D[Nhisto]->GetEntries()<<endl;
	QCDweight = QCDweight/ Histos2D[Nhisto]->GetEntries();
	cout<<QCDweight*Histos2D[Nhisto]->GetEntries()<<"  + "
	    <<nEWK+nGJ<<" = "<<QCDweight*Histos2D[Nhisto]->GetEntries()+nEWK
	    <<" <> "<<nData<<endl;
	Weight = QCDweight;
	if(!QCDAutoWeight) {
	  Weight = 1;
	  cout<<" Weight by lumi : no AutoWeight "<<endl;
	}
      }
      else 
	Weight =1;
      
      weights[ name[i] ] = Weight;
    }
  }

  }

  return weights;

}



void UseTree::PrintSelectedEvents() {

  cout<<" Event selected "<<endl;
  for(EventIter=Events.begin();EventIter!=Events.end();EventIter++)
    {
      cout<<"run "<<(*EventIter).first.first<<"  :event "<<(*EventIter).first.second
	  <<"  :Sample  "<<((*EventIter).second.first)
	  <<"  :name  "<<((*EventIter).second.second)<<endl;//"  --->  "<<EvtsInFile[i]<<endl;
    }

}

void UseTree::SaveSelectedEvents() {

  ofstream file("SavedEvents.txt", ios::out | ios::trunc);
  
  if(file) {

    string tmp="";

    for(EventIter=Events.begin();EventIter!=Events.end();EventIter++)
      {
	if(tmp != (*EventIter).second.first ) {
	  file << endl << (*EventIter).second.first <<endl;
	  tmp = (*EventIter).second.first;
	}
	file <<(*EventIter).first.first<<"  "
	     <<(*EventIter).first.second<<endl;
      }

  }
  else cout<<" No such file to save events "<<endl;

  file.close();

}


float UseTree::FitReweight(vector<TH1F*> histos, TH1F* dataH) {

  float al=0.5;
  float tmp=10000;
  float LL=0;

  int NBins = dataH->GetNbinsX();
  // TCanvas* c=new TCanvas("c","LogLikelihood",300,300,479,510);
  TGraph* LogLike=new TGraph(NBins);
  
  float NormEWK = histos[0]->Integral(0,histos[0]->GetNbinsX()+1) + histos[3]->Integral(0,histos[3]->GetNbinsX()+1) + histos[4]->Integral(0,histos[4]->GetNbinsX()+1);
  float NormQCD = histos[1]->Integral(0,histos[1]->GetNbinsX()+1) + histos[2]->Integral(0,histos[2]->GetNbinsX()+1);
  float NormData = dataH->Integral(0,dataH->GetNbinsX()+1);

  for(int k=0;k<2000;k++) {
    
    float alpha = k*0.0005;
    float S=0;
    float B=0;

    LL=0;
    for(int i=0;i<NBins;i++) {

      S=0;
      B=0;
      
      B += histos[1]->GetBinContent(i) + 
	histos[2]->GetBinContent(i) ;
	//histos[2]->GetBinContent(i);

      S +=histos[0]->GetBinContent(i) + histos[4]->GetBinContent(i) + histos[3]->GetBinContent(i) ;
      
      B = B/NormQCD;
      S = S/NormEWK;

      // cout<<S<<"     "<<B<<"  "<<data->GetBinContent(i)<<endl;
      if(S!=0 && B!=0 && dataH->GetBinContent(i)!=0)
	LL += UseTree::LnL(S, B, dataH->GetBinContent(i),NormData, alpha);
    }
    LogLike->SetPoint(k,alpha,LL);
    // cout<<LL<<"    "<<alpha<<endl;
    if(tmp>LL)
    {
      tmp = LL;
      al=alpha;
      //   cout<<" test "<<al<<endl;
    }
  
  }
  /*  LogLike->Draw("A*");
  LogLike->GetXaxis()->SetTitle("#alpha");
  LogLike->GetYaxis()->SetTitle("LL");
  c2->cd();*/
  return al;

}


float UseTree::LnL(float si, float bi, int ni,int NData, float a) {

  double L=0;
  double lmu = (a*si + (1-a)*bi )*NData;
  if(ni<20) {
    double lnmu= log( pow(lmu,(double)ni)/Fact(ni) );
    L = lmu - lnmu;
  }
  else {
    L = log(sqrt(2*acos(-1)*ni) ) + pow((ni- lmu),2)/(2*ni) ;
  }
    
  return L;

}

float UseTree::Fact(int x)
{

  if(x==0) return 1;

  float resultat = 1;
  for(int i = 1; i <= x; i++)
    resultat *= i;
  return resultat;
}


float UseTree::Chi2(float si, float bi, int ni,int NData , float a){

  float chi2=0;
  chi2 = pow( ((a*si + (1-a)*bi )*NData)-ni, 2);
  return chi2;
}


float UseTree::Kurt(vector<double> values, double mean, double sigma) {

  if(values.size()<2)
    return 0;

  float kurt=0;
  int N = values.size();
  float mom4=0;

  for(int i=0;i<N;i++) {
    mom4 += pow( values[i] - mean, 4);
  }
  mom4 /= N;

  kurt = (mom4- ( (N-3)/(N-1) * sigma ))/N;

  return kurt;

}


TPolyLine* UseTree::PolyLineFromGraph(TGraph* graph) {

  vector<double> xL,yL,xH,yH,X,Y;
  double x,y;
  for(int i=0;i<graph->GetN();i++) {
    graph->GetPoint(i,x,y);
    X.push_back(x);
    Y.push_back(y);
    xL.push_back( graph->GetErrorXlow(i) );
    xH.push_back( graph->GetErrorXhigh(i) );
    yL.push_back( graph->GetErrorYlow(i) );
    yH.push_back( graph->GetErrorYhigh(i) );
  }

  
  TPolyLine* tmp = new TPolyLine(yL.size()*2 );

  for(size_t i=0;i<yL.size();i++) {
    tmp->SetPoint(i,X[i],Y[i]-yL[i] );
  }
  
  for(size_t i=0;i<yH.size();i++) {
    if(Y[yH.size()-1-i]+yH[yH.size()-1-i] < 3.14)
      tmp->SetPoint(i+yH.size(),X[yH.size()-1-i],Y[yH.size()-1-i]+yH[yH.size()-1-i] );
    else
      tmp->SetPoint(i+yH.size(),X[yH.size()-1-i],3.14 );
  }

  return tmp;

}

TPolyLine* UseTree::PolyLineFromGraph(vector<TGraph*> graphs) {

  vector<double> x,y;
  
  double tmpx, tmpy;
  for(size_t i=0;i<graphs.size();i++) {

    for(int k=0;k<graphs[i]->GetN();k++) {
      graphs[i]->GetPoint(k,tmpx,tmpy);
      x.push_back(tmpx);
      y.push_back(tmpy);
    }
  }


  TPolyLine* tmp = new TPolyLine(x.size()-2);

  for(size_t i=0;i<x.size()/2-1;i++) {
    tmp->SetPoint(i,x[i],y[i]);

    cout<<i<<" x "<<x[i]<<"  "<<y[i]<<endl;
  }
  for(size_t i=1;i<x.size()/2;i++) {
    tmp->SetPoint(i-2+x.size()/2,x[x.size()-1-i],y[y.size()-1-i]);

    cout<<i-2+x.size()/2<<" x "<<x[x.size()-i-1]<<"  "<<y[x.size()-1-i]<<endl;
  }

  return tmp;

}

void UseTree::FillAddWeight(string dataset) {

  
  float Weight=1;
  
  map<string,float>::const_iterator iterLumi;

  float lum =1.;
 
  iterLumi=Luminosities.find(dataset);
  if(iterLumi!=Luminosities.end() )
    lum = (*iterLumi).second;
  else if(useXSect){
    float nevt =1.;
    float xsect =1.;
    map<string,float>::const_iterator iterNIE;
    map<string,float>::const_iterator iterXSect;
    iterNIE=internalNEvt.find(dataset);
    nevt = (*iterNIE).second;
    
    iterXSect=XSections.find(dataset);
    if(iterXSect==XSections.end()) {
      xsect =1.;
      cout<<" Be careful, no specified Xsect for "<<dataset<<" , 1 by default *****======== "<<endl;
    }
    else
      xsect = (*iterXSect).second; 
    
    cout<<" Number of events "<<nevt<<" / Xs "<<xsect<<" = "<<nevt/xsect<<"  :----> ";
    lum = nevt/xsect;
  }
  else
    cout<<" Be careful, no specified lumi for "<<dataset<<" , 1 by default *****======== "<<endl;

  map<string,float>::const_iterator iterKF;
  float KFactor =1.;
  iterKF=KFactors.find(dataset);
  if(iterKF!=KFactors.end() )
    KFactor = (*iterKF).second;
   
  Weight =(1./lum)*Lumi*KFactor;
  
  Weights[dataset]=Weight;

  cout<<" dataset "<<dataset<<" ----> "<<Weight<<endl;

}

void UseTree::FillEventMap(string dataset, int NevtB, int channel ) {

  std::pair<string, int> tmp(dataset,NevtB); 
  ChanNum[channel].push_back(tmp);

}


void UseTree::PrintChanNum() {

  for(size_t i=0;i<ChanNum.size();i++) {
    cout<<name[i]<<endl;
    for(size_t j=0;j<ChanNum[i].size();j++) {
      cout<<"   --> "<<ChanNum[i][j].first<<"    "<<ChanNum[i][j].second<<endl;
    }
  }
}


float UseTree::GetWeight(int channel, int evt) {

  if(evt==-1) return 0; //shield against initalization step

  //Get dataset
  for(size_t i=0;i<ChanNum[channel].size();i++) {
    
    if(evt >= (ChanNum[channel][i]).second && evt < (ChanNum[channel][i+1]).second ) {
      
      return Weights[ (ChanNum[channel][i]).first ];
    }
  }
  
  cout<<" ********************** Warning Weight "<<endl;
  return 1;

}

string UseTree::GetDataset(int channel, int evt) {

  if(evt==-1) return ""; //shield against initalization step
  if(channel==nt) return "data"; //shield against data

  //Get dataset
  for(size_t i=0;i<ChanNum[channel].size();i++) {
    
    if(evt >= (ChanNum[channel][i]).second && evt < (ChanNum[channel][i+1]).second ) {
      
      return (ChanNum[channel][i]).first;
    }
  }
  
  cout<<" ********************** Warning Dataset "<<endl;
  return "";

}



/*vector<vector<vector<float> > >  UseTree::GetFitWeight() {

}
*/

float
UseTree::GetWeightFromShape(int met, float var, string shref, string shvar) {

  if(shape[met]==NULL) {
  int Sref =  histoManager.FindNVar(shref);
  int Svar =  histoManager.FindNVar(shvar);

  int nhref= histoManager.access(Sref,nt);
  int nhvar= histoManager.access(Svar,nt);
  TH1F* href = (TH1F*)Histos[ nhref ]->Clone();
  TH1F* hvar = (TH1F*)Histos[ nhvar ]->Clone();
  
  int Obs = histoManager.FindNVar(observable);
  vector<float> templ = histoManager.GetTemplate(Obs);
  //float pas = (templ[2]-templ[1])/templ[0];
  

  href->Rebin( /*Bin*pas*/10./href->GetBinWidth(1) );
  hvar->Rebin( /*Bin*pas*/10./hvar->GetBinWidth(1) );
 
  float integref = href->Integral(0, href->GetNbinsX() +1 );
  float integvar = hvar->Integral(0, href->GetNbinsX() +1 );

  href->Scale( 1./integref );
  hvar->Scale( 1./integvar );

  shape[met] = HistoManager::ShapeWeight(href,hvar);
  delete href;
  delete hvar;
  }

  float w = shape[met]->GetBinContent( shape[met]->FindBin(var) );
  /* if(met==0)
     cout<<shape[met]->FindBin(var)<<"   "<<weight<<endl;*/
  int ii= shape[met]->FindBin(var);
  while(w==0 && ii!=0 )
    { w =  shape[met]->GetBinContent( ii ); ii--; }
  
  return w;
}



float
UseTree::GetMCWeightFromDataShape(int met, float var, string shref, string shvar, double bin, int ds) {

  if(shape[met]==NULL) {
    //cout<<" MET "<<met<<"  encore nul "<<"   "<<shref<<"   "<<shvar<<endl;
    int Sref =  histoManager.FindNVar(shref);
    int Svar =  histoManager.FindNVar(shvar);
    //cout<<"  "<<Sref<<"   "<<Svar<<endl;
    int nhref= histoManager.access(Sref,nt);
    int nhvar= histoManager.access(Svar,ds);
    TH1F* href = (TH1F*)Histos[ nhref ]->Clone();
    TH1F* hvar = (TH1F*)Histos[ nhvar ]->Clone();
    //cout<<" numbers "<<" ref  "<<nhref<<" var  "<<nhvar<<"   "<<observable<<endl;
    int Obs = histoManager.FindNVar(shref);
    if(Obs == -1) return 1;
    vector<float> templ = histoManager.GetTemplate(Obs);
    //float pas = (templ[2]-templ[1])/templ[0];
    //    cout<<templ.size()<<endl;

    href->Rebin( /*Bin*pas*/bin/href->GetBinWidth(1) );
    hvar->Rebin( /*Bin*pas*/bin/hvar->GetBinWidth(1) );
 
    float integref = href->Integral(0, href->GetNbinsX() +1 );
    float integvar = hvar->Integral(0, href->GetNbinsX() +1 );
    //cout<<integref<<"   "<<integvar<<endl;
    href->Scale( 1./integref );
    hvar->Scale( 1./integvar );

    shape[met] = HistoManager::ShapeWeight(href,hvar);
    delete href;
    delete hvar;
  }

  float w = shape[met]->GetBinContent( shape[met]->FindBin(var) );
  /*if(met==1)
    cout<<shape[met]->FindBin(var)<<"   "<<weight<<endl;*/
  int ii= shape[met]->FindBin(var);
  while(w==0 && ii!=0 )
    { w =  shape[met]->GetBinContent( ii ); ii--; }
  
  // cout << " Weight -->" <<w<<endl;

  return w;
}


float
UseTree::GetMCWeightFromDataShape2D( float var1, float var2, string shref, string shvar) {

  if(shape2D==NULL) {
    int Sref =  histoManager.FindNVar2D(shref);
    int Svar =  histoManager.FindNVar2D(shvar);
 
    int nhref= histoManager.access2D(Sref,nt);
    int nhvar= histoManager.access2D(Svar,nt-1);
    TH2F* href = (TH2F*)Histos2D[ nhref ]->Clone();
    TH2F* hvar = (TH2F*)Histos2D[ nhvar ]->Clone();
    cout<<" numbers "<<" ref  "<<nhref<<" var  "<<nhvar<<endl;
    //int Obs = histoManager.FindNVar(observable);
  
    float integref = href->Integral(0, href->GetNbinsX() +1, 0, href->GetNbinsY() +1 );
    float integvar = hvar->Integral(0, href->GetNbinsX() +1, 0, href->GetNbinsY() +1 );

    href->Scale( 1./integref );
    hvar->Scale( 1./integvar );

    shape2D = HistoManager::ShapeWeight2D(href,hvar);
    delete href;
    delete hvar;
  }

  float w = shape2D->GetBinContent( shape2D->GetXaxis()->FindBin(var1), shape2D->GetYaxis()->FindBin(var2) );
  
  int iix= shape2D->GetXaxis()->FindBin(var1);
  int iiy= shape2D->GetYaxis()->FindBin(var2);
  
  while(w==0 && iix!=0 )
    {
      while(w==0 && iiy!=0 )
	{
	  w =  shape2D->GetBinContent( iix, iiy ); iiy--; 
	}
      iix--; 
    }
  
  return w;
}

vector<float> UseTree::GetAlpha(TString obs) {

   vector<TH1F*> ClonesMT;
  string obsMT = (string)obs;
  
   int ObsMT = histoManager.FindNVar(obsMT);
 
   for(int i=0;i<nt;i++) {
    int NhistoMT = histoManager.access(ObsMT, i);
    ClonesMT.push_back( (TH1F*)Histos[NhistoMT]->Clone() );
    ClonesMT[i]->Rebin(2);
   
  }
  int NdataMT= histoManager.access(ObsMT,nt);
  TH1F* CloneDataMT = (TH1F*)Histos[ NdataMT ]->Clone();
  CloneDataMT->Rebin(2);

  if(Remove2V) {
    TH1F* V2Cor = GetVtxCorrection(obsMT,nt,2);
    CloneDataMT->Add(V2Cor,-1);
    for(int i=1;i<3;i++) {
      TH1F* V2Cor = GetVtxCorrection(obsMT, nt-i,Bin);
      ClonesMT[nt-i]->Add(V2Cor,-1);
    }
  }

  float a= FitReweight(ClonesMT, CloneDataMT);
  cout<<obsMT<<" Best a "<<a<<" ; 0.585 as reference  -> r="<<a/0.585<<endl;
  
  float NormEWK = ClonesMT[0]->Integral(0, ClonesMT[0]->GetNbinsX()+1) + ClonesMT[3]->Integral(0, ClonesMT[3]->GetNbinsX()+1) + ClonesMT[4]->Integral(0, ClonesMT[4]->GetNbinsX()+1);
  float NormQCD = ClonesMT[1]->Integral(0, ClonesMT[1]->GetNbinsX()+1) + ClonesMT[2]->Integral(0, ClonesMT[2]->GetNbinsX()+1);
  float NormData = CloneDataMT->Integral(0, CloneDataMT->GetNbinsX()+1);

  vector<float> Norms(2,0);
 
  Norms[0] = (a*NormData/NormEWK) ;
  Norms[1] = ((1-a)*NormData/NormQCD) ;
  return Norms;
  
}

void UseTree::ConfigureLumi(map<string,float> LumisXS, map<string,float> Kfac, float l,bool useXS) {
  
  KFactors = Kfac;
  Lumi =l;
  useXSect = useXS;

  if(useXSect) {
  XSections = LumisXS;
  } else {
    Luminosities = LumisXS;
  }
}



void UseTree::KolmogorovTest(TH1F* MC, TH1F* dataH, string obs) {

  double Kolm;
  Kolm= dataH->KolmogorovTest( MC,"UO");
  // cout<<" Kolmogorov Test on observable "<<observable<<" : "<<Kolm<<endl;

  double Chi2norm,Chi2prob;
  Chi2prob = dataH->Chi2Test(MC,"UW P");
  Chi2norm = dataH->Chi2Test(MC,"UW Chi2/NDF");
  //cout<<" Chi2/NDF : "<< Chi2norm<<endl;

  cout<<obs<<" ---> KS prob = "<<Kolm<<"   ; Chi2 prob = "<< Chi2prob<<"  ; Chi2/NDF = "<<Chi2norm<<" whithout weight "<<dataH->Chi2Test(MC,"P")<<endl;

}



float UseTree::Poisson(float mc, int dataN) {

  float poisson;

  if(dataN<20) {
    double lnmu= log( pow(mc,(float)dataN)/Fact(dataN) );
    poisson = mc - lnmu;
  }
   else //normal law approx
    poisson = log(sqrt(2*acos(-1)*dataN) ) + pow((dataN- mc),2)/(2*dataN) ;
  return poisson;
}



void UseTree::Chi2Test(TH1F* MC, TH1F* dataH) {
  
  float chi2=0;
  //  float ei=0;
  float ni=0;
  float li=0;
  int N=0;

  for(int i=0;i<MC->GetNbinsX();i++) {

    ni = dataH->GetBinContent(i);
    li = MC->GetBinContent(i);

    // ei=li;
    if( li!=0 ) 
      {
	chi2 += Poisson(li,ni);
	N++;
      }
    }
  cout<<"Chi2 = "<<chi2<<"   ;  ndf = "<<N<<"  chi2/ndof "<<chi2/(N)<<endl;

}



/*TH1F*
UseTree::GetVtxCorrection(string obs, int bin) {
  
}
*/


void 
UseTree::PrepareDatasets() {


}

void
UseTree::GetNProcessedEvent() {

}

void
UseTree::AddMCSample( string str, string sname, int col) {
  

 
   if(sname != "") //not the same sample
    {
      colors.push_back(col);
      name.push_back(sname);
      nt++;
      
      snametmp = sname;
    }
   //   cout<<str<<"   "<<nt<<endl;
   datasets.push_back( str );
   dt.push_back( nt );

}

void 
UseTree::ParserNorm(string str) {

  Norm=false;
  NormPlot=false;
  QCDAutoWeight = false;
  FitR = false;
  MultiWeight = false;
  Remove2V = false;
  ShapeWeight = false;

  if(str=="QCDAuto")
    QCDAutoWeight = true;
  if(str=="Fit")
    FitR = true;
  if(str=="Multi")
    MultiWeight = true;
  if(str=="Remove2V")
    {
      cout<<" Removing of 2 vertices  contamination "<<endl;
      NVert=1;
      Remove2V = true;
      FitR = true;
      AddSyst = true;
    }
  if(str=="ShapeWeight") 
    {
      FitR=true;
      ShapeWeight=true;
      //NVert=2;
    }
  if(str=="ShapeWeightNorm") 
    {
      Norm=true;
      ShapeWeight=true;
      AddSyst = true;
      //NVert=2;
    }
  if("ShapePlusPUClean"==str)
    {
      cout<<" Removing of 2 vertices  contamination "<<endl;
      NVert=1;
      Remove2V = true;
      FitR = true; 
      ShapeWeight=true;
      NormPlot=false;
    }
  if(str=="Norm")
    {
      Norm=true;
    }
  if(str=="NormClean")
    {
      cout<<" Removing of 2 vertices  contamination "<<endl;
      NVert=1;
      Norm=true;
      Remove2V = true;
      AddSyst = true;
    }
  if(str=="NormPlot")
    {
      Norm=true;
      NormPlot=true;
    }
  if(str=="NormPlotClean")
    {
      cout<<" Removing of 2 vertices  contamination "<<endl;
      NormPlot=true;
      NVert=1;
      Remove2V = true;
    }
  if(str=="NormSyst")
    {
      Norm=true;
      AddSyst = true;
    }

}




void
UseTree::AddSystematics(string obs, TH1F* MC, float W) {

  // TH1F* MCUncert =  (TH1F*)MC->Clone();
  
  // int Obs = histoManager.FindNVar(obs+"Unc");
  // int NS= histoManager.access(Obs,nt-1);
  // TH1F* CDatU = (TH1F*)Histos[ NS ]->Clone();
  // CDatU->Rebin(Bin);
  // cout<<" bin "<<CDatU->GetBinWidth(3)<<endl;
  // double syst;

  // for(int ip=0;ip<MC->GetNbinsX();ip++) {
    
  //   if(Remove2V)
  //   syst = sqrt( pow(W*fabs(CDatU->GetBinContent(ip+1) ),2) + pow(VtxCorSyst[ip],2) );
  //   else
  //     syst = sqrt( pow(W*fabs(CDatU->GetBinContent(ip+1) ),2));
  //   cout<<ip<<"   "<<W<<"   "<<syst<<"   "<<CDatU->GetBinContent(ip+1)<<endl;
  //   MCUncert->SetBinError(ip,syst);

  // }

  TGraphAsymmErrors* MCUncert=ComputeSystematics(obs,MC);

  MCUncert->SetMarkerSize(0); 
  MCUncert->SetMarkerStyle(1);
  MCUncert->SetMarkerColor(1);
  MCUncert->SetFillStyle(3001);
   /* MCUncert->SetLineColor(kGreen+2);
   MCUncert->SetLineWidth(Wline);
   MCUncert->SetLineStyle(1);*/
   
   //   MCUncert->SetFillColor(kBlack);
  //MCUncert->SetFillStyle(3254); //3003

  //gStyle->SetErrorX(0.);
    MCUncert->SetFillColor(1); //12
   //MCUncert->SetFillStyle(3002);
   //   MCUncert->Draw("E2same");
   MCUncert->Draw("p E2");
   
}


TGraphAsymmErrors*
UseTree::ComputeSystematics(string obs, TH1F* totMC) {

  systM systsUp = histoManager.FindSysts(obs,"Up");
  systM systsDo = histoManager.FindSysts(obs,"Do");
  systM systs = histoManager.FindSysts(obs);
  
  //First, rebin the histograms ==============
  for(itSystM itS=systs.begin();itS!=systs.end();itS++)
    ((*itS).second)->Rebin(Bin);
  for(itSystM itS=systsUp.begin();itS!=systsUp.end();itS++)
    ((*itS).second)->Rebin(Bin);
  for(itSystM itS=systsDo.begin();itS!=systsDo.end();itS++)
    ((*itS).second)->Rebin(Bin);
  //==========================================
  

  TGraphAsymmErrors* unc=new TGraphAsymmErrors(totMC->GetNbinsX()+2);  
  for(int ib=0;ib<totMC->GetNbinsX()+2;ib++) {
    
    float systU=0;
    float systD=0;

    //first symetric systs
    for(itSystM itS=systs.begin();itS!=systs.end();itS++)  {

      float s = fabs( ((*itS).second)->GetBinContent(ib) - totMC->GetBinContent(ib) );
      
      systU += s*s;
      systD += s*s;
   } //sym
    
    //now asymetric systs
    for(itSystM itS=systsUp.begin();itS!=systsUp.end();itS++)  {

      string n =  (*itS).first.substr( 0, (*itS).first.size() -2 );
      float sU = ((*itS).second)->GetBinContent(ib) - totMC->GetBinContent(ib);
      float sD = (systsDo[ n+"Do" ])->GetBinContent(ib) - totMC->GetBinContent(ib);
      // cout<<totMC->GetBinContent(ib)<<"   "<<((*itS).second)->GetBinContent(ib)<<"   "<<(systsDo[ n+"Do" ])->GetBinContent(ib)<<"  "<<sU<<"   "<<sD<<"   ----->  "<<((*itS).second)->GetEntries()<<"  --> "<<((*itS).second)->GetBinCenter(ib)<<"  "<<((*itS).second)->GetName()<<"   "<<totMC->GetBinCenter(ib)<<"   "<<((*itS).second)->GetXaxis()->GetXmax()<<endl;
      
      systU += sU>0?sU*sU:sD*sD;
      systD += sU<=0?sU*sU:sD*sD;
    } //asym
    
    unc->SetPoint(ib, totMC->GetXaxis()->GetBinCenter(ib), totMC->GetBinContent(ib) );
    
    float width=totMC->GetXaxis()->GetBinWidth(ib);

    unc->SetPointEXlow(ib, width/2. );
    unc->SetPointEXhigh(ib, width/2. );
    unc->SetPointEYhigh(ib, sqrt(systU) );
    unc->SetPointEYlow(ib, sqrt(systD) );
    
  }//bin

  return unc;

}



/*
void UseTree::PlotRMS(string Obs, int sample) {

  c2=new TCanvas("c2","Test",300,300,600,600);
  leg = new TLegend(0.22,0.71,0.49,0.88);
  leg->SetLineColor(0);
  leg->SetFillColor(0);

  vector<string> Observables;
  Observables.push_back("caloT2");
  Observables.push_back("tc");
  Observables.push_back("pf");

  vector<string> Names;
  Names.push_back("caloT2");
  Names.push_back("tc");
  Names.push_back("pf");
  
  vector<int> Colors;
  // Colors.push_back(kViolet+4);
  Colors.push_back(kCyan+1);
  Colors.push_back(kGreen+1);
  Colors.push_back(kOrange-2);
  Colors.push_back(kRed+1); //
  Colors.push_back(kOrange+7);
  Colors.push_back(896);
  Colors.push_back(38);

  
  vector<int> Markers;
  Markers.push_back(20);
  Markers.push_back(21);
  Markers.push_back(29);

  

  vector<TGraphErrors*> RMS;

  for(int unsigned i=0;i<Observables.size();i++) {

    int nvar = histoManager.FindNVar2D((Observables[i]+Obs));
    int Nhisto = histoManager.access2D(nvar, sample);
    YTitle = (histoManager.FindLeg2D( Observables[i]+Obs)).second;
    TGraph* dataG=  (TGraph*)Histos2D[Nhisto]->Clone();
  
    double VarBin2[11] = {0,1,2,3,4,5,6,7,8,9,10};
    TGraphErrors* test = histoManager.ConvertGraphToRMSGraphVarBin(dataG,10,VarBin2,"gregre",nt);
   
    test->SetLineWidth(2);
    test->SetMarkerSize(1.2);
    test->SetMarkerStyle(Markers[i]);
    test->SetLineColor(Colors[i]);
    test->SetMarkerColor(Colors[i]);
    test->GetXaxis()->SetTitle("number of vertices");
    test->GetYaxis()->SetTitle(YTitle.c_str());
    test->GetYaxis()->SetRangeUser(5,30);
    test->GetYaxis()->SetTitleOffset(1.15);
    //test->Draw("samePE");
    RMS.push_back(test);
    if(i==0)
      RMS[i]->Draw("AP");
    else
      RMS[i]->Draw("P");
    
    leg->AddEntry(RMS[i],(Names[i]+" RMS").c_str(),"pl");
    
  }
  leg->Draw("same");
  cmsPrel(Lumi);
}
*/


void UseTree::GetNumbers() {

  int unsigned MaxCut=0;
  vector<string> cNames;

  float globEff=1, globErr=0;

  int nd=nt+1;
  if(NoData)  nd =nt;

  for(int ids=0;ids<nd;ids++) { //datasets
    globEff=1; globErr=0;
    cout<< "Begin efficiencies for "<<name[ids]<<"   "<<effMap[ids].size()<<"  ********* "<<endl;

    if(effTotal[ids].size()>MaxCut)
      MaxCut = effMap[ids].size();

    for(int unsigned ic=0;ic<effMap[ids].size();ic++) { //cuts
  
      if(cNames.size()<ic+1)
	cNames.push_back(effTotal[ids][ic].first);

      //get total value
      float eff = effMap[ids][ic].second / effTotal[ids][ic].second;
      float error = HistoManager::BinomError( TotalN[ids][ic].second, eff);
    
      //get High and low values for systematics (if available)
      float effL = (systMap[ids][ic].second)[0] / effTotal[ids][ic].second;
      float effH = (systMap[ids][ic].second)[1] / effTotal[ids][ic].second;
   
      //FIXME
      if( effL ==0 )
	effL = eff;
      if( effH ==0 )
	effH = eff;
      
      float systL = fabs(eff-effH); //inverted for the good way
      float systH = fabs(eff-effL);


      //==========================================================

      //FIXME
      globEff = effMap[ids][ic].second / effTotal[ids][0].second ;//  * 0.7561*1.031334*0.9902; 
      globErr = HistoManager::BinomError( TotalN[ids][0].second, globEff);// *0.7561*1.0313*0.9902);
      
      cout<<" --> "<<setw(30)<<effMap[ids][ic].first<<"\t  = "<<eff<<" +- "<<error
	  <<" + "<<systH<<" - "<<systL
	  <<" \t\t "<<effMap[ids][ic].second<<"  / "<<effTotal[ids][ic].second<<"  ---> \t = "<<globEff*100<<" +- "<<globErr*100<<endl;
    }

  }

  //Numbers =======================================================

  //For Latex
  for(int unsigned ic=0;ic<MaxCut;ic++) { //cuts
    
    if(ic==0) {
      cout<<" Cut  "<<fixed<<setprecision(3);
      for(int ids=0;ids<nd+1;ids++) {
	if(ids<nt) {
	  cout<<" & "<<name[ids]<<"   ";
	}
	else if(ids==nt){
	  cout<<" & MC   ";
	}
	else
	  {
	    cout<<"  & "<<name[nt]<<"   ";
	  }
      }
      cout<<" \\\\ "<<endl;
    }

    cout<<cNames[ic]<<"    ";

    for(int ids=0;ids<nd+1;ids++) { //datasets
      cout<<" & ";
      if(ids<nt) {
	if(effTotal[ids].size()<ic+1)
	  { cout<<" - "; }
	else {
	  if(effMap[ids][ic].second>0.000001 ) {
	    cout<<effMap[ids][ic].second;
	    cout<<" $\\pm$ "<<sqrt(TotalNWeights[ids][ic].second);
	  }
	  else
	    cout<<" - "; 
	}
	cout<<"  ";
      }
      else if(ids==nt){
	if(TotalMC.size()<ic+1)
	  { cout<<" - "; }
	else {
	  if(TotalMC[ic].second>0.000001 && (TotalMC[ic].second<1000000000) )
	    cout<<TotalMC[ic].second<<" $\\pm$ "<<sqrt(TotalMCError[ic].second);
	  else
	    cout<<" - "; 
	}
      }
      else {
	if(effTotal[nt].size()<ic+1)
	  { cout<<" - "; }
	else {
	  if(effMap[nt][ic].second>0.000001 && (effMap[nt][ic].second<1000000000) )
	    cout<<effMap[nt][ic].second;
	  else
	    cout<<" - "; 
	}
	cout<<"  ";
      }
    }
    cout<<"\\\\"<<endl;
  }



  cout<<endl<<endl;



  //For Visu

  for(int unsigned ic=0;ic<MaxCut;ic++) { //cuts
    
    if(ic==0) {
      cout<<setw(30)<<" Cut"<< fixed<<setprecision(3)<<"\t";
      for(int ids=0;ids<nd+1;ids++) {
	if(ids<nt) {
	  cout<<name[ids]<<"\t";
	}
	else if(ids==nt){
	  cout<<"MC\t";
	}
	else
	  {
	    cout<<name[nt]<<"\t";
	  }
      }
      cout<<endl;
    }

    cout<<setw(30)<<cNames[ic]<<"\t";

    for(int ids=0;ids<nd+1;ids++) { //datasets
      
      if(ids<nt) {
	if(effTotal[ids].size()<ic+1)
	  { cout<<"-"; }
	else {
	  if(effMap[ids][ic].second>0.000001 )
	    cout<<effMap[ids][ic].second;
	  else
	    cout<<"-"; 
	}
	cout<<"\t";
      }
      else if(ids==nt){
	if(TotalMC.size()<ic+1)
	  { cout<<"-"; }
	else {
	  if(TotalMC[ic].second>0.000001 && (TotalMC[ic].second<100000000000) )
	    cout<<TotalMC[ic].second;
	  else
	    cout<<"-"; 
	}
	cout<<"\t";
      }
      else {
	if(effTotal[nt].size()<ic+1)
	  { cout<<"-"; }
	else {
	  if(effMap[nt][ic].second>0.000001 && (effMap[nt][ic].second<100000000000) )
	    cout<<effMap[nt][ic].second;
	  else
	    cout<<"-"; 
	}
	cout<<"\t";
      }
    }
    cout<<"\t"<<effMap[0][ic].second/sqrt(TotalMC[ic].second)<<endl;
  }

}


void UseTree::SetEfficiency(int i, string vkind, float w, bool acc ) {

  if(histoManager.GetInitStatus()) return;

  if((int)effTotal.size() != i+1)
    {
      vector<std::pair<string,float> > tmp;
      effTotal.push_back(tmp);
      TotalN.push_back(tmp);
      TotalNWeights.push_back(tmp);
    }

  int n=-1;
  for(int unsigned k=0;k<effTotal[i].size();k++) {
    
    if(effTotal[i][k].first  == vkind )
      { n=k; break; }
  }

  //  Acceptance
    if(useAccForEff)
      if(!inAcc) w=0;

  if(n==-1) {
    std::pair<string,float> tmp(vkind,w); //w
    effTotal[i].push_back(tmp);
    std::pair<string,float> tmp2(vkind,1); //w
    TotalN[i].push_back(tmp2);
    std::pair<string,float> tmp3(vkind,w*w); //w
    TotalNWeights[i].push_back(tmp3);
  }
  else {
    effTotal[i][n].second += w;  //w
    TotalN[i][n].second += 1;  //w
    TotalNWeights[i][n].second += w*w;  //w
  }

  //Only passing events
  
 
  
  if((int)effMap.size() != i+1)
    {
      vector<std::pair<string,float> > tmp;
      effMap.push_back(tmp);
    }

  n=-1;
  for(int unsigned k=0;k<effMap[i].size();k++) {
    
    if(effMap[i][k].first  == vkind )
      { n=k; break; }
  }
  if(acc) {
    if(n==-1) {
      std::pair<string,float> tmp(vkind,w); //w
      effMap[i].push_back(tmp);
    }
    else {
      effMap[i][n].second += w; //w
    }
  }
  else{
    if(n==-1) {
      std::pair<string,float> tmp(vkind,0); //not passed
      effMap[i].push_back(tmp);
    }
    else {
      effMap[i][n].second += 0;  //not passed
    }
  }

  if(i==nt) return;

  n=-1; 
  for(int unsigned k=0;k<TotalMC.size();k++) {
    
    if(TotalMC[k].first  == vkind )
      { n=k; break; }
  }
  int n2=-1;
  for(int unsigned k=0;k<TotalMCError.size();k++) {
    if(TotalMCError[k].first  == vkind )
      { n2=k; break; }
  }

  if(acc) {
    if(n==-1) {
      std::pair<string,float> tmp(vkind,w); //w
      TotalMC.push_back(tmp);
    }
    else {
      TotalMC[n].second += w; //w
    }
    if(n2==-1) {
      std::pair<string,float> tmp2(vkind,w*w); //w
      TotalMCError.push_back(tmp2);
    }
    else {
      TotalMCError[n2].second += w*w; //w
    }
  }
  else {
    if(n==-1) {
      std::pair<string,float> tmp(vkind,0); //w
      TotalMC.push_back(tmp);
    }
    else {
      TotalMC[n].second += 0; //w
    }
    if(n2==-1) {
      std::pair<string,float> tmp2(vkind,0); //w
      TotalMCError.push_back(tmp2);
    }
    else {
      TotalMCError[n2].second += 0; //w
    }
  }
  
}

void UseTree::SetSystematics(string type, int i, string cName, float w) {

  if(histoManager.GetInitStatus()) return;

  bool accept[2] = {false,false};
 
  if((int)systMap.size() != i+1 )
    {
      vector<std::pair<string, vector<float> > > tmp;
      systMap.push_back(tmp); 
      //if(i==2)    
      //cout<<systMap.size()<<endl;
    }
  
    int n=-1;
    for(int k=0;k<(int)systMap[i].size();k++) {
    
      if(systMap[i][k].first  == cName )
	{ n=k; break; }
    }
    
    if(n==-1) {
      vector<float> vtmp;
      for(int j=0;j<2;j++) {
	float val=w;
	if(!accept[j]) val =0;
	vtmp.push_back(val);
      }
      std::pair<string,vector<float> > tmp(cName,vtmp); //w
      systMap[i].push_back(tmp);
    }
    else {
      for(int j=0;j<2;j++) {
	float val=w;
	if(!accept[j]) val =0;
	(systMap[i][n].second)[j] += val; //w
      }
    }

}

// inline template < typename T >
// bool UseTree::MakeCut( T value, T valcut, string type, int i, string cName, float w, T seccut) {

//   bool accept;
  
//   if(!Cbool(skipCut, cName==Nm1Var) )
//     { return true; }

//   if( Cbool( invCut, cName==Nm1Var) ) {
//     type = InvCut(type);
//   }
 
//   if(type=="<") {
//     accept = (value < valcut);
//   }
//   else if(type=="<=") {
//     accept = (value <= valcut);
//   }
//   else if( type==">") {
//     accept = (value > valcut);
//   }
//   else if( type==">=") {
//     accept = (value >= valcut);
//   }
//   else if( type=="=") {
//     accept = (value == valcut);
//   }
//   else if(type=="!=") {
//     accept = (value != valcut);
//   }
//   else if(type=="[]") {
//     accept = (value >= valcut && value<= seccut );
//   }
//   else if(type=="][") {
//     accept = (value > valcut && value< seccut );
//   }
//   else if(type=="[!]") {
//     accept = !(value >= valcut && value <= seccut );
//   }
//   else if(type=="]![") {
//     accept = !(value > valcut && value < seccut );
//   }

//   else {
//     accept =false; cout<<" Warning cut "<<endl;
//   }  
  
//   SetSystematics( type, i, cName, w);
//   SetEfficiency(i, cName, w, accept);
 
//   return accept;
// }

string 
UseTree::InvCut(string type) {
  if(type=="<") {
    return ">=";
  }
  else if(type=="<=") {
   return ">";
  }
  else if( type==">") {
   return "<=";
  }
  else if( type==">=") {
    return "<";
  }
  else if( type=="=") {
   return "!=";
  }
  else if(type=="!=") {
   return "=";
  }
  else if(type=="[]") {
   return "]![";
  }
  else if(type=="][") {
    return "[!]";
  }
  else if(type=="[!]") {
     return "][";
  }
  else if(type=="]![") {
    return "[]";
  }
  else {
    cout<<" Warning cut "<<endl;
    return "=";
  }
} 


vector<TPad*>
UseTree::PreparePads(int npad) {

  float m=1.3;

  size_t n_ = npad;
  vector<TPad*> pad_;
  Double_t e_ = 15*m; //15
  Double_t t_ = 30*m; //30
  Double_t b_ = 80*m; //80
  Double_t g_ = 87*m; //87
  Double_t d_ = 30*m; //30
  Double_t h_ = 375*m; //400
  Double_t w_ = 350*m; //350
  Double_t ww_ = g_ + w_ + e_ ;
  Double_t W_ = g_ + n_*w_ + 2*(n_-1)*e_ + d_;
  Double_t H_ = t_ + h_ + b_ ;
  
  c2=new TCanvas("c2","Test",300,300,W_,H_);
  //  c2=new TCanvas("c2","Test",300,300,1.2*W_,1.2*H_);
  Double_t xlow_= 0;
  Double_t ylow_=0.;
  Double_t xup_=0;
  Double_t yup_=1.;
  for(int unsigned i=0;i<n_;i++) {
      
    TString padname_("pad_");
    padname_+=i;
    xup_ = xlow_ + ww_ /W_;
      
    TPad* padtmp = new TPad( padname_, padname_, 
			  xlow_, ylow_, xup_, yup_,
			  kWhite,0,0);
    xlow_ += (w_ + 2*e_)/W_;
    padtmp->SetLeftMargin(  g_/ww_ );
    padtmp->SetRightMargin( e_/ww_ );
    padtmp->SetTopMargin(  t_/H_ );
    padtmp->SetBottomMargin( b_/H_ );
    padtmp->SetFillColor(0);
    padtmp->SetTickx(1);
    padtmp->SetTicky(1);
    padtmp->SetFrameFillStyle(0);
    padtmp->SetFrameLineWidth(Wline);
    padtmp->SetFrameBorderMode(0);
      
    pad_.push_back(padtmp);
	
  }

  return pad_;

}




vector<vector<TPad*> >
UseTree::PreparePadsWithRatio(int npad) {

  float m=1.3;

  size_t n_ = npad;
  vector<TPad*> padhigh_;
  vector<TPad*> padlow_;
  Double_t e_ = 15*m; //15
  Double_t k_ = 7*m; //15
  Double_t t_ = 30*m; //30
  Double_t b_ = 50*m; //80
  Double_t g_ = 87*m; //87
  Double_t d_ = 30*m; //30
  Double_t h_ = 350*m; //400
  Double_t w_ = 350*m; //350
  Double_t hl_ = 70*m;

  Double_t ww_ = g_ + w_ + e_ ;
  Double_t W_ = g_ + n_*w_ + 2*(n_-1)*e_ + d_;
  Double_t H_ = t_ + h_ + b_ + hl_ +2*k_ ;
    
  c2=new TCanvas("c2","Test",300,300,W_,H_);
  //  c2=new TCanvas("c2","Test",300,300,1.2*W_,1.2*H_);
  Double_t xlow_= 0;
  Double_t ylow_=0.;
  Double_t xup_=0;
  Double_t yup_=1.;
  for(int unsigned i=0;i<n_;i++) {
      
    TString padname_("padhigh_");
    padname_+=i;
    xup_ = xlow_ + ww_ /W_;
    yup_ = 1.; 
    ylow_ = (hl_ + b_ + k_ ) /H_;
  
      
    TPad* padtmp = new TPad( padname_, padname_, 
			  xlow_, ylow_, xup_, yup_,
			  kWhite,0,0);
    xlow_ += (w_ + 2*e_)/W_;
    padtmp->SetLeftMargin(  g_/ww_ );
    padtmp->SetRightMargin( e_/ww_ );
    padtmp->SetTopMargin(  t_/H_ );
    padtmp->SetBottomMargin( k_/H_ );
    padtmp->SetFillColor(0);
    padtmp->SetTickx(1);
    padtmp->SetTicky(1);
    padtmp->SetFrameFillStyle(0);
    padtmp->SetFrameLineWidth(Wline);
    padtmp->SetFrameBorderMode(0);
      
    padhigh_.push_back(padtmp);
	
  }
  cout<<" second "<<endl;
  xlow_= 0; //reset
  c2->cd(); //reset
  for(int unsigned i=0;i<n_;i++) {

    TString padname_("padlow_");
    padname_+=i;
    xup_ = xlow_ + ww_ /W_;
    yup_ = (hl_ + b_ + k_ ) /H_; 
    ylow_ = 0;


    TPad* padtmp  = new TPad( padname_, padname_, 
			      xlow_, ylow_, xup_, yup_,
			      kWhite,0,0);
    xlow_ += (w_ + 2*e_)/W_;
    padtmp->SetLeftMargin(  g_/ww_ );
    padtmp->SetRightMargin( e_/ww_ );
    padtmp->SetTopMargin( k_/H_ );
    padtmp->SetBottomMargin( b_ /(hl_ + b_ + k_ ) );
    padtmp->SetFillColor(0);
    padtmp->SetTickx(1);
    padtmp->SetTicky(1);
    padtmp->SetFrameFillStyle(0);
    padtmp->SetFrameLineWidth(Wline);
    padtmp->SetFrameBorderMode(0);

    padlow_.push_back(padtmp);
  } 
  

 vector<vector<TPad*> > pad_;
 pad_.push_back(padhigh_);
 pad_.push_back(padlow_);
  
 return pad_;

}


void UseTree::Integral(string Obs, double xmin, double xmax) {

  //Get the MC
  TH1F* mc = (TH1F*)savedMC->Clone();
 
  //Get the data
  int Sdata = histoManager.FindNVar(Obs);
  int ndata= histoManager.access(Sdata,nt);
  TH1F* data = (TH1F*)Histos[ ndata ]->Clone();
  data->Rebin(Bin);

  int Binm = data->FindBin( xmin );
  int BinM = data->FindBin( xmax );

  double integMC = mc->Integral(Binm,BinM);
  double integdata = data->Integral(Binm,BinM);

  cout<<endl<<" *** Integral "<<xmin<<" -> "<<xmax<<" ---> "<<Binm<<"   "<<BinM<<" : MC -> "<<integMC<<"   ; data ->  "<<integdata<< " +- "<<sqrt(integdata)<<endl;

}


void
UseTree::PlotDataMCRatio(string Obs, bool removelabel) {
 
  //Get the data
  int Sdata = histoManager.FindNVar(Obs);
  int ndata= histoManager.access(Sdata,nt);
  TH1F* data = (TH1F*)Histos[ ndata ]->Clone();
 
  //The MC
  TH1F* MC = (TH1F*)savedMC->Clone();
 
  //The bidon histo
  vector<float> templates = histoManager.GetTemplate(Sdata);
  TH1F* emptyHisto= new TH1F(("bidonlow"+Obs).c_str(),"bidonlow",(int)templates[0],templates[1],templates[2]);
  
  TH1F* ratio = (TH1F*)data->Clone();
  ratio->SetName( ("ratio"+Obs).c_str());
  ratio->Rebin(Bin);
  ratio->Divide(MC);
 
  for(int ib=0;ib<emptyHisto->GetNbinsX()+2;ib++)
    { emptyHisto->SetBinContent(ib,-1000);   }

  // double max = ratio->GetMaximum();
  // double min = ratio->GetMaximum();

  // float dm = 1-min;
  // float dM = max - 1;
  // float rangemax = 1;
  // if(dm>dM)
  //     rangemax = 1.1*dm;
  //   else
  //     rangemax = 1.1*dM;
      
  float xmin = RangeX[0];
  float xmax = RangeX[1];

  

  if(removelabel)  {
    emptyHisto->GetYaxis()->SetTitleOffset(1000);
    ratio->GetYaxis()->SetTitleOffset(1000);
  }
  
  if(!OverFlowBin && !UnderFlowBin) {
    emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1] );
    ratio->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1] );
  }
  else {
    if(OverFlowBin) {
      emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1] + ratio->GetBinWidth(1) -0.001);
      ratio->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1] + ratio->GetBinWidth(1) -0.001);
      xmax = RangeX[1] + ratio->GetBinWidth(1) -0.001;
    }
    if(UnderFlowBin) {
      emptyHisto->GetXaxis()->SetRangeUser(RangeX[0] - ratio->GetBinWidth(1) ,RangeX[1] + ratio->GetBinWidth(1) -0.001);
      ratio->GetXaxis()->SetRangeUser(RangeX[0] - ratio->GetBinWidth(1) ,RangeX[1] + ratio->GetBinWidth(1) -0.001);
      xmin = RangeX[0] - ratio->GetBinWidth(1);
    }
  }
 
  string Xtitle = ((histoManager.FindLeg( Obs )).first).c_str();
  //cout<<" ---------> "<<emptyHisto->GetXaxis()->GetXmax()<<endl;
  if(xmax>=1000) {
    if(Xtitle.find("GeV") != (size_t)-1 ) {

      Xtitle.replace( Xtitle.find("GeV") , 3, "TeV" );
    }
  }

  if(xmin < emptyHisto->GetXaxis()->GetXmin() )
    xmin = emptyHisto->GetXaxis()->GetXmin();
  if(xmax >  emptyHisto->GetXaxis()->GetXmax() )
    xmax = emptyHisto->GetXaxis()->GetXmax();
  
  //if systematics exists, draw them!
  TPolyLine* sysBand;
  if(AddSyst) {
    TGraphAsymmErrors* MCUncert=ComputeSystematics(Obs,MC);
    double x,y/*,xl,xh*/,yl,yh;
  
    vector<float> xs;
    vector<float> yls;
    vector<float> yhs;
    
    for(int ib=0;ib<MCUncert->GetN();ib++) {
      MCUncert->GetPoint(ib,x,y);
      //   MCUncertRatio->SetPoint(ib, x, 1 );

      yl = MCUncert->GetErrorYlow(ib);
      yh = MCUncert->GetErrorYhigh(ib);
  
      if(x<xmin-ratio->GetXaxis()->GetBinWidth(1)/2. || x>xmax+ratio->GetXaxis()->GetBinWidth(1)/2. ) continue;
      
      xs.push_back(x);
      if(y!=0) {
	yls.push_back(yl/y);
	yhs.push_back(yh/y);
      }
      else {
	yls.push_back(0);
	yhs.push_back(0);
      }
    }
    
    sysBand = GetSystBand(xs,yls,yhs,xmin,xmax);
  }

  emptyHisto->GetYaxis()->SetRangeUser( 0.4, 1.6);// 1-rangemax, 1+rangemax );
  emptyHisto->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
  emptyHisto->GetYaxis()->SetNdivisions(3,Ydiv[1],Ydiv[2]);
  emptyHisto->GetXaxis()->SetTitle( Xtitle.c_str() );
  emptyHisto->GetYaxis()->SetTitle( "Data/MC" );
  emptyHisto->GetXaxis()->SetTitleSize(0.20);
  emptyHisto->GetXaxis()->SetTitleOffset(0.83);
  emptyHisto->GetXaxis()->SetLabelSize(0.165);
  emptyHisto->GetYaxis()->SetLabelSize(0.14);
  emptyHisto->GetYaxis()->SetLabelOffset(0.011);
  emptyHisto->GetYaxis()->SetTitleSize(0.15);
  emptyHisto->GetYaxis()->SetTitleOffset(0.54);
  emptyHisto->GetXaxis()->SetTickLength(0.09);
  emptyHisto->GetYaxis()->SetTickLength(0.05);
  
  ratio->GetYaxis()->SetRangeUser(  0.4, 1.6);//1-rangemax, 1+rangemax );
  ratio->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
  ratio->GetYaxis()->SetNdivisions(3,Ydiv[1],Ydiv[2]);
  ratio->GetXaxis()->SetTitle( Xtitle.c_str() );
  ratio->GetYaxis()->SetTitle( "Data/MC" );
  ratio->GetXaxis()->SetTitleSize(0.20);
  ratio->GetXaxis()->SetTitleOffset(0.83);
  ratio->GetXaxis()->SetLabelSize(0.165);
  ratio->GetYaxis()->SetLabelSize(0.14);
  ratio->GetYaxis()->SetLabelOffset(0.015);
  ratio->GetYaxis()->SetTitleSize(0.15);
  ratio->GetYaxis()->SetTitleOffset(0.54);
  ratio->GetXaxis()->SetTickLength(0.09);
  ratio->GetYaxis()->SetTickLength(0.05);

 if(removelabel)  {
   ratio->GetYaxis()->SetLabelOffset(1000);
  emptyHisto->GetYaxis()->SetLabelOffset(1000);
 }
   //  cout<<" *** Plot Ratio   "<< 1-rangemax<<"   " <<1+rangemax<<endl;

  emptyHisto->Draw();
  if(AddSyst)
    sysBand->Draw("F");
  ratio->Draw("same");


  TLine* line=new TLine(xmin,1,xmax,1);
  line->SetLineColor(kRed+1);
  line->SetLineStyle(7);
  line->SetLineWidth(2);

  line->Draw("same");

  //Blanco
  //TPave* pv = new TPave(0.9722,-0.126432,1.00676,0.378049,4,"NDC");
  TPave* pv = new TPave(1008.95,0.42,1053.69,0.67);
  pv->SetFillColor(0);
  pv->SetShadowColor(0);
  pv->SetLineColor(0);
  pv->Draw();
  

}




void UseTree::PlotRatio(string var1, string var2, string ds) {
  
  int nds=-1;
  for(int i=0;i<nt+1;i++) {
    if(ds==name[i])
      {nds=i; break;}
  }

  //Get the first var
  int Sh1 = histoManager.FindNVar(var1);
  int nh1= histoManager.access(Sh1,nds);
  TH1F* h1 = (TH1F*)Histos[ nh1 ]->Clone();
  
  int Sh2 = histoManager.FindNVar(var2);
  int nh2= histoManager.access(Sh2,nds);
  TH1F* h2 = (TH1F*)Histos[ nh2 ]->Clone();

  cout<<"Observables -> "<<Sh1<<"   "<<Sh2<<endl;

  //Rebinning
  h1->Rebin(Bin);
  h2->Rebin(Bin);

  string Xtitle =  (histoManager.FindLeg( var1 )).first;
  //Normalize plots
  float Int1 = h1->Integral(0, h1->GetNbinsX()+1 );
  float Int2 = h2->Integral(0, h2->GetNbinsX()+1 );

  h1->Scale(1./Int1);
  h2->Scale(1./Int2);
  
  TCanvas*  cRatio = new TCanvas("cRat","Ratio",300,300,485,485);
  cRatio->SetLeftMargin(  87./485 );  //87
  cRatio->SetRightMargin( 42./485 ); //42
  cRatio->SetTopMargin(  30./485 ); //30
  cRatio->SetBottomMargin( 80./485 ); //80
  cRatio->SetFillColor(0);
  cRatio->SetTickx(1);
  cRatio->SetTicky(1);
  cRatio->SetFrameFillStyle(0);
  cRatio->SetFrameLineWidth(4);
  cRatio->SetFrameBorderMode(0);

  TH1F* ratio = (TH1F*)h1->Clone();
  ratio->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);

  double xm = RangeX[0];
  double xM = RangeX[1];
  if((RangeX[0]<ratio->GetXaxis()->GetXmin()) )
    xm=ratio->GetXaxis()->GetXmin();
  if((RangeX[1]>ratio->GetXaxis()->GetXmax()) )
    xM=ratio->GetXaxis()->GetXmax();

  TLine* line = new TLine(xm,1,xM,1);
  line->SetLineColor(2);
  line->SetLineWidth(2);
  ratio->GetYaxis()->SetTitle(" ratio");
  ratio->GetXaxis()->SetTitle(Xtitle.c_str() );
  ratio->GetXaxis()->SetNdivisions(5,5,0);
  ratio->Divide(h2);
  ratio->Scale( 1 );
  ratio->SetMarkerSize(MarkerSize);
  ratio->Draw();
  line->Draw("same");

  cout<<" Ratio on "<<Xtitle<<" made "<<endl;

  float Bilan=0;
  for(int i=0;i<ratio->GetNbinsX()+2;i++)
    Bilan += ratio->GetBinContent(i);

  cout<<" Bilan : "<<Bilan<<endl;

}




void
UseTree::DrawNumbers() {

 
  
  int unsigned MaxCut=0;
 vector<string> cNames;
 
  for(int ids=0;ids<nt+1;ids++) { //datasets
    if(effTotal[ids].size()>MaxCut)
      MaxCut = effMap[ids].size();

    for(int unsigned ic=0;ic<effMap[ids].size();ic++) { //cuts
      
      if(cNames.size()<ic+1)
	cNames.push_back(effTotal[ids][ic].first);
    }
    
  }


  TLegend* legend = new TLegend(0.7,0.5,0.9,0.8);
  legend->SetFillColor(0);
  legend->SetShadowColor(0);
  legend->SetLineColor(0);

  vector<TH1F*> HyMC;
  TH1F* HyMCError= new TH1F("ErrorMCyield","yieldError",MaxCut,0,MaxCut);
  TH1F* HyData= new TH1F("yieldData","yieldD",MaxCut,0,MaxCut);

  //Fill MC yields
  for(int ids=0;ids<nt;ids++) { //datasets
    TH1F* tmp = new TH1F( ("yield"+name[ids]).c_str(),("yield"+name[ids]).c_str(),MaxCut,0,MaxCut);
    
    for(int unsigned ic=0;ic<MaxCut;ic++) { //cuts
      if(effTotal[ids].size()>=ic+1 && 
	 effMap[ids][ic].second>0.000001 && 
	 effMap[ids][ic].second<100000000000 )
	{
	  
	  tmp->SetBinContent(ic+1, effMap[ids][ic].second);
	}
      else {
	tmp->SetBinContent(ic+1, 0 );
      }
    }
    
    //Colors
    tmp->SetFillColor(colors[ids]);
    tmp->SetLineColor(colors[ids]);
    
    //sum
    if(ids==0) {
      HyMC.push_back(tmp);
    }
    else{
      tmp->Add(HyMC[ids-1]);
      HyMC.push_back(tmp);
    }

    legend->AddEntry(tmp,name[ids].c_str(),"f");

  }
 
  //Fill Data yields
  for(int unsigned ic=0;ic<MaxCut;ic++) { //cuts
    
    if(effTotal[nt].size()>=ic+1) {
	if(effMap[nt][ic].second>0.000001 && 
	   effMap[nt][ic].second<100000000000 )
	  {
	    HyData->SetBinContent(ic+1, effMap[nt][ic].second);
	  }
	else {
	  HyData->SetBinContent(ic+1, 0 );
	}
      }
      else {
	HyData->SetBinContent(ic+1, 0 );
      }

    //And uncertainties of MC
    HyMCError->SetBinContent(ic+1,HyMC[nt-1]->GetBinContent(ic+1));
    HyMCError->SetBinError(ic+1, sqrt(TotalMCError[ic].second) );

    cout<<TotalMCError[ic].first<<"    "<<TotalMCError[ic].second<<"    "<<HyMCError->GetBinError(ic+1)<<endl;
  }

  //Convert into a rooHist
  RooHist* roohist =new RooHist((*HyData));
  TGraphAsymmErrors* DataGraph(0);
  if(HyData->GetEntries()!=0) {
    int Nn0=0;
    vector<double> vY;
    vector<double> vX;
    vector<double > veY;
    vector<double > veX;
    vector<double> tmp(0,2);
  
    for(int ip=0;ip<roohist->GetN();ip++) {
      double Y,X;
      
      roohist->GetPoint(ip,X,Y);
    
      if(Y!=0) {
	Nn0++;
      
	vY.push_back(Y);
	vX.push_back(X);
	veX.push_back( roohist->GetErrorXlow(ip) );
	veX.push_back( roohist->GetErrorXhigh(ip) );
	veY.push_back( roohist->GetErrorYlow(ip) );
	veY.push_back( roohist->GetErrorYhigh(ip) );
      }
    }
  
    DataGraph=new TGraphAsymmErrors(Nn0);
    for(int ip=0;ip<Nn0;ip++) {
      DataGraph->SetPoint(ip,vX[ip],vY[ip]);
      DataGraph->SetPointError(ip,0,0/*veX[ip*2],veX[ip*2+1]*/,veY[ip*2],veY[ip*2+1]);
      //	cout<<" Errors "<<ip<<"  "<<vX[ip]<<"  --->  "<<veY[ip*2]<<"   "<<veY[ip*2+1]<<endl;
    }
    if(RangeX[1]> 110)
      DataGraph->SetMarkerSize(MarkerSize);//0.9
    else
      DataGraph->SetMarkerSize(MarkerSize); 
    DataGraph->SetMarkerStyle(20);
    DataGraph->SetMarkerColor(1);
    DataGraph->SetLineColor(1);
    DataGraph->SetLineWidth(LineWidth);
      
    double x,y;
    for(int k=0;k<DataGraph->GetN();k++)
      {
	  
	DataGraph->GetPoint(k,x,y);
	if(y <= RangeY[0])  {
	  DataGraph->SetPointError(k,0,0,0,0);
	}
      }

    legend->AddEntry(DataGraph,"data","pl");
  }
  
  //2 sigma sideband
  TGraphAsymmErrors* Sig2band=(TGraphAsymmErrors*) DataGraph->Clone();
  for(int ip=0;ip<DataGraph->GetN();ip++) {
    double Y,X;
    Sig2band->GetPoint(ip,X,Y);

    double yH = ErrorPH(Y,2);
    double yL = ErrorPL(Y,2);
    
    Sig2band->SetPointError(ip,0.2,0.2,yL,yH);
  }
  
  Sig2band->SetFillColor(kAzure+2);
  Sig2band->SetFillStyle(3001);
  Sig2band->SetLineColor(kAzure+2);
  Sig2band->SetLineWidth(2);

  TCanvas* cY = new TCanvas("Yields","Yields",1400,600);
  TH1F* emptyH= new TH1F("emptyY","emptyY",MaxCut,0,MaxCut);

  emptyH->SetFillColor(0);
  emptyH->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
  emptyH->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
  emptyH->GetYaxis()->SetNdivisions(Ydiv[0],Ydiv[1],Ydiv[2]);
  emptyH->GetYaxis()->SetTitle(YTitle.c_str());
  for(int unsigned ib=0;ib<MaxCut;ib++) {
    emptyH->GetXaxis()->SetBinLabel(ib+1, cNames[ib].c_str() );
  }

  emptyH->Draw();
  for(int ids=nt-1;ids>-1;ids--) {
    HyMC[ids]->Draw("same hist");
  }
  HyMCError->SetMarkerStyle(1);
  HyMCError->SetMarkerColor(1);
  HyMCError->SetLineColor(1);
  HyMCError->SetFillStyle(3002);
  HyMCError->SetFillColor(12);
  HyMCError->Draw("Sames:E2");
  // Sig2band->Draw("PZ 2");
  DataGraph->Draw("PZ"); 
  // HyData->Draw("same");
  
  cY->RedrawAxis();
  legend->Draw("same");

  
}




float 
UseTree::ErrorPH(int N, int sigma) {

  float alpha;

  switch (sigma)
    {
    case 1 : alpha = 0.3173; break;
    case 2 : alpha = 0.0455; break;
    case 3 : alpha = 0.0027; break;
    case 4 : alpha = 0.000063; break;
    case 5 : alpha = 0.00000057; break;
    default: {cout<<" By default, sigma=1"<<endl; alpha = 0.3173; break;}
    }
  
  float eh =  TMath::ChisquareQuantile(1- alpha/2,2*(N +1))*0.5 - N;

  if(N == 0)
    eh =  TMath::ChisquareQuantile(1- alpha,2*(N +1))*0.5  - N;
  
  return eh;
  
}

float 
UseTree::ErrorPL(int N, int sigma) {

  float alpha;

  switch (sigma)
    {
    case 1 : alpha = 0.3173; break;
    case 2 : alpha = 0.0455; break;
    case 3 : alpha = 0.0027; break;
    case 4 : alpha = 0.000063; break;
    case 5 : alpha = 0.00000057; break;
    default: {cout<<" By default, sigma=1"<<endl; alpha = 0.3173; break;}
    }
 
  float el = N - TMath::ChisquareQuantile(alpha/2,2*N)*0.5;

  return el;
  
}




float 
UseTree::dR(float e1, float e2, float p1, float p2, bool conv) {

  float p1r,p2r;
  //Convert phi from ∞ to rad
  if( conv ) {
    p1r = p1*TMath::Pi()/180;
    p2r = p2*TMath::Pi()/180;
  }
  else {
    p1r= p1;
    p2r= p2;
  }
    float dp2 = pow( dPhi(p1r,p2r) , 2);
  
  float de2 = pow(e1-e2,2);  

  return sqrt(de2 + dp2 );

}

float 
UseTree::phi( float x, float y )
{
  float phi_ =atan2(y, x);
  return (phi_>=0) ?  phi_ : phi_ + 2*TMath::Pi();
}

float 
UseTree::dPhi( float phi1, float phi2 )
{
  float phi1_= phi( cos(phi1), sin(phi1) );
  float phi2_= phi( cos(phi2), sin(phi2) );
  float dphi_= phi1_-phi2_;
  if( dphi_> TMath::Pi() ) dphi_-=2*TMath::Pi();
  if( dphi_<-TMath::Pi() ) dphi_+=2*TMath::Pi();

  return dphi_;
}



string
UseTree::FindProcess(int ievent, int ids) {

  if(ids==nt) return "   ";

  for(int unsigned j=1;j<ChanNum[ids].size();j++) {
    if( ievent  < ChanNum[ids][j].second  &&  ievent >= ChanNum[ids][j-1].second) {
      return ChanNum[ids][j-1].first;
    }
  }
  
  return "  ";
}




TPolyLine* 
UseTree::GetSystBand(vector<float> xs, vector<float> yl, vector<float> yh, float xmin, float xmax) {
  

  int nbin = xs.size()-1;
  



  vector< float > vecx_(8*nbin,-100000);
  vector< float > vecyl_(8*nbin,-100000);
  vector< float > vecy_(8*nbin,-100000);
  float lsyst;
  
  bool findM=false,findm=false;

  for( int ibin=0; ibin<nbin; ibin++ )
    {
      //  if(xmax > 1.0001*(xs[ibin]+xs[ibin+1])/2. && xmin < 0.9999*(xs[ibin]+xs[ibin+1])/2.) {
    
      //cout<<(xs[ibin]+xs[ibin+1])/2.<<"  "<<xs[ibin]<<endl;

      //if( (xs[ibin]*3+xs[ibin+1])/4.>xmin)
	vecx_[4*ibin]   = (xs[ibin]*3+xs[ibin+1])/4.;

      vecx_[4*ibin+1]   = (xs[ibin]+xs[ibin+1])/2.;

      //      if( (xs[ibin]+xs[ibin+1]*3)/4.>xmin && (xs[ibin]+xs[ibin+1]*3)/4.<xmax)
	vecx_[4*ibin+2]   = (xs[ibin]+xs[ibin+1]*3)/4.;
      
	//if (xs[ibin+1]>xmin && xs[ibin+1]<xmax)
	vecx_[4*ibin+3]   = xs[ibin+1];

      // if((xs[ibin]+xs[ibin+1])/2. == xmax )
      // 	findM=true;
      // if( (xs[ibin]+xs[ibin+1])/2. == xmin )
      // 	findm=true;

      if(findM) {
	vecx_[ibin] = vecx_[ibin-1];
	findM=false;
      }
      if(findm && ibin!=0) {
	vecx_[ibin] = vecx_[0];
	findm=false;
      }

      vecy_[4*ibin]   = yh[ibin];
      vecy_[4*ibin+1]  = (yh[ibin]+yh[ibin+1])/2.;
      vecy_[4*ibin+2]  = yh[ibin+1];//(ibin==(nbin-1))?yh[ibin+1]:(yh[ibin+1]+yh[ibin+2])/2.; //// (yh[ibin+1]+yh[ibin+2])/2.;
      vecy_[4*ibin+3]  = yh[ibin+1];
      vecyl_[4*ibin]  = yl[ibin];
      vecyl_[4*ibin+1]  = (yl[ibin]+yl[ibin+1])/2.;
      vecyl_[4*ibin+2]  = yl[ibin+1];// (ibin==(nbin-1))?yl[ibin+1]:(yl[ibin+1]+yl[ibin+2])/2.;
      vecyl_[4*ibin+3]  = yl[ibin+1];//yl[ibin+1];
      //lsyst = yh[ibin];
      // }
      // else {
      // 	vecx_[2*ibin]   = xmax;
      // 	vecy_[2*ibin]   = lsyst;
      // }

	//if(xmax > 0.9999*(xs[ibin+1]+xs[ibin+2])/2.) {
      // vecx_[2*ibin+1] = (xs[ibin+1]+xs[ibin+2])/2.;
      // vecy_[2*ibin+1] = yh[ibin];
      // vecyl_[2*ibin+1] = yl[ibin];
      //	lsyst = yl[ibin];
      // } else {
      // 	vecx_[2*ibin+1] = xmax;
      // 	vecy_[2*ibin+1] = lsyst;
      // }
	
	// cout<<ibin<<"   "<<1.0001*(xs[ibin]+xs[ibin+1])/2.
	//     <<"   "<<yh[ibin]<<" -> "<<0.9999*(xs[ibin+1]+xs[ibin+2])/2.<<"  "<<yh[ibin]<<endl;

      //vecx_[2*nbin-1-ibin]   = vecx_[ibin+1];
      vecx_[8*nbin-4*ibin-1] = vecx_[4*ibin];
      vecx_[8*nbin-4*ibin-2] = vecx_[4*ibin+1];
      vecx_[8*nbin-4*ibin-3] = vecx_[4*ibin+2];
      vecx_[8*nbin-4*ibin-4] = vecx_[4*ibin+3];
      //vecy_[2*nbin-1-ibin]   = -vecyl_[ibin+1];
      vecy_[8*nbin-4*ibin-1] = -vecyl_[4*ibin];
      vecy_[8*nbin-4*ibin-2] = -vecyl_[4*ibin+1];
      vecy_[8*nbin-4*ibin-3] = -vecyl_[4*ibin+2];
      vecy_[8*nbin-4*ibin-4] = -vecyl_[4*ibin+3];
    }
  
  // vecx_[2*nbin] = vecx_[0]; //vecy_[2*nbin] = vecy_[2*nbin-1];
  // vecx_[2*nbin+1] = vecx_[0]; 
  // vecy_[2*nbin+1] = vecyl_[0];
  
  int n=0;
  for( size_t ipt=0; ipt<vecx_.size(); ipt++ )
    {
      if(vecx_[ipt]>=xmin && vecx_[ipt]<=xmax)
	n++;
    }
  TPolyLine* band = new TPolyLine(n); 
  int cnt=0;
  for( size_t ipt=0; ipt<vecx_.size(); ipt++ )
    {
      if(vecx_[ipt]>=xmin && vecx_[ipt]<=xmax) {
	band->SetPoint(cnt, vecx_[ipt], (1+vecy_[ipt])>1.6?1.6:((1+vecy_[ipt])<0.4?0.4:1+vecy_[ipt]) );
	//cout<<ipt<<"   "<< vecx_[ipt]<<"   "<<1+vecy_[ipt]<<endl;
	cnt++;
      }
    }
  band->SetLineColor(kGray+1);
  band->SetFillColor(kGray+1);
  band->SetFillStyle(3002);

  return band;
}


void
UseTree::InitSkimming(int i, int ie, int nent) {
  
  _storeTuple=false;
  if(ie==-1) return;

  string tmpDs = GetDataset(i,ie);
  
  //  cout<<tmpDs<<"  "<<wDs<<endl;

  if(wDs!=tmpDs) {
    oFile = new TFile( ("Skimming/"+tmpDs+"_skim.root").c_str(),"RECREATE");
    tChains[i]->LoadTree(0);
    skimTree = (TTree*)tChains[i]->CloneTree(0);
    skimTree->SetDirectory( oFile );
    wDs = tmpDs;
  }

  if(ie+1<nent)
    tmpDs = GetDataset(i,ie+1);
  if(wDs!=tmpDs || ie+1==nent) {
    _storeTuple=true;

  }

//   cout<<i<<"  "<<ie<<"   "<<nent<<"   "<<wDs
//       <<"   "<<tmpDs<<"   "<<_storeTuple<<endl;
  
//   cout<<i<<"  "<<ie<<"   "<<nent<<"   "<<wDs
//       <<"   "<<tmpDs<<"   "<<_storeTuple<<"  "
//       <<skimTree<<endl;

}

void
UseTree::FinalizeSkimming() {

   if(!_storeTuple) return;

 //   cout<<skimTree<<"  "<<oFile<<endl;
//    cout<<" writing "<<skimTree->GetEntries()<<endl;
   
   oFile->cd();
   skimTree->Write();
   oFile->Write();
   oFile->Close();
    
}
