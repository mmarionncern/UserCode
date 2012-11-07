#include "HistoManager.hh"
#include <assert.h>

using namespace std;

ClassImp(HistoManager)

HistoManager::HistoManager() {
  TH1::SetDefaultSumw2(true);
  TH1::AddDirectory(kFALSE);
  //Histos = histos;

  _initsequence=true;
}

int HistoManager::FindNVar(string var) {
  itVar = Variables.find(var);
  if(itVar!=Variables.end())
    return Variables[ var ];
  else
    return -1;
}

int HistoManager::FindNUnc(string var) {
  itVar = Unc.find(var);
  if(itVar!=Unc.end())
    return Unc[ var ];
  else
    return -1;
}


systM
HistoManager::FindSysts(string var,string type) { //for uncertainties

  systM vars;
  //  cout<<" SYst=== "<<var<<"   "<<type<<endl;

  vector<string> names;
  for(itVar=Variables.begin();itVar!=Variables.end();itVar++) {

    string name = (*itVar).first;

    if( (name).find(var+"Unc")==(size_t)-1) continue;
      
    if(type=="Up" && (name).find("Up")!=(size_t)-1) {
      names.push_back(name);
    }
    if(type=="Do" && name.find("Do")!=(size_t)-1) {
      names.push_back(name);
    }
    if(type=="" && (name).find("Up")==(size_t)-1 
       &&  (name).find("Do")==(size_t)-1 ) {
      names.push_back(name);
    }
  }//loop
  
  //No idea of why it is needed to make that in two time...
  for(size_t in=0;in<names.size();in++) {
    // cout<<names[in]<<endl;
    vars[ names[in] ] = getHisto( names[in] ,nt,true);
  }
  
  return vars;
}

void HistoManager::PrintVariables() {
  
  std::map<string, int>::const_iterator itVarb = Variables.begin(); 
  std::map<string, int>::const_iterator itVare = Variables.end(); 

  for(itVar=itVarb;itVar!=itVare;itVar++) {
    cout<<(*itVar).first<<" --->   "<<(*itVar).second<<endl;
  }

}

vector<string> HistoManager::GetVariables() {

  vector<string> vars;

  std::map<string, int>::const_iterator itVarb = Variables.begin(); 
  std::map<string, int>::const_iterator itVare = Variables.end(); 

  for(itVar=itVarb;itVar!=itVare;itVar++) {
    vars.push_back((*itVar).first );
  }

  return vars;
}

int HistoManager::access(int nvar, int DS) {

  //  if(nvar<Nvar)
    return nvar+DS*Nvar;
  // else //uncertainty
  //   return nvar+nt*(nvar+1);
}

int HistoManager::accessUnc(int nvar, int nunc) {

  return nvar+nt*(Nvar)+nunc;

}

TH1F* HistoManager::getHisto(string obs, int DS, bool unc) {
  int Obs = FindNVar(obs);
  int NS= access(Obs,DS);

  if(unc) {
    int Nu = FindNUnc(obs);
    NS= NS+Nu;
  }

  //  cout<<" NS= "<<NS<<endl;
  return (TH1F*)Histos[ NS ]->Clone();
}


legend HistoManager::FindLeg(string var) {
  legend tmp;
  itVarLeg  = VarLegend.find(var);
  if(itVarLeg == VarLegend.end() )
    return tmp;
  else {
    tmp = (*itVarLeg).second;
    Name = (tmp).second;
    //VarLegend[ var ]
    return tmp;// VarLegend[ var ];
  }
}

void HistoManager::ConfigAnalysis(int N) {
  nt = N;
  
  Nunc =0;
  Nvar =0;
  Nvar2D =0;
}

void HistoManager::AddVariable(string var, int nBin, double min, double max,
		 string Xleg, string Yleg) {
  
  //Protection against overdeclaration
  itVar = Variables.find(var);
  if(itVar!=Variables.end()) {
    cout<<" Be careful : "<< var <<"already declared "<<endl;
    return;
  }

  vector<float> temp(3,0);
  for(int i=0;i<nt+1;i++) {
    TH1F* tmp =NULL; 
  
    Histos.push_back(tmp); 
  }

  temp[0] = nBin;  temp[1] = min; temp[2] = max;
  templates.push_back(temp);
  
  Variables[ var ] = Nvar;
  legend leg(Xleg,Yleg);
  VarLegend[ var ] = leg;

  Templates[ var ] = temp;
  if(Xleg!="Unc" && Xleg!="Unc")
  Nvar++;
}


void HistoManager::PrepareHistos(vector<string> Datasets) {
  
  for(int unsigned i=0;i<Variables.size();i++) {
    for(int j=0;j<nt+1;j++) {
    

      ostringstream os;
      os << i;
    
      string name = Datasets[j] +"_"+ os.str();
      string title = os.str();
      TH1F* tmp = new TH1F(name.c_str(), title.c_str(),
			   (int)templates[i][0], templates[i][1],
			   templates[i][2]  );
    
      Histos[i+Nvar*j] = tmp;

      //  if(name.substr(0,7)=="pfSumEt")
      //cout<<" the histo to be declared  "<<templates[i][2]<<"  "<<i<<"  "<<j<<"  "<<name<<"  "<<i+Nvar*j<<endl;
    }
  }
  
}

void HistoManager::PrepareUncHisto(string var, int nb, float min, float max) {
  
  int nv = Variables[ var ];
  //int nv = Variables[ var ];
  // nv+=Nunc;

  //for(int j=0;j<nt+1;j++) { //datasets
    ostringstream os,os2;
    os << nv+Nunc;
    // os2 << j;
    string name = "Unc_"+ os.str();
    string title = os.str();
    TH1F* tmp = new TH1F(name.c_str(), title.c_str(),
		         nb, min, max);
    
    //cout<<" ======> unc histo : "<<var<<"   "<<nv<<"  "<<Nvar<<"  "<<nt<<"  "<<nv+(Nvar)*nt+Nunc<<" :: "<<name
    //	<<endl;


    Histos[ nv+(Nvar)*nt+Nunc ] = tmp;
    Unc[ var ] = Nunc;
    Nunc++;
    //   cout<<" tapoué "<<endl;

    //}
}


void HistoManager::fill(string var, int DS , float value ) {

  if(_initsequence) return;

  itVar = Variables.find(var);
  
  if(itVar!=Variables.end()) {
    int nv = Variables[ var ];
    int nh = nv+DS*Nvar;
    
    Histos[nh]->Fill(value);
  }
  else {
    cout<<" Problème "<<var<<"  non déclarée "<<endl;
    abort();
  }
}

void HistoManager::fill(string var, int DS , float value, float weight ) {

  if(_initsequence) return;

  itVar = Variables.find(var);
  // if(var.substr(0,3)=="NVe")
  //   cout<<" var "<<var<<"  "<<DS<<"  "<<value<<"  "<<weight<<endl;
  if(itVar!=Variables.end()) {
    // if(var.substr(0,3)=="NVe")
      //  cout<<" yiyiyi "<<endl;
    int nv = Variables[ var ];
    int nh = nv+DS*Nvar;
    // if(var.substr(0,3)=="NVe")
    //   cout<<" the histo to fill  "<<var<<"  "<<nh<<"  "<<Histos.size()<<"  "<<Histos[nh]<<endl;
    Histos[nh]->Fill(value,weight);
  }
  else {
    cout<<" Problème "<<var<<"  non déclarée "<<endl;
    abort();
  }
  
}


void HistoManager::fillUnc(string var, string type , int DS , float value, float weight, string dir ) {

  if(_initsequence) return;

  if(DS==nt){ return;} //no uncertainty on data
  //cout<<"aui? "<<var+type+dir<<endl;
  itVar = Variables.find(var);
  
  if(itVar!=Variables.end()) {

    //Build default Unc
    string varUnc=var+"Unc"+type+dir;
    itVar = Variables.find(varUnc);

    if(itVar!=Variables.end()) {
      
      int nv = Variables[ varUnc ];
      int nu = Unc[ varUnc ];
      int nh = nv+nt*(Nvar)+nu;//+DS*Nvar;
      //cout<<nv<<"  "<<Nvar<<"  "<<nu<<" "<<nt<<"  "<<nh<<endl;
      /*  if(var.substr(0,7)=="pfSumEt")
	  cout<<" the histo to fill  "<<var<<"  "<<nh<<endl;*/
      Histos[nh]->Fill(value,weight);
  
      //      if(var.substr(0,3)=="pat")
      //cout<<" the histo to fill as Unc  "<<var<<"  "<<varUnc<<"  "<<nh<<endl;
  
    
    }//unc variable
    else {

      
      int nv = FindNVar( var );
      vector<float> tmp=GetTemplate(nv);
      //   cout<< " declaration "<<var<<"  "<<nv<<"  "<<varUnc<<"   "<<tmp[0]<<"  "<<tmp[1]<<"  "<<tmp[2]<<endl;
      AddVariable(varUnc, (int)tmp[0], tmp[1], tmp[2], "Unc", "Unc");
      PrepareUncHisto(varUnc, (int)tmp[0], tmp[1], tmp[2]);
    }

  } //variable
  else {
    
    cout<<" Problème "<<var<<"  non déclarée "<<endl;
    abort();
  }
  
}


vector<TH1F*> HistoManager::GetHistos() {
  return Histos;
}


vector<float> HistoManager::GetTemplate(int var) {

  vector<float> tmp(0,3);
  tmp = templates[var];
  return tmp;
}
vector<float> HistoManager::GetTemplate2D(int var) {

  vector<float> tmp(0,6);
  tmp = Templates2D[var];
  return tmp;
}

vector<float> HistoManager::GetTemplateP(int var) {

  vector<float> tmp(0,3);
  tmp = TemplatesP[var];
  return tmp;
}

double* HistoManager::GetVarTemplateP(string var) {
  //FIXME
  double* tmp = NULL;
  tmp = VarTemplates[ var ];
  if(tmp == NULL) {
    cout<<" Take care, no variable bin for observable "<<var<<endl; }
  return tmp;
}


void HistoManager::Add2DVariable(string var, string Xleg, string Yleg,string name, int NBinsX, float Xm,float XM, int NBinsY, float Ym, float YM) {
  
  //Protection against overdeclaration
  itVar2D = Variables2D.find(var);
  if(itVar2D!=Variables2D.end()) {
    cout<<" Be careful : "<< var <<"already declared "<<endl;
    return;
  }


  for(int i=0;i<nt+1;i++) {
    TH2F* tmp; 
    Histos2D.push_back(tmp); 
    
    //  vector<float> vtmp;
    //  Weight2D[ tmp ] = vtmp;
  }

  Variables2D[ var ] = Nvar2D;
  legend leg(Xleg,Yleg);
  VarLegend2D[ var ] = leg;

  VarName2D[ var ] = name;
  
  vector<float> temp(6,0);
  temp[0] = NBinsX; temp[1] = Xm; temp[2] = XM;
  temp[3] = NBinsY; temp[4] = Ym; temp[5] = YM;

  Templates2D.push_back(temp);
  
  Nvar2D++;
}

void HistoManager::Add2DVariable(string var, string Xleg, string Yleg,string name, int NBinsX, double VarBinsX[], int NBinsY, double VarBinsY[]) {
  
  for(int i=0;i<nt+1;i++) {
    TH2F* tmp; 
    Histos2D.push_back(tmp); 
  }

  Variables2D[ var ] = Nvar2D;
  legend leg(Xleg,Yleg);
  VarLegend2D[ var ] = leg;

  VarName2D[ var ] = name;
  
  vector<float> temp(2,0);
  temp[0] = NBinsX;
  temp[1] = NBinsY;

  VarTemplates2D[ var ].first  = VarBinsX;
  VarTemplates2D[ var ].second = VarBinsY;

  Templates2D.push_back(temp);
  
  Nvar2D++;
}


/*
void HistoManager::Add2DVariable(string var, string Xleg, string Yleg,string name, int NBinsX, float[] VarBin, int NBinsY, float Ym, float YM) {
  
  for(int i=0;i<nt+1;i++) {
    TGraph* tmp; 
    Histos2D.push_back(tmp); 
  }

  Variables2D[ var ] = Nvar2D;
  legend leg(Xleg,Yleg);
  VarLegend2D[ var ] = leg;

  VarName2D[ var ] = name;

  vector<float> temp(6,0);
  temp[0] = NBinsX; temp[1] = Xm; temp[2] = XM;
  temp[3] = NBinsY; temp[4] = Ym; temp[5] = YM;

  Templates2D[ var ] = temp;
  
  Nvar2D++;
}

*/

void HistoManager::AddProfVariable(string var, int nBin, double min, double max,
		 string Xleg, string Yleg) {
  
  //Protection against overdeclaration
  itVar = VariablesP.find(var);
  if(itVar!=VariablesP.end()) {
    cout<<" Be careful : "<< var <<"already declared "<<endl;
    return;
  }

  vector<float> temp(3,0);
  for(int i=0;i<nt+1;i++) {
    TProfile* tmp; 
  
    Profiles.push_back(tmp); 
  }

  temp[0] = nBin;  temp[1] = min; temp[2] = max;
  TemplatesP.push_back(temp);
  
  VariablesP[ var ] = NvarP;
  legend leg(Xleg,Yleg);
  VarLegendP[ var ] = leg;

  NvarP++;
}

void HistoManager::AddProfVariable(string var, int nBin, double varBins[],
		 string Xleg, string Yleg) {
  
  //Protection against overdeclaration
  itVar = VariablesP.find(var);
  if(itVar!=VariablesP.end()) {
    cout<<" Be careful : "<< var <<"already declared "<<endl;
    return;
  }

  for(int i=0;i<nt+1;i++) {
    TProfile* tmp; 
  
    Profiles.push_back(tmp); 
  }

  VarTemplates[ var ] = varBins;
  vector<float> temp(1,0);
  temp[0] = nBin;
  TemplatesP.push_back(temp);

  VariablesP[ var ] = NvarP;
  legend leg(Xleg,Yleg);
  VarLegendP[ var ] = leg;

  NvarP++;
}

void HistoManager::Prepare2DHistos(vector<string> Datasets) {
  
  //for(int unsigned i=0;i<Variables2D.size();i++) {
  int i=0;
  for(map<string, int>::const_iterator itvartmp= Variables2D.begin();
      itvartmp!=Variables2D.end();itvartmp++ ) {
    
    for(int j=0;j<nt+1;j++) {
      
      ostringstream os;
      os << i;
   
      string name = Datasets[j] +"_2D_"+ os.str();
      
      string title = os.str();
      TH2F* tmp;
    
      itVarT2D = VarTemplates2D.find( (*itvartmp).first );
      if(itVarT2D == VarTemplates2D.end() ) {
	tmp =new TH2F(name.c_str(), title.c_str(),(int)Templates2D[i][0], Templates2D[i][1],
		      Templates2D[i][2] , (int)Templates2D[i][3], Templates2D[i][4],
		      Templates2D[i][5] );
      }
      else {
	tmp=new TH2F(name.c_str(), title.c_str(),
			 (int)Templates2D[(*itvartmp).second][0],
			 (*itVarT2D).second.first,
			 (int)Templates2D[(*itvartmp).second][1], 
			 (*itVarT2D).second.second );
      }
    
      Histos2D[i+Nvar2D*j] = tmp;
      
      // delete tmp;

    }

    i++;
  }
}

void HistoManager::fill2D(string var, int DS , float valueX, float valueY ) {
 
  if(_initsequence) return;

  itVar = Variables2D.find(var);
  if(itVar!=Variables2D.end()) {
    int nv = Variables2D[ var ];
    int nh = nv+DS*Nvar2D;
  
    Histos2D[nh]->Fill(valueX,valueY);
  }
  else {
    cout<<" Error, 2D Histo Variable "<<var<<" not declared "<<endl;
  }


}

void HistoManager::fill2D(string var, int DS , float valueX, float valueY, float weight) {

  if(_initsequence) return;

  itVar = Variables2D.find(var);
  if(itVar!=Variables2D.end()) {
    int nv = Variables2D[ var ];
    int nh = nv+DS*Nvar2D;
    
    Histos2D[nh]->Fill(valueX,valueY,weight);

  }
  else{
    cout<<" Problème "<<var<<"  non déclarée "<<endl;
    abort();
  }
    
}

void HistoManager::PrepareProfiles(vector<string> Datasets) {
  
  int i=0;
  for(map<string, int>::const_iterator itvartmp= VariablesP.begin();
      itvartmp!=VariablesP.end();itvartmp++ ) {
    for(int j=0;j<nt+1;j++) {
     
      ostringstream os;
      os << i;
    
      string name = Datasets[j] +"_Prof_"+ os.str();
      string title = os.str();

      TProfile* tmp;

      itVarTP = VarTemplates.find( (*itvartmp).first );
      if(itVarTP == VarTemplates.end() ) {
	tmp=new TProfile(name.c_str(), title.c_str(),(int)TemplatesP[(*itvartmp).second][0],
			 TemplatesP[(*itvartmp).second][1], TemplatesP[(*itvartmp).second][2] );
      }
      else {
	tmp=new TProfile(name.c_str(), title.c_str(),(int)TemplatesP[(*itvartmp).second][0], (*itVarTP).second );
      }
      Profiles[(*itvartmp).second+NvarP*j] = tmp;
      
    }
    
    i++;

  }
}

void HistoManager::fillProf(string var, int DS , float valueX, float valueY ) {
 
  if(_initsequence) return;

  itVar = VariablesP.find(var);
  if(itVar!=VariablesP.end()) {
    int nv = VariablesP[ var ];
    int nh = nv+DS*NvarP;
    
    Profiles[nh]->Fill(valueX,valueY);
  }
}

void HistoManager::fillProf(string var, int DS , float valueX, float valueY, float weight) {

  if(_initsequence) return;

  itVar = VariablesP.find(var);
  if(itVar!=VariablesP.end()) {
    int nv = VariablesP[ var ];
    int nh = nv+DS*NvarP;

    Profiles[nh]->Fill(valueX,valueY,weight);
  }
}

vector<TH2F*> HistoManager::GetHistos2D() {
  return Histos2D;
}

vector<TProfile*> HistoManager::GetProfiles() {
  return Profiles;
}


int HistoManager::FindNVar2D(string var) {

  itVar = Variables2D.find(var);
  if(itVar!=Variables2D.end())
    return Variables2D[ var ];
  else
    return -1;
}

int HistoManager::FindNVarP(string var) {

  itVar = VariablesP.find(var);
  if(itVar!=VariablesP.end())
    return VariablesP[ var ];
  else {
    cout<<" Take Care, no such observable "<<endl;
    return -1;
  }

}

int HistoManager::access2D(int nvar, int DS) {
  return nvar+DS*Nvar2D;
}

int HistoManager::accessProf(int nvar, int DS) {
  return nvar+DS*NvarP;
}

std::pair<string,string> HistoManager::FindLeg2D(string var) {

  Name = VarName2D[ var ];

  return VarLegend2D[ var ];
}

std::pair<string,string> HistoManager::FindLegP(string var) {

  Name = VarNameP[ var ];

  return VarLegendP[ var ];
}
/*
void HistoManager::Cleaning2DHisto(vector<int> Nmax) {

  for(int unsigned i=0;i<Variables2D.size();i++) {
    for(int j=0;j<nt+1;j++) {
 

      double x,y;
      //int N = NEntries[j];
     
      for(int k=0;k<Nmax[j];k++) {

	Histos2D[i+Nvar2D*j]->GetPoint(k,x,y);
	if(x==0 && (y==0 || y<-1000) )
	  Histos2D[i+Nvar2D*j]->RemovePoint(k);
      }
    }
   }
   
}
*/

TH2F* HistoManager::ConvertGraphToHisto(TGraph* graph, int NbinX, double mX, double MX,
				      int NbinY, double mY, double MY,string name) {

  TH2F* tmp=new TH2F(name.c_str(),name.c_str(),NbinX,mX,MX,NbinY,mY,MY);

  double x,y;

  for(int i=0;i<graph->GetN();i++) {
    
    graph->GetPoint(i,x,y);
 
    if(x!=0 && y!=0) {
      tmp->Fill(x,y);
    }
  }

  return tmp;
  
}
 

TH2F* HistoManager::ConvertGraphToHisto(TGraph* graph, int NbinX, double BinX[],
				      int NbinY, double BinY[],string name) {

  TH2F* tmp=new TH2F(name.c_str(),name.c_str(),NbinX,BinX,NbinY,BinY);
  double x,y;

  for(int i=0;i<graph->GetN();i++) {
    
    graph->GetPoint(i,x,y);
    
    if(x!=0 && y!=0) {
      tmp->Fill(x,y);
    }
  }

  return tmp;
  
}

 
TProfile* HistoManager::ConvertGraphToProfile(TGraph* graph, int NbinX, double mX, double MX,
					    string name,string opt) {

  TProfile* tmp=new TProfile(name.c_str(),name.c_str(),NbinX,mX,MX);
  tmp->SetErrorOption(opt.c_str());
  double x,y;

  for(int i=0;i<graph->GetN();i++) {
    
    graph->GetPoint(i,x,y);
    
    if(x!=0 && y!=0) {
      tmp->Fill(x,y);
    }
  }

 
  return tmp;
  
}

vector<TGraph*> HistoManager::ConvertGraphToErrorProfile(TGraph* graph, int NbinX, double mX, double MX ) {

  double pas = (MX-mX)/NbinX;
  // TProfile* tmp=new TProfile((name+"_M").c_str(),name.c_str(),NbinX,mX,MX);
  // TProfile* tmpEm=new TProfile((name+"_Em").c_str(),name.c_str(),NbinX,mX,MX);
  // TProfile* tmpEM=new TProfile((name+"_EM").c_str(),name.c_str(),NbinX,mX,MX);
  // tmp->SetErrorOption(opt.c_str());
  TGraph* tmp=new TGraph(NbinX);
  TGraph* tmpEm=new TGraph(NbinX);
  TGraph* tmpEM=new TGraph(NbinX);

 double x,y;

  vector<TH1F*> tmp1h;
  for(int i=0;i<NbinX;i++) {
    TH1F* tmpf = new TH1F((TString)i, "", 2000, -100,100);
    tmp1h.push_back(tmpf);
  }

  for(int i=0;i<graph->GetN();i++) {
    
    graph->GetPoint(i,x,y);
    
    for(int k=0;k<NbinX;k++) {
      if(x!=0 && y!=0) {
	
	if(x>(k*pas) && x<( (k+1)*pas) )
	  tmp1h[k]->Fill(y);
      }
    }
  }
  

  double Em,EM,M;
  for(int i=0;i<NbinX;i++) {
    ComputeStatFromHisto(tmp1h[i], Em , EM, M);
    tmp->SetPoint(i,i*pas+pas/2.,M);
    tmpEm->SetPoint(i,i*pas+pas/2.,Em);
    tmpEM->SetPoint(i,i*pas+pas/2.,EM);
    cout<<i<<" Mean "<<M<<" ; Emin "<<Em<<" ; Emax "<<EM<<endl;
    delete tmp1h[i];
    // cout<<tmp->GetBinContent(i)<<endl;
  }
  
  
  vector<TGraph*> Graphs;
  Graphs.push_back(tmp);
  Graphs.push_back(tmpEm);
  Graphs.push_back(tmpEM);
 
  return Graphs;
  
}


vector<TGraph*> HistoManager::ConvertGraphToErrorProfile(TGraph* graph, int NbinX, int BinVar[]) {

  TGraph* tmp=new TGraph(NbinX);
  TGraph* tmpEm=new TGraph(NbinX);
  TGraph* tmpEM=new TGraph(NbinX);

 double x,y;

  vector<TH1F*> tmp1h;
  for(int i=0;i<NbinX;i++) {
    TH1F* tmpf = new TH1F((TString)i, "", 2000, -100,100);
    tmp1h.push_back(tmpf);
  }

  for(int i=0;i<graph->GetN();i++) {
    
    graph->GetPoint(i,x,y);
    
    for(int k=0;k<NbinX;k++) {
      if(x!=0 && y!=0) {
	
	if(x>=BinVar[k] && x<BinVar[k+1] )
	  tmp1h[k]->Fill(y);
      }
    }
  }
  

  double Em,EM,M;
  for(int i=0;i<NbinX;i++) {

    double pas = (BinVar[i+1]-BinVar[i]);
    ComputeStatFromHisto(tmp1h[i], Em , EM, M);
    if(i==0) {
      tmp->SetPoint(i,BinVar[i],M);
      tmpEm->SetPoint(i,BinVar[i],Em);
      tmpEM->SetPoint(i,BinVar[i],EM);
    }
    else {
      tmp->SetPoint(i,BinVar[i]+pas/2.,M);
      tmpEm->SetPoint(i,BinVar[i]+pas/2.,Em);
      tmpEM->SetPoint(i,BinVar[i]+pas/2.,EM);
    }
    cout<<i<<" Mean "<<M<<" ; Emin "<<Em<<" ; Emax "<<EM<<endl;
    delete tmp1h[i];
    // cout<<tmp->GetBinContent(i)<<endl;
  }
  
  
  vector<TGraph*> Graphs;
  Graphs.push_back(tmp);
  Graphs.push_back(tmpEm);
  Graphs.push_back(tmpEM);
 
  return Graphs;
  
}
/*
TH1F* HistoManager::ConvertGraphBinToHisto(TGraph* graph, int NbinX, double Xmin, double Xmax,
					 double Binm, double BinM,string name, int nh) {
  
  TH1F* tmp = new TH1F(name.c_str(), "", NbinX, Xmin, Xmax );
  double x,y;
   for(int i=0;i<graph->GetN();i++) {
    
     graph->GetPoint(i,x,y);
    
     if(x!=0 && y!=0) {
	if(x>=Binm && x<BinM )
	  {
	    tmp->Fill(y, (Weight2D[ Histos2D[nh] ])[i]);//FIXME
	  }
     }
   }
   
   return tmp;
   }*/
/*
TH1F* HistoManager::ConvertGraphBinToHisto(TGraph* graph, int NbinX, double Xmin, double Xmax,
					 double Binm, double BinM,string name, int& meanx, int nh) {
  
  TH1F* tmp = new TH1F(name.c_str(), "", NbinX, Xmin, Xmax );
  int N=0;
   int meanxtmp = 0;
  double x,y;
   for(int i=0;i<graph->GetN();i++) {
    
     graph->GetPoint(i,x,y);
    
     if(x!=0 && y!=0) {
	if(x>=Binm && x<BinM )
	  {
	    tmp->Fill(y, (Weight2D[ Histos2D[nh] ])[i] );//FIXME
	    meanxtmp +=x;
	    N++;
	  }
     }
   }
   meanxtmp /= N;
   meanx = meanxtmp;
   
   return tmp;
}
*/


TGraphAsymmErrors* HistoManager::GraphReduction(TGraph* graph, int NbinX, double mX, double MX,
					      string name,string opt) {

  double pas = (MX-mX)/NbinX;

  TGraphAsymmErrors* tmp=new TGraphAsymmErrors(NbinX);

  double x,y;

  vector<TH1F*> tmp1h;
  vector<TH1F*> tmp2h;
  for(int i=0;i<NbinX;i++) {
    TString s= "_"+name;
    s +=i;
    TH1F* tmpf = new TH1F(s, "", 2000, -100,100);
    tmp1h.push_back(tmpf);
    TH1F* tmpf2 = new TH1F(s+"_", "", 2000, -100,100);
    tmp2h.push_back(tmpf2);
  }

  
  for(int i=0;i<graph->GetN();i++) {
    
    graph->GetPoint(i,x,y);
    for(int k=0;k<NbinX;k++) {
      if(x!=0 && y!=0) {
	
	if(x>=(k*pas+mX) && x<( (k+1)*pas+mX) ) {
	  tmp1h[k]->Fill(y);
	  tmp2h[k]->Fill(x);
	}
      }
    }
  }
  /* cout<<"  "<<mX<<"  "<<MX<<endl;
  for(int k=0;k<NbinX;k++)
  cout<<"  "<<k*pas+mX<<"  "<<(k+1)*pas+mX<<endl;*/
     
  double Em,EM,M;
  double Em2,EM2,M2;
  for(int i=0;i<NbinX;i++) {
    ComputeStatFromHisto(tmp1h[i], Em , EM, M);
    ComputeStatFromHisto(tmp2h[i], Em2 , EM2, M2);

    if(opt!="rms") {
      tmp->SetPoint(i,M2,M);
      if(opt=="")
	tmp->SetPointError(i,fabs(Em2-M2),fabs(EM2-M2),fabs(Em-M)/tmp1h[i]->GetEntries(),fabs(EM-M)/tmp1h[i]->GetEntries() );
      else
	tmp->SetPointError(i,fabs(Em2-M2),fabs(EM2-M2),fabs(Em-M),fabs(EM-M) );
    }
    else {
      tmp->SetPoint(i,M2,fabs(EM-M) );
    }

    //  tmp->SetPointError(i,fabs(Em2-M2),fabs(EM2-M2),fabs(Em-M),fabs(EM-M) );
    //   cout<<"i "<<i<<" -> ("<<M2<<","<<M<<"  :  "<<Em2-M2<<" , "<<EM2-M2<<"  ;  "<<Em-M<<" , "<<EM-M<<endl;
 
    delete tmp1h[i];
    delete tmp2h[i];
  }
  
  // tmp2h[10]->Draw();
  return tmp;
  
}



TGraphAsymmErrors* HistoManager::GraphReduction(TGraph* graph, int NbinX, int BinVar[],
							string name,string opt) {

  TGraphAsymmErrors* tmp=new TGraphAsymmErrors(NbinX);
 
 double x,y;

  vector<TH1F*> tmp1h;
  vector<TH1F*> tmp2h;
  for(int i=0;i<NbinX;i++) {
    TString s= "_"+name;
    s +=i;
    TH1F* tmpf = new TH1F(s, "", 2000, -100,100);
    tmp1h.push_back(tmpf);
    TH1F* tmpf2 = new TH1F(s+"_", "", 2000, -100,100);
    tmp2h.push_back(tmpf2);
  }

  for(int i=0;i<graph->GetN();i++) {
    
    graph->GetPoint(i,x,y);
    
    for(int k=0;k<NbinX;k++) {
      if(x!=0 && y!=0) {
	
	if(x>=BinVar[k] && x<BinVar[k+1] ) {
	  tmp1h[k]->Fill(y);
	  tmp2h[k]->Fill(x);
	}
      }
    }
  }
  
  double Em,EM,M;
  double Em2,EM2,M2;
  for(int i=0;i<NbinX;i++) {

    ComputeStatFromHisto(tmp1h[i], Em , EM, M);
    ComputeStatFromHisto(tmp2h[i], Em2 , EM2, M2);

    if(opt=="")
      {
	Em /= tmp1h[i]->GetEntries();
	EM /= tmp1h[i]->GetEntries(); 
      }

    if(i==0) {
      tmp->SetPoint(i,M2,M);
      tmp->SetPointError(i,Em2,EM2,Em,EM);
    }
    else {
      tmp->SetPoint(i,M2,M);
      tmp->SetPointError(i,Em2,EM2,Em,EM);
    }
    // cout<<"i "<<i<<" -> ("<<M2<<","<<M<<"  :  "<<Em2<<" , "<<EM2<<"  ;  "<<Em<<" , "<<EM<<endl;
    delete tmp1h[i];
    delete tmp2h[i];
    // cout<<tmp->GetBinContent(i)<<endl;
  }
  
  
  return tmp;
  
}

/*
TH1D* HistoManager::ConvertHistoToResponsePlotVarBin(TH2F* histo2D, int n, double VarBin[], string namevar) {
  
  TH1D* RMSP = new TH1D( ("ResponseProf"+namevar).c_str(),"",n,VarBin); 
  for(int i=0;i<n;i++) {
    
    int Binm = histo2D->GetXaxis()->FindBin(VarBin[i]);
    int BinM = histo2D->GetXaxis()->FindBin(VarBin[i+1]);
    float sigma=0, sigerr=0;

    TH1D* tmp = (TH1D*)histo2D->ProjectionY("",Binm, BinM); 
    tmp->Rebin(4); cout<<tmp->GetEntries()<<endl;
    
    if(tmp->GetEntries()!=0) {
      cout<<tmp->GetRMS()<<"   "<<i<<"  "<<Binm<<"  "<<BinM<<"  "<<n<<" ---->  ";
      TF1* Gaus = new TF1("f1","gaus");
      histo2D->ProjectionY("",Binm, BinM)->Fit("f1","QN0");
      sigma = Gaus->GetParameter(2);
      sigerr = Gaus->GetParError(2);
      if(i==4) {Gaus->Draw("same");}

      if(rmConta) {

	TH1D* tmpConta = GetContamination2D(namevar, i, VarBin);
	TH1D* cleaned=(TH1D*)tmp->Clone();
	cleaned->Reset("ICEM");
	float xm=0;
	float xm2=0;
	float nevt=0;
	
	for(int ibin=1;ibin<cleaned->GetNbinsX()+1;ibin++) {
	  //  if(tmp->GetBinContent(ibin)-tmpConta->GetBinContent(ibin) > 0) {
	    nevt += tmp->GetBinContent(ibin)-tmpConta->GetBinContent(ibin);
	    cleaned->SetBinContent(ibin, tmp->GetBinContent(ibin)-tmpConta->GetBinContent(ibin) );
	    xm += (tmp->GetBinContent(ibin)-tmpConta->GetBinContent(ibin)) * tmp->GetBinCenter(ibin);
	    xm2 += (tmp->GetBinContent(ibin)-tmpConta->GetBinContent(ibin))*pow(  tmp->GetBinCenter(ibin),2);
	    // cout<<ibin<<"  "<<tmp->GetBinContent(ibin)<<"  "<<tmpConta->GetBinContent(ibin)<<"  "<<tmp->GetBinContent(ibin)-tmpConta->GetBinContent(ibin)<<"  "<<tmp->GetBinCenter(ibin)<<"  "<<xm<<"  "<<xm2 <<endl;
	    // }
	}
	//	cout<<" xm "<<xm/nevt<<"  "<<xm*xm/(nevt*nevt) <<"  "<<xm2/nevt<<endl;
	// sigma=sqrt( xm2/nevt - pow(xm/nevt,2) );
	cleaned->Fit("f1","QN0");
	sigma = Gaus->GetParameter(2);
	sigerr = Gaus->GetParError(2);
	cout<<"   "<<tmpConta->GetRMS()<<"   "<<cleaned->GetRMS()<<"   "<<sigma<<"    "<<tmp->GetEntries()<<endl;
	RMSP->SetBinContent(i+1,sigma  ); //FIXME
      }
      else {
	cout<<"   "<<sigma<<"    "<<tmp->GetEntries()<<endl;
	RMSP->SetBinContent(i+1,sigma );	
	
      }
      RMSP->SetBinError(i+1, sigerr );
    }
    delete tmp;
  }
  return RMSP;
}
*/



/*
TProfile* HistoManager::ConvertGraphToResponsePlotVarBin(TGraph* graph, int n, double VarBin[], string name, int nh ) {
  vector<int> Nentries(n,0);
  //  TProfile* Phisto = new TProfile(( "firstProf"+name).c_str(),"",500,0,500);
  TProfile* Phisto = new TProfile( ("EndProf"+name).c_str(),"",n,VarBin);
  double x,y;
  //  Phisto->SetErrorOption("s");
  for(int i=0;i<graph->GetN();i++) {
    graph->GetPoint(i,x,y);
    if(x!=0 && y!=0  && x>2 && x<500) {
  
      Phisto->Fill(x,-y/x, (Weight2D[ Histos2D[nh] ])[i] ); //FIXME
      for(int k=0;k<n;k++)
	if(VarBin[k]<=x && VarBin[k+1]>x)
	  { Nentries[k]++; break;}
    }
      
  }
  TProfile* Phisto2 = new TProfile( ("EndProf2"+name).c_str(),"",n,VarBin);
  float cont,error;
  // Phisto2->SetErrorOption("s");
  for(int i=0;i<Phisto->GetNbinsX();i++) {
    cont = Phisto->GetBinContent(i);
    error = Phisto->GetBinError(i);
    cont /= Phisto->GetBinCenter(i);
    if(cont!=0) {
      Phisto2->Fill(Phisto->GetBinCenter(i),(double)fabs(cont));
    }
  }

  for(int k=0;k<n;k++)
    { cout<<" Bin "<<VarBin[k]<<"-"<<VarBin[k+1]<<" Number of entries "<<Nentries[k]<<endl;
      // cout<<Phisto1->GetBinEntries(k)<<endl;
    }

  return Phisto;
}
*/

TH1D* HistoManager::ConvertHistoToRMSPlotVarBin(TH2F* histo2D, int n, double VarBin[], string namevar, bool rmConta=false, bool cplot=false) {
  TCanvas* cc = new TCanvas(("can"+(string)(histo2D->GetName()) ).c_str(),"testrms ");
  cc->Divide(4,4);
  TH1D* RMSP = new TH1D( ("RMSProf"+namevar).c_str(),"",n,VarBin); 
  for(int i=0;i<n;i++) {
    
    int Binm = histo2D->GetXaxis()->FindBin(VarBin[i]);
    int BinM = histo2D->GetXaxis()->FindBin(VarBin[i+1]);
    
    float sigma=0, sigerr=0;

    ostringstream os;
    os<<i+1;
    string nam = (string)(histo2D->GetName()) + os.str();

    TH1D* tmp = (TH1D*)histo2D->ProjectionY(nam.c_str() ,Binm, BinM); 
    
    //cout<<tmp->GetEntries()<<endl;
    // if(i==9) {
    //TH1D* tmp2 = (TH1D*)tmp->Clone();
    cc->cd(i+1);
      //cout<<"=========================== "<<tmp2->GetRMS()<<endl;
    
    histo2D->ProjectionY("grau",Binm, BinM)->Draw();
      // }
      if(tmp->GetEntries()>=11) {

	//TH1* tmp2 =  Get90HistoInterval(tmp);
	/*	tmp->Rebin(2); 
      if(tmp->GetEntries()<200)
	tmp->Rebin(2); 
      if(tmp->GetEntries()<100)
      tmp->Rebin(2); */
      if(cplot)
	tmp->Draw();
      //histo2D->ProjectionY("grau",Binm, BinM)->Draw();
    
      double range[2]; Get90PCInterval( tmp, range[0], range[1] );
      if(tmp->GetEntries() > 10000 )
	{ range[0] = -1000; range[1]=1000; }
      cout<<tmp->GetRMS()<<"   "<<i<<"  "<<Binm<<"  "<<BinM<<"  "<<"  ;  "<<range[0]<<"  "<<range[1]<<"  "<<n<<" ---->  ";
      TF1* Gaus = new TF1("f1","gaus",range[0], range[1]+20);
      cout<<range[0]<<"  "<<range[1]<<"  "<<n<<" ---->  ";
      if(cplot)
	tmp->Fit("f1","QR");
      else
	tmp->Fit("f1","QR");
      sigma = Gaus->GetParameter(2);
      sigerr = Gaus->GetParError(2);
      // if(i==9) {Gaus->Draw("same");}
      //  Gaus->Draw("same");

      float xm=0;
      float xm2=0;
      float nevt=0;
      sigma = tmp->GetRMS();
      sigerr = tmp->GetRMSError();
      
      /*      for(int ibin=1;ibin<tmp2->GetNbinsX()+1;ibin++) {
	nevt += tmp2->GetBinContent(ibin);
	xm += (tmp2->GetBinContent(ibin)) * tmp2->GetBinCenter(ibin);
	xm2 += (tmp2->GetBinContent(ibin))*pow(  tmp2->GetBinCenter(ibin),2);
      }
      sigma=sqrt( xm2/nevt - pow(xm/nevt,2) );*/

      if(rmConta) {

	TH1D* tmpConta = GetContamination2D(namevar, i, VarBin);
	if(tmp->GetEntries()<200)
	  tmpConta->Rebin(2);
	if(tmp->GetEntries()<70)
	  tmpConta->Rebin(2); 
	TH1D* cleaned=(TH1D*)tmp->Clone();
	cleaned->Reset("ICEM");

	xm=0;
	xm2=0;
	nevt=0;
	
	for(int ibin=1;ibin<cleaned->GetNbinsX()+1;ibin++) {
	  if(tmp->GetBinContent(ibin)-tmpConta->GetBinContent(ibin) > 0) {
	  nevt += tmp->GetBinContent(ibin)-tmpConta->GetBinContent(ibin);
	  cleaned->SetBinContent(ibin, tmp->GetBinContent(ibin)-tmpConta->GetBinContent(ibin) );
	  xm += (tmp->GetBinContent(ibin)-tmpConta->GetBinContent(ibin)) * tmp->GetBinCenter(ibin);
	  xm2 += (tmp->GetBinContent(ibin)-tmpConta->GetBinContent(ibin))*pow( tmp->GetBinCenter(ibin),2);
	    // cout<<ibin<<"  "<<tmp->GetBinContent(ibin)<<"  "<<tmpConta->GetBinContent(ibin)<<"  "<<tmp->GetBinContent(ibin)-tmpConta->GetBinContent(ibin)<<"  "<<tmp->GetBinCenter(ibin)<<"  "<<xm<<"  "<<xm2 <<endl;
	  }
	}
	//	cout<<" xm "<<xm/nevt<<"  "<<xm*xm/(nevt*nevt) <<"  "<<xm2/nevt<<endl;
	// sigma=sqrt( xm2/nevt - pow(xm/nevt,2) );
	cleaned->Fit("f1","QN0");
	sigma = Gaus->GetParameter(2);
	sigerr = Gaus->GetParError(2);
	cout<<"   "<<tmpConta->GetRMS()<<"   "<<cleaned->GetRMS()<<"   "<<sigma<<"    "<<tmp->GetEntries()<<endl;
	RMSP->SetBinContent(i+1, sigma  ); //FIXME
      }
      else {
	cout<<"   "<<sigma<<" +- "<<sigerr<<"    "<<tmp->GetEntries()<<endl;
	RMSP->SetBinContent(/*VarBin[i]+(VarBin[i+1]-VarBin[i])/2.*/i+1, /*tmp->GetRMS()*/ sigma );	
	
      }
      RMSP->SetBinError(i+1,/*tmp->GetRMSError()*/ sigerr );

    }

   
  }

  //  if(!cplot)
  // cc->Close();
  
  return RMSP;
}

/*
TGraphErrors* HistoManager::ConvertGraphToRMSGraphVarBin(TGraph* graph, int n, double VarBin[], string name, int nh) {

  TGraphErrors* RMSP = new TGraphErrors(); 
  int kk=0;
  // TCanvas* c=new TCanvas("d","d",300,300,600,600);
  for(int i=0;i<n;i++) {
    TH1F* tmp = ConvertGraphBinToHisto(graph,40,-100,100,VarBin[i],VarBin[i+1],name, nh);
    TF1* Gaus = new TF1("f1","gaus",-50,50);
    if(tmp->GetEntries()>20) {
      tmp->Fit("f1","QN0R");
      RMSP->SetPoint(kk,VarBin[i], Gaus->GetParameter(2) );
      RMSP->SetPointError(kk,(VarBin[i+1]-VarBin[i])/2., Gaus->GetParError(2) );
      kk++;
      cout<<" coin "<<VarBin[i]<<"  "<<Gaus->GetParameter(2)<<"  "<<Gaus->GetParError(2)<<endl;
      }
  }
  cout<<RMSP->GetN()<<endl;
  return RMSP;
}

*/

/*
TProfile* HistoManager::ConvertGraphToProfileVarBin(TGraph* graph, int n, int VarBin[],
						  string name,string opt) {
  
  TProfile* tmp=new TProfile(name.c_str(),name.c_str(),n,VarBin);
   tmp->SetErrorOption(opt.c_str());
  double x,y;

  for(int i=0;i<graph->GetN();i++) {
    
    graph->GetPoint(i,x,y);
    
    if(y<0)
      cout<<" x "<<x<<"  y  "<<y<<endl;
    
    if(x!=0 && y!=0) {
      tmp->Fill(x,y);
    }
  }
  
 
  return tmp;
  
}*/
/*
vector<TH2F*> HistoManager::ConvertAll(vector<TGraph*> graphs2d, string obs, vector<string> names,string k) {

  vector<TH2F*> HistosP;

  vector<float> temp = Templates2D[ obs ];

  for(int unsigned i=0;i<graphs2d.size();i++) {
    HistosP.push_back(  ConvertGraphToHisto( graphs2d[i], (int)temp[0], temp[1], temp[2],
					    (int)temp[3], temp[4], temp[5], names[i]+k ) );
  }
  
  return HistosP;

  }*/
/*
vector<TProfile*> HistoManager::ConvertAllProfile(vector<TGraph*> graphs2d, string obs,
						vector<string> names,string k,string opt) {

  vector<TProfile*> HistosP;

  vector<float> temp = Templates2D[ obs ];

  for(int unsigned i=0;i<graphs2d.size();i++) {
    HistosP.push_back( ConvertGraphToProfile( graphs2d[i],(int)temp[0], temp[1], temp[2],names[i]+k,opt )  );
  }
  return HistosP;
  
}
*/
vector<TProfile*> HistoManager::ConvertAllProfileVarBin(vector<TGraph*> graphs2d, int n, int VarBin[],
						      vector<string> names,string k,string opt) {

  vector<TProfile*> HistosP;

  for(int unsigned i=0;i<graphs2d.size();i++) {
    // HistosP.push_back( ConvertGraphToProfileVarBin( graphs2d[i],n, VarBin,names[i]+k,opt )  );
  }
  return HistosP;
  
}
/*
vector<TH1F*> HistoManager::ConvertAllGraphBinToHisto(vector<TGraph*> graphs2d, int NbinX, double Xmin, double Xmax, double Binm, double BinM, string name) {

  vector<TH1F*> HistosBin;
  
  for(int unsigned i=0;i<graphs2d.size();i++) {
    HistosBin.push_back(ConvertGraphBinToHisto(graphs2d[i], NbinX, Xmin, Xmax, Binm, BinM,name, nt) );
  }

  return HistosBin;
}
*/
vector<TH2F*> HistoManager::Stack2DHistos(vector<TH2F*> Histos2d, vector<string> datanames, map<string,float> weights) {

  vector<TH2F*> StackedHistos;
  
  //Initialisation
  for(int unsigned i=0;i<datanames.size();i++) {
    TH2F* tmp = (TH2F*)Histos2d[i]->Clone();
    tmp->Scale(weights[ datanames[i] ]);
    StackedHistos.push_back(tmp);
  }

  //Stacking
  for(int unsigned i=0;i<datanames.size();i++) {
    for(int unsigned j=1;j<datanames.size();j++) {
      if(j>i) {
	StackedHistos[i]->Add( Histos2d[j], weights[ datanames[j] ]);
      }
    }
  }

  return StackedHistos;

}


Double_t HistoManager::BinomError(float Nt, double eff) {
  
  double error=0;
  if(Nt==0) return 1;
  error = sqrt(eff*(1-eff)/Nt) ;
  return error;
}


TGraphErrors* HistoManager::EffRej(TH1F* signal, TH1F* bck) {

  signal->Rebin(2);
  bck->Rebin(2);

  int nBins = signal->GetNbinsX() +1 ;

  TGraphErrors* EffRejtmp=new TGraphErrors(nBins);
   
  double effrej[2]={0,0};
  double Errors[2]={0,0};
  cout<<signal->Integral()<<"    "<<bck->Integral()<<endl;
  for(int i=0;i<nBins;i++) {
    
    effrej[0] =  ((double)signal->Integral(i,nBins)/
		  signal->Integral(0,1000) );
    effrej[1] =  ((double)bck->Integral(i,nBins)/
		    bck->Integral(0,1000) );
    //Integral()/Lumi
    Errors[0] = BinomError(signal->GetEntries(), effrej[0]  );
    Errors[1] = BinomError(bck->GetEntries(), effrej[1]  );
    cout<<" i "<<i<<"  -> "<<signal->GetBinCenter(i)<<"  : S "<<effrej[0]<<" B "<<effrej[1]<<"   "<<Errors[0]<<"   "<<Errors[1]<<endl;
     if(bck->Integral()!=0 && signal->Integral()!=0 ) {
       EffRejtmp->SetPoint(i,effrej[0],effrej[1]);
       EffRejtmp->SetPointError(i, Errors[0], Errors[1]);	
     }
  }
  
  return EffRejtmp;
}


TGraphErrors* HistoManager::SigniErr(TH1F* signal, TH1F* bck) {
  
  signal->Rebin(2);
  bck->Rebin(2);

  int nBins = signal->GetNbinsX();

  TGraphErrors* EffRejtmp=new TGraphErrors(nBins);
   
  double effrej[2]={0,0};
  double Errors[2]={0,0};
  cout<<signal->Integral()<<"    "<<bck->Integral()<<"  "<<nBins<<endl;
  for(int i=0;i<nBins;i++) {
    
    effrej[1] =  ((double)signal->Integral(i,nBins)/
		  signal->Integral() );
    effrej[0] =  ((double)signal->Integral(i,nBins)/
		  (double)bck->Integral(i,nBins) );
    
    Errors[0] = BinomError(signal->GetEntries(), effrej[0]  );
    Errors[1] = BinomError(bck->GetEntries(), effrej[1]  );
    //  cout<<" i "<<i<<"  : S "<<effrej[0]<<" B "<<effrej[1]<<"   "<<Errors[0]<<"   "<<Errors[1]<<endl;
     if(bck->Integral(i,nBins)!=0 && signal->Integral()!=0 ) {
       EffRejtmp->SetPoint(i,effrej[1],effrej[0]);
       EffRejtmp->SetPointError(i, Errors[1], Errors[0]);	
     }
  }
  
  return EffRejtmp;
}


void HistoManager::ComputeStatFromHisto(TH1F* histo, double& Em, double& EM, double& M) {

  M = histo->GetMean();
  double integral = histo->Integral(0,histo->GetNbinsX() );
  double pas = histo->GetBinWidth(3);
  double em=0, eM=0;
  bool Limf[2]={false,false};
  for(int i=0;i<histo->GetNbinsX();i++) {
    em = histo->Integral(0,i)/integral;
    eM = histo->Integral(i,histo->GetNbinsX() )/integral;
   
    if(!Limf[0] && em >=0.16)
      {Limf[0]=true;  Em=i*pas-100;}
    if(!Limf[1] && eM <=0.16)
      {Limf[1]=true; EM=i*pas-100;}
  }

}



void HistoManager::RatioDistributionsNorm(TH1F* h1, TH1F* h2, string Xtitle, double rangeX[2]) {

  //Normalize plots
  float Int1 = h1->Integral(0, h1->GetNbinsX()+1 );
  float Int2 = h2->Integral(0, h2->GetNbinsX()+1 );

  h1->Scale(1./Int1);
  h2->Scale(1./Int2);

  cRatio = new TCanvas("cRatio","Ratio",300,300,485,485);
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
  ratio->GetXaxis()->SetRangeUser(rangeX[0],rangeX[1]);

  double xm = rangeX[0];
  double xM = rangeX[1];
  if((rangeX[0]<ratio->GetXaxis()->GetXmin()) )
    xm=ratio->GetXaxis()->GetXmin();
  if((rangeX[1]>ratio->GetXaxis()->GetXmax()) )
    xM=ratio->GetXaxis()->GetXmax();

  TLine* line = new TLine(xm,1,xM,1);
  line->SetLineColor(2);
  line->SetLineWidth(2);
  ratio->GetYaxis()->SetTitle(" ratio (1V/2V)");
  ratio->GetXaxis()->SetTitle(Xtitle.c_str() );
  ratio->GetXaxis()->SetNdivisions(5,5,0);
  ratio->Divide(h2);
  ratio->Scale( 1 );
  ratio->Draw();
  line->Draw("same");

  cout<<" Ratio on "<<Xtitle<<" made "<<endl;

  float Bilan=0;
  for(int i=0;i<ratio->GetNbinsX()+2;i++)
    Bilan += ratio->GetBinContent(i);

  cout<<" Bilan : "<<Bilan<<endl;


}


TH1F*
HistoManager::ShapeWeight(TH1F* href, TH1F* hvar) {
  
  TH1F* result = (TH1F*)href->Clone();
 
  result->Divide(hvar);

  // for(int i=0;i<href->GetNbinsX();i++)
  //   cout<<i<<"  "<<href->GetBinContent(i+1)
  // 	<<" / "<<hvar->GetBinContent(i+1)
  // 	<<" =  "<<result->GetBinContent(i+1)<<endl;
	
  return result;
}


TH2F*
HistoManager::ShapeWeight2D(TH2F* href, TH2F* hvar) {
  
  TH2F* result = (TH2F*)href->Clone();
 
  result->Divide(hvar);

  for(int i=0;i<href->GetNbinsX();i++) {
    for(int j=0;j<href->GetNbinsY();j++) {
      cout<<i<<"  "<<j<<"  "<<href->GetBinContent(i+1,j+1)
	  <<" / "<<hvar->GetBinContent(i+1,j+1)
	  <<" =  "<<result->GetBinContent(i+1,j+1)<<endl;
    }
  }
  
  return result;
}


vector<TH1F*>
HistoManager::ReduceHistoToNoEmptyBin(TH1F* source,TH1F* MC) {

  float pas = source->GetBinWidth(1);
  
  vector<vector<float> > coord;
  vector<float> tmpc(2,0);
  vector<float > errors;

  vector<vector<float> > coord2;
  vector<float > errors2;

  for(int i=1;i<source->GetNbinsX()+1;i++) {

    if(source->GetBinContent(i)!=0) {
      //Data
      tmpc[0] = source->GetBinCenter(i);
      tmpc[1] = source->GetBinContent(i);
      coord.push_back(tmpc);
      errors.push_back( source->GetBinError(i) );

      //MC
      tmpc[0] = MC->GetBinCenter(i);
      tmpc[1] = MC->GetBinContent(i);
      coord2.push_back(tmpc);
      errors2.push_back( MC->GetBinError(i) );
    }
  }

  int NBin = (int)((coord[ coord.size()-1 ][0]-coord[0][0])/pas +1);
  // cout<<coord[0][0]<<"    "<<coord[ coord.size()-1 ][0]<<"   "<<pas<<"   "<<NBin<<endl;
  TH1F* dataClean = new TH1F("d2","",NBin,coord[0][0],coord[ coord.size()-1 ][0]);
  TH1F* MCClean = new TH1F("d2","",NBin,coord[0][0],coord[ coord.size()-1 ][0]);

  for(int unsigned i=0;i<coord.size();i++) {
    dataClean->SetBinContent( (int)(coord[i][0]/pas)+1 , coord[i][1] );
    MCClean->SetBinContent( (int)(coord2[i][0]/pas)+1 , coord2[i][1] );
    
  }
 for(int unsigned i=0;i<coord.size();i++) {
   dataClean->SetBinError( (int)(coord[i][0]/pas)+1 , errors[i] );
   MCClean->SetBinError( (int)(coord2[i][0]/pas)+1 , errors2[i] );
   
   // cout<<(int)(coord[i][0]/pas)<<"   "<<errors[i]<<"   "<<dataClean->GetBinError( (int)(coord[i][0]/pas) )<<endl;
 }

  vector<TH1F*> vtmp;
  vtmp.push_back(dataClean);
  vtmp.push_back(MCClean);
  return vtmp;
}


void
HistoManager::Get90PCInterval(TH1* histo,double &range0,double &range1) {


  double integral = histo->Integral(0, histo->GetNbinsX() +1 );

  double tmp=0;
  int bins[2]={0,10000};
  int dbin=10000;

  for(int ib=0;ib<histo->GetNbinsX() +1;ib++) {

    tmp=0;

    for(int ib2=ib;ib2<histo->GetNbinsX() +1;ib2++) {

      tmp =  histo->Integral(ib,ib2);

      if(tmp/integral>0.95) {
	//cout<<ib2<<"   "<<ib<<"   "<<dbin<<"  "<<endl;

	if(dbin>(ib2-ib) )
	  {
	    //	    cout<<ib2<<"   "<<ib<<"   "<<dbin<<"  "<<" coin "<<tmp<<"   "<<integral<<endl;
	    dbin = ib2-ib;
	    bins[0] = ib -1;
	    bins[1] = ib2 +1;
	  }

      }
    }
  }

  // double range[2]={ histo->GetBinCenter(bins[0]) , histo->GetBinCenter(bins[1]) };
  range0 = histo->GetBinCenter(bins[0]);
  range1 = histo->GetBinCenter(bins[1]);

  //  return range;


}

TH1* 
HistoManager::Get90HistoInterval(TH1* histo) {


  double integral = histo->Integral(0, histo->GetNbinsX() +1 );

  double tmp=0;
  int bins[2]={0,10000};
  int dbin=10000;

  for(int ib=0;ib<histo->GetNbinsX() +1;ib++) {

    tmp=0;

    for(int ib2=ib;ib2<histo->GetNbinsX() +1;ib2++) {

      tmp =  histo->Integral(ib,ib2);
      //    cout<<" ib "<<ib<<"   "<<" ib2 "<<ib2<<"   "<<tmp/integral<<"    "<<tmp<<"   "<<histo->GetBinContent(ib2)<<"   "<<integral<<endl;
      if(tmp/integral>0.95) {
	if(dbin>(ib2-ib) )
	  {
	    dbin = ib2-ib;
	    bins[0] = ib;
	    bins[1] = ib2;
	  }
      }
    }
  }
  //cout<<" biiiiiiiiin "<<bins[0]<<"   "<<bins[1]<<endl;
  TH1* tmph = (TH1*)histo->Clone();
  tmph->Reset("ICEM");
  for(int ib=bins[0];ib<bins[1];ib++)
    tmph->SetBinContent(ib, histo->GetBinContent(ib) );

  //  cout<<" tmp "<<histo->Integral()<<"    "<<tmph->Integral()<<endl;

  return tmph;
}


TH1D*
HistoManager::GetContamination2D(string obs, int ibin, double VarBin[]) {
  
  int Obs = FindNVar2D(obs);
  int Ndata= access2D(Obs,nt);
  TH2F* dataTMP = (TH2F*)Histos2D[ Ndata ]->Clone();

  float Int = dataTMP->GetEntries();
  
  float Nevts = Int*contamination2D;
  Obs = FindNVar2D(obs+"Control");
  Ndata= access2D(Obs,nt); 
  TH2F* dataTMP2 = (TH2F*)Histos2D[ Ndata ]->Clone();
 
  Int = dataTMP2->GetEntries();
 
  int Binm = dataTMP2->GetXaxis()->FindBin(VarBin[ibin]);
  int BinM = dataTMP2->GetXaxis()->FindBin(VarBin[ibin+1]);
    
  TH1D* tmp = (TH1D*)dataTMP2->ProjectionY("",Binm, BinM); 

  tmp->Scale( Nevts/Int );
   
  return tmp;

}

void 
HistoManager::ConfigureContamination(float contamination) {

  contamination2D = contamination;

}

void 
HistoManager::SaveCanvas(string path,TCanvas* c) {

  string name = path + "/" + Name + ".pdf";
  string name2 = path + "/" + Name + ".eps";
  string name3 = path + "/" + Name + ".png";

  c->SaveAs(name.c_str() );
  c->SaveAs(name2.c_str() );
  c->SaveAs(name3.c_str() );

}

void 
HistoManager::SaveRoot(string path,TCanvas* c) {
  string name =  path + "/" + Name + ".root";

  c->SaveAs(name.c_str() );

}



void 
HistoManager::SavePlots(string path,TCanvas* c, string advname) {

  string extension[4]={"pdf","eps","png","root"};

  if(advname!="")
    Name=advname;

  string name = path + "/"+extension[0]+"/" + Name + "."+extension[0];
  string name2 = path + "/"+extension[1]+"/" + Name + "."+extension[1];
  string name3 = path + "/"+extension[2]+"/" + Name + "."+extension[2];
  string name4 =  path + "/"+extension[3]+"/" + Name + "."+extension[3];
 
  FILE *test_;
  for(int i=0;i<4;i++) {
    TString dirname_ = path +"/"+extension[i];
    if(i==1)
      cout<<" Dirname "<<dirname_<<endl;
    test_=fopen( dirname_.Data(), "r" );
    if( test_==0 )
      {
	TString command_ = TString("mkdir -p ") + dirname_; 
	system( command_.Data());
	//	assert( system( command_.Data() )==0 );
      }
    else
      fclose( test_ );

  }
  
  c->SaveAs(name.c_str() );
  c->SaveAs(name2.c_str() );
  c->SaveAs(name3.c_str() );
  c->SaveAs(name4.c_str() );
  
}
