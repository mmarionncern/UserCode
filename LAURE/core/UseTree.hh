#ifndef _UseTree_
#define _UseTree_


#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <map>
#include <math.h>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TObjString.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPad.h"
#include "TMath.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TEllipse.h"
#include "TArrow.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TPolyLine.h"
#include "TMarker.h"
#include "RooHistError.h"
#include "RooHist.h"

#ifndef _HistoMan_
#include "HistoManager.hh"
#endif

#include <time.h>

using namespace std;


class UseTree {

  //variables

protected:

string XTitle;
string YTitle;

  float _IsoCuts[2][4][3];
  float _IDCuts[2][4][4];
  float _IDcuts[2][4];
  float _Isocuts[2][3];
  float MTCut;
  float METcut;
  float PTcut;

int Xdiv[3];
int Ydiv[3];

  vector<string> SevObs;

  int METType;
  string MEtType;
  int Bin;
  int BinBkg;
  double RangeY[2];
  double RangeX[2];
  vector<int> Lines;
  bool Draw3on1;

  bool logYScale;

  int ISO; int ID;

  

  //
  int nt;
  string reposi;
  vector<string> datasets;
  string snametmp;
  vector<int> dt;
  vector<int> colors;
  vector<string> name;
  vector<int> NEntries;
  vector<TFile*> files;
  vector<TTree*> trees;
  vector<TChain*> tChains;
  vector<TH1F*> Histos;
  vector<TH2F*> Histos2D;
  vector<TProfile*> Profiles;
  TLegend* leg;
  TLegend* leg2D;
  
  bool saveFile;
  
  double Lumi;
  bool EventFilter;
  int EventNum;

  bool NoData;
  string QCDType;
  string WType;
  bool QCDAutoWeight;
  bool FitR;
  bool getRatio;

  vector<string> data;
  
  double QCDweight;
  double normQCDWeight;
  bool MultiWeight;
  bool Remove2V;
  bool ShapeWeight;
  bool Norm;
  bool NormPlot;

  bool ShowDMCRatio;

  bool AddSyst;

  int NVert;

  bool response;

  bool OverFlowBin;
  bool UnderFlowBin;

  bool defStyle;

  bool unbinned;
  bool SuperColor;
  bool switchRMS;
  string RMSOpt;
  int Wline;

  float LineWidth;
  float MarkerSize;

  bool FillProf;
  bool BasicProf;
  
  int lumi;
  vector<bool> Weighted;
  vector<string> Suffix;

  map< std::pair<int,int> , std::pair<string,string> > Events;
  map< std::pair<int,int> , std::pair<string,string> >::iterator EventIter;
  vector<int> EvtsInFile;
  vector<float> METs;

  TH1F* savedMC;

  string GetDataset(int channel, int evt);

  //Weighting
  float weight;
  void FillAddWeight(string dataset);
  void FillEventMap(string dataset, int NevtB, int Chan );
  float GetWeight(int channel, int evt);
  vector< vector<std::pair<string,int> > > ChanNum;
  map<string,float> KFactors;
  map<string,float> Luminosities;
  map<string,float> XSections;
  map<string,float> internalNEvt;
  map<string, float> Weights;
  bool useXSect;
  //std::pair<string,string> ShapeVar;
  vector<TH1F*> shape;
  TH2F* shape2D;
  
  //Errors
  vector<double> VtxCorSyst;

  //Efficiencies
  vector<std::pair<string,float> > TotalMC;
  vector<std::pair<string,float> > TotalMCError;
  vector< vector<std::pair<string,float> > > effTotal;
  vector< vector<std::pair<string,float> > > TotalN;
  vector< vector<std::pair<string,float> > > TotalNWeights;
  vector< vector<std::pair<string,float> > > effMap;
  
  //Acceptance
  bool inAcc;
  bool useAccForEff;

  //Systematics
  vector< vector<std::pair<string,vector<float> > > > systMap;
  vector< vector<std::pair<string,vector<float> > > >::const_iterator iterSMap;
  
  //
  string observable;
  int EcalP;
  bool convRej;
  bool vetoE;
  bool noDeltaCor;
  bool skipCut;
  bool invCut;
  string Nm1Var;
  
  vector<string> Sorder;


  //misc
  bool _isData;

public:
  vector<string> IDs;
  vector<string> Isos;

public:

  UseTree();
  HistoManager histoManager;
  static void cmsPrel(double lumi);
  bool Cbool(bool skip, bool namevar);
  void DrawLine(string variable, int WP, int EP,int RangeY);
  void SavePlot(bool sfile);
  void ConfigureIDIso(int iso, int id);
  void ConfigurePlots(int Binning, int AddBinBkg, bool overflow, bool underflow, double rangeY[2], double rangeX[2],string ytitle, int XDiv[3], int YDiv[3], float markerS, float lineS,vector<int> lines,vector<string> sevObs,  bool Draw3o1=false, bool sucolor=false,bool switchrms=false, string opt="", bool FillProfile=false,bool basicProf=true, bool ShowRatio=false, bool logY=false,bool addSyst=false);
  void ConfigureEndAnalysis(bool ForceCut,int iso, int id,float IdCuts[2][4], float IsoCuts[2][3],float PtCut, float MetCut, string MetType,float MTcut,int nvert, string Normalisation);
  void EndAnalysis( vector<string> data, string observable,
		    int EcalP, bool convRej, bool vetoE, bool Invc,
		    bool SkipCut, string SkipVar, bool noDeltaCor,
		    string QCDtype,string Wtype, string rep);
  
  void AddMCSample(string str, string sname, int col );
  void SetTreeName(std::string);

  // void AddSample(string str, int n);

  void PrintChanNum();
  std::string FindProcess(int event, int ds );

  void Plot2DHisto(string observable,bool fillLeg, bool removelabel);
  void PlotRMS(string observable,bool fillLeg, bool removelabel,int color);
  void Plot1DHisto(string observable,bool fillLeg, bool removelabel, int color, int color2);
  void PlotProfile(string observable, bool fillleg, bool removelabel, int color);
  void PlotDataMCRatio(string obs, bool removelabel);
  void PlotRatio(string var1, string var2, string ds);
  // void PlotHistoFromGraph(string observable, bool FillLeg, bool removelabel, int color, int color2);

  TTree* LoadChain(int dnum);

  void PlotDistribution(string observable, bool Fill2DHistos,bool FillProfile );
  virtual void PlotEffRej(string Obs, string ssig, string Obs2="", string Obs3="");
  virtual void PlotSigRej(string Obs, string ssig, string Obs2="", string Obs3="");
  // void PlotRMS(string Obs,int sample);
  map<string,float> Reweight(int Obs, bool OneDim);
   void Delete();
  void PrintSelectedEvents();
  void SaveSelectedEvents();
  TCanvas* c2;

  //Statistics and Errors
  float Kurt(vector<double>,double mean , double sigma);
  void KolmogorovTest(TH1F* MC, TH1F* data, string observable);
  float Poisson(float mc, int data);
  void Chi2Test(TH1F* MC, TH1F* data);

  static float ErrorPL(int N, int sigma=1);
  static float ErrorPH(int N, int sigma=1);

  //kinematics
  static float dR(float e1, float e2, float p1, float p2, bool conv=false);
  static float phi( float x, float y );
  static float dPhi( float phi1, float phi2 );

  void ParserNorm(string str);

  float FitReweight(vector<TH1F*> histos, TH1F* data);
  vector<float> GetAlpha(TString obs);
  virtual vector<vector<vector<float> > >  GetFitWeight() =0;
  float GetWeightFromShape(int met, float var, string shref, string shvar);
  virtual float GetMCWeightFromDataShape(int met, float var, string shref,
					 string shvar, double bin, int ds);
  virtual float GetMCWeightFromDataShape2D(float var, float var1, string shref, string shvar);
  virtual TH1F* GetVtxCorrection(string obs, int nds, int bin) =0;
  
  void AddSystematics(string observable, TH1F* Data, float W);
  TGraphAsymmErrors* ComputeSystematics(string obs, TH1F* totMC);
  TPolyLine* GetSystBand(vector<float> xs, vector<float> yl, vector<float> yh, float xmin, float xmax);

  static float LnL(float si, float bi, int ni, int NData, float a);
  static float Chi2(float si, float bi, int ni, int NData, float a);
  static float Fact(int x);
  TPolyLine* PolyLineFromGraph(vector<TGraph*> graphs);
  TPolyLine* PolyLineFromGraph(TGraph* graph);

  // void FillWGen(string observable="tcAbsRecoilvsGen", string observable2="pfAbsSpecRecoilvsGen");

  void ConfigureLumi(std::map<string,float> Lum ,std::map<string,float> Kfac, float l,bool useXS=false );
  void ConfigureData(bool, int,bool MCOnly);

  //Cuts
  template < typename T > inline
  bool MakeCut( T value, T valcut, string type, int i, string cName, float w, T seccut=0, bool noRegister=false) {

    bool accept;
  
    if(!Cbool(skipCut, cName==Nm1Var) )
      { return true; }

    if( !Cbool( invCut, cName==Nm1Var) ) {
      type = InvCut(type); }
 
    if(type=="<") {
      accept = (value < valcut);
    }
    else if(type=="<=") {
      accept = (value <= valcut);
    }
    else if( type==">") {
      accept = (value > valcut);
    }
    else if( type==">=") {
      accept = (value >= valcut);
    }
    else if( type=="=") {
      accept = (value == valcut);
    }
    else if(type=="!=") {
      accept = (value != valcut);
    }
    else if(type=="[]") {
      accept = (value >= valcut && value<= seccut );
    }
    else if(type=="][") {
      accept = (value > valcut && value< seccut );
    }
    else if(type=="[!]") {
      accept = !(value >= valcut && value <= seccut );
    }
    else if(type=="]![") {
      accept = !(value > valcut && value < seccut );
    }
    else {
      accept =false; cout<<" Warning cut "<<type<<endl;
    }  
   //  if(cName!="")
//     cout<< cName<<" "<<accept<<"  "<<noRegister<<endl;
    if(!noRegister) {
      SetSystematics( type, i, cName, w);
      SetEfficiency(i, cName, w, accept);
    } 

    if(histoManager.GetInitStatus()) {
      accept=true;
    }

    return accept;
  };

  template < typename T > inline
  bool SimpleCut( T value, T cut, string type, T seccut=0 ) {
    if(histoManager.GetInitStatus())
      return false;

    return MakeCut<T>( value, cut, type, -1, "", 0, seccut, true );
  };



  std::string InvCut(string s);

  // bool MakeGeneralCutF( float value, float valcut, string type,
  // 		       int i, string cName, float w);
  // bool MakeGeneralCutI( int value, int valcut, string type,
  // 		       int i, string cName, float w);
  // bool MakeGeneralCutB( bool value, bool valcut, string type,
  // 		       int i, string cName, float w);

  //Efficiencies and yields
  void SetEfficiency(int i, string dataset, float value, bool acc);
  void SetSystematics( string type, int i, string cName, float w);
  void GetNumbers();
  void DrawNumbers();

  //Integral
  void Integral(string Obs, double xmin, double xmax);

  //Name
  int GetNbyName(string sname);

  void SaveAllPlots(string path, bool erase=false);
  //  void SavePlots(string path,TCanvas* c);

protected:
  virtual void PrepareDatasets();
  virtual void GetNProcessedEvent();

protected:
  vector<string> treenames;


private:

  vector<TPad*> PreparePads(int npad);
  vector<vector<TPad*> > PreparePadsWithRatio(int npad);



protected: //skimming
  void InitSkimming(int i, int ie, int nent);
  void FinalizeSkimming();
  TFile* oFile;
  TTree* skimTree;
  string wDs;
  bool _storeTuple;
  

  ClassDef(UseTree,0)
};


#endif
