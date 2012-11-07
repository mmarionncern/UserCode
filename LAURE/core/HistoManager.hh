#ifndef _HistoMan_
#define _HistoMan_


#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <math.h>

#include "TF1.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

using namespace std;
typedef std::pair<string,string> legend;

typedef map<string,TH1F*> systM;
typedef map<string,TH1F*>::const_iterator itSystM;

class HistoManager {

private:

 
  
  std::map<string,int> Variables;
  std::map<string, int>::const_iterator itVar;
  std::map<string, legend> VarLegend;
  std::map<string, legend>::const_iterator itVarLeg;
  
  std::map<string,int> Variables2D;
  std::map<string, int>::const_iterator itVar2D;
  std::map<string, legend> VarLegend2D;
  std::map<string, string> VarName2D;
  

  std::map<string,int> VariablesP;
  std::map<string, legend> VarLegendP;
  std::map<string, string> VarNameP;
  //  std::map<string, vector<float> > TemplatesP;
  std::map<string, double* > VarTemplates;
  std::map<string, double* >::const_iterator itVarTP;

  std::map<string, std::pair<double*, double*> > VarTemplates2D;
  std::map<string, std::pair<double*, double*> >::const_iterator itVarT2D;

  std::map<string,int> Unc;

  std::map<string, vector<float> > Templates;
  std::map<string, vector<float> >::const_iterator itTemp;

  vector<vector<float> > templates;
  vector< vector<float> > Templates2D;
  vector< vector<float> > TemplatesP;

  vector<TH1F*> Histos;
  vector<TH2F*> Histos2D;
  vector<TProfile*> Profiles;

  vector<int> Nent2D;

  map<TGraph*, vector<float> > Weight2D;

  int nt;
  int Nunc;
  int Nvar;
  int Nvar2D;
  int NvarP;

  TCanvas* cRatio;

  float contamination2D;

  bool _initsequence;

public:
  string Name;

public:

  HistoManager();
  
  void StartFilling() {_initsequence=false;};
  bool GetInitStatus() {return _initsequence;};
  void InitStatus() {_initsequence=true;};

  void PrintVariables();
  vector<string> GetVariables();
  int FindNVar(string variable);
  int FindNUnc(string variable);
  int access(int nvar, int iDB);
  int accessUnc(int nvar,int Nu);
  legend FindLeg(string variable);

  int FindNVar2D(string variable);
  int access2D(int nvar, int iDB);
  legend FindLeg2D(string variable);
  
  int FindNVarP(string variable);
  int accessProf(int nvar, int DS);
  legend FindLegP(string variable);

  TH1F* getHisto(string obs, int DS,bool unc=false);

  systM FindSysts(string var,string type="");
  

  void ConfigAnalysis(int N);
  void AddVariable(string var, int nBin, double min, double max,
		   string Xleg, string Yleg);
  void PrepareHistos(vector<string> Datasets);
  void fill(string var, int DS , float value);
  void fill(string var, int DS , float value, float weight);

  void fill(string var, string type, int DS , float value, float weight=1.,string dir="");

  void Add2DVariable(string var, string Xleg, string Yleg,string name,int NBinsX, float Xm,float XM, int NBinY, float Ym, float YM);
  void Add2DVariable(string var, string Xleg, string Yleg,string name,int NBinsX, double varBinX[], int NBinY, double varBinY[]);
  void fill2D(string var, int DS , float valueX, float valueY );
  void fill2D(string var, int DS , float valueX, float valueY, float weight);
  void Prepare2DHistos(vector<string> Datasets);
  // void Cleaning2DHisto(vector<int> Nmax);

  void AddProfVariable(string var, int nBin, double min, double max,
		   string Xleg, string Yleg);
  void AddProfVariable(string var, int nBin, double varBin[],
		   string Xleg, string Yleg);
  void fillProf(string var, int DS , float valueX, float valueY );
  void fillProf(string var, int DS , float valueX, float valueY, float weight);
  void PrepareProfiles(vector<string> Datasets);
  
  void fillUnc(string var, string type , int DS , float value, float weight, string dir );
  void PrepareUncHisto(string var, int nb, float min, float max);



  static TH2F* ConvertGraphToHisto(TGraph* graph, int NbinX, double mX, double MX,
			    int NbinY, double mY, double MY,string name);
  static TH2F* ConvertGraphToHisto(TGraph* graph, int NbinX, double BinX[],
				   int NbinY, double BinY[],string name);
  /*  vector<TH2F*> ConvertAll(vector<TGraph*> graphs2d, string obs,
      vector<string> names,string k);*/
  static TProfile* ConvertGraphToProfile(TGraph* graph, int NbinX, double mX, double MX,
					 string name,string opt);
  /*  vector<TProfile*> ConvertAllProfile(vector<TGraph*> graphs2d, string obs,
      vector<string> names,string k,string opt);*/
 
  vector<TGraph*> ConvertGraphToErrorProfile(TGraph* graph, int NbinX, double mX, double MX );

  vector<TGraph*> ConvertGraphToErrorProfile(TGraph* graph, int NbinX, int BinVar[]);

  /*  TH1F* ConvertGraphBinToHisto(TGraph* graph, int NbinX, double Xmin, double Xmax,
			       double Binm, double BinM,string name, int nh);
  TH1F* ConvertGraphBinToHisto(TGraph* graph, int NbinX, double Xmin, double Xmax,
  double Binm, double BinM,string name,int& meanx, int nh);*/
  /* TProfile* ConvertGraphToResponsePlotVarBin(TGraph* graph, int n, double VarBin[], string name, int nh );*/
  // TH1D* ConvertGraphToRMSPlotVarBin(TGraph* graph, int n, double VarBin[], string name, int nh, bool rmConta);
  TH1D* ConvertHistoToRMSPlotVarBin(TH2F* histo2D, int n, double VarBin[], string name, bool rmConta,bool cplot);
  /* TGraphErrors* ConvertGraphToRMSGraphVarBin(TGraph* graph, int n, double VarBin[], string name, int nh);*/

  static TGraphAsymmErrors* GraphReduction(TGraph* graph, int NbinX, double mX, double MX,
					   string name,string opt);

  static TGraphAsymmErrors* GraphReduction(TGraph* graph, int NbinX, int BinVar[],
					   string name,string opt);

  static vector<TH1F*> ReduceHistoToNoEmptyBin(TH1F* source,TH1F* MC);

  // static TProfile* ConvertGraphToProfileVarBin(TGraph* graph, int n, int VarBin[],
  //					 string name,string opt);
  vector<TProfile*> ConvertAllProfileVarBin(vector<TGraph*> graphs2d, int n, int VarBin[],
				      vector<string> names,string k,string opt);
  /*  vector<TH1F*> ConvertAllGraphBinToHisto(vector<TGraph*> graphs2d, int NbinX, double Xmin,
					  double Xmax, double Binm, double BinM,
					  string name);*/

 vector<TH2F*> Stack2DHistos(vector<TH2F*> Histos2d,
					 vector<string> datanames,
					 map<string,float> weights);

  void RatioDistributionsNorm(TH1F* h1, TH1F* h2, string Xtitle, double rangex[2]);
  static TH1F* ShapeWeight(TH1F* href, TH1F* hvar);
  static TH2F* ShapeWeight2D(TH2F* href, TH2F* hvar);
  vector<TH1F*> GetHistos();
  vector<float> GetTemplate(int var);

  vector<float> GetTemplate2D(int var);
  vector<TH2F*> GetHistos2D();

  vector<float> GetTemplateP(int var);
  double* GetVarTemplateP(string var);
  vector<TProfile*> GetProfiles();

  static Double_t BinomError(float Nt, double eff);

  static void ComputeStatFromHisto(TH1F* histo, double& Em, double& EM, double& M);

  TGraphErrors* EffRej(TH1F* signal, TH1F* bck);
  TGraphErrors* SigniErr(TH1F* signal, TH1F* bck);

  TH1D* GetContamination2D(string obs, int ibin, double VarBin[]);
  void ConfigureContamination(float contamination);

  void Get90PCInterval(TH1* histo, double& range0, double& range1);
  TH1* Get90HistoInterval(TH1* histo);

  void SaveCanvas(string path,TCanvas* c);
  void SaveRoot(string path,TCanvas* c);
  void SavePlots(string path,TCanvas* c,string advname="");

  
  ClassDef(HistoManager,0)

};



#endif
