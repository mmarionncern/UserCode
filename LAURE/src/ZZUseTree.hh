#include "core/UseTree.hh"
#include <algorithm>


#include "plugins/NeuralNet.hh"
#include "plugins/METNN.hh"
#include <TF2.h>
#include <TVirtualFitter.h>
#include <TColor.h>
//#include <TGraphErrors.h>


class ZZUseTree : public UseTree
{
public :
  ZZUseTree();
  void FillZZTree();
  // void SpecialZZVariables(bool Response);

  // void PlotEffRej(string Obs,string addObs="");

  bool MakeSkim;

  void AddVariables(bool usepid, int pid,int mj,int mp,int ml, float dPhijmc, float NMcutH, float NMcutL, float nnmcut, bool fj, bool invj, bool btag, bool invmet,float nncut, bool useAFEff);

  void Get2DNumber(string obs, string chan, ostream& out);
  void GetAll2DNumbers(string var="ZMassVsSCPMET");
  void PrintGlobal2DNumbers(string var);

  void PlotAllVtx(string chan, string var="PF", int N=5);
  void DrawJetMultiMult(string chan);
  void CheckCharge(string chan);
  void CheckBtag();

  static double FuncProjMET(double* x, double* par);

  vector<vector<float> > GetReweightedYields(float f4M, float f5M, int nf4, int nf5);
  void  DrawContour(float lim95, float lim68=2, float limU=100000);
  void DrawATGCYields();

private :

  vector<vector<float> > DBWeightsZZ;
  vector<vector<float> > DBWeightsWZ;
  vector<vector<float> > DBWeightsWW;
  void LoadDBWeightZZ();
  void LoadDBWeightWZ();
  void LoadDBWeightWW();
  float SearchWeightZZ(float Zpt);
  float SearchWeightWZ(float Zpt);
  float SearchWeightWW(float Wpt);

  void LoadATGCWeights();
  float GetATGCWeight(float f4, float f5, int ibin);
  vector<map<vector<float>, vector<float> > > aTGCWeights;

  
  bool UsePdgId;
  int pdgID;
  int MaxJet;
  int MaxPhoton;
  int MaxLepton;
  float dPhiJMCut;
  float NormMETCutHigh;
  float NormMETCutLow;
  float NNoutCut;
  float NNoutMETCut;

  bool fixJet;
  bool invJet;
  bool invMet;
  bool bTagFlag;
  
  float Rapidity(float p, float mass, float pt, float eta);

  void PrepareDatasets();
  TH1F* GetVtxCorrection(string obs,int nds, int bin);
  
  vector<vector<vector<float> > >  GetFitWeight();


  //  vector<vector<float> > DBWeights;
  vector<string> Type;
  
  float GetNSigmaPara(float ZPt, float RPara, int Nv, bool isData);
  float GetNSigmaPerp(float ZPt, float RPerp, int Nv, bool isData);
  float Response(float Zpt);

  void StoreSelectedEvents(float mZZgen, float w);
  std::vector< std::pair<float, float> > selZZ;
  vector<TF2*> LoadDBaTGC();


  //Acceptance
  bool BasicAcc(float eta, float pt );
  bool isInAcc(float genlep[2][5], float genZnn[5], float lPart[5]);

 
  

  IClassifierReader* CreateNN();
  IClassifierReader* CreateMETNN();

ClassDef(ZZUseTree,0)
};
