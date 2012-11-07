#include "core/UseTree.hh"
#include <algorithm>


#include "plugins/NeuralNet.hh"
#include "plugins/METNN.hh"


class WWCCUseTree : public UseTree
{
public :
  WWCCUseTree();
  void FillZZTree();
  // void SpecialZZVariables(bool Response);

  // void PlotEffRej(string Obs,string addObs="");

  bool MakeSkim;

  void AddVariables(int pid,int mj,int mp,int ml, float dPhijmc, float NMcutH, float NMcutL, float nncut, bool fj);

  void Get2DNumber(string obs, string chan, ostream& out);
  void GetAll2DNumbers();
private :

  int pdgID;
  int MaxJet;
  int MaxPhoton;
  int MaxLepton;
  float dPhiJMCut;
  float NormMETCutHigh;
  float NormMETCutLow;
  float NNoutCut;

  bool fixJet;

  float Rapidity(float p, float mass, float pt);

  void PrepareDatasets();
  TH1F* GetVtxCorrection(string obs,int nds, int bin);
  float SearchWeight(float Wpt);

  vector<vector<vector<float> > >  GetFitWeight();


  vector<vector<float> > DBWeights;
  vector<string> Type;
  
  float GetNSigmaPara(float ZPt, float RPara, int Nv, bool isData);
  float GetNSigmaPerp(float ZPt, float RPerp, int Nv, bool isData);
  float Response(float Zpt);

  IClassifierReader* CreateNN();
  IClassifierReader* CreateMETNN();

ClassDef(WWCCUseTree,0)
};
