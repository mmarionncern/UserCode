#include "core/UseTree.hh"
#include <algorithm>

class WZUseTree : public UseTree
{
public :
  WZUseTree();
  void FillWZTree();
  // void SpecialZZVariables(bool Response);

  void PlotEffRej(string Obs,string addObs="");

  bool MakeSkim;

private :

  void PrepareDatasets();
  TH1F* GetVtxCorrection(string obs,int nds, int bin);
  float SearchWeight(float Wpt);

  vector<vector<vector<float> > >  GetFitWeight();


  vector<vector<float> > DBWeights;

  vector<string> Type;
  

  ClassDef(WZUseTree,1)
};
