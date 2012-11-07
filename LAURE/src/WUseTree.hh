#include "core/UseTree.hh"

using namespace std;

typedef std::pair<float,float> pff;

class WUseTree : public UseTree
{
public :
  WUseTree();
  void FillWTree();
  
 
private :

  TH1F* GetVtxCorrection(string obs, int nds, int bin);
  void PrepareDatasets();
  void LoadDBWeight();
  void LoadDBResponse();
  void LoadRecoilCorrections();
  float SearchWeight(float Wpt);
  float SearchResponse(float Wpt,int met);
  float SearchRecoilCor(float GenWpt);

  float GetParaRecoilCor(int met, int Nvtx, float GenWqT);
  float GetPerpRecoilCor(int met, int Nvtx, float GenWqT);
  float GetParaRecoilCorRes(int met, int Nvtx, float GenWqT);
  float GetPerpRecoilCorRes(int met, int Nvtx, float GenWqT);
  
  vector<vector<vector<float> > >  GetFitWeight();

  vector<vector<float> > DBWeights;
  //vector<vector<pff> > ResponseDB; 
  vector<vector<vector<float> > > MCResponseCor;
  vector<vector<vector<float> > > DataResponseCor;
  vector<vector<float> > RecoilCors;

  //Fraction 1vtx events
  float Frac1V;

  float GetMCWeightFromDataShape(int met, float var, string shref, string shvar, double bin);

  ClassDef(WUseTree,0)
};
