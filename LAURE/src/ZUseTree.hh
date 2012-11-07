#ifndef _ZUseTree_
#define _ZUseTree_

#include "core/UseTree.hh"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"

using namespace std;

class ZUseTree : public UseTree
{
public :
  
  ZUseTree();
  void FillZTree();
  void SpecialZVariables(bool Response);

  void PlotResolutionPerVtx(string met, string type, bool MC=false);
  void PlotGlobalResolutionPerVtx( string type,bool isMC=false);
  void PlotGlobalResponse(bool isMC=false);
  // void FillZResponse(string observable, string observable2 );
  void SaveResolutionPerVtx(string met, string type,bool isMC, bool cor);
  void PlotResolutionUsingRooFit(int Nvtx, string type,string proj,bool isMC);

  void PrintResponseSummary();

  void GetRecoilCorrections(string comp, string met);

private :

  void PrepareDatasets();
  TH1F* GetVtxCorrection(string obs, int nds, int bin);


  vector<vector<vector<float> > >  GetFitWeight();

  vector<vector<float> > DBWeights;
  vector<vector<vector<float> > > MCResponseCor;
  vector<vector<vector<float> > > DataResponseCor;

  //Spec kine computations
  float ComputeMinv(float Et1, float Et2, float eta1, float eta2, float phi1, float phi2);
  float ConversionEtaTheta(float eta);
  float ConversionThetaEta(float theta);
  float ConversionEpt(float energy,float eta );
  float ConversionPtE(float pt, float eta);
  void Conversion_REP_carte(float coord[3], float& x, float& y, float &z);

  //Fraction 1vtx events
  float Frac1V;

  void LoadDBWeight();
  float SearchWeight(float Zpt);
  
  void LoadResponseDB();
  float SearchESWeight(float Zpt, int met);

  float Rapidity(float p, float mass, float pt, float eta);

  ClassDef(ZUseTree,0)
};


#endif
