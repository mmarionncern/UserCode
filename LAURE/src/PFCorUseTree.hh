#ifndef _PFCorUseTree_
#define _PFCorUseTree_

#include "core/UseTree.hh"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"

using namespace std;

class PFCorUseTree : public UseTree
{
public :
  
  // ~PFCorUseTree()
  PFCorUseTree();
  //PFCorUseTree(int n);
  void FillTree();
  
 

private :

  void PrepareDatasets();
  
  //Spec kine computations
  float ComputeMinv(float Et1, float Et2, float eta1, float eta2, float phi1, float phi2);
  float ConversionEtaTheta(float eta);
  float ConversionThetaEta(float theta);
  float ConversionEpt(float energy,float eta );
  float ConversionPtE(float pt, float eta);
  void Conversion_REP_carte(float coord[3], float& x, float& y, float &z);

  float computeBR(float sigeta, float sigphi) {return sigphi/sigeta;};

  float GetCalibCorEnergy(float brLinear, float e, float eta);
  float GetCalibCorEB(float brLinear, float e, float eta);
  float GetCalibCorEE(float brLinear, float e, float eta);
 
  float GetBremCorEB(float brLinear, float et);
  float GetBremCorEE(float brLinear, float eta);
  
  float CalibCorEB(float et, float eta);
  float CalibCorEE(float et, float eta);

  TH1F* GetVtxCorrection(string obs, int nds, int bin);
  vector<vector<vector<float> > >  GetFitWeight();


  ClassDef(PFCorUseTree,0)
};

#endif
