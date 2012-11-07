#include "core/UseTree.hh"
#include <algorithm>

#include "VarClassLQ.hh"

#include <TF2.h>
#include <TVirtualFitter.h>
#include <TColor.h>
#include <TRandom3.h>
//#include <TGraphErrors.h>


class METUseTree : public UseTree
{
public :
  METUseTree();
  void FillMETTree();

  void AddVariables();


  void SetAdditionnalTrees(string s , vector<string> vs );
  
private :
  
  void PrepareDatasets();
  void GetNProcEvent(string dataset);
  
  TH1F* GetVtxCorrection(string obs,int nds, int bin);
  vector<vector<vector<float> > >  GetFitWeight();
  
  void PrepareHistograms();
  void FillPlots( VarClassLQ* vc, string anlvl, int i, float Weight );

  void LoadPUWeights();
  float SearchWeight(float trueNint);
  TH1F* puweights;

  // void StoreSelectedEvents(float mZZgen, float ZPt, float w);
  // std::vector< std::pair<float, float> > selZZ;
  // std::vector< float > selZPt;
  // vector<TF2*> LoadDBaTGC();

  TTree* GetAdditionnalTree(string type, int ds );
  void AddTreeLoading( string type, vector<string> AddTNames );
  std::map< string, vector<TChain*> > AddChains;

  TVector2 phiCorrection(TVector2 met);
  void fillMETHistos(string name, TVector2 met, TVector2 qT, int nvtx, float sumEt, int i, float Weight);
  
  void fillMETUnc(string name, TVector2 met, TVector2 qT, int i, float Weight, VarClassLQ* vc);
  map<string, pair<string,string> > uncMap;
  map<string, pair<string,string> >::const_iterator itUM;
  void addUnc(string tn, string n, string d);
  TVector2 mettmp,mettmpUnc;
  float metUnc,angleUnc;

  void prepareMETHistos(string name);
  void addMETname(string n, string p);
  vector<std::pair<string,string> > mets;
  TVector2 phiCorrection(TVector2 met, int Nvtx, bool isD,int run);

  TH1F* vtxW;
  void LoadVtxReweight();
  float GetVtwW(int nvtx);
  //  void DisableBranchs(TTree* tree);

  TLorentzVector l1,l2,l1U,l2U,Zv4;
  TRandom3 rnd;

  //Variables

  // vector<string> Type;

  // bool UsePdgId;
  // int pdgId;
  // bool applyBtag1;
  // bool applyBtag2;
  // float ptLep;
  // float ptTau;
  // int chrReq;

  //bool skimming;


ClassDef(METUseTree,0)
};
