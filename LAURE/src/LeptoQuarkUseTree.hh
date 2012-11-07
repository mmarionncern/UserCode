#include "core/UseTree.hh"
#include <algorithm>

#include "VarClassLQ.hh"

#include <TF2.h>
#include <TVirtualFitter.h>
#include <TColor.h>
//#include <TGraphErrors.h>


class LeptoQuarkUseTree : public UseTree
{
public :
  LeptoQuarkUseTree();
  void FillLQTree();

  void AddVariables(/*bool usepid, int pid*/int nbtag, float ptlep, float pttau, int chr, bool skim);


  void SetAdditionnalTrees(string s , vector<string> vs );
  
  void ComputeNFakes();


private :

    void PrepareDatasets();
  void GetNProcEvent(string dataset);

  TH1F* GetVtxCorrection(string obs,int nds, int bin);
  
  vector<vector<vector<float> > >  GetFitWeight();


  bool isMatch(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2);
  
  float Rapidity(float p, float mass, float pt, float eta);

  int MuonClassification( double pt, double eta, int isID );

  int ElectronClassification( double pt, double eta, double sieie, double deta, double dphi, double hoe );

  int TauClassification(double pt, double eta);

  int JetClassification( double pt, double eta, double phi,
			 double l1Eta, double l1Phi,
			 double l2Eta, double l2Phi );


  bool HLTRequirement(VarClassLQ* vc, int nds);
  bool ElectronID(VarClassLQ* vc);
  bool VtxID(VarClassLQ* vc);
  bool TauID(VarClassLQ* vc);
  bool JetID(VarClassLQ* vc);

  void PrepareHistograms();
  void FillPlots( VarClassLQ* vc, string anlvl, int i, float Weight );

 

  void StoreSelectedEvents(float mZZgen, float ZPt, float w);
  std::vector< std::pair<float, float> > selZZ;
  std::vector< float > selZPt;
  vector<TF2*> LoadDBaTGC();

  TTree* GetAdditionnalTree(string type, int ds );
  void AddTreeLoading( string type, vector<string> AddTNames );
  std::map< string, vector<TChain*> > AddChains;

  void DisableBranchs(TTree* tree);

  //Variables

  vector<string> Type;

  bool UsePdgId;
  int pdgId;
  bool applyBtag1;
  bool applyBtag2;
  float ptLep;
  float ptTau;
  int chrReq;

  bool skimming;


ClassDef(LeptoQuarkUseTree,0)
};
