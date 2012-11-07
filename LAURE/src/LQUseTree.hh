#include "core/UseTree.hh"
#include <algorithm>

#include "VarClassLQ2012.hh"


class LQUseTree : public UseTree
{
public :
  LQUseTree();
  void FillMETTree();

  void AddVariables(bool skim);


  void SetAdditionnalTrees(string s , vector<string> vs );
  
private :
  

  //Default functions needed for any analysis
  //(at least for back-compatibility)
  
  void PrepareDatasets();
  void GetNProcEvent(string dataset);
  
  TH1F* GetVtxCorrection(string obs,int nds, int bin);
  vector<vector<vector<float> > >  GetFitWeight();
  
  void PrepareHistograms();
  
  TTree* GetAdditionnalTree(string type, int ds );
  void AddTreeLoading( string type, vector<string> AddTNames );
  std::map< string, vector<TChain*> > AddChains;

  //====================================
  //user functions
  bool passLooseElectronID(size_t ie);
  bool passElectronID(size_t ie);
  float Aeff( float eta);
  bool passTauID(size_t it);
  bool passHLT();
  bool triggerDecision(size_t ihlt, string ref);
  bool passJetMETFilters();

  

  //user variables
  VarClassLQ2012 vc;

  bool skimming;


  ClassDef(LQUseTree,0)
};
