#include "core/UseTree.hh"
#include <algorithm>


#include "plugins/NeuralNet.hh"
#include "plugins/METNN.hh"
#include <TF2.h>
#include <TVirtualFitter.h>
#include <TColor.h>
//#include <TGraphErrors.h>


class DibosonUseTree : public UseTree
{
public :
  DibosonUseTree();
  void FillZZTree();
  // void SpecialZZVariables(bool Response);

  // void PlotEffRej(string Obs,string addObs="");

  bool MakeSkim;

  void AddVariables(bool usepid, int pid,int mj,int ml, float dPhijmc, 
		    float dPhizmc, float NMcutL, float NMcutH, bool btag,
		    float advMm, float advMM,bool useAFEff, string tjv);

  void Get2DNumber(string obs, string chan, ostream& out);
  void GetAll2DNumbers(string var="ZMassVsSCPMET");
  void PrintGlobal2DNumbers(string var);

  void PlotAllVtx(string chan, string var="PF", int N=5);
  void DrawJetMultiMult(string chan);
  void CheckCharge(string chan);
  void CheckBtag();

  vector<vector<float> > GetReweightedYields(float f4M, float f5M, int nf4, int nf5);
  void  DrawContour(float lim95, float lim68=2, float limU=100000);
  void DrawATGCYields();
  
  void DrawATGCComponent(float f4, float f5);


  void SetAdditionnalTrees(string s , vector<string> vs );
  
  void ComputeNFakes();


private :

  vector<vector<float> > DBWeightsZZ;
  vector<vector<float> > DBWeightsWZ;
  vector<vector<float> > DBWeightsWW;
  vector<vector<float> > DBWeightsZZJV;
  vector<vector<float> > DBWeightsWZJV;
  vector<vector<float> > DBWeightsWWJV;
  void LoadDBWeightZZ();
  void LoadDBWeightWZ();
  void LoadDBWeightWW();
  float SearchWeightZZ(float Zpt, bool veto=false);
  float SearchWeightWZ(float Zpt, bool veto=false);
  float SearchWeightWW(float Wpt, bool veto=false);

  void LoadATGCWeights();
  float GetATGCWeight(float f4, float f5, int ibin);
  vector<map<vector<float>, vector<float> > > aTGCWeights;

  
  bool UsePdgId;
  int pdgID;
  int MaxJet;
  int MaxLepton;
  float dPhiJMCut;
  float dPhiZMCut;
  float balMin;
  float balMax;

  float advMassMin;
  float advMassMax;

  string typeJetV;
  bool bTagFlag;
  
  bool isMatch(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2);

  float Rapidity(float p, float mass, float pt, float eta);

  void PrepareDatasets();
  TH1F* GetVtxCorrection(string obs,int nds, int bin);
  
  vector<vector<vector<float> > >  GetFitWeight();


  //  vector<vector<float> > DBWeights;
  vector<string> Type;
 

  void StoreSelectedEvents(float mZZgen, float ZPt, float w);
  std::vector< std::pair<float, float> > selZZ;
  std::vector< float > selZPt;
  vector<TF2*> LoadDBaTGC();


  //Acceptance
  bool BasicAcc(float pt, float eta, int pdg );
  bool isInAcc( float ptV1, float ptV2, float mZ,
		const vector<double>* gLPt, const vector<double>* gLEta,
		const vector<double>* gLPdg, int offset=0 );
  
  float METCorrectionData(float met, int nVertex);
  float METCorrectionMC(float met, int nVertex);
  float METCorrectionMin(float min, int nVertex, bool isMC, int pdg);

  //Additionnal trees;
  void AddTreeLoading( string type, vector<string> AddTNames );
  std::map< string, vector<TChain*> > AddChains;
  TTree* GetAdditionnalTree(string type, int ds );


  //Fake rate
  TH1F* FakeDB;
  void LoadFakeRate();
  std::vector<float> getFRate(float pt, float eta);
  std::vector<std::vector<std::vector<float> > > LooseElec;


ClassDef(DibosonUseTree,0)
};
