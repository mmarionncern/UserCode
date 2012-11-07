#include "DibosonUseTree.hh"

#include <iomanip>
#include <TLorentzVector.h>


TGraphErrors* contour=new TGraphErrors(0);;
//Epplis
void Ellips(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t);


using namespace std;

ClassImp(DibosonUseTree)


DibosonUseTree::DibosonUseTree():
UseTree()
{

 float ZeIsoCuts[2][4][3]= { { {0.05, 0.06, 0.03}, {0.09, 0.07, 0.1},
			       {0.12, 0.09, 0.1}, {0.15, 2.0, 0.12} },
			     { {0.025, 0.025, 0.02}, {0.04, 0.05, 0.025},
			       {0.05, 0.06, 0.03}, {0.08, 0.06, 0.05} } };
  
 float ZeIDCuts[2][4][4]= { { {0.01,0.004,0.03,0.025}, {0.01,0.004,0.06,0.04},
			      {0.01,0.007,0.8,0.12}, {0.01,0.007,0.8,0.15} },
			    { {0.03,0.005,0.02,0.025}, {0.03,0.007,0.03,0.025},
			      {0.03,0.009,0.7,0.05}, {0.03,0.01,0.7,0.07} } };

 
  //ID Iso initialization
  for(int j=0;j<2;j++){
    for(int i=0;i<4;i++) {
      if(i<3)
	_Isocuts[j][i]=10;
      _IDcuts[j][i]=10;
    }
  }

  //ID Iso Initialisation
  for(int j=0;j<2;j++){
    for(int i=0;i<4;i++) {
      for(int k=0;k<4;k++) {
	if(k<3)
	  _IsoCuts[j][i][k] = ZeIsoCuts[j][i][k];
	_IDCuts[j][i][k] = ZeIDCuts[j][i][k];
      }
    }
  }

  //Force Remove 2vertex cleaning
  Remove2V = false;

  MakeSkim = false;


  LoadDBWeightZZ();
  LoadDBWeightWZ();
  LoadDBWeightWW();

  LoadATGCWeights();

  FakeDB=NULL;

}

void
DibosonUseTree::AddVariables(bool usepid, int pid,int mj, int ml, float dphijmc, float dphizmc,  float NMcutL, float NMcutH, bool btag, float advMm, float advMM, bool useAFEff, string tjv="<") {

  UsePdgId = usepid;
  pdgID = pid;
  MaxJet = mj;
  MaxLepton = ml;
  dPhiJMCut = dphijmc;
  dPhiZMCut = dphizmc;
  balMin = NMcutL;
  balMax = NMcutH;
  advMassMin = advMm;
  advMassMax = advMM;

  typeJetV = tjv;
  bTagFlag = btag;
  

  useAccForEff = useAFEff;

}

void
DibosonUseTree::FillZZTree() {

  //Add Variables to be computed

  //Prepare suffix for outfiles
  string suffix = "Multivar";
  if(switchRMS && FillProf)
    suffix = "RMS_METs";

  if(Draw3on1)
    for(int i=0;i<5;i++)
      Suffix[i] = suffix;

  // string metT[5]={"pf","tc","calo","caloT1","caloT2"};
  //string lep[2]={"l1","l2"};

  //Vertex
  histoManager.AddVariable("NVertex",20,0,20,"number of vertices","Vertex");
  histoManager.AddVariable("NVertexControl",20,0,20,"number of vertices","Vertex");
  
  //Z Plots
  histoManager.AddVariable("ZMassPresel",1000,0,1000,"M_{ll} [GeV]","ZMassPresel");
  histoManager.AddVariable("ZMassMET",400,0,200,"M_{ll} [GeV]","ZMassMET");
  histoManager.AddVariable("ZMassBal",400,0,200,"M_{ll} [GeV]","ZMassBal");
  histoManager.AddVariable("ZMassJVeto",400,0,200,"M_{ll} [GeV]","ZMassJVeto");
  histoManager.AddVariable("ZMassDPhiJ",400,0,200,"M_{ll} [GeV]","ZMassDPhiJ");
  histoManager.AddVariable("ZMassDPhiZ",400,0,200,"M_{ll} [GeV]","ZMassDPhiZ");
  histoManager.AddVariable("ZMassLVeto",400,0,200,"M_{ll} [GeV]","ZMassLVeto");
  histoManager.AddVariable("ZMassBTag",400,0,200,"M_{ll} [GeV]","ZMassBTag");
  histoManager.AddVariable("ZMassMass",400,0,200,"M_{ll} [GeV]","ZMassMass");
 
  histoManager.AddVariable("ZPtPresel",400,0,400,"q_{T} [GeV]","ZPtPresel");
  histoManager.AddVariable("ZPtMET",400,0,400,"q_{T} [GeV]","ZPtMET");
  histoManager.AddVariable("ZPtBal",400,0,400,"q_{T} [GeV]","ZPtBal");
  histoManager.AddVariable("ZPtJVeto",400,0,400,"q_{T} [GeV]","ZPtJVeto");
  histoManager.AddVariable("ZPtDPhiJ",400,0,400,"q_{T} [GeV]","ZPtDPhiJ");
  histoManager.AddVariable("ZPtDPhiZ",400,0,400,"q_{T} [GeV]","ZPtDPhiZ");
  histoManager.AddVariable("ZPtLVeto",400,0,400,"q_{T} [GeV]","ZPtLVeto");
  histoManager.AddVariable("ZPtBTag",400,0,400,"q_{T} [GeV]","ZPtBTag");
  histoManager.AddVariable("ZPtMass",400,0,400,"q_{T} [GeV]","ZPtMass");

  histoManager.AddVariable("ZYPresel",100,-5,5,"Y_{ll}","ZYPresel");
  histoManager.AddVariable("ZYMET",100,-5,5,"Y_{ll}","ZYMET");
  histoManager.AddVariable("ZYBal",100,-5,5,"Y_{ll}","ZYBal");
  histoManager.AddVariable("ZYJVeto",100,-5,5,"Y_{ll}","ZYJVeto");
  histoManager.AddVariable("ZYDPhiJ",100,-5,5,"Y_{ll}","ZYDPhiJ");
  histoManager.AddVariable("ZYDPhiZ",100,-5,5,"Y_{ll}","ZYDPhiZ");
  histoManager.AddVariable("ZYLVeto",100,-5,5,"Y_{ll}","ZYLVeto");
  histoManager.AddVariable("ZYBTag",100,-5,5,"Y_{ll}","ZYBTag");
  histoManager.AddVariable("ZYMass",100,-5,5,"Y_{ll}","ZYMass");

  //MET plots
  histoManager.AddVariable("METPresel",400,0,400,"#slash{E}_{T} [GeV]","METPresel");
  histoManager.AddVariable("METMET",400,0,400,"#slash{E}_{T} [GeV]","METMET");
  histoManager.AddVariable("METBal",400,0,400,"#slash{E}_{T} [GeV]","METBal");
  histoManager.AddVariable("METJVeto",400,0,400,"#slash{E}_{T} [GeV]","METJVeto");
  histoManager.AddVariable("METDPhiJ",400,0,400,"#slash{E}_{T} [GeV]","METDPhiJ");
  histoManager.AddVariable("METDPhiZ",400,0,400,"#slash{E}_{T} [GeV]","METDPhiZ");
  histoManager.AddVariable("METLVeto",400,0,400,"#slash{E}_{T} [GeV]","METLVeto");
  histoManager.AddVariable("METBTag",400,0,400,"#slash{E}_{T} [GeV]","METBTag");
  histoManager.AddVariable("METMass",400,0,400,"#slash{E}_{T} [GeV]","METMass");

  histoManager.AddVariable("CorMETPresel",400,0,400,"cor #slash{E}_{T} [GeV]","CorMETPresel");
  histoManager.AddVariable("CorMETMET",400,0,400,"cor #slash{E}_{T} [GeV]","CorMETCorMET");
  histoManager.AddVariable("CorMETBal",400,0,400,"cor #slash{E}_{T} [GeV]","CorMETBal");
  histoManager.AddVariable("CorMETJVeto",400,0,400,"cor #slash{E}_{T} [GeV]","CorMETJVeto");
  histoManager.AddVariable("CorMETDPhiJ",400,0,400,"cor #slash{E}_{T} [GeV]","CorMETDPhiJ");
  histoManager.AddVariable("CorMETDPhiZ",400,0,400,"cor #slash{E}_{T} [GeV]","CorMETDPhiZ");
  histoManager.AddVariable("CorMETLVeto",400,0,400,"cor #slash{E}_{T} [GeV]","CorMETLVeto");
  histoManager.AddVariable("CorMETBTag",400,0,400,"cor #slash{E}_{T} [GeV]","CorMETBTag");
  histoManager.AddVariable("CorMETMass",400,0,400,"cor #slash{E}_{T} [GeV]","CorMETMass");

  histoManager.AddVariable("VtxCleanMETPresel",400,0,400,"#slash{E}_{T} [GeV]","VtxCleanMETPresel");
  histoManager.AddVariable("VtxCleanMETMET",400,0,400,"#slash{E}_{T} [GeV]","VtxCleanMETMET");
  histoManager.AddVariable("VtxCleanMETBal",400,0,400,"#slash{E}_{T} [GeV]","VtxCleanMETBal");
  histoManager.AddVariable("VtxCleanMETJVeto",400,0,400,"#slash{E}_{T} [GeV]","VtxCleanMETJVeto");
  histoManager.AddVariable("VtxCleanMETDPhiJ",400,0,400,"#slash{E}_{T} [GeV]","VtxCleanMETDPhiJ");
  histoManager.AddVariable("VtxCleanMETDPhiZ",400,0,400,"#slash{E}_{T} [GeV]","VtxCleanMETDPhiZ");
  histoManager.AddVariable("VtxCleanMETLVeto",400,0,400,"#slash{E}_{T} [GeV]","VtxCleanMETLVeto");
  histoManager.AddVariable("VtxCleanMETBTag",400,0,400,"#slash{E}_{T} [GeV]","VtxCleanMETBTag");
  histoManager.AddVariable("VtxCleanMETMass",400,0,400,"#slash{E}_{T} [GeV]","VtxCleanMETMass");

  histoManager.AddVariable("SumCleanMETPresel",400,0,400,"#slash{E}_{T} [GeV]","SumCleanMETPresel");
  histoManager.AddVariable("SumCleanMETMET",400,0,400,"#slash{E}_{T} [GeV]","SumCleanMETMET");
  histoManager.AddVariable("SumCleanMETBal",400,0,400,"#slash{E}_{T} [GeV]","SumCleanMETBal");
  histoManager.AddVariable("SumCleanMETJVeto",400,0,400,"#slash{E}_{T} [GeV]","SumCleanMETJVeto");
  histoManager.AddVariable("SumCleanMETDPhiJ",400,0,400,"#slash{E}_{T} [GeV]","SumCleanMETDPhiJ");
  histoManager.AddVariable("SumCleanMETDPhiZ",400,0,400,"#slash{E}_{T} [GeV]","SumCleanMETDPhiZ");
  histoManager.AddVariable("SumCleanMETLVeto",400,0,400,"#slash{E}_{T} [GeV]","SumCleanMETLVeto");
  histoManager.AddVariable("SumCleanMETBTag",400,0,400,"#slash{E}_{T} [GeV]","SumCleanMETBTag");
  histoManager.AddVariable("SumCleanMETMass",400,0,400,"#slash{E}_{T} [GeV]","SumCleanMETMass");

  histoManager.AddVariable("SqrtCleanMETPresel",400,0,400,"#slash{E}_{T} [GeV]","SqrtCleanMETPresel");
  histoManager.AddVariable("SqrtCleanMETMET",400,0,400,"#slash{E}_{T} [GeV]","SqrtCleanMETMET");
  histoManager.AddVariable("SqrtCleanMETBal",400,0,400,"#slash{E}_{T} [GeV]","SqrtCleanMETBal");
  histoManager.AddVariable("SqrtCleanMETJVeto",400,0,400,"#slash{E}_{T} [GeV]","SqrtCleanMETJVeto");
  histoManager.AddVariable("SqrtCleanMETDPhiJ",400,0,400,"#slash{E}_{T} [GeV]","SqrtCleanMETDPhiJ");
  histoManager.AddVariable("SqrtCleanMETDPhiZ",400,0,400,"#slash{E}_{T} [GeV]","SqrtCleanMETDPhiZ");
  histoManager.AddVariable("SqrtCleanMETLVeto",400,0,400,"#slash{E}_{T} [GeV]","SqrtCleanMETLVeto");
  histoManager.AddVariable("SqrtCleanMETBTag",400,0,400,"#slash{E}_{T} [GeV]","SqrtCleanMETBTag");
  histoManager.AddVariable("SqrtCleanMETMass",400,0,400,"#slash{E}_{T} [GeV]","SqrtCleanMETMass");

  histoManager.AddVariable("CorMinMETPresel",400,0,400,"cor clean #slash{E}_{T} [GeV]","CorMinMETPresel");
  histoManager.AddVariable("CorMinMETMET",400,0,400,"cor clean #slash{E}_{T} [GeV]","CorMinMETCorMET");
  histoManager.AddVariable("CorMinMETBal",400,0,400,"cor clean #slash{E}_{T} [GeV]","CorMinMETBal");
  histoManager.AddVariable("CorMinMETJVeto",400,0,400,"cor clean #slash{E}_{T} [GeV]","CorMinMETJVeto");
  histoManager.AddVariable("CorMinMETDPhiJ",400,0,400,"cor clean #slash{E}_{T} [GeV]","CorMinMETDPhiJ");
  histoManager.AddVariable("CorMinMETDPhiZ",400,0,400,"cor clean #slash{E}_{T} [GeV]","CorMinMETDPhiZ");
  histoManager.AddVariable("CorMinMETLVeto",400,0,400,"cor clean #slash{E}_{T} [GeV]","CorMinMETLVeto");
  histoManager.AddVariable("CorMinMETBTag",400,0,400,"cor clean #slash{E}_{T} [GeV]","CorMinMETBTag");
  histoManager.AddVariable("CorMinMETMass",400,0,400,"cor clean #slash{E}_{T} [GeV]","CorMinMETMass");

  histoManager.AddVariable("BalPresel",200,0,10,"q_{T}/#slash{E}_{T}","BalPresel");
  histoManager.AddVariable("BalMET",200,0,10,"q_{T}/#slash{E}_{T}","BalMET");
  histoManager.AddVariable("BalBal",200,0,10,"q_{T}/#slash{E}_{T}","BalBal");
  histoManager.AddVariable("BalJVeto",200,0,10,"q_{T}/#slash{E}_{T}","BalJVeto");
  histoManager.AddVariable("BalDPhiJ",200,0,10,"q_{T}/#slash{E}_{T}","BalDPhiJ");
  histoManager.AddVariable("BalDPhiZ",200,0,10,"q_{T}/#slash{E}_{T}","BalDPhiZ");
  histoManager.AddVariable("BalLVeto",200,0,10,"q_{T}/#slash{E}_{T}","BalLVeto");
  histoManager.AddVariable("BalBTag",200,0,10,"q_{T}/#slash{E}_{T}","BalBTag");
  histoManager.AddVariable("BalMass",200,0,10,"q_{T}/#slash{E}_{T}","BalMass");

  //Z + MET plots
  histoManager.AddVariable("ZMETMassPresel",500,0,1000,"M_{T}(ZZ) [GeV]","ZMETMassPresel");
  histoManager.AddVariable("ZMETMassMET",500,0,1000,"M_{T}(ZZ) [GeV]","ZMETMassMET");
  histoManager.AddVariable("ZMETMassBal",500,0,1000,"M_{T}(ZZ) [GeV]","ZMETMassBal");
  histoManager.AddVariable("ZMETMassJVeto",500,0,1000,"M_{T}(ZZ) [GeV]","ZMETMassJVeto");
  histoManager.AddVariable("ZMETMassDPhiJ",500,0,1000,"M_{T}(ZZ) [GeV]","ZMETMassDPhiJ");
  histoManager.AddVariable("ZMETMassDPhiZ",500,0,1000,"M_{T}(ZZ) [GeV]","ZMETMassDPhiZ");
  histoManager.AddVariable("ZMETMassLVeto",500,0,1000,"M_{T}(ZZ) [GeV]","ZMETMassLVeto");
  histoManager.AddVariable("ZMETMassBTag",500,0,1000,"M_{T}(ZZ) [GeV]","ZMETMassBTag");
  histoManager.AddVariable("ZMETMassMass",500,0,1000,"M_{T}(ZZ) [GeV]","ZMETMassMass");

  //Angles Plots
  histoManager.AddVariable("dPhiJPresel",128,-180,180,"#Delta#phi(jet, #slash{E}_{T}) [#circ]","dPhiJPresel");
  histoManager.AddVariable("dPhiJMET",128,-180,180,"#Delta#phi(jet, #slash{E}_{T}) [#circ]","dPhiJMET");
  histoManager.AddVariable("dPhiJBal",128,-180,180,"#Delta#phi(jet, #slash{E}_{T}) [#circ]","dPhiJBal");
  histoManager.AddVariable("dPhiJJVeto",128,-180,180,"#Delta#phi(jet, #slash{E}_{T}) [#circ]","dPhiJJVeto");
  histoManager.AddVariable("dPhiJDPhiJ",128,-180,180,"#Delta#phi(jet, #slash{E}_{T}) [#circ]","dPhiJDPhiJ");
  histoManager.AddVariable("dPhiJDPhiZ",128,-180,180,"#Delta#phi(jet, #slash{E}_{T}) [#circ]","dPhiJDPhiZ");
  histoManager.AddVariable("dPhiJLVeto",128,-180,180,"#Delta#phi(jet, #slash{E}_{T}) [#circ]","dPhiJLVeto");
  histoManager.AddVariable("dPhiJBTag",128,-180,180,"#Delta#phi(Z, #slash{E}_{T}) [#circ]","dPhiZBTag");
  histoManager.AddVariable("dPhiJMass",128,-180,180,"#Delta#phi(jet, #slash{E}_{T}) [#circ]","dPhiJMass");

  histoManager.AddVariable("dPhiZPresel",128,-180,180,"#Delta#phi(Z, #slash{E}_{T}) [#circ]","dPhiZPresel");
  histoManager.AddVariable("dPhiZMET",128,-180,180,"#Delta#phi(Z, #slash{E}_{T}) [#circ]","dPhiZMET");
  histoManager.AddVariable("dPhiZBal",128,-180,180,"#Delta#phi(Z, #slash{E}_{T}) [#circ]","dPhiZBal");
  histoManager.AddVariable("dPhiZJVeto",128,-180,180,"#Delta#phi(Z, #slash{E}_{T}) [#circ]","dPhiZJVeto");
  histoManager.AddVariable("dPhiZDPhiJ",128,-180,180,"#Delta#phi(Z, #slash{E}_{T}) [#circ]","dPhiZDPhiJ");
  histoManager.AddVariable("dPhiZDPhiZ",128,-180,180,"#Delta#phi(Z, #slash{E}_{T}) [#circ]","dPhiZDPhiZ");
  histoManager.AddVariable("dPhiZLVeto",128,-180,180,"#Delta#phi(Z, #slash{E}_{T}) [#circ]","dPhiZLVeto");
  histoManager.AddVariable("dPhiZBTag",128,-180,180,"#Delta#phi(Z, #slash{E}_{T}) [#circ]","dPhiZBTag");
  histoManager.AddVariable("dPhiZMass",128,-180,180,"#Delta#phi(Z, #slash{E}_{T}) [#circ]","dPhiZMass");


  //Jet plots
  histoManager.AddVariable("FJetPresel",400,0,400,"p_{T} first jet [GeV]","FJetPresel");
  histoManager.AddVariable("FJetMET",400,0,400,"p_{T} first jet [GeV]","FJetMET");
  histoManager.AddVariable("FJetBal",400,0,400,"p_{T} first jet [GeV]","FJetBal");
  histoManager.AddVariable("FJetJVeto",400,0,400,"p_{T} first jet [GeV]","FJetJVeto");
  histoManager.AddVariable("FJetDPhiJ",400,0,400,"p_{T} first jet [GeV]","FJetDPhiJ");
  histoManager.AddVariable("FJetDPhiZ",400,0,400,"p_{T} first jet [GeV]","FJetDPhiZ");
  histoManager.AddVariable("FJetLVeto",400,0,400,"p_{T} first jet [GeV]","FJetLVeto");
  histoManager.AddVariable("FJetBTag",400,0,400,"p_{T} first jet [GeV]","FJetBTag");
  histoManager.AddVariable("FJetMass",400,0,400,"p_{T} first jet [GeV]","FJetMass");
  
  histoManager.AddVariable("SJetPresel",400,0,400,"p_{T} second jet [GeV]","SJetPresel");
  histoManager.AddVariable("SJetMET",400,0,400,"p_{T} second jet [GeV]","SJetMET");
  histoManager.AddVariable("SJetBal",400,0,400,"p_{T} second jet [GeV]","SJetBal");
  histoManager.AddVariable("SJetJVeto",400,0,400,"p_{T} second jet [GeV]","SJetJVeto");
  histoManager.AddVariable("SJetDPhiJ",400,0,400,"p_{T} second jet [GeV]","SJetDPhiJ");
  histoManager.AddVariable("SJetDPhiZ",400,0,400,"p_{T} second jet [GeV]","SJetDPhiZ");
  histoManager.AddVariable("SJetLVeto",400,0,400,"p_{T} second jet [GeV]","SJetLVeto");
  histoManager.AddVariable("SJetBTag",400,0,400,"p_{T} second jet [GeV]","SJetBTag");
  histoManager.AddVariable("SJetMass",400,0,400,"p_{T} second jet [GeV]","SJetMass");

  //b-Tag plots
  histoManager.AddVariable("BTagTKFJetPresel",200,-10,40,"track counting b-tag","BTagTKFJetPresel");
  histoManager.AddVariable("BTagTKFJetMET",200,-10,40,"track counting b-tag","BTagTKFJetMET");
  histoManager.AddVariable("BTagTKFJetBal",200,-10,40,"track counting b-tag","BTagTKFJetBal");
  histoManager.AddVariable("BTagTKFJetJVeto",200,-10,40,"track counting b-tag","BTagTKFJetJVeto");
  histoManager.AddVariable("BTagTKFJetDPhiJ",200,-10,40,"track counting b-tag","BTagTKFJetDPhiJ");
  histoManager.AddVariable("BTagTKFJetDPhiZ",200,-10,40,"track counting b-tag","BTagTKFJetDPhiZ");
  histoManager.AddVariable("BTagTKFJetLVeto",200,-10,40,"track counting b-tag","BTagTKFJetLVeto");
  histoManager.AddVariable("BTagTKFJetBTag",200,-10,40,"track counting b-tag","BTagTKFJetBTag");
  histoManager.AddVariable("BTagTKFJetMass",200,-10,40,"track counting b-tag","BTagTKFJetMass");
  
  histoManager.AddVariable("BTagSMFJetPresel",55,-0.1,1,"soft muon b-tag","BTagSMFJetPresel");
  histoManager.AddVariable("BTagSMFJetMET",55,-0.1,1,"soft muon b-tag","BTagSMFJetMET");
  histoManager.AddVariable("BTagSMFJetBal",55,-0.1,1,"soft muon b-tag","BTagSMFJetBal");
  histoManager.AddVariable("BTagSMFJetJVeto",55,-0.1,1,"soft muon b-tag","BTagSMFJetJVeto");
  histoManager.AddVariable("BTagSMFJetDPhiJ",55,-0.1,1,"soft muon b-tag","BTagSMFJetDPhiJ");
  histoManager.AddVariable("BTagSMFJetDPhiZ",55,-0.1,1,"soft muon b-tag","BTagSMFJetDPhiZ");
  histoManager.AddVariable("BTagSMFJetLVeto",55,-0.1,1,"soft muon b-tag","BTagSMFJetLVeto");
  histoManager.AddVariable("BTagSMFJetBTag",55,-0.1,1,"soft muon b-tag","BTagSMFJetBTag");
  histoManager.AddVariable("BTagSMFJetMass",55,-0.1,1,"soft muon b-tag","BTagSMFJetMass");

  histoManager.AddVariable("BTagTKSJetPresel",200,-10,40,"track counting b-tag","BTagTKSJetPresel");
  histoManager.AddVariable("BTagTKSJetMET",200,-10,40,"track counting b-tag","BTagTKSJetMET");
  histoManager.AddVariable("BTagTKSJetBal",200,-10,40,"track counting b-tag","BTagTKSJetBal");
  histoManager.AddVariable("BTagTKSJetJVeto",200,-10,40,"track counting b-tag","BTagTKSJetJVeto");
  histoManager.AddVariable("BTagTKSJetDPhiJ",200,-10,40,"track counting b-tag","BTagTKSJetDPhiJ");
  histoManager.AddVariable("BTagTKSJetDPhiZ",200,-10,40,"track counting b-tag","BTagTKSJetDPhiZ");
  histoManager.AddVariable("BTagTKSJetLVeto",200,-10,40,"track counting b-tag","BTagTKSJetLVeto");
  histoManager.AddVariable("BTagTKSJetBTag",200,-10,40,"track counting b-tag","BTagTKSJetBTag");
  histoManager.AddVariable("BTagTKSJetMass",200,-10,40,"track counting b-tag","BTagTKSJetMass");
  
  histoManager.AddVariable("BTagSMSJetPresel",55,-0.1,1,"soft muon b-tag","BTagSMSJetPresel");
  histoManager.AddVariable("BTagSMSJetMET",55,-0.1,1,"soft muon b-tag","BTagSMSJetMET");
  histoManager.AddVariable("BTagSMSJetBal",55,-0.1,1,"soft muon b-tag","BTagSMSJetBal");
  histoManager.AddVariable("BTagSMSJetJVeto",55,-0.1,1,"soft muon b-tag","BTagSMSJetJVeto");
  histoManager.AddVariable("BTagSMSJetDPhiJ",55,-0.1,1,"soft muon b-tag","BTagSMSJetDPhiJ");
  histoManager.AddVariable("BTagSMSJetDPhiZ",55,-0.1,1,"soft muon b-tag","BTagSMSJetDPhiZ");
  histoManager.AddVariable("BTagSMSJetLVeto",55,-0.1,1,"soft muon b-tag","BTagSMSJetLVeto");
  histoManager.AddVariable("BTagSMSJetBTag",55,-0.1,1,"soft muon b-tag","BTagSMSJetBTag");
  histoManager.AddVariable("BTagSMSJetMass",55,-0.1,1,"soft muon b-tag","BTagSMSJetMass");

  //Lepton Mult
  histoManager.AddVariable("LepMultPresel",10,0,10,"lepmult","LepMultPresel");
  histoManager.AddVariable("LepMultMET",10,0,10,"lepmult","LepMultMET");
  histoManager.AddVariable("LepMultBal",10,0,10,"lepmult","LepMultBal");
  histoManager.AddVariable("LepMultJVeto",10,0,10,"lepmult","LepMultJVeto");
  histoManager.AddVariable("LepMultDPhiJ",10,0,10,"lepmult","LepMultDPhiJ");
  histoManager.AddVariable("LepMultDPhiZ",10,0,10,"lepmult","LepMultDPhiZ");
  histoManager.AddVariable("LepMultLVeto",10,0,10,"lepmult","LepMultLVeto");
  histoManager.AddVariable("LepMultBTag",10,0,10,"lepmult","LepMultBTag");
  histoManager.AddVariable("LepMultMass",10,0,10,"lepmult","LepMultMass");

  histoManager.AddVariable("LMTPresel",200,0,200,"MT [GeV]","LMTPresel");
  histoManager.AddVariable("LMTMET",200,0,200,"MT [GeV]","LMTMET");
  histoManager.AddVariable("LMTBal",200,0,200,"MT [GeV]","LMTBal");
  histoManager.AddVariable("LMTJVeto",200,0,200,"MT [GeV]","LMTJVeto");
  histoManager.AddVariable("LMTDPhiJ",200,0,200,"MT [GeV]","LMTDPhiJ");
  histoManager.AddVariable("LMTDPhiZ",200,0,200,"MT [GeV]","LMTDPhiZ");
  histoManager.AddVariable("LMTLVeto",200,0,200,"MT [GeV]","LMTLVeto");
  histoManager.AddVariable("LMTBTag",200,0,200,"MT [GeV]","LMTBTag");
  histoManager.AddVariable("LMTMass",200,0,200,"MT [GeV]","LMTMass");

  histoManager.AddVariable("LepPdgPresel",40,-20,20,"lepPdg","LepPdgPresel");
  histoManager.AddVariable("LepPdgMET",40,-20,20,"lepPdg","LepPdgMET");
  histoManager.AddVariable("LepPdgBal",40,-20,20,"lepPdg","LepPdgBal");
  histoManager.AddVariable("LepPdgJVeto",40,-20,20,"lepPdg","LepPdgJVeto");
  histoManager.AddVariable("LepPdgDPhiJ",40,-20,20,"lepPdg","LepPdgDPhiJ");
  histoManager.AddVariable("LepPdgDPhiZ",40,-20,20,"lepPdg","LepPdgDPhiZ");
  histoManager.AddVariable("LepPdgLVeto",40,-20,20,"lepPdg","LepPdgLVeto");
  histoManager.AddVariable("LepPdgBTag",40,-20,20,"lepPdg","LepPdgBTag");
  histoManager.AddVariable("LepPdgMass",40,-20,20,"lepPdg","LepPdgMass");

  //Jet Mult
  histoManager.AddVariable("JetMultPresel",10,0,10,"jetmult","JetMultPresel");
  histoManager.AddVariable("JetMultMET",10,0,10,"jetmult","JetMultMET");
  histoManager.AddVariable("JetMultBal",10,0,10,"jetmult","JetMultBal");
  histoManager.AddVariable("JetMultJVeto",10,0,10,"jetmult","JetMultJVeto");
  histoManager.AddVariable("JetMultDPhiJ",10,0,10,"jetmult","JetMultDPhiJ");
  histoManager.AddVariable("JetMultDPhiZ",10,0,10,"jetmult","JetMultDPhiZ");
  histoManager.AddVariable("JetMultLVeto",10,0,10,"jetmult","JetMultLVeto");
  histoManager.AddVariable("JetMultBTag",10,0,10,"jetmult","JetMultBTag");
  histoManager.AddVariable("JetMultMass",10,0,10,"jetmult","JetMultMass");

  //Photon plots
  histoManager.AddVariable("FPhotPresel",400,0,400,"p_{T} first phot [GeV]","FPhotPresel");
  histoManager.AddVariable("FPhotMET",400,0,400,"p_{T} first phot [GeV]","FPhotMET");
  histoManager.AddVariable("FPhotBal",400,0,400,"p_{T} first phot [GeV]","FPhotBal");
  histoManager.AddVariable("FPhotJVeto",400,0,400,"p_{T} first phot [GeV]","FPhotJVeto");
  histoManager.AddVariable("FPhotDPhiJ",400,0,400,"p_{T} first phot [GeV]","FPhotDPhiJ");
  histoManager.AddVariable("FPhotDPhiZ",400,0,400,"p_{T} first phot [GeV]","FPhotDPhiZ");
  histoManager.AddVariable("FPhotLVeto",400,0,400,"p_{T} first phot [GeV]","FPhotLVeto");
  histoManager.AddVariable("FPhotBTag",400,0,400,"p_{T} first phot [GeV]","FPhotBTag");
  histoManager.AddVariable("FPhotMass",400,0,400,"p_{T} first phot [GeV]","FPhotMass");



  //Misc
  histoManager.AddVariable("dRJL",200,-2,2,"D pt","DPT");
  histoManager.AddVariable("nL",10,0,10,"Pos Lepton","DPT");

  histoManager.AddVariable("ZCateg", 4, 0, 4, "Categ Z", "ZCateg");
  histoManager.AddVariable("LepAcc", 2, 0, 2, "Categ lep", "ZCateg");

  histoManager.AddVariable("AddLepCategE", 2, 0, 2, "categ Lepton","categ");
  histoManager.AddVariable("AddLepCategM", 2, 0, 2, "categ Lepton","categ");

  histoManager.AddVariable("J1DR1", 200, -3, 3, "J1DR1", "J1DR1");
  histoManager.AddVariable("J1DR2", 200, -3, 3, "J1DR2", "J1DR2");
  histoManager.AddVariable("J2DR1", 200, -3, 3, "J2DR1", "J2DR1");
  histoManager.AddVariable("J2DR2", 200, -3, 3, "J2DR2", "J2DR2");

  //End Declaration, now preparation
  histoManager.PrepareHistos(name);
  histoManager.Prepare2DHistos(name);
 

  vector<vector<vector<float> > > Wghts;
  // if(ShapeWeight)
    Wghts = GetFitWeight();
  
  string Epart="EE";
  if(EcalP==0)
    Epart = "EB";
  if(EcalP==2)
    Epart="Combined";
  if(EcalP==3)
    Epart="EBEE for Z studies";

  cout<<" Ecal Part "<<Epart<<endl;

  cout<<" check cuts EB : id "<<_IDcuts[0][0]<<"  "
      <<_IDcuts[0][1]<<"  "<<_IDcuts[0][2]<<"  "<<_IDcuts[0][3]<<endl;
  cout<<" check cuts EB : iso "<<_Isocuts[0][0]<<"  "
      <<_Isocuts[0][1]<<"  "<<_Isocuts[0][2]<<endl;

 cout<<" check cuts EE : id "<<_IDcuts[1][0]<<"  "
     <<_IDcuts[1][1]<<"  "<<_IDcuts[1][2]<<"  "<<_IDcuts[1][3]<<endl;
  cout<<" check cuts EE : iso "<<_Isocuts[1][0]<<"  "
      <<_Isocuts[1][1]<<"  "<<_Isocuts[1][2]<<endl;

  int S=0,B=0;
 
  vector<int> NumberEntries(nt+1,0);


  //sélection
  for(int i=0;i<nt+1;i++) {
   
    //Déclaration des variables

    int run;
    int event;
    string* categ=0;
    string* c_str=0;
    TBits* c_bits=0;
    int nZCand;
    bool isTight;
    float mll;
    float qTll;
    float phill;
    float yll;
    float hll;
    float mTZZ;
    int nLoose;
    int nVertex;
    vector<double>* LPt=NULL;
    vector<double>* LEta=NULL;
    vector<double>* LPhi=NULL;
    vector<int>* LID=NULL;
    vector<int>* LPdg=NULL;
    vector<double>* LMT=NULL;
    vector<double>* LTrkIso=NULL;
    vector<double>* LEcalIso=NULL;
    vector<double>* LHcalIso=NULL;
    vector<double>* LCombIso=NULL;
    float MET;
    float METPhi;
    float sumEt;
    float mEtSig;
    float projMET;
    float corProjMET;
    float sigMET;
    float dPhiMin;
    float balZMET;
    float cpuMETVtx;
    float cpuMETSum;
    float cpuMETSqrt;
    vector<double>* jEt=NULL;
    vector<double>* jEta=NULL;
    vector<double>* jPhi=NULL;
    vector<double>* jDR1=NULL;
    vector<double>* jDR2=NULL;
    vector<double>* jDPhiMET=NULL;
    vector<double>* jBtagTkC=NULL;
    vector<double>* jBtagSoftM=NULL;
    float JZB15;
    float JZB20;
    float JZB25;
    int jMult30;
    int jMult25;
    int jMult20;
    int jMult15;
    int jMult10;
    float dPhiJMET;

    vector<double>* pEt=NULL;
    vector<double>* pEta=NULL;
    vector<double>* pPhi=NULL;
    vector<double>* pDR1=NULL;
    vector<double>* pDR2=NULL;
    vector<double>* pDPhiMET=NULL;

    tChains[i]->SetBranchAddress("run", &run);
    tChains[i]->SetBranchAddress("event", &event);
    tChains[i]->SetBranchAddress("categ", &categ);
    tChains[i]->SetBranchAddress("c_str", &c_str);
    tChains[i]->SetBranchAddress("c_bits", &c_bits);
    tChains[i]->SetBranchAddress("nZCand", &nZCand);
    tChains[i]->SetBranchAddress("isTight", &isTight);
    tChains[i]->SetBranchAddress("mll", &mll);
    tChains[i]->SetBranchAddress("qTll", &qTll);
    tChains[i]->SetBranchAddress("phill", &phill);
    tChains[i]->SetBranchAddress("yll", &yll);
    tChains[i]->SetBranchAddress("hll", &hll);
    tChains[i]->SetBranchAddress("mTZZ", &mTZZ);
    tChains[i]->SetBranchAddress("nLoose", &nLoose);
    tChains[i]->SetBranchAddress("nVertex", &nVertex);
    tChains[i]->SetBranchAddress("LPt", &LPt);
    tChains[i]->SetBranchAddress("LEta", &LEta);
    tChains[i]->SetBranchAddress("LPhi", &LPhi);
    tChains[i]->SetBranchAddress("LID", &LID);
    tChains[i]->SetBranchAddress("LPdg", &LPdg);
    tChains[i]->SetBranchAddress("LMT", &LMT);
    tChains[i]->SetBranchAddress("LTrkIso", &LTrkIso);
    tChains[i]->SetBranchAddress("LEcalIso", &LEcalIso);
    tChains[i]->SetBranchAddress("LHcalIso", &LHcalIso);
    tChains[i]->SetBranchAddress("LCombIso", &LCombIso);
    tChains[i]->SetBranchAddress("MET", &MET);
    tChains[i]->SetBranchAddress("METPhi", &METPhi);
    tChains[i]->SetBranchAddress("sumEt", &sumEt);
    tChains[i]->SetBranchAddress("mEtSig", &mEtSig);
    tChains[i]->SetBranchAddress("projMET", &projMET);
    tChains[i]->SetBranchAddress("corProjMET", &corProjMET);
    tChains[i]->SetBranchAddress("sigMET", &sigMET);
    tChains[i]->SetBranchAddress("dPhiMin", &dPhiMin);
    tChains[i]->SetBranchAddress("balZMET", &balZMET);
    tChains[i]->SetBranchAddress("cpuMETVtx", &cpuMETVtx);
    tChains[i]->SetBranchAddress("cpuMETSum", &cpuMETSum);
    tChains[i]->SetBranchAddress("cpuMETSqrt", &cpuMETSqrt);
    tChains[i]->SetBranchAddress("jEt", &jEt);
    tChains[i]->SetBranchAddress("jEta", &jEta);
    tChains[i]->SetBranchAddress("jPhi", &jPhi);
    tChains[i]->SetBranchAddress("jDR1", &jDR1);
    tChains[i]->SetBranchAddress("jDR2", &jDR2);
    tChains[i]->SetBranchAddress("jDPhiMET", &jDPhiMET);
    tChains[i]->SetBranchAddress("jBtagTkC", &jBtagTkC);
    tChains[i]->SetBranchAddress("jBtagSoftM", &jBtagSoftM);
    tChains[i]->SetBranchAddress("JZB15", &JZB15);
    tChains[i]->SetBranchAddress("JZB20", &JZB20);
    tChains[i]->SetBranchAddress("JZB25", &JZB25);
    tChains[i]->SetBranchAddress("jMult30", &jMult30);
    tChains[i]->SetBranchAddress("jMult25", &jMult25);
    tChains[i]->SetBranchAddress("jMult20", &jMult20);
    tChains[i]->SetBranchAddress("jMult15", &jMult15);
    tChains[i]->SetBranchAddress("jMult10", &jMult10);
    tChains[i]->SetBranchAddress("dPhiJMET", &dPhiJMET);
    
    tChains[i]->SetBranchAddress("pEt", &pEt);
    tChains[i]->SetBranchAddress("pEta", &pEta);
    tChains[i]->SetBranchAddress("pPhi", &pPhi);
    tChains[i]->SetBranchAddress("pDR1", &pDR1);
    tChains[i]->SetBranchAddress("pDR2", &pDR2);
    tChains[i]->SetBranchAddress("pDPhiMET", &pDPhiMET);


    //Addtionnal trees ===================================================
    //Monte carlo Truth

    float mVV;
    float mV1;
    float ptV2;
    float ptV1;
    vector<double>* gLPt=NULL;
    vector<double>* gLEta=NULL;
    vector<double>* gLPdg=NULL;

    TTree* MCTree = GetAdditionnalTree("MCTruth", i);
    MCTree->SetBranchAddress("mVV", &mVV);
    MCTree->SetBranchAddress("mV1", &mV1);
    MCTree->SetBranchAddress("ptV1", &ptV1);
    MCTree->SetBranchAddress("ptV2", &ptV2);
    MCTree->SetBranchAddress("ptL", &gLPt);
    MCTree->SetBranchAddress("etaL", &gLEta);
    MCTree->SetBranchAddress("pdgL", &gLPdg);

    //Lepton Tree ====

    int expInnHits;
    bool isConv;
    bool isFPH;
    float LPT;
    TTree* LepTree = GetAdditionnalTree("Leptons", i);
    LepTree->SetBranchAddress("expInnHits", &expInnHits );
    LepTree->SetBranchAddress("isConv", &isConv );
    LepTree->SetBranchAddress("isFPHitTrk", &isFPH );
    LepTree->SetBranchAddress("pt", &LPT);
    int nEL =0;
    bool expI[2]={true,true};

    //====================================================================




    //bool Sel[2]={true,true};

    float Weight=1;
    float Weight2=1;

    cout<<" beginning tree  " <<name[i]<<"   "<<tChains[i]->GetEntries()<<endl;
    //Boucle et sélection, remplissage histos
    //int EP=0;
    int ent = tChains[i]->GetEntries();
    for(int ie=0;ie<ent;ie++) {
      // cout<<i<<"  "<<nt<<"   "<<ie<<"  "<<ent<<"   "<<tChains[i]<<endl;

      // if(ie%1000==0)
      // 	fprintf(stdout, "Done : %d / %d\r", ie, ent );

      if(i!=nt)
      	Weight = GetWeight(i, ie);
      else
      	Weight = 1;

      tChains[i]->GetEntry(ie);  
      MCTree->GetEntry(ie);  
      
      if(i!=nt)
      	Weight *= GetMCWeightFromDataShape(i, nVertex, 
      					   "NVertexControl", 
      					   "NVertexControl",1, i);
   
      Weight2 = Weight;


      //Lepton counter //FIXME not ready
      if( abs( (*LPdg)[0] ) ==11) {
	for(int il=0;il<2;il++) {
	  LepTree->GetEntry(nEL);
	  // cout<<" n--> "<<ie<<"  ==>  "<<nEL<<"     "<<expInnHits<<"    "<<isConv<<"     "
	  //     <<(*LPt)[0]<<"   "<<(*LPt)[1]<<"   "<<LPT<<"   "<<LPT<<endl;

	  // expInnHits

	  if( expInnHits != 0) expI[il]=false; 
	  nEL++;
	}
      }

      //FIXME
      //if(qTll< 30) continue;
     

      //Acceptance
      inAcc=true;
   
      if(name[i].substr(0,2)=="di" || name[i].substr(0,2)=="ZZ" ||
      	 name[i].substr(0,2)=="ZV" || name[i].substr(0,2)=="WZ" ||
      	 name[i].substr(0,2)=="WW") {
	//	cout<<Weight<<"   ";
	string proc = FindProcess( ie, i );

	if( proc.substr(0,2)=="ZZ") {
	  Weight *= SearchWeightZZ(ptV1); 
	  Weight2 *= SearchWeightZZ(ptV1 , true );
	  //Acceptance
	  inAcc = isInAcc( ptV1, ptV2, mV1, gLPt, gLEta, gLPdg);
	  //FIXME
	  //	  if(!inAcc) continue;
	}
	if( proc.substr(0,2)=="WZ") {
	  Weight *= SearchWeightWZ(ptV2);
	  Weight2 *= SearchWeightWZ(ptV2 , true );
	  inAcc = isInAcc( ptV1, ptV2, mV1, gLPt, gLEta, gLPdg, 2);
	}
	if( proc.substr(0,2)=="WW") {
	  Weight *= SearchWeightWW(ptV1); //based on W+
	  Weight2 *= SearchWeightZZ(ptV1 , true );
	}
	//	cout<<ptV1<<"     "<<GetWeight(i, ie)<<"      "<<Weight<<"    "<<Weight2<<endl;
      }
      //      cout<<inAcc<<endl;
     
      if(i==nt)//name[i].substr(0,4)=="data"
      	{
      	  //cout<<run<<"   "<<event<<endl;
      	  bool doubleCount=false;
      	  for(size_t ik=0;ik<Events.size();ik++)
      	    {
      	      if( ((Events[ik]).first.first) == run && ((Events[ik]).first.second) == event)
      	      	{doubleCount=true; cout<<" ============>>> found  double event  "<<run<<"   "<<event<<endl; break;}
      	    }
      	  if(doubleCount || (EventFilter && run>EventNum ) )
      	    {  continue; }
      	}

      //All events


      if(nVertex!=NVert && NVert!=0 ) continue; //NVertex
      if(nVertex==0) continue;

      //=====================================================================================
 
      if(UsePdgId && ( abs( (*LPdg)[0]) !=pdgID || abs( (*LPdg)[1]) !=pdgID) )
       	continue;

      MakeCut<bool>(true,true,"=", i, "Beginning", Weight);

      // if( (*LPt)[0] *  (*LPt)[1] >0 ) continue;  //MM

      //      cout<<(*LID)[0]%2<<"    "<<(*LID)[1]%4<<"    "<<(*LID)[1]%2<<"    "<<(*LID)[0]%4<<endl;

      if(!MakeCut<float>(mll,60.,"[]", i, "window mass", Weight, 120. )) continue;

      if( abs( (*LPdg)[0]) ==11) {
	if( !( ( ((*LID)[0])>=7 && ((*LID)[1])>=7) ||
	       ( ((*LID)[1])>=7 && ((*LID)[0])>=7) ) ) continue;
      }
      else {
	if( !( ( ((*LID)[0])>=7 && ((*LID)[1])>=3) ||
	       ( ((*LID)[1])>=7 && ((*LID)[0])>=3) ) ) continue;
      }

      // if( UsePdgId && pdgID==13 ) {
      // 	  if( ((*LCombIso)[0])>0.20 || ((*LCombIso)[1])>0.20 ) 
      // 	    continue;
      // 	}

      //cout<<(*LID)[0]<<"   "<<(*LID)[1]<<"    "<<mll<<"    "<<isTight<<endl;
      MakeCut<bool>(true,true,"=", i, "preselection", Weight);
      if( nZCand!=1 ) continue;
      
      //  if( !isTight ) continue;
      
      //if(!(*c_bits)[0] /* || !(*c_bits)[2]*/) continue;

      //CorMET Computation
      float CorMET;
      if(i!=nt) {
      	CorMET = METCorrectionMC(MET, nVertex );
      }
      else {
      	CorMET = METCorrectionData(MET, nVertex );
      }
      //and dPhi ZMET
      float dphimetz = dPhi( METPhi , phill+ TMath::Pi() );


      //And with minimum
      float CorMMET = METCorrectionMin( min(cpuMETVtx,MET), nVertex, (i!=nt), fabs( (int)((*LPdg)[0] ) ) );
      
      //Preselection
      //MakeCut<bool>(true,true,"=",i,"preselection", Weight);
      

      if(!MakeCut<float>(qTll,30.,">=", i, "qT", Weight )) continue;

      //top control
      //if( mll < 30 || (mll>80 && mll<100) ) continue;
      //if(!MakeCut<float>(mll,80.,"][", i, "preselection", Weight, 100. )) continue;

      histoManager.fill("ZMassPresel",i,mll,Weight);
      histoManager.fill("ZPtPresel",i,qTll,Weight);
      histoManager.fill("ZYPresel",i,yll,Weight);
      histoManager.fill("LepMultPresel",i,nLoose,Weight);
      //histoManager.fill("jetMultPresel",i,nLoose,Weight);
      histoManager.fill("METPresel",i,MET,Weight);
      histoManager.fill("CorMETPresel",i,CorMET,Weight);
      histoManager.fill("VtxCleanMETPresel",i,min(cpuMETVtx,MET),Weight);
      histoManager.fill("SumCleanMETPresel",i,min(cpuMETSum,MET),Weight);
      histoManager.fill("SqrtCleanMETPresel",i,min(cpuMETSqrt,MET),Weight);
      histoManager.fill("CorMinMETPresel",i,CorMMET,Weight);
      histoManager.fill("BalPresel",i,MET/qTll,Weight);
      histoManager.fill("ZMETMassPresel",i,mTZZ,Weight);
      histoManager.fill("dPhiJPresel",i,dPhiJMET,Weight);
      histoManager.fill("dPhiZPresel",i,dphimetz*180/TMath::Pi(),Weight);
      if( (*jEt).size() != 0) {
      	histoManager.fill("FJetPresel",i,(*jEt)[0],Weight);
      	histoManager.fill("BTagTKFJetPresel",i,(*jBtagTkC)[0],Weight);
      	histoManager.fill("BTagSMFJetPresel",i,(*jBtagSoftM)[0],Weight);
    	histoManager.fill("SJetPresel",i,(*jEt)[1],Weight);
      	histoManager.fill("BTagTKSJetPresel",i,(*jBtagTkC)[1],Weight);
      	histoManager.fill("BTagSMSJetPresel",i,(*jBtagSoftM)[1],Weight);
      }
      else {
      	histoManager.fill("FJetPresel",i,-1000,Weight);
      	histoManager.fill("BTagTKFJetPresel",i,-1000,Weight);
      	histoManager.fill("BTagSMFJetPresel",i,-1000,Weight);
 	histoManager.fill("SJetPresel",i,-1000,Weight);
      	histoManager.fill("BTagTKSJetPresel",i,-1000,Weight);
      	histoManager.fill("BTagSMSJetPresel",i,-1000,Weight);
      }
      if( (*pEt).size() != 0) {
	histoManager.fill("FPhotPresel",i,(*pEt)[0],Weight);
      }
      else {
	histoManager.fill("FPhotPresel",i, -1000,Weight);
      }
      

      if(!MakeCut<int>(jMult30,MaxJet,typeJetV, i, "jet veto", Weight2 ) ) continue;
      
      // // //Recompute multiplicity for SYst
      // int nJ=0;
      // for(size_t ij=0;ij<(*jEt).size();ij++) {
	
      // 	if( (*jEt)[ij]*0.97 >= 30 && fabs((*jEta)[ij]) <= 5 )
      // 	  nJ++;
      // }
      
      // if(!MakeCut<int>(nJ,MaxJet,typeJetV, i, "jet veto", Weight2 ) ) continue;

      Weight = Weight2;

      histoManager.fill("ZMassJVeto",i,mll,Weight);
      histoManager.fill("ZPtJVeto",i,qTll,Weight);
      histoManager.fill("ZYJVeto",i,yll,Weight);
      histoManager.fill("LepMultJVeto",i,nLoose,Weight);
      histoManager.fill("METJVeto",i,MET,Weight);
      histoManager.fill("CorMETJVeto",i,CorMET,Weight);
      histoManager.fill("VtxCleanMETJVeto",i,min(cpuMETVtx,MET),Weight);
      histoManager.fill("SumCleanMETJVeto",i,min(cpuMETSum,MET),Weight);
      histoManager.fill("SqrtCleanMETJVeto",i,min(cpuMETSqrt,MET),Weight);
      histoManager.fill("CorMinMETJVeto",i,CorMMET,Weight);
      histoManager.fill("BalJVeto",i,MET/qTll,Weight);
      histoManager.fill("ZMETMassJVeto",i,mTZZ,Weight);
      histoManager.fill("dPhiJJVeto",i,dPhiJMET,Weight);
      histoManager.fill("dPhiZJVeto",i,dphimetz*180/TMath::Pi(),Weight);
      if( (*jEt).size() != 0) {
      	histoManager.fill("FJetJVeto",i,(*jEt)[0],Weight);
      	histoManager.fill("BTagTKFJetJVeto",i,(*jBtagTkC)[0],Weight);
      	histoManager.fill("BTagSMFJetJVeto",i,(*jBtagSoftM)[0],Weight);
      }
      else {
      	histoManager.fill("FJetJVeto",i,-1000,Weight);
      	histoManager.fill("BTagTKFJetJVeto",i,-1000,Weight);
      	histoManager.fill("BTagSMFJetJVeto",i,-1000,Weight);
      }
      if( (*jEt).size() > 0) {
      histoManager.fill("SJetJVeto",i,(*jEt)[1],Weight);
      histoManager.fill("BTagTKSJetJVeto",i,(*jBtagTkC)[1],Weight);
      histoManager.fill("BTagSMSJetJVeto",i,(*jBtagSoftM)[1],Weight);
      }      
      else {
      	histoManager.fill("SJetJVeto",i,-1000,Weight);
      	histoManager.fill("BTagTKSJetJVeto",i,-1000,Weight);
      	histoManager.fill("BTagSMSJetJVeto",i,-1000,Weight);
      }

      if(nLoose == 1) {
	histoManager.fill("LepPdgJVeto",i,(*LPdg)[2],Weight);
	histoManager.fill("LMTJVeto",i,(*LMT)[2],Weight);
      }
      if( (*pEt).size() != 0) {
	histoManager.fill("FPhotJVeto",i,(*pEt)[0],Weight);
      }
      else {
	histoManager.fill("FPhotJVeto",i, -1000,Weight);
      }


      if(!MakeCut<float>(CorMMET,METcut,">", i, "met", Weight ) ) continue;
      
      histoManager.fill("ZMassMET",i,mll,Weight);
      histoManager.fill("ZPtMET",i,qTll,Weight);
      histoManager.fill("ZYMET",i,yll,Weight);
      histoManager.fill("LepMultMET",i,nLoose,Weight);
      histoManager.fill("METMET",i,MET,Weight);
      histoManager.fill("CorMETMET",i,CorMET,Weight);
      histoManager.fill("VtxCleanMETMET",i,min(cpuMETVtx,MET),Weight);
      histoManager.fill("SumCleanMETMET",i,min(cpuMETSum,MET),Weight);
      histoManager.fill("SqrtCleanMETMET",i,min(cpuMETSqrt,MET),Weight);
      histoManager.fill("CorMinMETMET",i,CorMMET,Weight);
      histoManager.fill("BalMET",i,MET/qTll,Weight);
      histoManager.fill("ZMETMassMET",i,mTZZ,Weight);
      histoManager.fill("dPhiJMET",i,dPhiJMET,Weight);
      histoManager.fill("dPhiZMET",i,dphimetz*180/TMath::Pi(),Weight);
      if( (*jEt).size() != 0) {
      	histoManager.fill("FJetMET",i,(*jEt)[0],Weight);
      	histoManager.fill("BTagTKFJetMET",i,(*jBtagTkC)[0],Weight);
      	histoManager.fill("BTagSMFJetMET",i,(*jBtagSoftM)[0],Weight);
      }
      else {
      	histoManager.fill("FJetMET",i,-1000,Weight);
      	histoManager.fill("BTagTKFJetMET",i,-1000,Weight);
      	histoManager.fill("BTagSMFJetMET",i,-1000,Weight);
      }
      if( (*jEt).size() > 0) {
      	histoManager.fill("SJetMET",i,(*jEt)[1],Weight);
      	histoManager.fill("BTagTKSJetMET",i,(*jBtagTkC)[1],Weight);
      	histoManager.fill("BTagSMSJetMET",i,(*jBtagSoftM)[1],Weight);
      }      
      else {
      	histoManager.fill("SJetMET",i,-1000,Weight);
      	histoManager.fill("BTagTKSJetMET",i,-1000,Weight);
      	histoManager.fill("BTagSMSJetMET",i,-1000,Weight);
      }
      if( (*pEt).size() != 0) {
	histoManager.fill("FPhotMET",i,(*pEt)[0],Weight);
      }
      else {
	histoManager.fill("FPhotMET",i, -1000,Weight);
      }
      
      
      if(!MakeCut<float>(MET/qTll, balMin ,"[]", i, "balance", Weight, balMax ) ) continue;

      histoManager.fill("ZMassBal",i,mll,Weight);
      histoManager.fill("ZPtBal",i,qTll,Weight);
      histoManager.fill("ZYBal",i,yll,Weight);
      histoManager.fill("LepMultBal",i,nLoose,Weight);
      histoManager.fill("METBal",i,MET,Weight);
      histoManager.fill("CorMETBal",i,CorMET,Weight);
      histoManager.fill("VtxCleanMETBal",i,min(cpuMETVtx,MET),Weight);
      histoManager.fill("SumCleanMETBal",i,min(cpuMETSum,MET),Weight);
      histoManager.fill("SqrtCleanMETBal",i,min(cpuMETSqrt,MET),Weight);
      histoManager.fill("CorMinMETBal",i,CorMMET,Weight);
      histoManager.fill("BalBal",i,MET/qTll,Weight);
      histoManager.fill("ZMETMassBal",i,mTZZ,Weight);
      histoManager.fill("dPhiJBal",i,dPhiJMET,Weight);
      histoManager.fill("dPhiZBal",i,dphimetz*180/TMath::Pi(),Weight);
      if( (*jEt).size() != 0) {
      	histoManager.fill("FJetBal",i,(*jEt)[0],Weight);
      	histoManager.fill("BTagTKFJetBal",i,(*jBtagTkC)[0],Weight);
      	histoManager.fill("BTagSMFJetBal",i,(*jBtagSoftM)[0],Weight);
      }
      else {
      	histoManager.fill("FJetBal",i,-1000,Weight);
      	histoManager.fill("BTagTKFJetBal",i,-1000,Weight);
      	histoManager.fill("BTagSMFJetBal",i,-1000,Weight);
      }
      if( (*jEt).size() > 0) {
      	histoManager.fill("SJetBal",i,(*jEt)[1],Weight);
      	histoManager.fill("BTagTKSJetBal",i,(*jBtagTkC)[1],Weight);
      	histoManager.fill("BTagSMSJetBal",i,(*jBtagSoftM)[1],Weight);
      }      
      else {
      	histoManager.fill("SJetBal",i,-1000,Weight);
      	histoManager.fill("BTagTKSJetBal",i,-1000,Weight);
      	histoManager.fill("BTagSMSJetBal",i,-1000,Weight);
      }
      if( (*pEt).size() != 0) {
	histoManager.fill("FPhotBal",i,(*pEt)[0],Weight);
      }
      else {
	histoManager.fill("FPhotBal",i, -1000,Weight);
      }
      
      // float dpjm=100;
      // for(size_t ij=0;ij<(*jEt).size();ij++) {
      // 	bool lmatch=false;
      // 	for(size_t il=0;il<(*LPt).size();il++) {
      // 	  if( isMatch( (*jEt)[ij], (*jEta)[ij], (*jPhi)[ij],
      // 		       (*LPt)[il], (*LEta)[il], (*LPhi)[il] ) )
      // 	    { lmatch=true; 
      // 	      break;}
      // 	}

      // 	if(!lmatch)
      // 	  { dpjm = (*jDPhiMET)[ij]; break; }
      // }


      if(!MakeCut<float>( fabs(dPhiJMET),dPhiJMCut,">=", i, "#Delta#Phi(jet,met)", Weight ) ) continue;
      
      //if(!MakeCut<float>( fabs(dpjm),dPhiJMCut,">", i, "#Delta#Phi(jet,met)", Weight ) ) continue;

      histoManager.fill("ZMassDPhiJ",i,mll,Weight);
      histoManager.fill("ZPtDPhiJ",i,qTll,Weight);
      histoManager.fill("ZYDPhiJ",i,yll,Weight);
      histoManager.fill("LepMultDPhiJ",i,nLoose,Weight);
      histoManager.fill("METDPhiJ",i,MET,Weight);
      histoManager.fill("CorMETDPhiJ",i,CorMET,Weight);
      histoManager.fill("VtxCleanMETDPhiJ",i,min(cpuMETVtx,MET),Weight);
      histoManager.fill("SumCleanMETDPhiJ",i,min(cpuMETSum,MET),Weight);
      histoManager.fill("SqrtCleanMETDPhiJ",i,min(cpuMETSqrt,MET),Weight);
      histoManager.fill("CorMinMETDPhiJ",i,CorMMET,Weight);
      histoManager.fill("BalDPhiJ",i,MET/qTll,Weight);
      histoManager.fill("ZMETMassDPhiJ",i,mTZZ,Weight);
      histoManager.fill("dPhiJDPhiJ",i,dPhiJMET,Weight);
      histoManager.fill("dPhiZDPhiJ",i,dphimetz*180/TMath::Pi(),Weight);
      if( (*jEt).size() != 0) {
      	histoManager.fill("FJetDPhiJ",i,(*jEt)[0],Weight);   
      	histoManager.fill("BTagTKFJetDPhiJ",i,(*jBtagTkC)[0],Weight);
      	histoManager.fill("BTagSMFJetDPhiJ",i,(*jBtagSoftM)[0],Weight);
      }
      else {
      	histoManager.fill("FJetDPhiJ",i,-1000,Weight);
      	histoManager.fill("BTagTKFJetDPhiJ",i,-1000,Weight);
      	histoManager.fill("BTagSMFJetDPhiJ",i,-1000,Weight);
      }
      if( (*jEt).size() > 0) {
      	histoManager.fill("SJetDPhiJ",i,(*jEt)[1],Weight);
      	histoManager.fill("BTagTKSJetDPhiJ",i,(*jBtagTkC)[1],Weight);
      	histoManager.fill("BTagSMSJetDPhiJ",i,(*jBtagSoftM)[1],Weight);
      } 
      else {
      	histoManager.fill("SJetDPhiJ",i,-1000,Weight);
      	histoManager.fill("BTagTKSJetDPhiJ",i,-1000,Weight);
      	histoManager.fill("BTagSMSJetDPhiJ",i,-1000,Weight);
      }
      if( (*pEt).size() != 0) {
	histoManager.fill("FPhotDPhiJ",i,(*pEt)[0],Weight);
      }
      else {
	histoManager.fill("FPhotDPhiJ",i, -1000,Weight);
      }
      
      if(nLoose == 1) {
	histoManager.fill("LepPdgDPhiJ",i,(*LPdg)[2],Weight);
	histoManager.fill("LMTDPhiJ",i,(*LMT)[2],Weight);
      }

      
      if(!MakeCut<double>( fabs(dphimetz*180/TMath::Pi()),(double)dPhiZMCut,"<=", i, "#Delta#Phi(Z,met)", Weight ) ) continue;
      
      histoManager.fill("ZMassDPhiZ",i,mll,Weight);
      histoManager.fill("ZPtDPhiZ",i,qTll,Weight);
      histoManager.fill("ZYDPhiZ",i,yll,Weight);
      histoManager.fill("LepMultDPhiZ",i,nLoose,Weight);
      histoManager.fill("METDPhiZ",i,MET,Weight);
      histoManager.fill("CorMETDPhiZ",i,CorMET,Weight);
      histoManager.fill("VtxCleanMETDPhiZ",i,min(cpuMETVtx,MET),Weight);
      histoManager.fill("SumCleanMETDPhiZ",i,min(cpuMETSum,MET),Weight);
      histoManager.fill("SqrtCleanMETDPhiZ",i,min(cpuMETSqrt,MET),Weight);
      histoManager.fill("CorMinMETDPhiZ",i,CorMMET,Weight);
      histoManager.fill("BalDPhiZ",i,MET/qTll,Weight);
      histoManager.fill("ZMETMassDPhiZ",i,mTZZ,Weight);
      histoManager.fill("dPhiJDPhiZ",i,dPhiJMET,Weight);
      histoManager.fill("dPhiZDPhiZ",i,dphimetz*180/TMath::Pi(),Weight);
      if( (*jEt).size() != 0) {
      	histoManager.fill("FJetDPhiZ",i,(*jEt)[0],Weight);
      	histoManager.fill("BTagTKFJetDPhiZ",i,(*jBtagTkC)[0],Weight);
      	histoManager.fill("BTagSMFJetDPhiZ",i,(*jBtagSoftM)[0],Weight);
      }
      else {
      	histoManager.fill("FJetDPhiZ",i,-1000,Weight);
      	histoManager.fill("BTagTKFJetDPhiZ",i,-1000,Weight);
      	histoManager.fill("BTagSMFJetDPhiZ",i,-1000,Weight);
      }
      if( (*jEt).size() > 0) {
      	histoManager.fill("SJetDPhiZ",i,(*jEt)[1],Weight);
      	histoManager.fill("BTagTKSJetDPhiZ",i,(*jBtagTkC)[1],Weight);
      	histoManager.fill("BTagSMSJetDPhiZ",i,(*jBtagSoftM)[1],Weight);
      }      
      else {
      	histoManager.fill("SJetDPhiZ",i,-1000,Weight);
      	histoManager.fill("BTagTKSJetDPhiZ",i,-1000,Weight);
      	histoManager.fill("BTagSMSJetDPhiZ",i,-1000,Weight);
      }
      if( (*pEt).size() != 0) {
	histoManager.fill("FPhotDPhiZ",i,(*pEt)[0],Weight);
      }
      else {
	histoManager.fill("FPhotDPhiZ",i, -1000,Weight);
      }
      
      if(nLoose == 1) {
	histoManager.fill("LepPdgDPhiZ",i,(*LPdg)[2],Weight);
	histoManager.fill("LMTDPhiZ",i,(*LMT)[2],Weight);
      }


      //Btag
      bool tag=false;
      for(size_t ij=0;ij<(*jEt).size();ij++) {
      	if( (*jEt)[ij] < 10 || fabs((*jEta)[ij]) > 5 ) continue;
      	if( (*jBtagTkC)[ij] > 2.5 || (*jBtagSoftM)[ij] > 0.3)
      	  { tag = true;}
      }
      if(!MakeCut<bool>( tag, bTagFlag,"=", i, "b tagging", Weight ) ) continue;
      
      histoManager.fill("ZMassBTag",i,mll,Weight);
      histoManager.fill("ZPtBTag",i,qTll,Weight);
      histoManager.fill("ZYBTag",i,yll,Weight);
      histoManager.fill("LepMultBTag",i,nLoose,Weight);
      histoManager.fill("METBTag",i,MET,Weight);
      histoManager.fill("CorMETBTag",i,CorMET,Weight);
      histoManager.fill("VtxCleanMETBTag",i,min(cpuMETVtx,MET),Weight);
      histoManager.fill("SumCleanMETBTag",i,min(cpuMETSum,MET),Weight);
      histoManager.fill("SqrtCleanMETBTag",i,min(cpuMETSqrt,MET),Weight);
      histoManager.fill("CorMinMETBTag",i,CorMMET,Weight);
      histoManager.fill("BalBTag",i,MET/qTll,Weight);
      histoManager.fill("ZMETMassBTag",i,mTZZ,Weight);
      histoManager.fill("dPhiJBTag",i,dPhiJMET,Weight);
      histoManager.fill("dPhiZBTag",i,dphimetz*180/TMath::Pi(),Weight);
      if( (*jEt).size() != 0) {
      	histoManager.fill("FJetBTag",i,(*jEt)[0],Weight);
      	histoManager.fill("BTagTKFJetBTag",i,(*jBtagTkC)[0],Weight);
      	histoManager.fill("BTagSMFJetBTag",i,(*jBtagSoftM)[0],Weight);
      }
      else {
      	histoManager.fill("FJetBTag",i,-1000,Weight);
      	histoManager.fill("BTagTKFJetBTag",i,-1000,Weight);
      	histoManager.fill("BTagSMFJetBTag",i,-1000,Weight);
      }
      if( (*jEt).size() > 0) {
      	histoManager.fill("SJetBTag",i,(*jEt)[1],Weight);
      	histoManager.fill("BTagTKSJetBTag",i,(*jBtagTkC)[1],Weight);
      	histoManager.fill("BTagSMSJetBTag",i,(*jBtagSoftM)[1],Weight);
      }      
      else {
      	histoManager.fill("SJetBTag",i,-1000,Weight);
      	histoManager.fill("BTagTKSJetBTag",i,-1000,Weight);
      	histoManager.fill("BTagSMSJetBTag",i,-1000,Weight);
      }
      if( (*pEt).size() != 0) {
	histoManager.fill("FPhotBTag",i,(*pEt)[0],Weight);
      }
      else {
	histoManager.fill("FPhotBTag",i, -1000,Weight);
      }

      if(nLoose == 1) {
	histoManager.fill("LepPdgBTag",i,(*LPdg)[2],Weight);
	histoManager.fill("LMTBTag",i,(*LMT)[2],Weight);

	if( fabs( (*LPdg)[2]) ==11) {
	  int cat=0; if( fabs( (*LEta)[2] ) > 1.5 ) cat=1;
	  histoManager.fill("AddLepCategE",i,cat,Weight);
	}
	else {
	  int cat=0; if( fabs( (*LEta)[2] ) > 1.3 ) cat=1;
	  histoManager.fill("AddLepCategM",i,cat,Weight);

	}

      }

      if(name[i].substr(0,2)=="WZ") 
	{
	  //  cout<<(*gLPdg)[0]<<"   "<< (*gLPdg)[1]<<endl;
	  if( (int)(fabs((*gLPdg)[0]))%2 != 0) {

	    if( BasicAcc( (*gLPt)[0], (*gLEta)[0],  (*gLPdg)[0] ) ) {
	      histoManager.fill("LepAcc",i,1,Weight);
	    }
	    else {
	      histoManager.fill("LepAcc",i,0,Weight);
	    }
	  } 
	  if( (int)(fabs((*gLPdg)[1]))%2 != 0) {
	    if( BasicAcc( (*gLPt)[1], (*gLEta)[1],  (*gLPdg)[1] ) ) {
	      histoManager.fill("LepAcc",i,1,Weight);
	    }
	    else {
	      histoManager.fill("LepAcc",i,0,Weight);
	    }

	  } 
	  
	}
      

      
      if(!MakeCut<int>( nLoose, MaxLepton,"<", i, "lepton veto", Weight ) ) continue;
      
      histoManager.fill("ZMassLVeto",i,mll,Weight);
      histoManager.fill("ZPtLVeto",i,qTll,Weight);
      histoManager.fill("ZYLVeto",i,yll,Weight);
      histoManager.fill("LepMultLVeto",i,nLoose,Weight);
      histoManager.fill("METLVeto",i,MET,Weight);
      histoManager.fill("CorMETLVeto",i,CorMET,Weight);
      histoManager.fill("VtxCleanMETLVeto",i,min(cpuMETVtx,MET),Weight);
      histoManager.fill("SumCleanMETLVeto",i,min(cpuMETSum,MET),Weight);
      histoManager.fill("SqrtCleanMETLVeto",i,min(cpuMETSqrt,MET),Weight);
      histoManager.fill("CorMinMETLVeto",i,CorMMET,Weight);
      histoManager.fill("BalLVeto",i,MET/qTll,Weight);
      histoManager.fill("ZMETMassLVeto",i,mTZZ,Weight);
      histoManager.fill("dPhiJLVeto",i,dPhiJMET,Weight);
      histoManager.fill("dPhiZLVeto",i,dphimetz*180/TMath::Pi(),Weight);
      if( (*jEt).size() != 0) {
      	histoManager.fill("FJetLVeto",i,(*jEt)[0],Weight);
      	histoManager.fill("BTagTKFJetLVeto",i,(*jBtagTkC)[0],Weight);
      	histoManager.fill("BTagSMFJetLVeto",i,(*jBtagSoftM)[0],Weight);
      }
      else {
      	histoManager.fill("FJetLVeto",i,-1000,Weight);
      	histoManager.fill("BTagTKFJetLVeto",i,-1000,Weight);
      	histoManager.fill("BTagSMFJetLVeto",i,-1000,Weight);
      }
      if( (*jEt).size() > 0) {
      	histoManager.fill("SJetLVeto",i,(*jEt)[1],Weight);
      	histoManager.fill("BTagTKSJetLVeto",i,(*jBtagTkC)[1],Weight);
      	histoManager.fill("BTagSMSJetLVeto",i,(*jBtagSoftM)[1],Weight);
      }      
      else {
      	histoManager.fill("SJetLVeto",i,-1000,Weight);
      	histoManager.fill("BTagTKSJetLVeto",i,-1000,Weight);
      	histoManager.fill("BTagSMSJetLVeto",i,-1000,Weight);
      }
      if( (*pEt).size() != 0) {
	histoManager.fill("FPhotLVeto",i,(*pEt)[0],Weight);
      }
      else {
	histoManager.fill("FPhotLVeto",i, -1000,Weight);
      }

      //Category
      int cat=10;
      if( fabs( (*LEta)[0]) < 1.5 &&  fabs( (*LEta)[1]) < 1.5 )
	cat =0;
      if( ( fabs( (*LEta)[0]) < 1.5 &&  fabs( (*LEta)[1]) > 1.5 ) || 
	  ( fabs( (*LEta)[0]) > 1.5 &&  fabs( (*LEta)[1]) < 1.5 ) )
	{
	  if( fabs( (*LEta)[0]) < 1.5 ) {
	    if( ((*LID)[0])>=7 )
	      cat =1; // tight in barrel
	    else
	      cat =2; // tight in endcap
	  }
	  else {
	    if( ((*LID)[0])>=7 )
	      cat =2; // tight in endcap
	    else
	      cat =1; // tight in barrel
	  }
	}
      if( fabs( (*LEta)[0]) > 1.5 &&  fabs( (*LEta)[1]) > 1.5 )
	cat =3;

      histoManager.fill("ZCateg",i, cat, Weight);

      //FIXME
    
      // bool failConv=false;
      // if( expI[0] || expI[1])
      // 	failConv=true;
      //cout<<expI[0]<<"  "<< expI[1]<<"    "<<failConv<<endl;
      //if(failConv) continue;
      //if(!MakeCut<bool>( failConv, false ,"=", i, "conv", Weight ) ) continue;

       if(!MakeCut<float>( mll, advMassMin,"[]", i, "mass", Weight, advMassMax ) ) continue;

      histoManager.fill("ZMassMass",i,mll,Weight);
      histoManager.fill("ZPtMass",i,qTll,Weight);
      histoManager.fill("ZYMass",i,yll,Weight);
      histoManager.fill("LepMultMass",i,nLoose,Weight);
      histoManager.fill("METMass",i,MET,Weight);
      histoManager.fill("CorMETMass",i,CorMET,Weight);
      histoManager.fill("VtxCleanMETMass",i,min(cpuMETVtx,MET),Weight);
      histoManager.fill("SumCleanMETMass",i,min(cpuMETSum,MET),Weight);
      histoManager.fill("SqrtCleanMETMass",i,min(cpuMETSqrt,MET),Weight);
      histoManager.fill("CorMinMETMass",i,CorMMET,Weight);
      histoManager.fill("BalMass",i,MET/qTll,Weight);
      histoManager.fill("ZMETMassMass",i,mTZZ,Weight);
      histoManager.fill("dPhiJMass",i,dPhiJMET,Weight);
      histoManager.fill("dPhiZMass",i,dphimetz*180/TMath::Pi(),Weight);
      if( (*jEt).size() != 0) {
      	histoManager.fill("FJetMass",i,(*jEt)[0],Weight);
      	histoManager.fill("BTagTKFJetMass",i,(*jBtagTkC)[0],Weight);
      	histoManager.fill("BTagSMFJetMass",i,(*jBtagSoftM)[0],Weight);
      }
      else {
      	histoManager.fill("FJetMass",i,-1000,Weight);
      	histoManager.fill("BTagTKFJetMass",i,-1000,Weight);
      	histoManager.fill("BTagSMFJetMass",i,-1000,Weight);
      }
      if( (*jEt).size() > 0) {
      	histoManager.fill("SJetMass",i,(*jEt)[1],Weight);
      	histoManager.fill("BTagTKSJetMass",i,(*jBtagTkC)[1],Weight);
      	histoManager.fill("BTagSMSJetMass",i,(*jBtagSoftM)[1],Weight);
      }      
      else {
      	histoManager.fill("SJetMass",i,-1000,Weight);
      	histoManager.fill("BTagTKSJetMass",i,-1000,Weight);
      	histoManager.fill("BTagSMSJetMass",i,-1000,Weight);
      }
      if( (*pEt).size() != 0) {
	histoManager.fill("FPhotMass",i,(*pEt)[0],Weight);
      }
      else {
	histoManager.fill("FPhotMass",i, -1000,Weight);
      }
   
      if(MaxJet==2) {
	
	  histoManager.fill("J1DR1",i, (*jDR1)[0], Weight);
	  histoManager.fill("J1DR2",i, (*jDR2)[0], Weight);

	  histoManager.fill("J2DR1",i, (*jDR1)[1], Weight);
	  histoManager.fill("J2DR2",i, (*jDR2)[1], Weight);
      }


     if(i==nt) {
	vector<vector<float> > lep(2,vector<float>(3,0));
	lep[0][0] = (*LPt)[0];       lep[1][0] = (*LPt)[1];
	lep[0][1] = (*LEta)[0];      lep[1][1] = (*LEta)[1];
	lep[0][2] = (*LID)[0];       lep[1][2] = (*LID)[1];
	LooseElec.push_back(lep);
      }


      //=====================================================================================

      if(name[i].substr(0,2)=="di" || name[i].substr(0,2)=="ZZ"
      	 ||  name[i].substr(0,2)=="ZV") {
      	string proc = FindProcess( ie, i );
      	if( proc.substr(0,2)=="ZZ") {
      	  StoreSelectedEvents( mVV , qTll, Weight );

	  // histoManager.fill("McZZMass",i,mVV, Weight );
      	}
      }
      
      //End Filling*********************
      // if(name[i]=="Zll" && METcut>4.) {
      // 	cout<<fileName<<"   "<<Event<<"   "<<JetMultiplicity<<"   "<<LJet[0]<<"   "<<NNmetOut<<"    "<<dPhiJetMET<<"   "<<TrkMet[0]<<"    "<<SigniCorProjMET<<"   "<<CorProjMET<<"   "<<Z[0]<<endl;
      // }
      if(name[i].substr(0,4)=="data")
	{
	  std::pair<int,int> tmp(run,event);
	  string t1(""),t2("");
	  std::pair<string,string> tmp2( t1, t2 );
	  std::pair< std::pair<int,int> , std::pair<string,string> > tmp3(tmp,tmp2);
	  Events.push_back(tmp3);
	  EvtsInFile.push_back(event);

	  // if((*LID)[0]<7 && (*LID)[1]<7) {
	  //   cout<<" L1 "<<(*LPt)[0]<<"   "<<(*LEta)[0]<<"   "<<(*LID)[0]<<endl;
	  //   cout<<" L2 "<<(*LPt)[1]<<"   "<<(*LEta)[1]<<"   "<<(*LID)[1]<<endl;
	  // }
	  // else {
	  //   cout<<(*LID)[0]<<"      "<<(*LID)[1]<<endl;
	  // }

	   
	}
      if(name[i]=="ZZ #rightarrow 2l2#nu")
	{ S++; 
	}
      else if(name[i]!="data")
	{B++; }
     
      NumberEntries[i]++;

    }//End events

    //Pointer deletion =====
    delete  categ;
    delete  c_str;
    delete  c_bits;
    delete  LPt;
    delete  LEta;
    delete  LPhi;
    delete  LID;
    delete  LPdg;
    delete  LMT;
    delete  LTrkIso;
    delete  LEcalIso;
    delete  LHcalIso;
    delete  jEt;
    delete  jEta;
    delete  jPhi;
    delete  jDR1;
    delete  jDR2;
    delete  jDPhiMET;
    delete  jBtagTkC;
    delete  jBtagSoftM;

    // delete categ;
    // delete c_str;
    // delete c_bits;

    //Additionnal tree deletion
    delete gLPt;
    delete gLEta;
    delete gLPdg;
    delete MCTree;
    
    delete LepTree;

    //=====================


  } //End datasets
 
  cout<<" End Filling , S="<<S<<" ;  B= "<<B<<endl;
  cout<<" detail "<<endl;
  for(int i=0;i<nt+1;i++) {
    cout<<" --->  "<<name[i]<<"  "<< NumberEntries[i]<<endl;
  }
  
  Histos = histoManager.GetHistos();
  Histos2D = histoManager.GetHistos2D();
  Profiles = histoManager.GetProfiles();
  
  for(size_t i=0;i<Histos.size();i++)
    Weighted.push_back(false);
  
  cout<<"Histos "<<Histos.size()<<"   "<<Histos2D.size()<<endl;
  
  cout<<" test "<<endl;

}



void
DibosonUseTree::PrepareDatasets() {

  //We
  //  reposi="ZZAnalysis";
 

  colors.push_back(kBlack); //End Line
  for(int unsigned i=0;i<datasets.size();i++) {
    FillAddWeight(datasets[i]);
    Type.push_back(datasets[i]);
  }
  name.push_back("data");

  for(int i=0;i<nt;i++) {
    cout<<name[nt-1-i]<<endl;
    Sorder.push_back(name[nt-1-i]);
  }
  
  cout<<" End dataset preparation : nt="<<nt<<endl;

}


TH1F* 
DibosonUseTree::GetVtxCorrection(string obs,int nds, int bin) {

  //Protection
  if(NVert!=1) {
    cout<<" Warning, no requirement of only 1 Vtx event !!!!!"<<endl;
    abort();
  }
  
  float contamination = 0.069;
  float Systcontamination = 0.02;

  int Obs = histoManager.FindNVar(obs);
  int Ndata= histoManager.access(Obs,nds);
  TH1F* dataTMP = (TH1F*)Histos[ Ndata ]->Clone();

  float Int = dataTMP->Integral(0, dataTMP->GetNbinsX()+1);
  
  float Nevts = Int*contamination;
  Obs = histoManager.FindNVar(obs+"Control");
  Ndata= histoManager.access(Obs,nds);
  TH1F* dataTMP2 = (TH1F*)Histos[ Ndata ]->Clone();
 
  TH1F* Syst = (TH1F*)Histos[ Ndata ]->Clone();
  float NevtsSyst = Int*Systcontamination;

  Int = dataTMP2->Integral(0, dataTMP2->GetNbinsX()+1);

  dataTMP2->Scale( Nevts/Int );
  dataTMP2->Rebin(bin);
  dataTMP->Rebin(bin);

  Syst->Scale( NevtsSyst/Int );
  Syst->Rebin(bin);
 
  VtxCorSyst.clear();
  if(VtxCorSyst.size()==0)
    for(int i=1;i< Syst->GetNbinsX();i++)
      VtxCorSyst.push_back(Syst->GetBinContent(i) );
  
  cout<<" Contamination  "<<contamination<<" ;  "<<Nevts<<"  ;  "<<Int<<"   ;  "<<Nevts/Int<<"   "<<endl;
  if(Draw3on1)
    getRatio =true;
  if(!getRatio) {
    string nameO =  (histoManager.FindLeg( observable )).first;
    histoManager.RatioDistributionsNorm( (TH1F*)dataTMP->Clone(),(TH1F*)dataTMP2->Clone(), nameO, RangeX);
    c2->cd();
    getRatio=true;
  }
  return dataTMP2;

}



vector<vector<vector<float> > >  
DibosonUseTree::GetFitWeight() {

  vector<int> NumberEntries(nt+1,0);
  
  cout<<" ==================================================== "<<endl;
  cout<<" ================= En cours de weight =============== "<<endl;

  //Needed by the NN
  vector<double> inputVec(8,0);

  //sÃ©lection
  for(int i=0;i<nt+1;i++) {

    //Déclaration des variables
    
    int nVertex;
    tChains[i]->SetBranchAddress("nVertex", &nVertex);

    cout<<" beginning tree " <<name[i]<<"   "<<tChains[i]->GetEntries()<<endl;
    //Boucle et sélection, remplissage histos
    int ent = tChains[i]->GetEntries();
    
    float Weight;

    for(int ie=0;ie<ent;ie++) {
      
      if(ie%1000==0)
	fprintf(stdout, "Done : %d / %d\r", ie, ent );

      if(i!=nt)
	Weight = GetWeight(i, ie);
      else
	Weight = 1;
      
      tChains[i]->GetEntry(ie);      

      histoManager.fill("NVertexControl",i,nVertex, Weight );

      
    }//End events
  }// End datasets
 
    
  Histos = histoManager.GetHistos();
  Histos2D = histoManager.GetHistos2D();
    
  vector<vector<vector<float> > > Wghts(5,vector<vector<float> >(2,vector<float>(2,1)));
 
  return Wghts;
    
}
 
float DibosonUseTree::Rapidity(float p, float mass, float pt, float eta) {

  //compute pz
  float pz = sqrt(p*p-pt*pt);

  //compute E
  float E = sqrt(p*p + mass*mass);

  //and now rapidity
  float y = 0.5 * log( (E+pz) / (E-pz) );

  return y;

}



void DibosonUseTree::Get2DNumber(string obs, string chan, ostream& os) {

  int Obs = histoManager.FindNVar2D(obs);
  if(Obs==-1) {
    cout<<" Be careful, no such observable "<<obs<<endl;
    return;
  }

  int NDS= GetNbyName(chan);
  int Nhisto = histoManager.access2D(Obs, NDS);

  TH2F* histo = (TH2F*)Histos2D[ Nhisto ]->Clone();

  vector<float> temp = histoManager.GetTemplate2D(Obs);
  XTitle = (histoManager.FindLeg2D( obs )).first;
  YTitle = (histoManager.FindLeg2D( obs )).second;

  float total=0;

  if(NDS!=nt)
    os<<chan<<"   "<<GetWeight(NDS,1)<<endl;
  else
    os<<chan<<"   "<<1<<endl; 
  for(int ix=0;ix<(int)temp[0];ix++) {

      for(int iy=0;iy<(int)temp[1];iy++) {
      
	os<< setw(8)<< histo->GetBinContent(ix+1,iy+1)<<"   ";

	if(ix==1)
	  total += histo->GetBinContent(ix+1,iy+1);

    }
      os<<endl;
  
    
  }


  cout<<" Total -> "<<total<<endl<<endl;

}


void
DibosonUseTree::LoadDBWeightZZ() {

  ifstream in("/home/mmarionn/Documents/CMS/CMSDatabase/kFactorDatabase", ios::in );  

  vector<vector<float> > tmpw;

  if(in) {
    cout<<" Loading Database "<<endl;
    while(!in.eof()) {
      vector<float> tmpv(2,0);
      in >> tmpv[0] >> tmpv[1];
      tmpw.push_back(tmpv);
    }
  }
  else 
    {    cout<<" No DB loaded !! "<<endl; abort();}
  
  DBWeightsZZ = tmpw;


  ifstream injv("/home/mmarionn/Documents/CMS/CMSDatabase/kFactorDatabaseZZ_jv", ios::in );  

  vector<vector<float> > tmpwjv;

  if(injv) {
    cout<<" Loading Database Jet Veto "<<endl;
    while(!injv.eof()) {
      vector<float> tmpv(2,0);
      injv >> tmpv[0] >> tmpv[1];
      tmpwjv.push_back(tmpv);
    }
  }
  else 
    {    cout<<" No DB loaded !! "<<endl; abort();}
  
  DBWeightsZZJV = tmpwjv;


}

void
DibosonUseTree::LoadDBWeightWZ() {

  ifstream in("/home/mmarionn/Documents/CMS/CMSDatabase/kFactorDatabaseWZ", ios::in );  

  vector<vector<float> > tmpw;

  if(in) {
    cout<<" Loading Database "<<endl;
    while(!in.eof()) {
      vector<float> tmpv(2,0);
      in >> tmpv[0] >> tmpv[1];
      tmpw.push_back(tmpv);
    }
  }
  else 
    {    cout<<" No DB loaded !! "<<endl; abort();}
  
  DBWeightsWZ = tmpw;


  ifstream injv("/home/mmarionn/Documents/CMS/CMSDatabase/kFactorDatabaseWZ_jv", ios::in );  

  vector<vector<float> > tmpwjv;

  if(injv) {
    cout<<" Loading Database Jet Veto "<<endl;
    while(!injv.eof()) {
      vector<float> tmpv(2,0);
      injv >> tmpv[0] >> tmpv[1];
      tmpwjv.push_back(tmpv);
    }
  }
  else 
    {    cout<<" No DB loaded !! "<<endl; abort();}
  
  DBWeightsWZJV = tmpwjv;


}

void
DibosonUseTree::LoadDBWeightWW() {

  ifstream in("/home/mmarionn/Documents/CMS/CMSDatabase/kFactorDatabaseWW", ios::in );  

  vector<vector<float> > tmpw;

  if(in) {
    cout<<" Loading Database "<<endl;
    while(!in.eof()) {
      vector<float> tmpv(2,0);
      in >> tmpv[0] >> tmpv[1];
      tmpw.push_back(tmpv);
    }
  }
  else 
    {    cout<<" No DB loaded !! "<<endl; abort();}
  
  DBWeightsWW = tmpw;



  ifstream injv("/home/mmarionn/Documents/CMS/CMSDatabase/kFactorDatabaseWW_jv", ios::in );  

  vector<vector<float> > tmpwjv;

  if(injv) {
    cout<<" Loading Database Jet Veto "<<endl;
    while(!injv.eof()) {
      vector<float> tmpv(2,0);
      injv >> tmpv[0] >> tmpv[1];
      tmpwjv.push_back(tmpv);
    }
  }
  else 
    {    cout<<" No DB loaded !! "<<endl; abort();}
  
  DBWeightsWWJV = tmpwjv;

}



float 
DibosonUseTree::SearchWeightZZ(float Zpt, bool veto) {
  
  float w=1;
  
  vector<vector<float> > tmp;
  int n;

  if(!veto) {
    tmp = DBWeightsZZ;
    n = DBWeightsZZ.size();
  }
  else {
    tmp = DBWeightsZZJV;
    n = DBWeightsZZJV.size();
  }

  for(int i=0;i<n;i++) {
    if(i!=n-1) {
      if(Zpt >= tmp[i][0] && Zpt < tmp[i+1][0])
	{ w = tmp[i][1]; break;}
    }
    else
      { w = tmp[i][1]; }
  }
  
  return w;
}


float 
DibosonUseTree::SearchWeightWZ(float Zpt, bool veto) {
  
  float w=1;
  
 
  vector<vector<float> > tmp;
  int n;

  if(!veto) {
    tmp = DBWeightsWZ;
    n = DBWeightsWZ.size();
  }
  else {
    tmp = DBWeightsWZJV;
    n = DBWeightsWZJV.size();
  }

  for(int i=0;i<n;i++) {
    if(i!=n-1) {
      if(Zpt >= tmp[i][0] && Zpt < tmp[i+1][0])
	{ w = tmp[i][1]; break;}
    }
    else
      { w = tmp[i][1]; }
  }
  
  return w;
}


float 
DibosonUseTree::SearchWeightWW(float Wpt, bool veto) {
  
  float w=1;
  

  vector<vector<float> > tmp;
  int n;

  if(!veto) {
    tmp = DBWeightsWW;
    n = DBWeightsWW.size();
  }
  else {
    tmp = DBWeightsWWJV;
    n = DBWeightsWWJV.size();
  }

  for(int i=0;i<n;i++) {
    if(i!=n-1) {
      if(Wpt >= tmp[i][0] && Wpt < tmp[i+1][0])
	{ w = tmp[i][1]; break;}
    }
    else
      { w = tmp[i][1]; }
  }
  
  return w;
}




void 
DibosonUseTree::StoreSelectedEvents(float mZZgen, float ZPt, float w) {

  std::pair<float, float> p(mZZgen, w);
  selZZ.push_back(p);
  selZPt.push_back(ZPt);
} 




vector<vector<float> >
DibosonUseTree::GetReweightedYields(float f4M, float f5M, int nf4, int nf5) {

  vector<TF2*> curves= LoadDBaTGC();
  if(nf4%2 ==1)
    nf4++;
  if(nf5%2 ==1)
    nf5++;

  float f4v, f5v;
  vector<vector<float> > vRatio;
  // cout<<" Size "<<selZZ.size()<<endl;

  for(int i4=0;i4<nf4+1;i4++) {
    f4v = -f4M + i4*(2*f4M/nf4);
    if( f4v > 0 )
      f4v  = -f4M + 0.001 + i4*(2*f4M/nf4);

    //for histomap
    cout<<" -> "<<f4v<<"    ";//<<(float)f4M/nf4<<"    "<<f4M/nf4<<endl;
    f4v = f4v*100000; //0.01*100000 = 1000
    f4v = f4v - (int)f4v%1000;
    f4v = (int)f4v/1000;
    f4v = f4v/100.;
    cout<<f4v<<endl;

    for(int i5=0;i5<nf5+1;i5++) {
      f5v = -f5M + i5*(2*f5M/nf5);
      if( f5v > 0 )
	f5v  = -f5M + 0.001 + i5*(2*f5M/nf5);

      //for histomap
      // cout<<i4<<"  "<<i5<<" -> "<<f4v<<"   "<<f5v<<"     ";
      f5v = f5v*100000; //0.01*100000 = 1000
      f5v = f5v - (int)f5v%1000;
      f5v = (int)f5v/1000;
      f5v = f5v/100.;
     

      float sumW=0; int nbmzz;
      float sumNLOW=0;
      for(int unsigned ie=0;ie<selZZ.size();ie++) {
	
	nbmzz = (int)(selZZ[ie].first)/20;
	if(nbmzz >= 100 )
	  { nbmzz=99; }

	sumW +=	GetATGCWeight(f4v, f5v, nbmzz)*selZZ[ie].second;
	sumNLOW += selZZ[ie].second;

      }
      vector<float> vtmp(3,0);
      vtmp[0] = f4v;
      vtmp[1] = f5v;
      vtmp[2] = sumW/(sumNLOW);
      vRatio.push_back(vtmp);
      //cout<<f4v<<"\t"<<f5v<<"\t  --->  "<<sumW<<"\t "<<sumW/selZZ.size()<<endl;
    }
  }

  return vRatio;

}


void 
DibosonUseTree::LoadATGCWeights() {


  map<vector<float>, vector<float> > tmp2;
  vector<map<vector<float>, vector<float> > > mapWeight(100,tmp2);
  
 ifstream file("/home/mmarionn/Documents/CMS/ZZAnalysis/ZZaTGCWeight/aTGCDB", ios::in);

 if(file) {
   
    map<vector<float>, vector<float> >::const_iterator iter;
    
    float f4,f5,w,we;
    string s, stmp;
    vector<float> tmp(2,0);
    vector<float> vtmp(2,0);
  
    while(!file.eof()) {

      file >> f4 >> f5 >> s;
      
      tmp[0] = f4;
      tmp[1] = f5;

      if(file.eof()) break;

      for(int i=0;i<100;i++) {
	
	file >> w >> we >> stmp;
	iter = mapWeight[i].find( tmp );
	if( iter == mapWeight[i].end() ) {
	  vtmp[0] = w;
	  vtmp[1] = we;
	  //cout<<f4<<"   "<<f5<<"   "<<i<<"   "<<w<<endl;
	  if(w==0)
	    w=1;
	  
	  mapWeight[i][ tmp ] = vtmp;
	}
	else 
	  cout<<" soucis, "<<f4<<"   "<<f5<<"  en double "<<endl;
	
      }
      
    }
  }
 else { cout<<" grou "<<endl; abort(); }
  aTGCWeights = mapWeight;

}


vector<TF2*>
DibosonUseTree::LoadDBaTGC() {

  ifstream ifile("/home/mmarionn/Documents/CMS/CMSDatabase/paraboloiddatabase" , ios::in );

  vector<TF2*> curves;

  if(ifile) {

    string tmp;
    float bi,x0,y0,a,b,c;
    //header
    ifile >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp;
    
    while(!ifile.eof()) {
      
      ifile >> bi >> x0 >> y0 >> a >> b >> c;
      
      if(ifile.eof()) break;
	
      ostringstream os;
      os << b;

      TF2* tmp = new TF2( ("f"+os.str()).c_str(),"1 + (x - [0])*(x- [0])*([2]) + (y - [1])*(y- [1])*([3]) + 2*(y - [1])*(x- [0])*([4])",-0.3,0.3,-0.3,0.3);

      tmp->FixParameter(0,x0);
      tmp->FixParameter(1,y0);
      tmp->FixParameter(2,a);
      tmp->FixParameter(3,b);
      tmp->FixParameter(4,c);

      curves.push_back( tmp );

    }
     

  }
  else {cout<<" no aTGC database !!!!!!!!!!!!!!"<<endl; abort();}

  return curves;

}


float 
DibosonUseTree::GetATGCWeight(float f4, float f5, int ibin) {

  map<vector<float>, vector<float> >::const_iterator iter;

  vector<float> v;
  v.push_back(f4);
  v.push_back(f5);

  iter = aTGCWeights[ ibin ].find(v);
  if( iter == aTGCWeights[ ibin ].end() )
    cout<<" youpla "<<endl;

  // cout<<(*iter)<<endl;
  return ((*iter).second)[0];

}


bool
DibosonUseTree::BasicAcc(float pt, float eta, int pdgId ) {

  if(pt < 20 ) return false; 
  if(abs(pdgId)==11) 
    {
      if( fabs(eta)>2.5 || (fabs(eta)>1.44 && fabs(eta)<1.56) ) return false;
    }
  else {
    if( fabs(eta)>2.4) return false;
  }

  
  return true;
}

bool 
DibosonUseTree::isInAcc(float pt1, float pt2, float mZ,
			const vector<double>* LPt, const vector<double>* LEta,
			const vector<double>* LPdg, int offset) {

  //FIXME
   if( !BasicAcc( (*LPt)[0+offset], (*LEta)[0+offset] , (*LPdg)[0+offset] ) ||
       !BasicAcc( (*LPt)[1+offset], (*LEta)[1+offset] , (*LPdg)[1+offset] ) ) return false;
   
  if( mZ > 120 || mZ < 60 ) return false;  
  if( pt1 < 30 ) return false;
  
  //  if( pt2 < 30 ) return false;

  return true;
}



float
DibosonUseTree::METCorrectionData(float met, int nVertex) {
  float alpha_0  = 6.54;
  float alpha_PU = 0.78;
  float sigma_0  = 4.50;
  float sigma_PU = 2.37;
  
  float sigma_T  = sqrt( pow(sigma_0,2) + (nVertex-1)*pow(sigma_PU,2) );
  float peak_T   = alpha_0 + (nVertex-1)*alpha_PU;
  
  float metSig = ( met - peak_T ) / sigma_T;
  return sigma_0 * metSig + alpha_0;
}

float
DibosonUseTree::METCorrectionMC(float met, int nVertex) {

  float   alpha_data_0  = 6.54;
  // float   alpha_data_PU = 0.78;
  float   sigma_data_0  = 4.50;
  //  float   sigma_data_PU = 2.37;

  float   alpha_0  = 4.82;
  float   alpha_PU = 0.87;
  float   sigma_0  = 3.08;
  float   sigma_PU = 2.33;


  float sigma_T  = sqrt( pow(sigma_0,2) + (nVertex-1)*pow(sigma_PU,2) );
  float peak_T   = alpha_0 + (nVertex-1)*alpha_PU;

  float metSig = ( met - peak_T ) / sigma_T;
  return sigma_data_0 * metSig + alpha_data_0;
  
}


float 
DibosonUseTree::METCorrectionMin(float min, int nVertex, bool isMC, int pdg) {

  // float  alpha_data_0  = 4.946;
  // float  alpha_data_PU = 0.240;
  // float  sigma_data_0  = 3.896;
  // float  sigma_data_PU = 1.103;
 
  float  alpha_data_0  = 4.997;
  float  alpha_data_PU = 0.234;
  float  sigma_data_0  = 3.962;
  float  sigma_data_PU = 1.069;
  if(pdg == 13 ) {
    alpha_data_0  = 4.900;
    alpha_data_PU = 0.245;
    sigma_data_0  = 3.773;
    sigma_data_PU = 1.170;
  }

  float  alpha_0 = alpha_data_0;
  float  alpha_PU = alpha_data_PU;
  float  sigma_0 = sigma_data_0;
  float  sigma_PU = sigma_data_PU;


  if(isMC) {
    // alpha_0  = 3.797;
    // alpha_PU = 0.310;
    // sigma_0  = 2.905;
    // sigma_PU = 1.275;
    
    alpha_0  = 3.807;
    alpha_PU = 0.305;
    sigma_0  = 2.924;
    sigma_PU = 1.267;

    if(pdg == 13 ) {
      alpha_0  = 3.776;
      alpha_PU = 0.315;
      sigma_0  = 2.908;
      sigma_PU = 1.262;
    }

  }





  float sigma_T  = sqrt( pow(sigma_0,2) +
			 (nVertex-1)*pow(sigma_PU,2) );
  float peak_T   = alpha_0 + (nVertex-1)*alpha_PU;

  //  theMET = cpuMETVtxMin;  // min de MET et cpuMETVtxMin

  float ghmSig = ( min - peak_T ) / sigma_T;
  return sigma_data_0 * ghmSig + alpha_data_0;

}


void
DibosonUseTree::AddTreeLoading( string type, vector<string> AddTNames ) {

  vector<TChain*> tmpch;
  
 //Preparation des TChains
   for(int i=0;i<nt+1;i++) {
     TChain* ctmp = new TChain(AddTNames[0].c_str());
     tmpch.push_back(ctmp);
   }
   
   cout<<" debut lecture "<<reposi<<endl;
   TFile* datafile;
   if(!NoData) {
     for(size_t i=0;i<data.size();i++) {
       if(data[i]!="") {
	 string NameF = "../root/"+reposi+"/"+data[i]+".root"; 
	 datafile = new TFile(NameF.c_str(), "READ");
	 if(datafile==NULL) { cout<<" No such file "<<histoManager.Name<<endl; return;}
	 if(AddTNames.size()!=1) {
	   for(size_t it=0; it<AddTNames.size();it++) {
	     TTree* tmptree = (TTree*)datafile->Get( (AddTNames[it]).c_str() );
	     if(tmptree != NULL ) {
	       tmpch[nt]->Add( (NameF+"/"+AddTNames[it]).c_str()); }
	     delete tmptree;
	   }
	 }
	 else{
	   tmpch[nt]->Add(NameF.c_str());
	 }
	 
	 datafile->Close();
       }
     }
   }
   
  int nTmp=0; int ndt=0;
  for(size_t i=0;i<datasets.size();i++) {
    string NameF = "../root/"+reposi+"/"+datasets[i]+".root";

    if(ndt!=dt[i])
      nTmp=0;
    ndt=dt[i];
    
    datafile = new TFile(NameF.c_str(), "READ");
    if(AddTNames.size()!=1) {
      for(size_t it=0; it<AddTNames.size();it++) {
	TTree* tmptree = (TTree*)datafile->Get( (AddTNames[it]).c_str() );
	if(tmptree != NULL ) {
	  tmpch[ dt[i] ]->Add( (NameF+"/"+AddTNames[it]).c_str()); }
	delete tmptree;
      }
    }
    else{
      tmpch[ dt[i] ]->Add(NameF.c_str());
    }
    datafile->Close();
  }

  AddChains[ type ] = tmpch;
  cout<<" Additionnal tree loaded : "<<type<<endl;

}



void 
DibosonUseTree::SetAdditionnalTrees(string s, vector<string> vs ) {
  AddTreeLoading( s, vs );
}


TTree*
DibosonUseTree::GetAdditionnalTree(string type, int ds ) {

  map<string, vector<TChain*> >::const_iterator iter;

  iter = AddChains.find(type);

  return (TTree*)(((*iter).second)[ds]);

}



bool 
DibosonUseTree::isMatch(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) {

  
  if( dR(eta1,eta2,phi1,phi2) < 0.2 && 
      ( (fabs(pt1)-fabs(pt2))/fabs(pt2) >0 ) && (fabs(pt1)-fabs(pt2))/fabs(pt2) <0.6 )
    return true;
  else
    return false;

}

// bool 
// DibosonUseTree::isInAcc(const vector<float>* gLPt, const vector<float>* gLEta, 
// 			const vector<float>* gLPdg,
// 			float ptV1, float ptV2) {

 
//   // lP.SetPtEtaPhiE(lPart[0], lPart[1], lPart[2], lPart[3]);
 
//   if( !BasicAcc( (*gLPt)[0], (*gLEta)[0], (*gLPdg)[0] ) ||
//       !BasicAcc( (*gLPt)[1], (*gLEta)[1], (*gLPdg)[1] ) ) return false;
//   if( ptV1 < 30 ) return false;
//   if( ptV2 < 30 ) return false;
  
//   //  if( lP.Pt() > 30 ) return false;

//   return true;
// }



void
DibosonUseTree::ComputeNFakes() {

  if(FakeDB==NULL)
    LoadFakeRate();


  float Ndf=0;
  float Nf=0;
  
  float Ndferror=0;
  float Nferror=0;

  float Ndfp=0;
  float Nfp=0;
  float Ndfm=0;
  float Nfm=0;


  int nLLep = LooseElec.size();

  cout<<" Size ---> "<<nLLep<<endl;

  vector<float> fr1, fr2;
  int id1, id2;

  for(int i=0;i<nLLep;i++) {

    id1 = LooseElec[i][0][2];
    id2 = LooseElec[i][1][2];
    
    fr1 = getFRate( fabs(LooseElec[i][0][0]), fabs(LooseElec[i][0][1])  );
    fr2 = getFRate( fabs(LooseElec[i][1][0]), fabs(LooseElec[i][1][1])  );

    if( id1<7 && id2<7) {

      float coef1 = fr1[0]/(1-fr1[0]);
      float coef2 = fr2[0]/(1-fr2[0]);

      Ndf += coef1*coef2;
      Ndferror += coef1*coef2 *  coef1*coef2;

      float coef1m = (fr1[0]-fr1[1])/( 1-(fr1[0]-fr1[1]) );
      float coef1p = (fr1[0]+fr1[1])/( 1-(fr1[0]+fr1[1]) );
      float coef2m = (fr2[0]-fr2[1])/( 1-(fr2[0]-fr2[1]) );
      float coef2p = (fr2[0]+fr2[1])/( 1-(fr2[0]+fr2[1]) );
      
      Ndfp += coef1p*coef2p;
      Ndfm += coef1m*coef2m;
    }
    else if( (id1<7 && id2>=7 ) || (id1>=7 && id2<7) ) {

      float coef = fr1[0]/(1-fr1[0]);
      if(id2<7)
	coef = fr2[0]/(1-fr2[0]);

      cout<<(id2<7)<<"    "<<coef<<"    "
	  <<LooseElec[i][1][0]<<"   "
	  <<LooseElec[i][1][1]<<endl;

      Nf += coef;
      Nferror += coef*coef;

      float coefm = (fr1[0]-fr1[1])/( 1-(fr1[0]-fr1[1]) );
      float coefp = (fr1[0]+fr1[1])/( 1-(fr1[0]+fr1[1]) );
      if(id2<7) {
	coefm = (fr2[0]-fr2[1])/( 1-(fr2[0]-fr2[1]) );
	coefp = (fr2[0]+fr2[1])/( 1-(fr2[0]+fr2[1]) );
      }

      Nfp += coefp;
      Nfm += coefm;
    }
    
  }
  cout<< "Number of fakes    -> "<<Ndf<<"  "<<Nf
      <<" sum    "<<Ndf+Nf<<" +- "<<sqrt( Ndferror + Nferror)/*sqrt(Ndf+Nf)*/
      <<" + "<<fabs(Ndf-Ndfm) + fabs(Nf-Nfm) 
      <<" - "<<fabs(Ndf-Ndfp) + fabs(Nf-Nfp) <<endl;
  
}


void
DibosonUseTree::LoadFakeRate() {

  TFile* file= new TFile("/home/mmarionn/Documents/CMS/ZZAnalysis/FakeElectron/FakeDBTightEle.root");

   FakeDB = (TH1F*)file->Get("FR");
   
}

vector<float> 
DibosonUseTree::getFRate(float pt, float eta) {
  

  if(pt > 70)
    pt = 70;

  int binx = FakeDB->GetXaxis()->FindBin(eta);
  int biny = FakeDB->GetYaxis()->FindBin(pt);

  float fr =  FakeDB->GetBinContent( binx , biny );
  float fre =  FakeDB->GetBinError( binx , biny );
  
  vector<float> VFR;
  VFR.push_back(fr);
  VFR.push_back(fre);

  return VFR;

}



void 
Ellips(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
  //minimisation function computing the sum of squares of residuals
  Int_t np = contour->GetN();
  f = 0;
  //   cout<<" new "<<f<<endl;

  // Double_t *x = contour->GetX();
  // Double_t *y = contour->GetY();
  // for (Int_t i=0;i<np;i++) {
  //   Double_t u = (x[i] - par[0])*cos(par[4]) - (y[i] - par[1])*sin(par[4]);
  //   Double_t v = (y[i] - par[1])*cos(par[4]) + (x[i] - par[0])*sin(par[4]);
  //   Double_t r = sqrt(x[i]*x[i] + y[i]*y[i] );
  //   Double_t el = u*u/(par[2]*par[2]) + v*v/(par[3]*par[3]);
  //   Double_t dr = par[5] - sqrt(el);
  //   f += dr*dr;
  //   //  cout<<"\t"<<r<<"    "<<sqrt(el)<<"   "<<dr<<"   "<<f<<"   "<<par[2]<<"   "<<par[3]<<endl;
  // }
  // cout<<par[2]<<"   "<<par[3]<<"  "<<f<<endl;

  Double_t *x = contour->GetX();
  Double_t *y = contour->GetY();
  for (Int_t i=0;i<np;i++) {

    Double_t u = (x[i] - par[0])*cos(par[4]) - (y[i] - par[1])*sin(par[4]);
    Double_t v = (y[i] - par[1])*cos(par[4]) + (x[i] - par[0])*sin(par[4]);
    TVector2 p( u, v );
     
    double phi = p.Phi();
    double rP = p.Mod();

    double e = sqrt( par[2]*par[2] - par[3]*par[3])/par[2];
    double rE = sqrt( par[3]*par[3]/(1- e*e*pow( cos( phi ), 2 ) ) );

    double dR = rP-rE;

    double ex = contour->GetErrorX(i);
     
    f += (dR*dR)/(ex*ex);
  }
   
}



void
DibosonUseTree::DrawContour(float lim95, float lim68, float limU) {

  vector<vector<float> > ratios = GetReweightedYields(0.1,0.1,20,20);

  vector< std::pair<TVector2, float> > bV;
  vector< std::pair<TVector2, float> > bV68;
  //init
  for(int i=0;i<180;i++) {
    TVector2 t(0,0);
    bV.push_back( std::pair<TVector2, float>(t,0));
    bV68.push_back( std::pair<TVector2, float>(t,0));
  }
  

  for(int unsigned i=0;i<ratios.size(); i++) {

    TVector2 pos(ratios[i][0] , ratios[i][1]);
    float rat= ratios[i][2];
 
    for(int it=0;it<180;it+=1) {
      
      float p = it*2*3.1415/180;

      if( fabs( phi( pos.X(), pos.Y() ) -p) < 3*3.1415/180 ) {
	// cout<<p<<"    "<<ratios[i][0] <<"   "<<ratios[i][1]<<" -> "
	//     <<phi( pos.X(), pos.Y() )<<"   "<<rat<<"   "<<bV[it].first.Mod()
	//     <<"     "<<(fabs(lim95-rat) < fabs(lim95-bV[it].second ))<<endl;
	if( fabs(lim95-rat) < fabs(lim95-bV[it].second ) ) {
	  bV[it] = std::pair<TVector2, float>(pos, rat);
	} 
	break;
      } 
    }
  }

  for(int unsigned i=0;i<ratios.size(); i++) {

    TVector2 pos(ratios[i][0] , ratios[i][1]);
    float rat= ratios[i][2];
 
    for(int it=0;it<180;it+=1) {
      
      float p = it*2*3.1415/180;

      if( fabs( phi( pos.X(), pos.Y() ) -p) < 3*3.1415/180 ) {
	// cout<<p<<"    "<<ratios[i][0] <<"   "<<ratios[i][1]<<" -> "
	//     <<phi( pos.X(), pos.Y() )<<"   "<<rat<<"   "<<bV[it].first.Mod()
	//     <<"     "<<(fabs(lim95-rat) < fabs(lim95-bV[it].second ))<<endl;
	if( fabs(lim68-rat) < fabs(lim68-bV68[it].second ) ) {
	  bV68[it] = std::pair<TVector2, float>(pos, rat);
	} 
	break;
      } 
    }
  }

  cout<<" grou1 "<<endl;
  //End  pointing

  contour=new TGraphErrors(0);
  int n=0;
  for(int i=0;i<180;i++) {
    //cout<<bV[i].first.Mod()/lim95<<endl;
    if(bV68[i].second/lim68 > 0.9 && bV68[i].second/lim68 < 1.1) {
      contour->SetPoint(n, bV68[i].first.X(), bV68[i].first.Y());
      contour->SetPointError(n, fabs(bV68[i].second-lim68)/50, fabs(bV68[i].second-lim68)/50 );
      n++;
    }
  }
  cout<<" grou2 "<<endl;
 
  TCanvas* c2=new TCanvas("contour","contour",600,600);
  
  TEllipse *ell68;
  {
    //Fit a circle to the graph points
    TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 5);
    fitter->SetFCN(Ellips);
   
    fitter->SetParameter(0, "x0",       0, 0.1, 0,0);
    fitter->SetParameter(1, "y0",       0, 0.1, 0,0);
    fitter->SetParameter(2, "a",      0.05, 0.12, 0,0);
    fitter->SetParameter(3, "b",      0.05, 0.12, 0,0);
    fitter->SetParameter(4, "theta", 1.7, 0.1, 0,3.30);
    // fitter->FixParameter(0);
    // fitter->FixParameter(1);

    Double_t arglist[1] = {0};
    fitter->ExecuteCommand("MIGRAD", arglist, 0);
   
    //Getparameters
    //float a = fitter->GetParameter(2)/(1-fitter->GetParameter(3)*fitter->GetParameter(3));

    //Draw the circle on top of the points
    ell68 = new TEllipse(fitter->GetParameter(0),
			 fitter->GetParameter(1),
			 fitter->GetParameter(2),
			 fitter->GetParameter(3), 
			 0, 360,
			 fitter->GetParameter(4) );
    ell68->SetLineColor(kBlue+1);
    ell68->SetLineWidth(3);
    ell68->SetLineStyle(7);
    ell68->SetFillStyle(0);
   
 
  }
  cout<<" grou3 "<<endl;
  
  TEllipse *ell95;

  {
   
    contour=new TGraphErrors(0);
    int n=0;
    for(int i=0;i<180;i++) {
      //cout<<bV[i].first.Mod()/lim95<<endl;
      if(bV[i].second/lim95 > 0.9 && bV[i].second/lim95 < 1.1) {
	contour->SetPoint(n, bV[i].first.X(), bV[i].first.Y());
	contour->SetPointError(n, fabs(bV[i].second-lim95)/50, fabs(bV[i].second-lim95)/50 );
	n++;
      }
    }
    

    //Fit a circle to the graph points
    TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 5);
    fitter->SetFCN(Ellips);
   
    fitter->SetParameter(0, "x0",       0, 0.1, 0,0);
    fitter->SetParameter(1, "y0",       0, 0.1, 0,0);
    fitter->SetParameter(2, "a",      0.08, 0.12, 0,0);
    fitter->SetParameter(3, "b",      0.08, 0.12, 0,0);
    fitter->SetParameter(4, "theta", 1.7, 0.1, 0,3.30);
    // fitter->FixParameter(0);
    // fitter->FixParameter(1);

    Double_t arglist[1] = {0};
    fitter->ExecuteCommand("MIGRAD", arglist, 0);
   
    //Getparameters
    //float a = fitter->GetParameter(2)/(1-fitter->GetParameter(3)*fitter->GetParameter(3));

    //Draw the circle on top of the points
    ell95 = new TEllipse(fitter->GetParameter(0),
			 fitter->GetParameter(1),
			 fitter->GetParameter(2),
			 fitter->GetParameter(3), 
			 0, 360,
			 fitter->GetParameter(4) );
    ell95->SetLineColor(kRed);
    ell95->SetLineWidth(3);
    ell95->SetLineStyle(1);
    ell95->SetFillStyle(0);
   
 
  }
  cout<<" grou4 "<<endl;

  TH1F* bidon= new TH1F("bidu","bidi",2000,-1,1);
  bidon->SetBinContent(1,-100);

  // bidon->GetXaxis()->SetLimits(-0.10,0.10);
  // bidon->GetYaxis()->SetLimits(-0.10,0.10);
  bidon->GetXaxis()->SetTitle("f_{4}^{Z}   ");
  bidon->GetYaxis()->SetTitle("f_{5}^{Z}   ");
  bidon->GetXaxis()->SetTitleSize(0.05);
  bidon->GetYaxis()->SetTitleSize(0.05);
  bidon->GetXaxis()->SetNdivisions(9,5,0);
  bidon->GetYaxis()->SetNdivisions(9,5,0);
  bidon->GetXaxis()->SetTitleOffset(1.1);
  bidon->GetYaxis()->SetTitleOffset(1.35);
  bidon->Draw();
  contour->Draw("p");
  bidon->GetYaxis()->SetRangeUser(-0.105,0.105);
  bidon->GetXaxis()->SetRangeUser(-0.105,0.105);
  ell95->Draw("same");
  ell68->Draw("same");

  TLegend* leg = new TLegend(0.713087,0.770979,0.91443,0.928322);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.034965);
  leg->AddEntry(ell95," 95 % C.L.","l");
  leg->AddEntry(ell68," 68 % C.L.","l");
  leg->Draw("same");

  TLine* line= new TLine(0,-0.105,0,0.105);
  line->SetLineColor(1);
  line->Draw("same");

}


void
DibosonUseTree::DrawATGCYields() {


  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);


  const Int_t NRGBs = 5;
  const Int_t NCont = 256;
    
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);



  vector<vector<float> > ratios = GetReweightedYields(0.2,0.2,40,40);

  TH2F* atgc = new TH2F("atgc","atgc",41,-0.205,0.205,41,-0.205,0.205);

  for(int unsigned i=0;i<ratios.size(); i++) {
    
    // TVector2 pos(ratios[i][0] , ratios[i][1]);
    // float rat= ratios[i][2];
    
    atgc->Fill( ratios[i][0], ratios[i][1], ratios[i][2] );
    
  }


  atgc->GetXaxis()->SetTitle("f^{4}_{Z}");
  atgc->GetYaxis()->SetTitle("f^{5}_{Z}");
  atgc->GetXaxis()->SetNdivisions(9,5,0);
  atgc->GetXaxis()->SetNdivisions(9,5,0);

  TCanvas* c = new TCanvas("atgcY","atgcY",600,600);
  atgc->Draw("colz");


}


void
DibosonUseTree::DrawATGCComponent(float f4v, float f5v) {

  //First, produce the ATGC histo ========================
  int Obs = histoManager.FindNVar("ZPtMass");

  vector<float> temp = histoManager.GetTemplate(Obs);

  TH1F* atgcH=new TH1F("atgcH","atgcH",(int)temp[0], temp[1], temp[2] );
  int nbmzz;
  for(int unsigned ie=0;ie<selZZ.size();ie++) {
    nbmzz = (int)(selZZ[ie].first)/20;
    if(nbmzz >= 100 )
      { nbmzz=99; }

    atgcH->Fill(selZPt[ie],  GetATGCWeight(f4v, f5v, nbmzz)*selZZ[ie].second );
  }
  atgcH->SetLineColor(kAzure-3);
  atgcH->SetLineStyle(2);
  atgcH->SetLineWidth(2);
  atgcH->Rebin(Bin);

  //Now get the Pad
  ( (TCanvas*)(gROOT->FindObject("c2")))->cd();
  ( (TPad*)(gROOT->FindObject("pad_0")))->cd();


  //Get The Plots
  ostringstream val;
  val<<Obs;
  TH1F* ZZh = (TH1F*)(gROOT->FindObject( (name[0]+"_"+val.str()).c_str() ));
  TH1F* All = (TH1F*)(gROOT->FindObject( (name[nt-1]+"_"+val.str()).c_str() ));
  cout<<name[0]+val.str()<<"   "<<name[nt-1]+val.str()<<endl;
  cout<<ZZh<<"   "<<All<<endl;

  All->Add(ZZh,-1);
  atgcH->Add(All);
  
  atgcH->Draw("same hist");

  //THe legend
  TLegend* leg= (TLegend*)(gROOT->FindObject("TPave"));
 ostringstream os,os2;
 os << f4v;
 os2 << f5v;
 string l="aTGC : f_{4}^{Z} = "+os.str() + ", f_{5}^{Z} = "+os2.str() ;
 leg->AddEntry(atgcH, l.c_str() , "l");
 
 leg->Draw("same");
 
}
