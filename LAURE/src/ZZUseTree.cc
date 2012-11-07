#include "ZZUseTree.hh"

#include <iomanip>
#include <TLorentzVector.h>


TGraphErrors* contour=new TGraphErrors(0);;
//Epplis
void Ellips(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t);

using namespace std;

ClassImp(ZZUseTree)


ZZUseTree::ZZUseTree():
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

}

void
ZZUseTree::AddVariables(bool usepid, int pid,int mj,int mp,int ml, float dphijmc,
			float NMcutH, float NMcutL, float nnmcut, bool fj, bool invj, 
			bool btag, bool invmet, float nncut, bool useAFEff ) {

  UsePdgId = usepid;
  pdgID = pid;
  MaxJet = mj;
  MaxPhoton = mp;
  MaxLepton = ml;
  dPhiJMCut = dphijmc;
  NormMETCutHigh = NMcutH;
  NormMETCutLow = NMcutL;
  NNoutMETCut =  nnmcut;
  fixJet = fj;
  invJet = invj;
  invMet = invmet;
  bTagFlag = btag;
  NNoutCut = nncut;

  useAccForEff = useAFEff;

}

void
ZZUseTree::FillZZTree() {

  //Add Variables to be computed

  //Prepare suffix for outfiles
  string suffix = "Multivar";
  if(switchRMS && FillProf)
    suffix = "RMS_METs";

  if(Draw3on1)
    for(int i=0;i<5;i++)
      Suffix[i] = suffix;

  // string metT[5]={"pf","tc","calo","caloT1","caloT2"};
  string lep[2]={"l1","l2"};

  //Vertex
  histoManager.AddVariable("NVertex",10,0,10,"number of vertices","Vertex");
  histoManager.AddVariable("NVertex1",10,0,10,"number of vertices","Vertex");
  histoManager.AddVariable("NVertexControl",10,0,10,"number of vertices","Vertex");

  histoManager.Add2DVariable("SumETNVertexControl","#sum E_{T} [GeV]","N vertices","SumEtNVtxCtr",150,0,1500,10,0,10);

  //Lepton
  for(int ii=0;ii<2;ii++) {
    histoManager.AddVariable("Pt"+lep[ii],400,0,200,"p_{T} [GeV]","Pt_"+lep[ii]);
    histoManager.AddVariable("Eta"+lep[ii],100,-3,3,"#eta","Eta_"+lep[ii]);
    histoManager.AddVariable("Phi"+lep[ii],128,0,8,"#phi","Phi_"+lep[ii]);
    histoManager.AddVariable("P"+lep[ii],500,0,500,"p [GeV]","P_"+lep[ii]);
    histoManager.AddVariable("PdgId"+lep[ii],50,-25,25,"PdgId","PdgId_"+lep[ii]);
  }
  histoManager.AddVariable("Charge",4,-2,2,"#Pi charge","Charge");
  histoManager.AddVariable("ChargeEnd",4,-2,2,"#Pi charge","Charge");

  histoManager.AddVariable("TrackIsoMuNoRel",400,0,40,"Trk IsoMu","IsoMu");
  histoManager.AddVariable("TrackIsoMu",400,0,4,"Trk IsoMu","IsoMu");
  histoManager.AddVariable("EcalIsoMu",400,0,4,"Ecal IsoMu","IsoMu");
  histoManager.AddVariable("HcalIsoMu",400,0,4,"Hcal IsoMu","IsoMu");
  histoManager.AddVariable("TrackIso",400,0,4,"TrackIso ","TrackIso");
  histoManager.AddVariable("EcalIso",400,0,4,"EcalIso ","EcalIso");
  histoManager.AddVariable("HcalIso",400,0,4,"HcalIso ","HcalIso");
    


  histoManager.Add2DVariable("ChargeVsPt","Z q_{T} [GeV]","#Pi charge","METPCor_Vtx",10,0,200,3,-1,2);
      
  histoManager.AddVariable("ElecMult",10,0,10,"Electron multiplicity","ElecMult");
  histoManager.AddVariable("MuonMult",10,0,10,"Muon multiplicity","MuonMult");
  histoManager.AddVariable("LeptonMult",10,0,10,"Lepton multiplicity","LeptonMult");
  histoManager.AddVariable("LeptonMultEnd",10,0,10,"Lepton multiplicity","LeptonMult");

  histoManager.AddVariable("expInHit1",2,0,2,"expInHit1","expInHit1");
  histoManager.AddVariable("expInHit2",2,0,2,"expInHit2","expInHit2");
  histoManager.AddVariable("inHit1",2,0,2,"inHit1","inHit1");
  histoManager.AddVariable("inHit2",2,0,2,"inHit2","inHit2");
  histoManager.AddVariable("convRej1",2,0,2,"convRej1","convRej1");
  histoManager.AddVariable("convRej2",2,0,2,"convRej2","convRej2");

  histoManager.AddVariable("LElecPt",200,0,100,"p_{T} [GeV]","Pt_LElec");
  histoManager.AddVariable("LElecEta",100,-3,3,"#eta Le","Eta_LElec");
  histoManager.AddVariable("LElecPhi",128,-180,180,"#phi Le","Phi_LElec");
  histoManager.AddVariable("LMuonPt",200,0,100,"p_{T} [GeV]","Pt_LMuon");
  histoManager.AddVariable("LMuonEta",100,-3,3,"#eta Le","Eta_LMuon");
  histoManager.AddVariable("LMuonPhi",128,-180,180,"#phi Le","Phi_LMuon");

  histoManager.AddVariable("LMuonID",2,0,2,"IsIDMu","IsIDMu");
  histoManager.AddVariable("LMuonIso",2,0,2,"IsIsoMu","IsIsoMu");
  histoManager.AddVariable("LElecID",2,0,2,"IsIDElec","IsIDElec");
  histoManager.AddVariable("LElecIso",2,0,2,"IsIsoElec","IsIsoElec");

  histoManager.AddVariable("dRLMuonJet",110,-1,10,"dRLMuonJet","dRLMuonJet");
  histoManager.AddVariable("dRLElecJet",110,-1,10,"dRLElecJet","dRLElecJet");
  
  histoManager.AddVariable("RatioPt",400,0,4,"Ratio Pt lepton","RatioPt");

  //Z
  histoManager.AddVariable("ZMass",400,0,200,"M_{ll} [GeV]","Z_mass");
  histoManager.AddVariable("ZPt",1000,0,500,"p_{T} Z [GeV]","Z_pt");
  histoManager.AddVariable("ZEta",100,-3,3,"#eta Z","Z_eta");
  histoManager.AddVariable("ZPhi",128,-180,180,"#phi Z","Z_phi");
  histoManager.AddVariable("ZP",500,0,500,"p Z [GeV]","Z_p");
  histoManager.AddVariable("ZCateg",3,0,3,"Category Z","Z_categ");
  histoManager.AddVariable("ZMult",10,0,10,"Z multiplicity","Z_multiplicity");

  histoManager.AddVariable("ZY",100,-5,5,"Y_{Z}","Z_rapidity");

  histoManager.AddVariable("PfZMass",400,0,200,"PF M_{ll} [GeV]","pfZ_mass");
  
  histoManager.AddVariable("ZPtMN",1000,0,500,"p_{T} Z [GeV]","Z_pt");
  histoManager.AddVariable("ZPtNoMN",1000,0,500,"p_{T} Z [GeV]","Z_pt");
  
  histoManager.AddVariable("ZPtNoNN",1000,0,500,"p_{T} Z [GeV]","Z_pt");

  //Met
  histoManager.AddVariable("MET",400,0,200," #slash{E}_{T} [GeV]","MET_pt");
  histoManager.AddVariable("PhiMET",128,0,8,"#Phi #slash{E}_{T}","MET_phi");

  histoManager.AddVariable("METNoNN",400,0,200," #slash{E}_{T} [GeV]","MET_pt");

  histoManager.AddVariable("SpecMET",400,0,200,"Spec #slash{E}_{T} [GeV]","SpecMET_pt");
  histoManager.AddVariable("METNorm",200,0,2," #slash{E}_{T} /q_{T}","METNorm_pt");
  histoManager.AddVariable("TcMET",400,0,200," tc #slash{E}_{T} [GeV]","TcMET_pt");
  histoManager.AddVariable("METPrime",400,0,200," #slash{E}_{T} Prime [GeV]","METPrime_pt");
  histoManager.AddVariable("METPrimeT",400,-200,200," #slash{E}_{T} PrimeT [GeV]","METPrime_pt");
  histoManager.AddVariable("METPrimeL",400,-200,200," #slash{E}_{T} PrimeL [GeV]","METPrime_pt");
  histoManager.AddVariable("METPrimeCor",400,0,200," #slash{E}_{T} PrimeCor [GeV]","METPrime_pt");
  histoManager.AddVariable("METPrimeCorJet",400,0,200," #slash{E}_{T} PrimeCorJet [GeV]","METPrime_pt");
  histoManager.AddVariable("METPrimeCorNorm",500,0,10," #slash{E}_{T} Cor /q_{T}","METNorm_pt");
  histoManager.AddVariable("ProjMET",400,0,200," proj#slash{E}_{T} [GeV] ","METProj");
  histoManager.AddVariable("CorProjMET",400,0,200,"cor#slash{E}_{T} [GeV] ","METProj");
  histoManager.AddVariable("SigniCorProjMET",400,0,20," S_{#slash{E}_{T}}  ","METProj");
  histoManager.AddVariable("SigniCorProjMETNorm",400,0,10," S cor projected #slash{E}_{T}/q_{T}","METProj");
  histoManager.AddVariable("CorProjMETNorm",500,0,10,"B_{#slash{E}_{T}}  ","METNorm_pt");
  histoManager.AddVariable("ProjMETNorm",500,0,10," proj #slash{E}_{T} /q_{T}","METNorm_pt");

  histoManager.AddVariable("CorProjMETNormNoMC",500,0,10," cor proj #slash{E}_{T} /q_{T}","METNorm_pt");
  histoManager.AddVariable("CorProjMETNormMC",500,0,10," cor proj #slash{E}_{T} /q_{T}","METNorm_pt");

  for(int k=1;k<14;k++) {
    ostringstream os;
    os<<k;
    string n="PFMET_V"+os.str();
    histoManager.AddVariable(n,400,0,200," pf#slash{E}_{T}[GeV]","pfMET");
    string prjn="ProjMET_V"+os.str();
    histoManager.AddVariable(prjn,400,0,200," projected #slash{E}_{T}[GeV]","METProj");
    string cprjn="CorProjMET_V"+os.str();
    histoManager.AddVariable(cprjn,400,0,200," cor proj #slash{E}_{T}[GeV]","METProj");
  }

  histoManager.AddVariable("MetPrimeTrans",400,-200,200," #slash{E}_{T} Prime Trans [GeV]","METPrimeT_pt");
  histoManager.AddVariable("MetPrimeLong",400,-200,200," #slash{E}_{T} Prime Long [GeV]","METPrimeL_pt");
  histoManager.AddVariable("MetPrimeRecoilTrans",400,-200,200," #slash{E}_{T} Prime Recoil Trans [GeV]","METPrimeRecoilT_pt");
  histoManager.AddVariable("MetPrimeRecoilLong",400,-200,200," #slash{E}_{T} Prime Recoil Long [GeV]","METPrimeRecoilL_pt");
  histoManager.AddVariable("MetPrimeCorRecoilTrans",400,-200,200," #slash{E}_{T} Prime Cor Recoil Trans [GeV]","METPrimeRecoilT_pt");
  histoManager.AddVariable("MetPrimeCorRecoilLong",400,-200,200," #slash{E}_{T} Prime Cor Recoil Long [GeV]","METPrimeRecoilL_pt");
  histoManager.AddVariable("MetPrimeCorTrans",400,-5,5," #slash{E}_{T} PrimeCor Trans [GeV]","METPrimeCorT_pt");
  histoManager.AddVariable("MetPrimeCorLong",400,-5,5," #slash{E}_{T} PrimeCor Long [GeV]","METPrimeCorL_pt");
  histoManager.AddVariable("MetPrimeCorCaloTrans",400,-100,100," #slash{E}_{T} Prime CorCalo Trans [GeV]","METPrimeT_pt");
  histoManager.AddVariable("MetPrimeCorCaloLong",400,-100,100," #slash{E}_{T} Prime CorCalo Long [GeV]","METPrimeL_pt");
   
  histoManager.Add2DVariable("METPCorVsVtx","#slash{E}_{T} prime cor [GeV]","N vertices","METPCor_Vtx",100,0,200,10,0,10);
  histoManager.Add2DVariable("ProjMETVsVtx","proj #slash{E}_{T} [GeV]","N vertices","METPCor_Vtx",100,0,200,10,0,10);

  //Jet
  histoManager.AddVariable("JetMult",10,0,10,"Jet multiplicity","Jet_multiplicity");
  histoManager.AddVariable("TcJetMult",10,0,10,"TcJet multiplicity","TcJet_multiplicity");
  histoManager.AddVariable("LJetPt",400,0,200,"p_{T} leading jet [GeV]","LJet_pt");
  histoManager.AddVariable("LJetEta",200,-6,6,"#eta leading jet","LJet_eta");
  histoManager.AddVariable("LJetPhi",128,-180,180,"#phi leading jet","LJet_phi");
  histoManager.AddVariable("LJetFrac",50,0,1,"hadronic fraction leading jet","LJet_hadfrac");
  histoManager.AddVariable("LJetSoftMuon",200,0,1,"soft muon btag","LJet_SMuBTag");
  histoManager.AddVariable("LJetTrkCount",400,-100,100,"Tracking count bTag","LJet_TrkCountBtag");

  for(int im=0;im<10;im++) {
    ostringstream os;
    os << im;
    string n = "JetMult_Pt_"+os.str();
    histoManager.AddVariable(n,20,0,20," Jet Multiplicity","JetMult");
 
    string m = "BTag_Jet_"+os.str();
    histoManager.AddVariable(m,2,0,2,"is B tag","Btag");
  }

  //Photon
  histoManager.AddVariable("PhotonMult",10,0,10,"Photon multiplicity","Photon_multiplicity");
  histoManager.AddVariable("LPhotonPt",400,0,200,"p_{T} leading photon [GeV]","LPhoton_pt");
  histoManager.AddVariable("LPhotonEta",100,-3,3,"#eta leading photon","LPhoton_eta");
  histoManager.AddVariable("LPhotonPhi",128,-180,180,"#phi leading photon","LPhoton_phi");
  
  histoManager.AddVariable("dRPhotJet",110,-1,10,"dRPhotJet","dRPhotJet");
  histoManager.AddVariable("dRPhotElec",110,-1,10,"dRPhotElec","dRPhotElec");

  //Lepton related variables
  histoManager.AddVariable("dPhiLepton",128,-180,180,"#Delta#phi(l,l)","Lepton_dphi");
  histoManager.AddVariable("dPhiLeptonZ",128,-180,180,"#Delta#phi(Z,l_{1})","Lepton_Z_dPhi");
  histoManager.AddVariable("dEtaLepton",100,-3,3,"#Delta#eta(l,l)","Lepton_deta");
  histoManager.AddVariable("dPtLepton",200,0,200,"#Delta p_{T}(l,l) [GeV]","Lepton_dpt");
  histoManager.AddVariable("dPLepton",500,-250,250,"#Delta p(l,l) [GeV]","Lepton_dp");
  histoManager.AddVariable("dRLepton",100,0,10,"#DeltaR(l,l)","Lepton_dr");

  //Boosted Leptons variables
  histoManager.AddVariable("dPhiLeptonS",128,-180,180,"#Delta#Phi(l,l)^{*}","Lepton_dphi_star");
  histoManager.AddVariable("dEtaLeptonS",100,-3,3,"#Delta#eta(l,l)^{*}","Lepton_deta_star");
  histoManager.AddVariable("dPtLeptonS",200,0,200,"#Delta p_{T}(l,l)^{*} [GeV]","Lepton_dpt_star");
  histoManager.AddVariable("dPLeptonS",500,-250,250,"#Delta p(l,l)^{*} [GeV]","Lepton_dp_star");
  histoManager.AddVariable("dRLeptonS",100,0,10,"#Delta R(l,l)^{*}","Lepton_dr_star");
  histoManager.AddVariable("CosThetaM",50,-1,1,"cos(#theta^{*}_{l^{-}})","Lepton_cosM_star");
  histoManager.AddVariable("CosThetaP",50,-1,1,"cos(#theta^{*}_{l^{+}})","Lepton_cosP_star");
  histoManager.AddVariable("CosThetaCSM",50,-1,1,"cos(#theta^{*}_{CS})","Lepton_cosM_star");
  histoManager.AddVariable("CosThetaCSP",50,-1,1,"cos(#theta^{*}_{CS})","Lepton_cosP_star");

  histoManager.Add2DVariable("CSAvsY","cos(#theta^{*}_{CS})","|Y_{Z}|","CSAvsY",25,-1,1,25,0,2.5);

  //MET related variables
  histoManager.AddVariable("ZMETMass",500,0,500,"M_{T}_{(Z,#slash{E}_{T})} [GeV]","ZMET_mass");
  histoManager.AddVariable("ZMassMET",500,0,500,"M_{Z} + #slash{E}_{T} [GeV]","Z_MET");
  histoManager.AddVariable("ZPtMET",500,0,500,"p_{T} Z + #slash{E}_{T} [GeV]","Z_MET_pt");
  histoManager.AddVariable("ZPMET",500,0,500,"p Z + #slash{E}_{T} [GeV]","Z_MET_p");
  histoManager.AddVariable("dPhiMETZ",128,-180,180,"#Delta#phi(Z,#slash{E}_{T})","Z_MET_dphi");
  histoManager.AddVariable("dPhiMETL1",128,-180,180,"#Delta#phi(#slash{E}_{T},L1)","MET_L1_dphi");
  histoManager.AddVariable("dPhiMETL2",128,-180,180,"#Delta#phi(#slash{E}_{T},L2)","MET_L2_dphi");
  histoManager.AddVariable("METPara",400,-200,200,"#slash{E}_{T,||} [GeV]","MET_para");
  histoManager.AddVariable("METPerp",400,-200,200,"#slash{E}_{T,#perp}  [GeV]","MET_perp");
  histoManager.AddVariable("L1METPara",400,-200,200,"#slash{E}_{T,||} L1 [GeV]","MET_L1_para");
  histoManager.AddVariable("L1METPerp",400,-200,200,"#slash{E}_{T,#perp}  L1 [GeV]","MET_L1_perp");
  histoManager.AddVariable("L2METPara",400,-200,200,"#slash{E}_{T,||} L2 [GeV]","MET_L2_para");
  histoManager.AddVariable("L2METPerp",400,-200,200,"#slash{E}_{T,#perp}  L2 [GeV]","MET_L2_perp");
  histoManager.AddVariable("L1METParaNorm",400,-2,2,"#slash{E}_{T,||}^{l1} / q_{T}","MET_L1_para");
  histoManager.AddVariable("L1METPerpNorm",400,-2,2,"#slash{E}_{T,#perp}^{l1}  / q_{T}","MET_L1_perp");
  histoManager.AddVariable("L2METParaNorm",400,-2,2,"#slash{E}_{T,||}^{l2}  / q_{T}","MET_L2_para");
  histoManager.AddVariable("L2METPerpNorm",400,-2,2,"#slash{E}_{T,#perp}^{l2}  / q_{T}","MET_L2_perp");
  histoManager.AddVariable("Recoil",400,0,200,"u_{T} [GeV]","Recoil_pt");
  histoManager.AddVariable("PhiRecoil",128,-180,180,"#Phi u_{T}","Recoil_phi");
  histoManager.AddVariable("ZRecoil",400,0,200,"p_{T} Z + u_{T}","Z_Recoil");
  histoManager.AddVariable("RecoilPara",400,-200,200,"u_{||} [GeV]","Recoil_para");
  histoManager.AddVariable("RecoilPerp",400,-200,200,"u_{#perp}  [GeV]","Recoil_perp");
  histoManager.AddVariable("L1RecoilPara",400,-200,200,"u_{||} L1 [GeV]","Recoil_L1_para");
  histoManager.AddVariable("L1RecoilPerp",400,-200,200,"u_{#perp} L1  [GeV]","Recoil_L1_perp");
  histoManager.AddVariable("L2RecoilPara",400,-200,200,"u_{||} L2 [GeV]","Recoil_L2_para");
  histoManager.AddVariable("L2RecoilPerp",400,-200,200,"u_{#perp} L2  [GeV]","Recoil_L2_perp");
  histoManager.AddVariable("dPhiRecoilZ",128,-180,180,"#Delta#Phi (Z,u_{T}) [GeV]","Z_recoil_dphi");
  histoManager.AddVariable("dPhiRecoilL1",128,-180,180,"#Delta#Phi (L1,u_{T}) [GeV]","L1_recoil_dphi");
  histoManager.AddVariable("dPhiRecoilL2",128,-180,180,"#Delta#Phi (L2,u_{T}) [GeV]","L2_recoil_dphi");
  histoManager.AddVariable("ZPtmMet",200,-100,100,"p_{T} Z - #slash{E}_{T} [GeV]","Z_ptmMET");

  histoManager.AddVariable("dPhiRecoilMET",180,0,180,"#Delta #Phi (u_{T},#slash{E}_{T}) [°]","dPhi_RMET");
  histoManager.AddVariable("HT",500,0,250,"H_{T} [GeV]","HT");

  histoManager.AddVariable("SumET",500,0,1000,"#sum E_{T} [GeV]","sumEt");

  //Resolution related variables
  histoManager.AddVariable("ResoPara",200,-10,10," RPara [#sigma( u_{||} )]","ResoPara");
  histoManager.AddVariable("ResoPerp",200,-10,10," RPerp [#sigma( u_{#perp}  )]","ResoPerp");
  
  //Jet related variables
  histoManager.AddVariable("dPhiJetMET",128,-180,180,"#Delta#Phi (jet,#slash{E}_{T}) [#circ] ","Jet_MET_dphi");
  histoManager.AddVariable("dPhiJetZ",128,-180,180,"#Delta#Phi (Z,jet)","Z_jet_dphi");
  histoManager.AddVariable("dEtaJetZ",100,-3,3,"#Delta#eta (Z,jet)","Z_jet_deta");
  histoManager.AddVariable("ZPtJet",500,0,500,"p_{T} Z + jet [GeV]","Z_jet_pt");
  histoManager.AddVariable("ZJetMass",500,0,500,"M_{(Z,jet)} [GeV]","Z_jet_mass");
  histoManager.AddVariable("dPhiJetL1",128,-180,180,"Delta#Phi (L1,jet)","L1_jet_dphi");
  histoManager.AddVariable("dPhiJetL2",128,-180,180,"Delta#Phi (L2,jet)","L2_jet_dphi");
  histoManager.AddVariable("dEtaJetL1",100,-3,3,"Delta#eta (L1,jet)","L1_jet_deta");
  histoManager.AddVariable("dEtaJetL2",100,-3,3,"Delta#eta (L2,jet)","L2_jet_deta");
      
  //Photon related variables
  histoManager.AddVariable("dPhiPhotonMET",128,-180,180,"#Delta#Phi (#gamma,#slash{E}_{T})","Photon_MET_dphi");
  histoManager.AddVariable("dPhiPhotonZ",128,-180,180,"#Delta#Phi (Z,#gamma)","Z_photon_dphi");
  histoManager.AddVariable("dEtaPhotonZ",100,-3,3,"#Delta#eta (Z,#gamma)","Z_photon_deta");
  histoManager.AddVariable("ZPtPhoton",500,0,500,"p_{T} Z + #gamma [GeV]","Z_photon_pt");
  histoManager.AddVariable("ZPhotonMass",500,0,500,"M_{(Z,#gamma)} [GeV]","Z_photon_mass");
  histoManager.AddVariable("dPhiPhotonL1",128,-180,180,"Delta#Phi (L1,#gamma)","L1_photon_dphi");
  histoManager.AddVariable("dPhiPhotonL2",128,-180,180,"Delta#Phi (L2,#gamma)","L2_photon_dphi");
  histoManager.AddVariable("dEtaPhotonL1",100,-3,3,"Delta#eta (L1,#gamma)","L1_photon_deta");
  histoManager.AddVariable("dEtaPhotonL2",100,-3,3,"Delta#eta (L2,#gamma)","L2_photon_deta");

  //NeuralNet output
  histoManager.AddVariable("NNout",100,0,1,"NN output","NNoutput");
  histoManager.AddVariable("NNoutMET",120,0,1.2,"#slash{E}_{T} NN output","NNoutputMET");
  histoManager.AddVariable("NNoutMETControl",120,0,1.2,"#slash{E}_{T} NN output","NNoutputMET");
  //Boosted variables


  //Control variables after selection
  histoManager.AddVariable("dPhiJetMETEnd",128,-180,180,"#Delta#Phi (jet,#slash{E}_{T}) [#circ]","Jet_MET_dphi");
  histoManager.AddVariable("NNoutMETEnd",120,0,1.2,"#slash{E}_{T} NN output","NNoutputMET");
 histoManager.AddVariable("METEnd",400,0,200," #slash{E}_{T} [GeV]","MET_pt");
 histoManager.AddVariable("ProjMETEnd",400,0,200," projected #slash{E}_{T}[GeV]","METProj");
  histoManager.AddVariable("CorProjMETEnd",400,0,200,"cor#slash{E}_{T}} [GeV]   ","METProj");
  histoManager.AddVariable("SigniCorProjMETEnd",400,0,20," S_{#slash{E}_{T}}   ","METProj");
  histoManager.AddVariable("SigniCorProjMETNormEnd",400,0,10," S_{#slash{E}_{T}}/q_{T}","METProj");
  histoManager.AddVariable("CorProjMETNormEnd",500,0,10," B_{#slash{E}_{T}}   ","METNorm_pt");
  histoManager.AddVariable("JetMultEnd",10,0,10,"Jet multiplicity","Jet_multiplicity");
  histoManager.AddVariable("NVertexEnd",10,0,10,"number of vertices","Vertex");

  //MonteCarlo Truth
  histoManager.AddVariable("McZZMass",400,0,2000,"#sqrt{#hat{s}} [GeV]","shat");

  //ZJet Estimation variable
  double binsX[4]={40,70,110,200};
  double binsY[5]={1.5,2.,3.,4.25,10};
  histoManager.Add2DVariable("ZMassVsSCPMET","M_{ll} [GeV]","Signi CPMET","ZMassVsSCPMET",3,binsX,4,binsY);

  double binsX2[3]={-1000,-30,1000};
  double binsY2[5]={1.5,2.,3.,4.25,1000};
  histoManager.Add2DVariable("ZMassVsUpara","u_{||} [GeV]","Signi CPMET","ZMassVsSCPMET",2,binsX2,4,binsY2);


  //End Declaration, now preparation
  histoManager.PrepareHistos(name);
  histoManager.Prepare2DHistos(name);
 

  vector<vector<vector<float> > > Wghts;
  if(ShapeWeight)
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

  //string metT[5]={"pf","tc","calo","caloT1","caloT2"};

  //Create Neural Net
  IClassifierReader* NN = CreateNN();
  vector<double> inputVec(7,0);

  //Create METNN
  IClassifierReader* MetNN = CreateMETNN();
  vector<double> inputVecMET(9,0);

  //sélection
  for(int i=0;i<nt+1;i++) {
   
    //Déclaration des variables

    //Leptons
    float IsoVar[2][3];
    float IDVar[2][4];
    float Lepton[2][4];
    int Charge[2];
    int PdgId[2];
    float EOP[2];
    
    bool convRej1;
    bool inHit1;
    bool expInHit1;
    bool convRej2;
    bool inHit2;
    bool expInHit2;
    
    int ElecMult;
    int MuonMult;
    float LeadAddElec[3];
    float LeadAddMuon[3];

    bool IsIDMu;
    bool IsIsoMu;
    bool IsIDElec;
    bool IsIsoElec;

    //Z
    float Z[4];
    float ZMass;
    int ZCateg;
    int ZMultiplicity;

    float PfZMass;
    float PfZ[4];
    
    //MET
    float Met[2];
    float Recoil[2];

    float SpecMET;
    float TrkMet[2];

    float METPrime[2];
    float METPrimeCor[2];
    float METPrimeRecoil[2];
    
    float METPrimeMETM[2];
    float METPrimeJetM[2];
    float METPrimeJetP[2];
    float METPrimeUnclusM[2];
    float METPrimeUnclusP[2];

    //Jet
    int JetMultiplicity;
    int TcJetMultiplicity;
    float LJet[4];
    float bTag[2];
    
    int MultiJetMultiplicity[10];

    //Photon
    int PhotonMultiplicity;
    float LPhoton[5];

     //Lepton related variables
    float dPhiLepton;
    float dPhiLeptonZ;
    float dEtaLepton;
    float dPtLepton;
    float dPLepton;
    float dRLepton;

    //Boosted leptons variables
    float dPhiLeptonS;
    float dEtaLeptonS;
    float dPtLeptonS;
    float dPLeptonS;
    float dRLeptonS;
    float CosThetaM;
    float CosThetaP;
    float CosThetaCSM;
    float CosThetaCSP;

    //MET related variables
    float ZMETMass;
    float ZMassMET;
    float ZPtMET;
    float ZPMET;
    float dPhiMETZ;
    float dPhiMETL1;
    float dPhiMETL2;
    float METProjZ[2];
    float METProjL1[2];
    float METProjL2[2];
    float ZRecoil;
    float RecoilProjZ[2];
    float RecoilProjL1[2];
    float RecoilProjL2[2];
    float dPhiRecoilZ;
    float dPhiRecoilL1;
    float dPhiRecoilL2;
    float ProjMET;

    float SumET;

    float HT[2];

    float RValMC[2];
    float RValData[2];

    //Jet related variables
    float dPhiJetMET;
    float dPhiJetZ;
    float dEtaJetZ;
    float ZPtJet;
    float ZJetMass;
    float dPhiJetL1;
    float dPhiJetL2;
    float dEtaJetL1;
    float dEtaJetL2;
    
    //Photon related variables
    float dPhiPhotonMET;
    float dPhiPhotonZ;
    float dEtaPhotonZ;
    float ZPtPhoton;
    float ZPhotonMass;
    float dPhiPhotonL1;
    float dPhiPhotonL2;
    float dEtaPhotonL1;
    float dEtaPhotonL2;
    
    //MonteCarlo truth
    float McLep[2][5];
    float McZl[5];
    float McZn[5];
    float McLP[5];

    //Vertices
    int NVertex;

    //Miscellaneous
    int Event;
    int AbsEvent;
    int Run;
    char sample[50];
    char fileName[50];

    float Weight=0;

    bool Sel[2]={true,true};

    tChains[i]->SetBranchAddress("NVertex",&NVertex);
    
    // Leptons 
    tChains[i]->SetBranchAddress("IDVar",IDVar);
    tChains[i]->SetBranchAddress("IsoVar",IsoVar);
    tChains[i]->SetBranchAddress("Lepton",Lepton);
    tChains[i]->SetBranchAddress("Charge",Charge);
    tChains[i]->SetBranchAddress("PdgId",PdgId);
    tChains[i]->SetBranchAddress("EOP",EOP);

    tChains[i]->SetBranchAddress("InHit1",&inHit1);
    tChains[i]->SetBranchAddress("ExpInHit1",&expInHit1);
    tChains[i]->SetBranchAddress("ConvRej1",&convRej1);
    tChains[i]->SetBranchAddress("InHit2",&inHit2);
    tChains[i]->SetBranchAddress("ExpInHit2",&expInHit2);
    tChains[i]->SetBranchAddress("ConvRej2",&convRej2);

    tChains[i]->SetBranchAddress("ElecMult",&ElecMult);
    tChains[i]->SetBranchAddress("MuonMult",&MuonMult);

    tChains[i]->SetBranchAddress("LeadAddElec",LeadAddElec);
    tChains[i]->SetBranchAddress("LeadAddMuon",LeadAddMuon);

    tChains[i]->SetBranchAddress("IsIDMu",&IsIDMu);
    tChains[i]->SetBranchAddress("IsIsoMu",&IsIsoMu);

    tChains[i]->SetBranchAddress("IsIDElec",&IsIDElec);
    tChains[i]->SetBranchAddress("IsIsoElec",&IsIsoElec);

  //Z
    tChains[i]->SetBranchAddress("Z",Z);
    tChains[i]->SetBranchAddress("ZMass",&ZMass);
    tChains[i]->SetBranchAddress("ZCateg",&ZCateg);
    tChains[i]->SetBranchAddress("ZMultiplicity",&ZMultiplicity);

    tChains[i]->SetBranchAddress("PfZMass",&PfZMass);
    tChains[i]->SetBranchAddress("PfZ",PfZ);
  
  //MET
    tChains[i]->SetBranchAddress("MET",Met);
    tChains[i]->SetBranchAddress("Recoil", Recoil);

    tChains[i]->SetBranchAddress("SpecMET",&SpecMET);
    tChains[i]->SetBranchAddress("TrkMET",TrkMet);

    tChains[i]->SetBranchAddress("METPrime",METPrime);
    tChains[i]->SetBranchAddress("METPrimeCor",METPrimeCor);
    tChains[i]->SetBranchAddress("METPrimeRecoil",METPrimeRecoil);
    //  tChains[i]->SetBranchAddress("METPrimeCorRecoil",METPrimeCorRecoil);
    //  tChains[i]->SetBranchAddress("METPrimeCorCalo",METPrimeCorCalo);

    tChains[i]->SetBranchAddress("METPrimeMETM",METPrimeMETM);
    tChains[i]->SetBranchAddress("METPrimeJetM",METPrimeJetM);
    tChains[i]->SetBranchAddress("METPrimeJetP",METPrimeJetP);
    tChains[i]->SetBranchAddress("METPrimeUnclusM",METPrimeUnclusM);
    tChains[i]->SetBranchAddress("METPrimeUnclusP",METPrimeUnclusP);

    tChains[i]->SetBranchAddress("ProjMET",&ProjMET);
    
  //Jet
    tChains[i]->SetBranchAddress("JetMultiplicity",&JetMultiplicity);
    tChains[i]->SetBranchAddress("TcJetMultiplicity",&TcJetMultiplicity);
    tChains[i]->SetBranchAddress("LeadingJet",LJet);
    tChains[i]->SetBranchAddress("bTag",bTag);
    
    tChains[i]->SetBranchAddress("MultiJetMultiplicity",MultiJetMultiplicity);

  //Photon
    tChains[i]->SetBranchAddress("PhotonMultiplicity",&PhotonMultiplicity);
    tChains[i]->SetBranchAddress("LeadingPhoton",LPhoton);
  
    //Lepton related variables
    tChains[i]->SetBranchAddress("dPhiLepton",&dPhiLepton);
    tChains[i]->SetBranchAddress("dPhiLeptonZ",&dPhiLeptonZ);
    tChains[i]->SetBranchAddress("dEtaLepton",&dEtaLepton);
    tChains[i]->SetBranchAddress("dPtLepton",&dPtLepton);
    tChains[i]->SetBranchAddress("dPLepton",&dPLepton);
    tChains[i]->SetBranchAddress("dRLepton",&dRLepton);

    //Boosted leptons variables
    tChains[i]->SetBranchAddress("dPhiLeptonS",&dPhiLeptonS);
    tChains[i]->SetBranchAddress("dEtaLeptonS",&dEtaLeptonS);
    tChains[i]->SetBranchAddress("dPtLeptonS",&dPtLeptonS);
    tChains[i]->SetBranchAddress("dPLeptonS",&dPLeptonS);
    tChains[i]->SetBranchAddress("dRLeptonS",&dRLeptonS);
    tChains[i]->SetBranchAddress("CosThetaM",&CosThetaM);
    tChains[i]->SetBranchAddress("CosThetaP",&CosThetaP);
    tChains[i]->SetBranchAddress("CosThetaCSM",&CosThetaCSM);
    tChains[i]->SetBranchAddress("CosThetaCSP",&CosThetaCSP);


    //MET related variables
    tChains[i]->SetBranchAddress("ZMETMass",&ZMETMass);
    tChains[i]->SetBranchAddress("ZMassMET",&ZMassMET);
    tChains[i]->SetBranchAddress("ZPtMET",&ZPtMET);
    tChains[i]->SetBranchAddress("ZPMET",&ZPMET);
    tChains[i]->SetBranchAddress("dPhiMETZ",&dPhiMETZ);
    tChains[i]->SetBranchAddress("dPhiMETL1",&dPhiMETL1);
    tChains[i]->SetBranchAddress("dPhiMETL2",&dPhiMETL2);
    tChains[i]->SetBranchAddress("METProjZ",METProjZ);
    tChains[i]->SetBranchAddress("METProjL1",METProjL1);
    tChains[i]->SetBranchAddress("METProjL2",METProjL2);
    tChains[i]->SetBranchAddress("ZRecoil",&ZRecoil);
    tChains[i]->SetBranchAddress("RecoilProjZ",RecoilProjZ);
    tChains[i]->SetBranchAddress("RecoilProjL1",RecoilProjL1);
    tChains[i]->SetBranchAddress("RecoilProjL2",RecoilProjL2);
    tChains[i]->SetBranchAddress("dPhiRecoilZ",&dPhiRecoilZ);
    tChains[i]->SetBranchAddress("dPhiRecoilL1",&dPhiRecoilL1);
    tChains[i]->SetBranchAddress("dPhiRecoilL2",&dPhiRecoilL2);

    tChains[i]->SetBranchAddress("HT",HT);

    tChains[i]->SetBranchAddress("SumEt",&SumET);
    
    //Resolutino related variables
    tChains[i]->SetBranchAddress("RValMC", RValMC);
    tChains[i]->SetBranchAddress("RValData", RValData);

    //Jet related variables
    tChains[i]->SetBranchAddress("dPhiJetMET",&dPhiJetMET);
    tChains[i]->SetBranchAddress("dPhiJetZ",&dPhiJetZ);
    tChains[i]->SetBranchAddress("dEtaJetZ",&dEtaJetZ);
    tChains[i]->SetBranchAddress("ZPtJet",&ZPtJet);
    tChains[i]->SetBranchAddress("ZJetMass",&ZJetMass);
    tChains[i]->SetBranchAddress("dPhiJetL1",&dPhiJetL1);
    tChains[i]->SetBranchAddress("dPhiJetL2",&dPhiJetL2);
    tChains[i]->SetBranchAddress("dEtaJetL1",&dEtaJetL1);
    tChains[i]->SetBranchAddress("dEtaJetL2",&dEtaJetL2);
  
    //Photon related variables
    tChains[i]->SetBranchAddress("dPhiPhotonMET",&dPhiPhotonMET);
    tChains[i]->SetBranchAddress("dPhiPhotonZ",&dPhiPhotonZ);
    tChains[i]->SetBranchAddress("dEtaPhotonZ",&dEtaPhotonZ);
    tChains[i]->SetBranchAddress("ZPtPhoton",&ZPtPhoton);
    tChains[i]->SetBranchAddress("ZPhotonMass",&ZPhotonMass);
    tChains[i]->SetBranchAddress("dPhiPhotonL1",&dPhiPhotonL1);
    tChains[i]->SetBranchAddress("dPhiPhotonL2",&dPhiPhotonL2);
    tChains[i]->SetBranchAddress("dEtaPhotonL1",&dEtaPhotonL1);
    tChains[i]->SetBranchAddress("dEtaPhotonL2",&dEtaPhotonL2);
    
    //Miscellaneous
    tChains[i]->SetBranchAddress("Event",&Event);
    tChains[i]->SetBranchAddress("AbsEvent",&AbsEvent);
    tChains[i]->SetBranchAddress("Run",&Run);
    tChains[i]->SetBranchAddress("sample",&sample);
    tChains[i]->SetBranchAddress("fileName",&fileName);

    //Montecarlotruth
    if( name[i].substr(0,2) == "di" || name[i].substr(0,2) == "ZZ" ||
	name[i].substr(0,2) == "WZ" || name[i].substr(0,2) == "WW" ||
	name[i].substr(0,2) == "ZV" ) {
    
      tChains[i]->SetBranchAddress("McZl",McZl);
      tChains[i]->SetBranchAddress("McZn",McZn);
      tChains[i]->SetBranchAddress("McLP",McLP);
      tChains[i]->SetBranchAddress("McLepton",McLep);
    }    

    cout<<" beginning tree  " <<name[i]<<"   "<<tChains[i]->GetEntries()<<endl;
    //Boucle et sélection, remplissage histos
    int EP=0;
    int ent = tChains[i]->GetEntries();
    for(int ie=0;ie<ent;ie++) {
   

      if(ie%1000==0)
	fprintf(stdout, "Done : %d / %d\r", ie, ent );

      if(i!=nt)
	Weight = GetWeight(i, ie);
      else
	Weight = 1;
     
      tChains[i]->GetEntry(ie);      
      
      //Acceptance
      inAcc=true;
   
      if(name[i].substr(0,2)=="di" || name[i].substr(0,2)=="ZZ" ||
	 name[i].substr(0,2)=="ZV" || name[i].substr(0,2)=="WZ" ||
	 name[i].substr(0,2)=="WW") {
	//	cout<<Weight<<"   ";
	string proc = FindProcess( ie, i );

	if( proc.substr(0,2)=="ZZ") {
	  Weight *= SearchWeightZZ(McZl[0]); 
	  //Acceptance
	  inAcc = isInAcc( McLep, McZn, McLP);
	  // cout<<McZl[0]<<"   "<<McZn[0]<<"   "<<inAcc<<endl;
	}
	if( proc.substr(0,2)=="WZ") {
	  Weight *= SearchWeightWZ(McZl[0]);
	}
	if( proc.substr(0,2)=="WW") {
	  Weight *= SearchWeightWW(McZn[0]); //based on W+
	}
	//	cout<<Weight<<endl;
      }
      // cout<<inAcc<<endl;
     

      if(name[i].substr(0,4)=="data")
	{
	  bool doubleCount=false;
	  for(size_t ik=0;ik<Events.size();ik++)
	    {
	      if( ((Events[ik]).first.first) == Run && ((Events[ik]).first.second) == AbsEvent)
	      	{doubleCount=true; break;}
	    }
	  if(doubleCount || (EventFilter && Run>EventNum ) )
	    { continue; }
	}

      if(NVertex!=NVert && NVert!=0 ) continue; //NVertex
      // if(i<nt-1 &&  NVert!=0 && pdgID==13) {Weight *= 0.29;}
      // if(i<nt-1 &&  NVert!=0 && pdgID==11) {Weight *= 0.29;}

    
	
      for(int id=0;id<2;id++) {
	Sel[id]=true;

	if(abs(PdgId[id])==13) {
	  histoManager.fill("TrackIsoMuNoRel",i,IsoVar[id][0],Weight);
	  histoManager.fill("TrackIsoMu",i,IsoVar[id][0]/Lepton[id][0],Weight);
	  histoManager.fill("EcalIsoMu",i,IsoVar[id][1]/Lepton[id][0],Weight);
	  histoManager.fill("HcalIsoMu",i,IsoVar[id][2]/Lepton[id][0],Weight);
	}

	else {
	  histoManager.fill("TrackIso",i,IsoVar[id][0]/Lepton[id][0],Weight);
	  histoManager.fill("EcalIso",i,IsoVar[id][1]/Lepton[id][0],Weight);
	  histoManager.fill("HcalIso",i,IsoVar[id][2]/Lepton[id][0],Weight);
	}
      }	
      for(int id=0;id<2;id++) {
	if(abs(PdgId[id])==11) {

	  if(fabs(Lepton[id][1])>1.479)
	    {EP=1; }
	  else
	    {EP=0; }
	  
	  if(EcalP!=EP && EcalP!=2 && EcalP!=3) Sel[id]=false;

	  if(IDVar[id][0] > _IDcuts[EP][0] && Cbool(skipCut, "sigieie"==observable.substr(0,7)) )  Sel[id]=false;
	  if( fabs(IDVar[id][1]) > _IDcuts[EP][1] && Cbool(skipCut, "Deta"==observable.substr(0,4)) )  Sel[id]=false; 
	  if( fabs(IDVar[id][2]) > _IDcuts[EP][2] && Cbool(skipCut, "Dphi"==observable.substr(0,4)) )  Sel[id]=false;
	  if(IDVar[id][3] > _IDcuts[EP][3] && Cbool(skipCut, "HoE"==observable.substr(0,3)) )  Sel[id]=false;
	  if(IsoVar[id][0]/Lepton[id][0] > _Isocuts[EP][0] && Cbool(skipCut, "TrackIso"==observable.substr(0,8)) )  Sel[id]=false;
	  if(IsoVar[id][1]/Lepton[id][0] > _Isocuts[EP][1] && Cbool(skipCut, "EcalIso"==observable.substr(0,7)) )  Sel[id]=false;

	  if(IsoVar[id][2]/Lepton[id][0] > _Isocuts[EP][2] && Cbool(skipCut, "HcalIso"==observable.substr(0,7)) )  Sel[id]=false;
	  
	  if(!Sel[id]) break;
	}
	else {
	  if( IsoVar[id][0] > 3 && Cbool(skipCut, "isoSum"==observable.substr(0,7)) )  Sel[id]=false;
	}
	  
      }
     

      if(!Sel[0] || !Sel[1] ) continue;

       //EOP cut
      //      if( (EOP[0]<0.2 && PdgId[0]) || (EOP[1]<0.2 && PdgId[1]) ) continue;

      //REweighting with vertex
       if(ShapeWeight && i==nt-1) {
	  Weight = Weight*GetMCWeightFromDataShape(0, NVertex,
						  "NVertexControl", 
						   "NVertexControl",1, i);
	  /*	 Weight = Weight*GetMCWeightFromDataShape2D( SumET, NVertex ,
						     "SumETNVertexControl" ,
						     "SumETNVertexControl");*/
       }
      //Filling Histos******************
    
       if(UsePdgId && ( abs(PdgId[0]) != pdgID || abs(PdgId[1]) != pdgID ) ) continue; 

       //Define CorProjMET, Gaussian Def ===============================/**0.2041/*/
       /*   float CorProjMET =  (ProjMET - 0.9072*(NVertex-1))/sqrt( 1 + (NVertex-1)*0.1774  );
       if(i==nt) {
	 CorProjMET =  (ProjMET - 0.9072*(NVertex-1))/sqrt( 1 + (NVertex-1)*0.295  ); //0.2261
       }
       float SigniCorProjMET = CorProjMET/7.36; //7.36
       */
       float CorProjMET =  (ProjMET - 0.8555*(NVertex-1))/sqrt( 1 + (NVertex-1)*0.178  ); //3.1
       if(i==nt) {
	 CorProjMET =  (ProjMET - 0.8555*(NVertex-1))/sqrt( 1 + (NVertex-1)*0.282  ); //3.9
       }
       float SigniCorProjMET = CorProjMET/7.4;
       //values taken on Z_2e_pu_powheg for MaxJet = 2 (1 leading jet accepted)

       /*   if(pdgID==13) { //test muon

	 CorProjMET =  (ProjMET - 0.60*(NVertex-1))/sqrt( 1 + (NVertex-1)*0.242  );
	 if(i==nt) {
	   CorProjMET =  (ProjMET - 0.60*(NVertex-1))/sqrt( 1 + (NVertex-1)*0.356  ); //3.9
	 }
	 SigniCorProjMET = CorProjMET/7.54;
	 }*/

       //=================================================

	//Define CorProjMET, Novosibirsk Def ===============================
       // float sPU = 3.024;
       // float s0  = 3.647;
       // float p0 = 0.3137;
       // float p1 = -0.0472;
       // float A = sqrt(log(4));
       // float tau = p0 + ((NVertex-1)*p1);
      
       // float sT = sqrt( s0*s0 + (NVertex-1)*(sPU*sPU) );
       // float se = sT*A/sinh(tau*A);
       // float se0 = s0*A/sinh(p0*A);
       
       // float CorProjMET = (ProjMET)*se0/se;
       // float SigniCorProjMET=CorProjMET/se0;

       // if(name[i]=="Z" && ie<100)
       // 	 cout<<ProjMET<<"    "<<CorProjMET<<"    " <<NVertex<<endl;

       //=================================================


      
      //Define NNout
      inputVecMET[0] = (double)METPrime[0];
      inputVecMET[1] = (double)METPrime[1];
      //    inputVecMET[2] = (double)METPrimeCor[0];
      //    inputVecMET[3] = (double)METPrimeCor[1];
      inputVecMET[2] = (double)METPrimeMETM[0];
      inputVecMET[3] = (double)METPrimeMETM[1];
      inputVecMET[4] = (double)METPrimeJetM[0]/NVertex;
      inputVecMET[5] = (double)METPrimeJetM[1]/NVertex;
      //  inputVecMET[8] = (double)METPrimeJetP[0]/NVertex;
      //    inputVecMET[9] = (double)METPrimeJetP[1]/NVertex;
      inputVecMET[6] = (double)METPrimeUnclusM[0]/NVertex;
      inputVecMET[7] = (double)METPrimeUnclusM[1]/NVertex;
      inputVecMET[8] = (double)METPrimeUnclusP[0]/NVertex;
      //  inputVecMET[13] = (double)METPrimeUnclusP[1]/NVertex;
      float NNmetOut = (float)(MetNN->GetMvaValue( inputVecMET ));
      //   if(NNmetOut<=NNoutCut) continue;
      //==================================================

      //estimation fake MET from Z+jet =====================
      {
	bool keep=true;
	if( Charge[0]*Charge[1] != -1 ) keep=false;
	if(convRej && ( !expInHit1 || !expInHit2 ) ) keep= false;
	if(JetMultiplicity >= MaxJet ) keep=false;
	if(PhotonMultiplicity >= MaxPhoton ) keep=false;
	if(ElecMult + MuonMult >= MaxLepton ) keep=false;
	
	if(CorProjMET/Z[0]>NormMETCutHigh) keep=false;
	if(CorProjMET/Z[0]<NormMETCutLow) keep=false;
	if( LJet[0]>25  && fabs(dPhiJetMET)< dPhiJMCut ) keep=false;
	//if( NNmetOut < NNoutCut ) keep =false;
	
	if(keep) {
	  //if(i==0)
	  // cout<<JetMultiplicity<<"   "<<PhotonMultiplicity<<"  "<<ElecMult + MuonMult<<endl;
	  histoManager.fill2D("ZMassVsSCPMET",i, ZMass, SigniCorProjMET , Weight );
	  histoManager.fill2D("ZMassVsUpara",i,RecoilProjZ[0], SigniCorProjMET , Weight );
	}
    }
  //estimation fake MET from Z+jet =====================
      
      //MakeGeneralCutB( true, true, "=", i, "candidat Z iso/ID", Weight );
      MakeCut( true, true, "=", i, "candidat Z iso/ID", Weight );
      
      //if( !  ) continue;

      if( (!invCut && (ZMass>120 || ZMass<MTCut) ) ||
	  (invCut && (ZMass<120 && ZMass>MTCut) ) ) 
	{  
	  // MakeGeneralCutB( false, true, "=", i, " 60 < M_{ll} < 120 GeV", Weight );
	  MakeCut( false, true, "=", i, " 60 < M_{ll} < 120 GeV", Weight );
	  continue;
	}
      else
	{ MakeCut( true, true, "=", i, " 60 < M_{ll} < 120 GeV", Weight ); }

      histoManager.fill("NVertex",i,NVertex, Weight);
      
    

      //(ProjMET -1.099*(NVertex-1))/( 0.896 + 0.1056*NVertex );
      histoManager.fill("ZCateg",i,ZCateg,Weight);
      histoManager.fill("Charge",i,Charge[0]*Charge[1],Weight);
      histoManager.fill2D("ChargeVsPt",i,Z[0],Charge[0]*Charge[1],Weight);

      if( ! MakeCut( Charge[0]*Charge[1] , -1, "=", i, "charge opposee", Weight ) ) continue;

      if( /*abs(PdgId[0])==11 &&  abs(PdgId[1])==11 &&*/ convRej) {


	histoManager.fill("expInHit1",i,expInHit1,Weight);
	histoManager.fill("inHit1",i,inHit1,Weight);
	histoManager.fill("convRej1",i,convRej1,Weight);

	histoManager.fill("expInHit2",i,expInHit2,Weight);
	histoManager.fill("inHit2",i,inHit2,Weight);
	histoManager.fill("convRej2",i,convRej2,Weight);

	bool unconv=true;

	if( abs(PdgId[0])==11 ) {
	  if(!expInHit1)
	    unconv =false;

	    //	  if( ! MakeGeneralCutB( expInHit1, true, "=", i, "ExpInHit1", Weight ) ) continue;
	}
	/*	else {
	  if( ! MakeGeneralCutB( true, true, "=", i, "ExpInHit1", Weight ) ) continue;
	  }*/

	if( abs(PdgId[1])==11 ) {
	  if(!expInHit2)
	    unconv =false;
	
	  //  if( ! MakeGeneralCutB( expInHit2, true, "=", i, "ExpInHit2", Weight ) ) continue;
	} 
	/*else {
	  if( ! MakeGeneralCutB( true, true, "=", i, "ExpInHit2", Weight ) ) continue;
	  }*/
	
	if( ! MakeCut( unconv, true, "=", i, "rejection #gamma", Weight ) ) continue;
      }
      
      //Matching between leading jet and leading lepton
      float dRJEl = dR(LeadAddElec[1], LJet[1],
		       LeadAddElec[2]*180/3.1415, LJet[2] );
      float dRJMu = dR(LeadAddMuon[1], LJet[1],
		       LeadAddMuon[2]*180/3.1415, LJet[2] );
      float dRJPh = dR(LPhoton[1], LJet[1],
		        LPhoton[2], LJet[2] );
      float dRPhEl = dR(LeadAddElec[1], LPhoton[1],
			LeadAddElec[2]*180/3.1415, LPhoton[2] );
      
      histoManager.fill("dRLMuonJet",i,dRJMu,Weight);
      histoManager.fill("dRLElecJet",i,dRJEl,Weight);
      histoManager.fill("dRPhotJet",i,dRJPh,Weight);
      histoManager.fill("dRPhotElec",i,dRPhEl,Weight);
      
      histoManager.fill("JetMult",i,JetMultiplicity,Weight);
      
      for(int im=0;im<10;im++) {
	ostringstream os;
	os << im;
	string n = "JetMult_Pt_"+os.str();
	histoManager.fill(n,i,MultiJetMultiplicity[im],Weight);
  }

      if( !fixJet && !invJet && ! MakeCut(JetMultiplicity, MaxJet, "<", i, "veto jet", Weight ) ) continue;
      if( !fixJet && invJet && ! MakeCut(JetMultiplicity, MaxJet, ">", i, "veto jet", Weight ) ) continue;
      if( fixJet && ! MakeCut(JetMultiplicity, MaxJet, "=", i, "veto jet", Weight ) ) continue;


      histoManager.fill("NVertex1",i,NVertex, Weight);

      histoManager.fill("PhotonMult",i,PhotonMultiplicity,Weight);
  
      // if(MaxJet!=1 && !fixJet && !invJet)
      //	if( ! MakeGeneralCutI(PhotonMultiplicity, MaxPhoton, "<", i, "PhotMult", Weight ) ) continue;
   
     

       histoManager.fill("MET",i,Met[0],Weight);   
       histoManager.fill("ProjMET",i,ProjMET,Weight);
       histoManager.fill("CorProjMET",i,CorProjMET,Weight);
      

      ostringstream os;
      if(NVertex<13)
	os<<NVertex;
      else
	os<<13;
      string mn="PFMET_V"+os.str();
      histoManager.fill(mn,i,Met[0],Weight);
      string prjn="ProjMET_V"+os.str();
      histoManager.fill(prjn,i,ProjMET,Weight);
      string cprjn="CorProjMET_V"+os.str();
      histoManager.fill(cprjn,i,CorProjMET,Weight);

      //Control Sample NN
      if(MaxJet!=1 && !fixJet && !invJet) {
	if( (SigniCorProjMET<=METcut && SigniCorProjMET>=15/7.36) &&
	    (CorProjMET/Z[0]<NormMETCutHigh) ) {
	  histoManager.fill("NNoutMETControl",i,NNmetOut,Weight);
	} 
      }
      else {
	if( (SigniCorProjMET<=METcut && SigniCorProjMET>=15/7.36) &&
	    (CorProjMET/Z[0]>NormMETCutLow && CorProjMET/Z[0]<NormMETCutHigh) &&
	    ( LJet[0]<=25 || ( LJet[0]>25 && fabs(dPhiJetMET)>dPhiJMCut) ) ) {
	
	  histoManager.fill("NNoutMETControl",i,NNmetOut,Weight);
	} 
      }
   
      histoManager.fill("SigniCorProjMETNorm",i,SigniCorProjMET/Z[0],Weight);
      histoManager.fill("CorProjMETNorm",i,CorProjMET/Z[0],Weight);
      histoManager.fill("ProjMETNorm",i,ProjMET/Z[0],Weight);


      {
	histoManager.fill("ZPtNoMN",i,Z[0],Weight);
	if( ! MakeCut(CorProjMET/Z[0], NormMETCutHigh, "<", i, "balance", Weight ) ) continue;
	histoManager.fill("ZPtMN",i,Z[0],Weight);
      }
      //   if(MaxJet!=1 && !fixJet && !invJet)
      //	if( ! MakeGeneralCutF(CorProjMET/Z[0], NormMETCutLow, ">", i, "METNL", Weight ) ) continue;
	  

      histoManager.fill("SigniCorProjMET",i,SigniCorProjMET,Weight);
      /*   if( ! MakeGeneralCutF(SigniCorProjMET, 0.25, ">", i, "  S_{#slash{E}_{T}}0.25", Weight ) ) continue;
      if( ! MakeGeneralCutF(SigniCorProjMET, 0.5, ">", i, "  S_{#slash{E}_{T}}0.5", Weight ) ) continue;
      if( ! MakeGeneralCutF(SigniCorProjMET, 0.75, ">", i, "  S_{#slash{E}_{T}}0.75", Weight ) ) continue;
      if( ! MakeGeneralCutF(SigniCorProjMET, 1, ">", i, "  S_{#slash{E}_{T}}1", Weight ) ) continue;
      */ //   if( ! MakeGeneralCutF(SigniCorProjMET, 1.5, ">", i, "  S_{#slash{E}_{T}}1.5", Weight ) ) continue;
      /*   if( ! MakeGeneralCutF(SigniCorProjMET, 2, ">", i, "  S_{#slash{E}_{T}}2", Weight ) ) continue;
	   if( ! MakeGeneralCutF(SigniCorProjMET, 2.5, ">", i, "  S_{#slash{E}_{T}}2.5", Weight ) ) continue;*/
      //    if( ! MakeGeneralCutF(SigniCorProjMET, 3, ">", i, "  S_{#slash{E}_{T}}3", Weight ) ) continue;
      
      histoManager.fill("CorProjMETNormNoMC",i,CorProjMET/Z[0],Weight);
      if(fabs(PdgId[0])==13 && fabs(PdgId[1])==13) {
	if(!invMet && !MakeCut(SigniCorProjMET, (float)(METcut+0.3), ">", i, "  S_{#slash{E}_{T}} finale", Weight ) ) continue;
	if(invMet && !MakeCut(SigniCorProjMET,  (float)(METcut+0.3), "<", i, "  S_{#slash{E}_{T}} finale", Weight ) ) continue;
      }
      else {
	if(!invMet && !MakeCut(SigniCorProjMET, METcut, ">", i, "  S_{#slash{E}_{T}} finale", Weight ) ) continue;
	if(invMet && !MakeCut(SigniCorProjMET, METcut, "<", i, "  S_{#slash{E}_{T}} finale", Weight ) ) continue;
      }
	histoManager.fill("CorProjMETNormMC",i,CorProjMET/Z[0],Weight);


      // if( SigniCorProjMET > 4 ) continue; //fixme, control sample for mistag top 

      //if(CorProjMET<20 || CorProjMET>27.5) continue;
    

      //Cuts and efficiencies ***********************************
      if( LJet[0]>=20 ) {
	histoManager.fill("dPhiJetMET",i,dPhiJetMET,Weight);
	if( ! MakeCut( fabs(dPhiJetMET), dPhiJMCut, ">", i, "#Delta#Phi(j,#slash{E}_{T})", Weight ) ) continue;
      }
      else{
	MakeCut( true, true, "=", i, "#Delta#Phi(j,#slash{E}_{T})", Weight );
      }
      if(MaxJet!=1 && !fixJet && !invJet) {
	if( LJet[0]>25 ) {
	  // histoManager.fill("dPhiJetMET",i,dPhiJetMET,Weight);
	  if( ! MakeCut( fabs(dPhiJetMET), dPhiJMCut, ">", i, "#Delta#Phi(j,#slash{E}_{T})", Weight ) ) continue;
	}
	else{
	  MakeCut( true, true, "=", i, "dPhiJM", Weight );
	}
      }      

      //temporaire
      // histoManager.fill("NNoutMET",i,NNmetOut,Weight);
      // if(MaxJet!=1 && !fixJet && !invJet) {
      // 	if( ! MakeGeneralCutF(NNmetOut, NNoutMETCut, ">", i, "NNMET", Weight ) ) continue;
      // }
      //if( ! MakeGeneralCutF( fabs(dPhiMETZ), 20 , ">", i, "dPhiMZ", Weight ) ) continue;

      
      histoManager.fill("LJetTrkCount",i,bTag[0],Weight);
      histoManager.fill("LJetSoftMuon",i,bTag[1],Weight);
      if(bTagFlag) {
	bool btagged=false;
	if(bTag[0]>2.5 || bTag[1]>0.23) 
	  btagged=true;
	if( ! MakeCut(btagged, false ,"=", i, "etiquetage b", Weight) ) continue;
	//	if( ! MakeGeneralCutF( bTag[0], 2.5, "<", i, "TrkCount", Weight ) ) continue;
	//	if( ! MakeGeneralCutF( bTag[1], 0.23, "<", i, "SoftMuon", Weight ) ) continue;
      }	  
	

	
      //BTag Control ===================================================
      for(int im=0;im<10;im++) {
	ostringstream os;
	os << im;
	string n = "BTag_Jet_"+os.str();
	if(MultiJetMultiplicity[im] >1 ) {
	  
	  if( bTag[0] < 2.5 && bTag[1]< 0.23)
	    histoManager.fill(n,i,true,Weight);
	  else
	    histoManager.fill(n,i,false,Weight);
	}
      }
      //================================================================

      //Resovar ================
      bool isData=false;
      if(name[i].substr(0,4)=="data")
	isData=true;
      
      float ResoPara =  GetNSigmaPara(Z[0], RecoilProjZ[0], NVertex, isData);
      float ResoPerp =  GetNSigmaPerp(Z[0], RecoilProjZ[1], NVertex, isData);

      histoManager.fill("ResoPara",i,ResoPara,Weight);
      histoManager.fill("ResoPerp",i,ResoPerp,Weight);
      //Resovar ================


      //NeuralNet
      inputVec[0] = (double)(METProjL1[1]/Z[0]);
      inputVec[1] = (double)(METProjL1[0]/Z[0]);
      inputVec[2] = (double)(METProjL2[0]/Z[0]);
      inputVec[3] = (double)(METProjL2[1]/Z[0]);
     
      inputVec[4] = (double)(Lepton[0][0]);
      inputVec[5] = (double)(dPhiLeptonZ);
      inputVec[6] = (double)(CosThetaP);
      /* inputVec[5] = (double)(dRLeptonS);
      inputVec[6] = (double)(METPrimeMETM[1]/Z[0]);
      inputVec[7] = (double)(METPrimeMETM[0]/Z[0]);*/
      float NNout = (float)(NN->GetMvaValue( inputVec ));
         
      int lepMult=0;
      // if(IsIDElec || IsIsoElec || ElecMult>1)
      lepMult +=ElecMult;
      // if(IsIsoMu || MuonMult>1)
      lepMult +=MuonMult;
      histoManager.fill("LeptonMult",i, ElecMult + MuonMult,Weight);
       if( ! MakeCut( lepMult, MaxLepton, "<", i, "veto lepton", Weight ) ) continue;
 

       histoManager.fill("ZPtNoNN",i,Z[0],Weight);
       histoManager.fill("METNoNN",i,Met[0],Weight);  
       
       histoManager.fill("NNout",i,NNout,Weight);
       if( ! MakeCut( NNout, NNoutCut, ">", i, "sortie NN", Weight ) ) continue;
       
      //lepton
      for(int id=0;id<2;id++) {
	histoManager.fill("Pt"+lep[id],i,Lepton[id][0],Weight);
	histoManager.fill("Eta"+lep[id],i,Lepton[id][1],Weight);
	histoManager.fill("Phi"+lep[id],i,Lepton[id][2],Weight);
	histoManager.fill("P"+lep[id],i,Lepton[id][3],Weight);
	histoManager.fill("PdgId"+lep[id],i,PdgId[id],Weight);
      }
      histoManager.fill("ElecMult",i,ElecMult,Weight);
      histoManager.fill("MuonMult",i,MuonMult,Weight);
      histoManager.fill("LeptonMultEnd",i,ElecMult+MuonMult,Weight);
  

      if(LeadAddElec[0]!=-100) {
	histoManager.fill("LElecPt",i,LeadAddElec[0],Weight);
	histoManager.fill("LElecEta",i,LeadAddElec[1],Weight);
	histoManager.fill("LElecPhi",i,LeadAddElec[2]*180/3.14,Weight);
       	histoManager.fill("LElecID",i,IsIDElec,Weight);
	histoManager.fill("LElecIso",i,IsIsoElec,Weight);
      }
      if(LeadAddMuon[0]!=-100) {
	histoManager.fill("LMuonPt",i,LeadAddMuon[0],Weight);
	histoManager.fill("LMuonEta",i,LeadAddMuon[1],Weight);
	histoManager.fill("LMuonPhi",i,LeadAddMuon[2]*180/3.14,Weight);
	histoManager.fill("LMuonID",i,IsIDMu,Weight);
	histoManager.fill("LMuonIso",i,IsIsoMu,Weight);
      }
      histoManager.fill("RatioPt",i,Lepton[0][0]/Lepton[1][0],Weight);

    
      //Variables After Selection
      histoManager.fill("ChargeEnd",i,Charge[0]*Charge[1],Weight);
      histoManager.fill("JetMultEnd",i,JetMultiplicity,Weight);
      histoManager.fill("NVertexEnd",i,NVertex, Weight);
      histoManager.fill("METEnd",i,Met[0],Weight);   
      histoManager.fill("ProjMETEnd",i,ProjMET,Weight);
      histoManager.fill("CorProjMETEnd",i,CorProjMET,Weight);
      histoManager.fill("SigniCorProjMETEnd",i,SigniCorProjMET,Weight);
      histoManager.fill("SigniCorProjMETNormEnd",i,SigniCorProjMET/Z[0],Weight);
      histoManager.fill("CorProjMETNormEnd",i,CorProjMET/Z[0],Weight);
      histoManager.fill("NNoutMETEnd",i,NNmetOut,Weight);
      histoManager.fill("dPhiJetMETEnd",i,dPhiJetMET,Weight);

      //Z
      histoManager.fill("ZMass",i,ZMass,Weight);
      histoManager.fill("PfZMass",i,PfZMass,Weight);
      histoManager.fill("ZPt",i,Z[0],Weight);
      histoManager.fill("ZEta",i,Z[1],Weight);
      histoManager.fill("ZPhi",i,Z[2],Weight);
      histoManager.fill("ZP",i,Z[3],Weight);
     
      histoManager.fill("ZMult",i,ZMultiplicity,Weight);
    
      histoManager.fill("ZY",i,Rapidity(Z[3], ZMass, Z[0], Z[1]), Weight);

      //Met
  
      histoManager.fill("PhiMET",i,Met[1],Weight);
      histoManager.fill("Recoil",i,Recoil[0],Weight);
      histoManager.fill("PhiRecoil",i,Recoil[1]*180/3.14,Weight);
    
      histoManager.fill("SpecMET",i,fabs(SpecMET),Weight);
      histoManager.fill("METNorm",i,Met[0]/Z[0],Weight);
      histoManager.fill("TcMET",i,TrkMet[0],Weight);
      
      histoManager.fill("METPrime",i,sqrt(METPrime[1]*METPrime[1]+(1.5*1.5*METPrime[0]*METPrime[0])),Weight);
      
      histoManager.fill("MetPrimeTrans",i,METPrime[0],Weight);
      histoManager.fill("MetPrimeLong",i,METPrime[1],Weight);
      histoManager.fill("MetPrimeRecoilTrans",i,METPrimeRecoil[0],Weight);
      histoManager.fill("MetPrimeRecoilLong",i,METPrimeRecoil[1],Weight);
      //histoManager.fill("MetPrimeCorRecoilTrans",i,METPrimeCorRecoil[0],Weight);
      // histoManager.fill("MetPrimeCorRecoilLong",i,METPrimeCorRecoil[1],Weight);
      histoManager.fill("MetPrimeCorTrans",i,METPrimeCor[0],Weight);
      histoManager.fill("MetPrimeCorLong",i,METPrimeCor[1],Weight);
      //histoManager.fill("MetPrimeCorCaloTrans",i,METPrimeCorCalo[0],Weight);
      //histoManager.fill("MetPrimeCorCaloLong",i,METPrimeCorCalo[1],Weight);
    

    
      //  cout<<ProjMET<<endl;
      
      //Jet
     
      histoManager.fill("TcJetMult",i,TcJetMultiplicity,Weight);
      histoManager.fill("LJetPt",i,LJet[0],Weight);
      if(LJet[0]>25) {
	histoManager.fill("LJetEta",i,LJet[1],Weight);
	histoManager.fill("LJetPhi",i,LJet[2],Weight);
	histoManager.fill("LJetFrac",i,LJet[3],Weight);

      }

      //Photon
      //    histoManager.fill("PhotonMult",i,PhotonMultiplicity,Weight);
      histoManager.fill("LPhotonPt",i,LPhoton[0],Weight);
      histoManager.fill("LPhotonEta",i,LPhoton[1],Weight);
      histoManager.fill("LPhotonPhi",i,LPhoton[2],Weight);
   
      //Lepton related variables
      histoManager.fill("dPhiLepton",i,dPhiLepton,Weight);
      histoManager.fill("dPhiLeptonZ",i,dPhiLeptonZ,Weight);
      histoManager.fill("dEtaLepton",i,dEtaLepton,Weight);
      histoManager.fill("dPtLepton",i,dPtLepton,Weight);
      histoManager.fill("dPLepton",i,dPLepton,Weight);
      histoManager.fill("dRLepton",i,dRLepton,Weight);
 
      //boosted leptons variables
      histoManager.fill("dPhiLeptonS",i,dPhiLeptonS,Weight);
      histoManager.fill("dEtaLeptonS",i,dEtaLeptonS,Weight);
      histoManager.fill("dPtLeptonS",i,dPtLeptonS,Weight);
      histoManager.fill("dPLeptonS",i,dPLeptonS,Weight);
      histoManager.fill("dRLeptonS",i,dRLeptonS,Weight);
      histoManager.fill("CosThetaM",i,CosThetaM,Weight);
      histoManager.fill("CosThetaP",i,CosThetaP,Weight);

      histoManager.fill("CosThetaCSM",i,CosThetaCSM,Weight);
      if(Rapidity(Z[3], ZMass, Z[0], Z[1]) > 1.7 ) {
	histoManager.fill("CosThetaCSP",i,CosThetaCSP,Weight);
      }
      histoManager.fill2D("CSAvsY",i,CosThetaCSP,Rapidity(Z[3], ZMass, Z[0], Z[1]),Weight);

      //MET related variables 
      histoManager.fill("ZMETMass",i,ZMETMass,Weight); //ZMETMass //sqrt(2*Z[0]*Met[0]*(1-cos(dPhiMETZ*3.1415/180)))
      histoManager.fill("ZMassMET",i,ZMassMET,Weight);
      histoManager.fill("ZPtMET",i,ZPtMET,Weight);
      histoManager.fill("ZPMET",i,ZPMET,Weight);
      histoManager.fill("dPhiMETZ",i,dPhiMETZ,Weight);
      histoManager.fill("dPhiMETL1",i,dPhiMETL1,Weight);
      histoManager.fill("dPhiMETL2",i,dPhiMETL2,Weight);
      histoManager.fill("METPara",i,METProjZ[0],Weight);
      histoManager.fill("METPerp",i,METProjZ[1],Weight);
      histoManager.fill("L1METPara",i,METProjL1[0],Weight);
      histoManager.fill("L1METPerp",i,METProjL1[1],Weight);
      histoManager.fill("L2METPara",i,METProjL2[0],Weight);
      histoManager.fill("L2METPerp",i,METProjL2[1],Weight);
      histoManager.fill("L1METParaNorm",i,METProjL1[0]/Z[0],Weight);
      histoManager.fill("L1METPerpNorm",i,METProjL1[1]/Z[0],Weight);
      histoManager.fill("L2METParaNorm",i,METProjL2[0]/Z[0],Weight);
      histoManager.fill("L2METPerpNorm",i,METProjL2[1]/Z[0],Weight);
      histoManager.fill("ZRecoil",i,ZRecoil,Weight);
      histoManager.fill("RecoilPara",i,RecoilProjZ[0],Weight);
      histoManager.fill("RecoilPerp",i,RecoilProjZ[1],Weight);
      histoManager.fill("L1RecoilPara",i,RecoilProjL1[0],Weight);
      histoManager.fill("L1RecoilPerp",i,RecoilProjL1[1],Weight);
      histoManager.fill("L2RecoilPara",i,RecoilProjL2[0],Weight);
      histoManager.fill("L2RecoilPerp",i,RecoilProjL2[1],Weight);
      histoManager.fill("dPhiRecoilZ",i,dPhiRecoilZ,Weight);
      histoManager.fill("dPhiRecoilL1",i,dPhiRecoilL1,Weight);
      histoManager.fill("dPhiRecoilL2",i,dPhiRecoilL2,Weight);
      histoManager.fill("ZPtmMet",i,Z[0]-Met[0],Weight);
    
      histoManager.fill("dPhiRecoilMET",i,360-fabs(dPhiRecoilZ)-fabs(dPhiMETZ),Weight);
      
      histoManager.fill("HT",i,HT[0],Weight);

      histoManager.fill("SumET",i,SumET,Weight);

      //Jet related variables
      histoManager.fill("dPhiJetZ",i,dPhiJetZ,Weight);
      histoManager.fill("dEtaJetZ",i,dEtaJetZ,Weight);
      histoManager.fill("ZPtJet",i,ZPtJet,Weight);
      histoManager.fill("ZJetMass",i,ZJetMass,Weight);
      histoManager.fill("dPhiJetL1",i,dPhiJetL1,Weight);
      histoManager.fill("dPhiJetL2",i,dPhiJetL2,Weight);
      histoManager.fill("dEtaJetL1",i,dEtaJetL1,Weight);
      histoManager.fill("dEtaJetL2",i,dEtaJetL2,Weight);
     
      //Photon related variables
      histoManager.fill("dPhiPhotonMET",i,dPhiPhotonMET,Weight);
      histoManager.fill("dPhiPhotonZ",i,dPhiPhotonZ,Weight);
      histoManager.fill("dEtaPhotonZ",i,dEtaPhotonZ,Weight);
      histoManager.fill("ZPtPhoton",i,ZPtPhoton,Weight);
      histoManager.fill("ZPhotonMass",i,ZPhotonMass,Weight);
      histoManager.fill("dPhiPhotonL1",i,dPhiPhotonL1,Weight);
      histoManager.fill("dPhiPhotonL2",i,dPhiPhotonL2,Weight);
      histoManager.fill("dEtaPhotonL1",i,dEtaPhotonL1,Weight);
      histoManager.fill("dEtaPhotonL2",i,dEtaPhotonL2,Weight);
     
      if(name[i].substr(0,2)=="di" || name[i].substr(0,2)=="ZZ"
	 ||  name[i].substr(0,2)=="ZV") {
	string proc = FindProcess( ie, i );
	if( proc.substr(0,2)=="ZZ") {
	  TLorentzVector Zll(0,0,0,0); Zll.SetPtEtaPhiM(McZl[0],McZl[1],McZl[2],McZl[3]);
	  TLorentzVector Znn(0,0,0,0); Znn.SetPtEtaPhiM(McZn[0],McZn[1],McZn[2],McZn[3]);
	  StoreSelectedEvents( (Zll+Znn).M(), Weight );

	  histoManager.fill("McZZMass",i,(Zll+Znn).M(), Weight );
	}
      }


      //if(METcut>25) {

	

	// }

      //End Filling*********************
	if(name[i]=="Zll" && METcut>4.) {
	  cout<<fileName<<"   "<<Event<<"   "<<JetMultiplicity<<"   "<<LJet[0]<<"   "<<NNmetOut<<"    "<<dPhiJetMET<<"   "<<TrkMet[0]<<"    "<<SigniCorProjMET<<"   "<<CorProjMET<<"   "<<Z[0]<<endl;
	}
	if(name[i].substr(0,4)=="data")
	  {
	    std::pair<int,int> tmp(Run,AbsEvent);
	    string t1(sample),t2(fileName);
	    std::pair<string,string> tmp2( t1, t2 );
	    std::pair< std::pair<int,int> , std::pair<string,string> > tmp3(tmp,tmp2);
	    Events.push_back(tmp3);
	    EvtsInFile.push_back(Event);
	   

	 
	    /* if(name[i]=="data") {
	      cout<<ProjMET<<"  "<<NVertex<<"  --> "<<CorProjMET<<"   "<<SigniCorProjMET<<"   "<<Run<<"   "<<AbsEvent<<endl;
	      }*/
	  }
      if(name[i]=="ZZ #rightarrow 2l2#nu")
	{ S++; 
	}
      else if(name[i]!="data")
	{B++; }
      
      NumberEntries[i]++;
	    
      //   if(MakeSkim)
      //	newtree->Fill();
    }//End events
    
    /*  if(MakeSkim) {
      newtree->Write();
      skimfile->Close();
    }*/
    //  cout<<" end dataset "<<name[i]<<endl;
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
ZZUseTree::PrepareDatasets() {

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

  /* Sorder.push_back("ZZ #rightarrow 2l2#nu");
  Sorder.push_back("ZZ #rightarrow 4l");
  Sorder.push_back("ZZ");
  Sorder.push_back("WZ");
  Sorder.push_back("WW");
  Sorder.push_back("Z #rightarrow 2#tau");
  Sorder.push_back("Z #rightarrow 2l");
  Sorder.push_back("t#bar{t}");*/

  cout<<" End dataset preparation : nt="<<nt<<endl;

}


TH1F* 
ZZUseTree::GetVtxCorrection(string obs,int nds, int bin) {

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
ZZUseTree::GetFitWeight() {

  vector<int> NumberEntries(nt+1,0);
  
  cout<<" ==================================================== "<<endl;
  cout<<" ================= En cours de weight =============== "<<endl;

  //Needed by the NN
  vector<double> inputVec(8,0);

  //sÃ©lection
  for(int i=0;i<nt+1;i++) {

    //Déclaration des variables
   
    //Vertices
    int NVertex;
   
    //Lepton
    float IsoVar[2][3];
    float IDVar[2][4];
    float Lepton[2][4];
    int PdgId[2];

    bool expInHit1;
    bool expInHit2;
    
    int ElecMult;
    int MuonMult;

    //MET
    float METPrime[2];
    float METPrimeRecoil[2];
  
    float METProjL1[2];
    float METProjL2[2];

    float SumET;

    //Jet
    int JetMultiplicity;
    
    //Photon
    int PhotonMultiplicity;
    
    //Boosted leptons variables
    float dPtLeptonS;
    float dRLeptonS;
 
    //Z
    float Z[4];

    float Weight=1;

   
 

    bool Sel[2]={true,true};

    tChains[i]->SetBranchAddress("NVertex",&NVertex);

    // Leptons 
    tChains[i]->SetBranchAddress("IDVar",IDVar);
    tChains[i]->SetBranchAddress("IsoVar",IsoVar);
    tChains[i]->SetBranchAddress("Lepton",Lepton);
    tChains[i]->SetBranchAddress("PdgId",PdgId);

    tChains[i]->SetBranchAddress("ExpInHit1",&expInHit1);
    tChains[i]->SetBranchAddress("ExpInHit2",&expInHit2);
	
     tChains[i]->SetBranchAddress("ElecMult",&ElecMult);
     tChains[i]->SetBranchAddress("MuonMult",&MuonMult);

     //Z
     tChains[i]->SetBranchAddress("Z",Z);

     //MET
     tChains[i]->SetBranchAddress("METPrime",METPrime);
     //     tChains[i]->SetBranchAddress("METPrimeCor",METPrimeCor);
     tChains[i]->SetBranchAddress("METPrimeRecoil",METPrimeRecoil);
     //  tChains[i]->SetBranchAddress("METPrimeCorRecoil",METPrimeCorRecoil);
     // tChains[i]->SetBranchAddress("METPrimeCorCalo",METPrimeCorCalo);

     tChains[i]->SetBranchAddress("METProjL1",METProjL1);
     tChains[i]->SetBranchAddress("METProjL2",METProjL2);

     //Jet
     tChains[i]->SetBranchAddress("JetMultiplicity",&JetMultiplicity);

     //Photon
     tChains[i]->SetBranchAddress("PhotonMultiplicity",&PhotonMultiplicity);

     //Boosted leptons variables
     tChains[i]->SetBranchAddress("dPtLeptonS",&dPtLeptonS);
     tChains[i]->SetBranchAddress("dRLeptonS",&dRLeptonS);

     tChains[i]->SetBranchAddress("SumEt",&SumET);
     
    int EP=0;
    cout<<" beginning tree " <<name[i]<<"   "<<tChains[i]->GetEntries()<<endl;
    //Boucle et sélection, remplissage histos
    int ent = tChains[i]->GetEntries();
    
    for(int ie=0;ie<ent;ie++) {
	
      if(ie%1000==0)
	fprintf(stdout, "Done : %d / %d\r", ie, ent );

      if(i!=nt)
	Weight = GetWeight(i, ie);
      else
	Weight = 1;

      tChains[i]->GetEntry(ie);      
     
      if(pdgID==11) {
	
	for(int id=0;id<2;id++) {
	
	  if(fabs(Lepton[id][1])>1.479)
	    {EP=1; }
	  else
	    {EP=0; }

	  Sel[id]=true;

	  if(IDVar[id][0] > _IDcuts[EP][0] && Cbool(skipCut, "sigieie"==Nm1Var.substr(0,7)) )  Sel[id]=false;
	  if( fabs(IDVar[id][1]) > _IDcuts[EP][1] && Cbool(skipCut, "Deta"==Nm1Var.substr(0,4)) )  Sel[id]=false; 
	  if( fabs(IDVar[id][2]) > _IDcuts[EP][2] && Cbool(skipCut, "Dphi"==Nm1Var.substr(0,4)) )  Sel[id]=false;
	  if(IDVar[id][3] > _IDcuts[EP][3] && Cbool(skipCut, "HoE"==Nm1Var.substr(0,3)) )  Sel[id]=false;
	  if(IsoVar[id][0]/Lepton[id][0] > _Isocuts[EP][0] && Cbool(skipCut, "TrackIso"==Nm1Var.substr(0,8)) )  Sel[id]=false;
	  if(IsoVar[id][1]/Lepton[id][0] > _Isocuts[EP][1] && Cbool(skipCut, "EcalIso"==Nm1Var.substr(0,7)) )  Sel[id]=false;

	  if(IsoVar[id][2]/Lepton[id][0] > _Isocuts[EP][2] && Cbool(skipCut, "HcalIso"==Nm1Var.substr(0,7)) )  Sel[id]=false;

	  if(!Sel[id]) break;
	}

	if(!Sel[0] || !Sel[1] ) continue;
      }

      //Filling Histos******************
    
      //      float METPCor;
      /*
      //1.45, 0.4
      double Tt = METPrime[0] + METPrimeRecoil[0] + 1.4*METPrimeCorCalo[0] + 1.4*METPrimeCor[0];
      double Tl = METPrime[1] + METPrimeRecoil[1] + 1.4*METPrimeCorCalo[1] + 1.4*METPrimeCor[1];
      Tt = max(Tt,0.);
      Tl = max(Tl,0.);

      METPCor = sqrt(1.5*1.5*Tt*Tt + Tl*Tl);   
     */
      //Cuts and efficiencies ***********************************
      if(abs(PdgId[0]) != pdgID || abs(PdgId[1]) != pdgID ) continue; 
    
      if(JetMultiplicity > MaxJet ) continue; 
      if(PhotonMultiplicity > MaxPhoton ) continue;
      if(ElecMult + MuonMult > MaxLepton ) continue;
      if(pdgID==11 && convRej) {
	if(!expInHit1) continue;
	if(!expInHit2) continue;
      } 
      
      //  if(METPCor<METcut ) continue; 
      // if(METPCor/Z[0] > NormMETCut ) continue;
     
      /*  //NeuralNet
      inputVec[0] = (double)(METProjL1[1]/Z[0]);
      inputVec[1] = (double)(METProjL1[0]/Z[0]);
      inputVec[2] = (double)(METProjL2[0]/Z[0]);
      inputVec[3] = (double)(METProjL2[1]/Z[0]);
     
      inputVec[4] = (double)(Lepton[0][0]);
      // inputVec[7] = (double)(dPtLeptonS);
      inputVec[5] = (double)(dRLeptonS);
      inputVec[6] = (double)(METPrimeMETM[1]/Z[0]);
      inputVec[7] = (double)(METPrimeMETM[1]/Z[0]);
      float NNout = (float)(NN->GetMvaValue( inputVec ));
      */
      //      if(NNout<=NNoutCut) continue;

      histoManager.fill("NVertexControl",i,NVertex, Weight);

      histoManager.fill2D("SumETNVertexControl",i,SumET,NVertex, Weight);

    }//End events
 }// End datasets
 
    
 Histos = histoManager.GetHistos();
 Histos2D = histoManager.GetHistos2D();
    
 vector<vector<vector<float> > > Wghts(5,vector<vector<float> >(2,vector<float>(2,1)));
 
 return Wghts;
    
    
}
 

/*
void 
ZZUseTree::PlotEffRej(string Obs,string addObs) {
  
  c2=new TCanvas("c2","Test",300,300,600,600);
  leg = new TLegend(0.70,0.62,0.9,0.82);
  leg->SetLineColor(0);
  leg->SetFillColor(0);

  TGraphErrors* EffRej;
  
  TH1F* signal(0);
  TH1F* totalObs(0);
  
  int nvar = histoManager.FindNVar(Obs);
  
  //MM FIXME
  Reweight(nvar, true);
  
  for(int j=0;j<nt;j++) {
    
    int Nhisto = histoManager.access(nvar, j);
      
    if(j==0) {
	signal = (TH1F*)Histos[Nhisto]->Clone();
      }
      else {
	if(j==1)
	  totalObs = (TH1F*)Histos[Nhisto]->Clone();
	else {
	  totalObs->Add(Histos[Nhisto]);
	}
      }
    }
    
  EffRej = histoManager.EffRej( signal, totalObs);
  EffRej->GetXaxis()->SetTitle("signal efficiency" );
  EffRej->GetYaxis()->SetTitle("bckgrd efficiency" );
  EffRej->SetMarkerStyle(20);
  EffRej->SetMarkerColor(kRed-2);
  EffRej->SetLineColor(kRed-2);
  leg->AddEntry(EffRej, (Obs).c_str(),"pl" );

  EffRej->Draw("APL");
  leg->Draw("same");


  if(addObs=="") return;

 TGraphErrors* EffRej2;
  
 TH1F* signal2(0);
 TH1F* totalObs2(0);

  int nvar2 = histoManager.FindNVar(addObs);
  
  //MM FIXME
  Reweight(nvar2, true);
  
  for(int j=0;j<nt;j++) {
    
    int Nhisto = histoManager.access(nvar2, j);
       if(j==0) {
	signal2 = (TH1F*)Histos[Nhisto]->Clone();
      }
      else {
	if(j==1)
	  totalObs2 = (TH1F*)Histos[Nhisto]->Clone();
	else {
	  totalObs2->Add(Histos[Nhisto]);
	}
      }
  }
    
  EffRej2 = histoManager.EffRej( signal2,totalObs2);
  EffRej2->GetXaxis()->SetTitle("signal efficiency" );
  EffRej2->GetYaxis()->SetTitle("bckgrd efficiency" );
  EffRej2->SetMarkerStyle(20);
  EffRej2->SetMarkerColor(kBlue-2);
  EffRej2->SetLineColor(kBlue-2);
  leg->AddEntry(EffRej2, (addObs).c_str(),"pl" );

  EffRej2->Draw("PL");
  leg->Draw("same");
  
  delete signal;
  delete signal2;
  delete totalObs;
  delete totalObs2;


}
*/

IClassifierReader* 
ZZUseTree::CreateNN() {

 //Définition des variables (comme pour le training)
  vector<string> inputVars;

  inputVars.push_back("METProjL1[1]/Z[0]");
  inputVars.push_back("METProjL1[0]/Z[0]");
  inputVars.push_back("METProjL2[0]/Z[0]");
  inputVars.push_back("METProjL2[1]/Z[0]");
 
  inputVars.push_back("Lepton[0][0]");
  inputVars.push_back("dPhiLeptonZ");
  inputVars.push_back("CosThetaCSP");
  /* inputVars.push_back("dRLeptonS");
  inputVars.push_back("METPrimeMETM[1]/Z[0]");
  inputVars.push_back("METPrimeMETM[0]/Z[0]");*/

  // on crée un neural net du type NeuralNetWWZZ
 
  IClassifierReader* NNZZWW = new NeuralNet( inputVars );

  return NNZZWW;

}


IClassifierReader* 
ZZUseTree::CreateMETNN() {

 //Définition des variables (comme pour le training)
  vector<string> inputVars;

  inputVars.push_back("METPrime[0]");
  inputVars.push_back("METPrime[1]");
  //  inputVars.push_back("METPrimeCor[0]");
  //  inputVars.push_back("METPrimeCor[1]");
  inputVars.push_back("METPrimeMETM[0]");
  inputVars.push_back("METPrimeMETM[1]");
  inputVars.push_back("METPrimeJetM[0]/NVertex");
  inputVars.push_back("METPrimeJetM[1]/NVertex");
  //  inputVars.push_back("METPrimeJetP[0]/NVertex");
  //  inputVars.push_back("METPrimeJetP[1]/NVertex");
  inputVars.push_back("METPrimeUnclusM[0]/NVertex");
  inputVars.push_back("METPrimeUnclusM[1]/NVertex");
  inputVars.push_back("METPrimeUnclusP[0]/NVertex");
  //  inputVars.push_back("METPrimeUnclusP[1]/NVertex");
  // on crée un neural net du type METNN
 
  IClassifierReader* mNN = new METNN( inputVars );

  return mNN;

}



float ZZUseTree::Rapidity(float p, float mass, float pt, float eta) {

  //compute pz
  float pz = sqrt(p*p-pt*pt);

  //compute E
  float E = sqrt(p*p + mass*mass);

  //and now rapidity
  float y = 0.5 * log( (E+pz) / (E-pz) );

  return y;

}




float ZZUseTree::GetNSigmaPara(float ZPt, float RPara, int Nv, bool isData) {

  float R = Response(ZPt);
  float sigma;

  RPara*= 1./R;
  
  // cout<<RPara*R<<"  ---> "<< RPara<<"   "<<ZPt<<endl;

  if(isData)
    sigma = sqrt( pow( sqrt(ZPt)*0.919644 + 1.71437  ,2 )   + pow(5.08886 ,2 )*pow( R, 2 ) + ( pow(3.41446,2)*(Nv-1) ) *pow( R, 2 ) );
  else
    sigma = sqrt( pow( sqrt(ZPt)*1.25721 + 0  ,2 )   + pow(4.9347 ,2 )*pow( R, 2 ) + ( pow(3.52464,2)*(Nv-1)) *pow( R, 2 ) );
 
  return (RPara+ZPt)/sigma;
}



float ZZUseTree::GetNSigmaPerp(float ZPt, float RPerp, int Nv, bool isData) {

  float R = Response(ZPt);// cout<<ZPt<<"   "<<R<<endl;
  float sigma;

  RPerp*= 1./R;

  if(isData)
    sigma = sqrt( pow( sqrt(ZPt)*0.626908 + 1.42442  ,2 )   + pow(4.83213 ,2 )*pow( R, 2 ) + ( pow(3.67182,2)*(Nv-1)) *pow( R, 2 ) );
  else
    sigma = sqrt( pow( sqrt(ZPt)*0.66288 + 0  ,2 )   + pow(4.90741 ,2 )*pow( R, 2 ) + ( pow(3.50873,2)*(Nv-1)) *pow( R, 2 ) );
 
  return (RPerp)/sigma;
}


float ZZUseTree::Response(float pt) {

  //For pfT1
  float r;
  r = 1.01379 + (-43.3708) / (125.404 + pt*pt ) + (0.562124) / (13.5828 + pt );
  return r;
}




void ZZUseTree::Get2DNumber(string obs, string chan, ostream& os) {

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

	//	if(ix==0 && iy==(int)temp[1]+1)
	//  cout<<"    "<<XTitle;
    }
      os<<endl;
      //  if(ix==(int)temp[0]-1)
      //    cout<<YTitle<<endl;
    
  }


  cout<<" Total -> "<<total<<endl<<endl;

}


void ZZUseTree::GetAll2DNumbers(string var) {

  ofstream out("/home/mmarionn/Documents/CMS/ZZAnalysis/ZJetEstimation/ZMassSCPMETNumbers.txt", ios::out | ios::trunc );

  if(!out)
    { cout<<" Error writing file, abort"<<endl; abort(); }

  for(int unsigned i=0;i<name.size();i++) {
    Get2DNumber(var,name[i], cout );

    if(i>=(size_t)(nt-1)) {
      Get2DNumber(var,name[i], out );
    }
    else {
      ofstream out2( ("/home/mmarionn/Documents/CMS/ZZAnalysis/ZJetEstimation/ZMassSCPMETNumbers_"+name[i]+".txt").c_str() , ios::out | ios::trunc );
      Get2DNumber(var,name[i], out2 );
    }
  }
  
}



void ZZUseTree::PrintGlobal2DNumbers(string var) {
  int Obs = histoManager.FindNVar2D(var);
  vector<float> temp = histoManager.GetTemplate2D(Obs);
  float total;

  vector<vector<float> > values((int)temp[0], vector<float>((int)temp[1],0));

  for(int unsigned i=0;i<name.size()-1;i++) {

  
    if(Obs==-1) {
      cout<<" Be careful, no such observable "<<var<<endl;
      return;
    }
    XTitle = (histoManager.FindLeg2D( var )).first;
    YTitle = (histoManager.FindLeg2D( var )).second;
    cout<<name[i]<<endl;
    int NDS= GetNbyName(var);
    int Nhisto = histoManager.access2D(Obs,i);

    TH2F* histo = (TH2F*)Histos2D[ Nhisto ]->Clone();

    float total=0;
    
    for(int ix=0;ix<(int)temp[0];ix++) {

      for(int iy=0;iy<(int)temp[1];iy++) {
	
	values[ix][iy] += histo->GetBinContent(ix+1,iy+1);

	if(ix==1)
	  total += histo->GetBinContent(ix+1,iy+1);

      }
    }
  }

  cout<<" Total -> "<<total<<endl<<endl;

  
 for(int ix=0;ix<(int)temp[0];ix++) {
   for(int iy=0;iy<(int)temp[1];iy++) {
	cout<<"   "<<values[ix][iy];
      }
      cout<<endl;
 }

}



void 
ZZUseTree::PlotAllVtx(string chan, string var, int N) {


  TCanvas* c4=new TCanvas("Mets","Mets",600,600);

  int colors[10]={kBlue+1,/*kBlue-3,*/
		  /* kCyan+1,*/kTeal+4,
		  /* kGreen+1,*/kOrange+2,
		  /* kOrange+7,*/kRed+1,
		  /*  kViolet-2,kMagenta+2*/};
  
  int NDS= GetNbyName(chan);

  TH1F* emptyHisto= new TH1F("bidon","bidon",200,0,200);

  emptyHisto->Draw();
  emptyHisto->GetYaxis()->SetRangeUser(0.00001,1);
  emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);

  TLegend* legend=new TLegend(0.589,0.537,0.916,0.928);
  legend->SetFillColor(0);
  legend->SetLineColor(0);
  legend->SetShadowColor(0);

  for(int i=1;i<N;i++) {

    ostringstream os;
    os << i;


    string obsnam = var + "MET_V" + os.str() ;
    int Obs = histoManager.FindNVar(obsnam);
    int Nhisto = histoManager.access(Obs, NDS);
    
    TH1F* histo = (TH1F*)Histos[ Nhisto ]->Clone();

    histo->Rebin(Bin);

    histo->SetLineColor(colors[i-1]);
    histo->SetLineWidth(2);
    if(i==0)
      histo->DrawNormalized("hist");
    else
      histo->DrawNormalized("same hist");

    ostringstream o;
    o<<i;
    string leglab = os.str() + " vertex";
    if(i==1)
      leglab += " (pas de PU)";

    legend->AddEntry(histo,leglab.c_str(), "l");

    if(N-1==1)
      {
	TF1* f=new TF1("f",ZZUseTree::FuncProjMET,0,40,5);
	histo->Fit("f","R");

	f->SetLineWidth(2);
	f->SetLineStyle(7);
	f->SetLineColor(kRed+1);

	f->Draw("same");
      }

  }

  legend->Draw();

}



double
ZZUseTree::FuncProjMET(double* x, double* par) {

  double val;

  if(x[0]>15)
    val = par[3]*exp(x[0]*par[4]);
  else {
    val = par[0]*TMath::Gaus(x[0],par[1],par[2]);
  }

  return val;

}



void 
ZZUseTree::DrawJetMultiMult(string chan) {

  TCanvas* c4=new TCanvas("JetMult","JetMult",600,600);

  int NDS= GetNbyName(chan);

  float ThrJetMult[10]={15,20,25,30,35,40,50,60,80,100};

  TH1F* Efficiencies = new TH1F("Efficiencies","Efficiencies",9,ThrJetMult);
  TH1F* EffData = new TH1F("EffData","EffData",9,ThrJetMult);
  
  for(int i=0;i<10;i++) {

    ostringstream os;
    os << i;

    string obsnam = "JetMult_Pt_" + os.str() ;
    int Obs = histoManager.FindNVar(obsnam);
    int Nhisto = histoManager.access(Obs, NDS);
    TH1F* histo = (TH1F*)Histos[ Nhisto ]->Clone();

    float integ = histo->Integral(0,50);
    float pass = histo->Integral(2,50);
    Efficiencies->SetBinContent(i, 1. - pass/integ );


    int Ndata = histoManager.access(Obs, nt);
    TH1F* hdata = (TH1F*)Histos[ Ndata ]->Clone();
    
    float integD = hdata->Integral(0,50);
    float passD = hdata->Integral(2,50);
    EffData->SetBinContent(i, 1. - passD/integD );
    EffData->SetBinError(i, HistoManager::BinomError(integD, (double)(1.-passD/integD)) );
  }


  Efficiencies->SetLineColor(kRed+1);
  Efficiencies->SetLineWidth(2);
  EffData->SetLineWidth(2);


  TH1F* emptyHisto= new TH1F("bidon","bidon",100,15,100);
  emptyHisto->Draw();
  emptyHisto->GetYaxis()->SetRangeUser(0.7,1.02);
  emptyHisto->GetYaxis()->SetTitle("Jet Veto Efficiency");
  emptyHisto->GetXaxis()->SetTitle("Jet Pt Threshold [GeV]");

  Efficiencies->Draw("same");
  EffData->Draw("same ");

  TCanvas* c5=new TCanvas("JetMultRatio","JetMultRatio",600,300);
  TH1F* emptyHisto2= new TH1F("bidon2","bidon",100,15,100);
  emptyHisto2->Draw();
  emptyHisto2->GetYaxis()->SetRangeUser(0.95,1.05);
  emptyHisto2->GetXaxis()->SetTitle("Jet Pt Threshold [GeV]");
  emptyHisto2->GetYaxis()->SetTitle("Data/MC");

  TLine* line = new TLine(15,1,100,1);
  line->SetLineColor(kRed+1);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw();

  TH1F* ratio = (TH1F*)EffData->Clone();
  ratio->Divide(Efficiencies);
  ratio->Draw("same");

  

}



void 
ZZUseTree::CheckCharge(string chan) {

  TCanvas* c4=new TCanvas("Charge","Charge",600,600);

  int NDS= GetNbyName(chan);

  TH1F* Efficiencies = new TH1F("Efficiencies","Efficiencies",10,0,200);
  TH1F* EffData = new TH1F("EffData","EffData",10,0,200);
  
  
  string obsnam = "ChargeVsPt" ;
  int Obs = histoManager.FindNVar2D(obsnam);
  int Nhisto = histoManager.access2D(Obs, NDS);
  TH2F* histo = (TH2F*)Histos2D[ Nhisto ]->Clone();

  int Ndata = histoManager.access2D(Obs, nt);
  TH2F* hdata = (TH2F*)Histos2D[ Ndata ]->Clone();

  for(int ix=0;ix<hdata->GetNbinsX();ix++) {

    float integ = histo->Integral(ix,ix+1, 0,3);
    float pass = histo->Integral(ix,ix+1,0,1);
    Efficiencies->SetBinContent(ix, pass/integ );

    float integD = hdata->Integral(ix,ix+1, 0,3);
    float passD = hdata->Integral(ix,ix+1,0,1);
    EffData->SetBinContent(ix, passD/integD );
    EffData->SetBinError(ix, HistoManager::BinomError(integD, (double)(1.-passD/integD)) );
    
  }
  
  Efficiencies->SetLineColor(kRed+1);
  Efficiencies->SetLineWidth(2);
  EffData->SetLineWidth(2);


  TH1F* emptyHisto= new TH1F("bidon","bidon",100,0,200);
  emptyHisto->Draw();
  emptyHisto->GetYaxis()->SetRangeUser(0.7,1.02);
  emptyHisto->GetYaxis()->SetTitle("Jet Veto Efficiency");
  emptyHisto->GetXaxis()->SetTitle("Jet Pt Threshold [GeV]");

  Efficiencies->Draw("same");
  EffData->Draw("same ");

  TCanvas* c5=new TCanvas("JetMultRatio","JetMultRatio",600,300);
  TH1F* emptyHisto2= new TH1F("bidon2","bidon",100,0,200);
  emptyHisto2->Draw();
  emptyHisto2->GetYaxis()->SetRangeUser(0.95,1.05);
  emptyHisto2->GetXaxis()->SetTitle("Z Pt [GeV]");
  emptyHisto2->GetYaxis()->SetTitle("Data/MC");

  TLine* line = new TLine(15,1,100,1);
  line->SetLineColor(kRed+1);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw();

  TH1F* ratio = (TH1F*)EffData->Clone();
  ratio->Divide(Efficiencies);
  ratio->Draw("same");

  

}





void 
ZZUseTree::CheckBtag() {

  TCanvas* c4=new TCanvas("Btag","Btag",600,600);

  int NDS= GetNbyName("t#bar{t}");

  float ThrJetMult[10]={15,20,25,30,35,40,50,60,80,100};

  TH1F* Efficiencies = new TH1F("Efficiencies","Efficiencies",9,ThrJetMult);
  TH1F* EffData = new TH1F("EffData","EffData",9,ThrJetMult);
  
  float nP=0,nF=0;
  float nPd=0,nFd=0;

  for(int i=9;i>-1;i--) {

    ostringstream os;
    os << i;

    string obsnam = "BTag_Jet_" + os.str() ;
    int Obs = histoManager.FindNVar(obsnam);
    int Nhisto = histoManager.access(Obs, NDS);
    TH1F* histo = (TH1F*)Histos[ Nhisto ]->Clone();

    float integ = histo->Integral(0,3);
    float pass = histo->Integral(2,3);

    //   cout<<i<<"  "<<integ<<"   "<<nP+nF<<"   "<<pass<<"   "<<nP<<endl;

    integ -= nP+nF;
    pass -= nP;

    Efficiencies->SetBinContent(9-i, pass/integ );
    Efficiencies->SetBinError(9-i, HistoManager::BinomError(integ/0.00485779, pass/integ ) );

    nP += pass;
    nF += integ -pass;
    
    int Ndata = histoManager.access(Obs, nt);
    TH1F* hdata = (TH1F*)Histos[ Ndata ]->Clone();
    
    float integD = hdata->Integral(0,3);
    float passD = hdata->Integral(2,3);

    // cout<<i<<"  "<<integD<<"   "<<nPd+nFd<<"   "<<passD<<"   "<<nPd<<"   --->  "<<integD-nPd-nFd<<"    "<<passD-nPd<<endl;

    integD -= nPd+nFd;
    passD -= nPd;
    
    // cout<<passD/integD<<endl;

    if(integD!=0) {
      EffData->SetBinContent(9-i, passD/integD );
      EffData->SetBinError(9-i, HistoManager::BinomError(integD, (double)(1.-passD/integD)) );
    }
    else {
      EffData->SetBinContent(9-i, 1 );
    }

    nPd += passD;
    nFd += integD -passD;


  }


  Efficiencies->SetLineColor(kRed+1);
  Efficiencies->SetLineWidth(2);
  EffData->SetLineWidth(2);


  TH1F* emptyHisto= new TH1F("bidon","bidon",100,15,100);
  emptyHisto->Draw();
  emptyHisto->GetYaxis()->SetRangeUser(0.,1.02);
  emptyHisto->GetYaxis()->SetTitle("BTag Efficiency");
  emptyHisto->GetXaxis()->SetTitle("Jet Pt Threshold [GeV]");

  Efficiencies->Draw("same");
  EffData->Draw("same ");

  TCanvas* c5=new TCanvas("JetMultRatio","JetMultRatio",600,300);
  TH1F* emptyHisto2= new TH1F("bidon2","bidon",100,15,100);
  emptyHisto2->Draw();
  emptyHisto2->GetYaxis()->SetRangeUser(0.5,2.0);
  emptyHisto2->GetXaxis()->SetTitle("Jet Pt Threshold [GeV]");
  emptyHisto2->GetYaxis()->SetTitle("Data/MC");

  TLine* line = new TLine(15,1,100,1);
  line->SetLineColor(kRed+1);
  line->SetLineWidth(2);
  line->SetLineStyle(7);
  line->Draw();

  TH1F* ratio = (TH1F*)EffData->Clone();
  ratio->Divide(Efficiencies);
  ratio->Draw("same");

  

}


void
ZZUseTree::LoadDBWeightZZ() {

  ifstream in("/home/mmarionn/Documents/CMS/CMSDatabase/kFactorDatabase", ios::in );  

  vector<vector<float> > tmpw;//(18,vector<float>(3,0));

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

}

void
ZZUseTree::LoadDBWeightWZ() {

  ifstream in("/home/mmarionn/Documents/CMS/CMSDatabase/kFactorDatabaseWZ", ios::in );  

  vector<vector<float> > tmpw;//(18,vector<float>(3,0));

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

}

void
ZZUseTree::LoadDBWeightWW() {

  ifstream in("/home/mmarionn/Documents/CMS/CMSDatabase/kFactorDatabaseWW", ios::in );  

  vector<vector<float> > tmpw;//(18,vector<float>(3,0));

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

}



float 
ZZUseTree::SearchWeightZZ(float Zpt) {
  
  float w=1;
  
  for(int unsigned i=0;i<DBWeightsZZ.size();i++) {
    if(i!=DBWeightsZZ.size()-1) {
      if(Zpt >= DBWeightsZZ[i][0] && Zpt < DBWeightsZZ[i+1][0])
	{ w = DBWeightsZZ[i][1]; break;}
    }
    else
      { w = DBWeightsZZ[i][1]; }
  }
  
  return w;
}


float 
ZZUseTree::SearchWeightWZ(float Zpt) {
  
  float w=1;
  
  for(int unsigned i=0;i<DBWeightsZZ.size();i++) {
    if(i!=DBWeightsWZ.size()-1) {
      if(Zpt >= DBWeightsWZ[i][0] && Zpt < DBWeightsWZ[i+1][0])
	{ w = DBWeightsWZ[i][1]; break;}
    }
    else
      { w = DBWeightsWZ[i][1]; }
  }
  
  return w;
}


float 
ZZUseTree::SearchWeightWW(float Wpt) {
  
  float w=1;
  
  for(int unsigned i=0;i<DBWeightsWW.size();i++) {
    if(i!=DBWeightsWW.size()-1) {
      if(Wpt >= DBWeightsWW[i][0] && Wpt < DBWeightsWW[i+1][0])
	{ w = DBWeightsWW[i][1]; break;}
    }
    else
      { w = DBWeightsWW[i][1]; }
  }
  
  return w;
}




void 
ZZUseTree::StoreSelectedEvents(float mZZgen, float w) {

  std::pair<float, float> p(mZZgen, w);
  selZZ.push_back(p);
} 




vector<vector<float> >
ZZUseTree::GetReweightedYields(float f4M, float f5M, int nf4, int nf5) {

  vector<TF2*> curves= LoadDBaTGC();
  LoadATGCWeights();
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
    // cout<<" -> "<<f4v<<"    "<<(float)f4M/nf4<<"    "<<f4M/nf4<<endl;
    f4v = f4v*100000; //0.01*100000 = 1000
    f4v = f4v - (int)f4v%1000;
    f4v = (int)f4v/1000;
    f4v = f4v/100.;
   

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
      for(int unsigned ie=0;ie<selZZ.size();ie++) {
	
	nbmzz = (int)(selZZ[ie].first)/20;
	if(nbmzz >= 100 )
	  { nbmzz=99; }
	//	cout<<"   "<<nbmzz<<endl;
	sumW +=	GetATGCWeight(f4v, f5v, nbmzz);
	//cout<<GetATGCWeight(f4v, f5v, nbmzz)<<endl;
	/*nbmzz = (int)(selZZ[ie].first)/20;
	if(nbmzz < 6 )
	  sumW += selZZ[ie].second;
	else {
	  if(nbmzz > 100 )
	    { nbmzz=99; }
	  
	  nbmzz -=6;
	  
	  cout<<curves[ nbmzz ]->Eval(f4v,f5v)<<endl;
	  sumW += (selZZ[ie].second)* curves[ nbmzz ]->Eval(f4v,f5v);
	}*/	
      }
      vector<float> vtmp(3,0);
      vtmp[0] = f4v;
      vtmp[1] = f5v;
      vtmp[2] = sumW/selZZ.size();
      vRatio.push_back(vtmp);
      //cout<<f4v<<"\t"<<f5v<<"\t  --->  "<<sumW<<"\t "<<sumW/selZZ.size()<<endl;
    }
  }

  return vRatio;

}

vector<TF2*>
ZZUseTree::LoadDBaTGC() {

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


bool
ZZUseTree::BasicAcc(float eta, float pt ) {

  if(pt < 20 ) return false;
  if(pdgID==11) 
    {
      if( fabs(eta)>2.5 || (fabs(eta)>1.44 && fabs(eta)<1.56) ) return false;
    }
  else {
    if( fabs(eta)>2.1) return false;
  }

  
  return true;
}

bool 
ZZUseTree::isInAcc(float genLep[2][5], float genZnn[5], float lPart[5]) {


  TLorentzVector l1, l2, Znn, lP;
  //Zll.SetPtEtaPhiM(genZll[0], genZll[1], genZll[2], genZll[3]);
  Znn.SetPtEtaPhiM(genZnn[0], genZnn[1], genZnn[2], genZnn[3]);
  lP.SetPtEtaPhiE(lPart[0], lPart[1], lPart[2], lPart[3]);
 
  if( !BasicAcc( genLep[0][1], genLep[0][0]) ||
      !BasicAcc( genLep[1][1], genLep[1][0]) ) return false;
  if( Znn.Pt() < 20 ) return false;
  
  if( lP.Pt() > 30 ) return false;

  return true;
}



void 
ZZUseTree::LoadATGCWeights() {


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
  else { cout<<" grou "<<endl; abort();}
  aTGCWeights = mapWeight;

}


float 
ZZUseTree::GetATGCWeight(float f4, float f5, int ibin) {

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
ZZUseTree::DrawContour(float lim95, float lim68, float limU) {

  vector<vector<float> > ratios = GetReweightedYields(0.1,0.1,20,20);

  vector< std::pair<TVector2, float> > bV;
  vector< std::pair<TVector2, float> > bV68;
  //init
  for(int i=0;i<180;i++) {
    TVector2 t(0,0);
    bV.push_back( std::pair<TVector2, float>(t,0));
    bV68.push_back( std::pair<TVector2, float>(t,0));
  }
  

  for(int i=0;i<ratios.size(); i++) {

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

 for(int i=0;i<ratios.size(); i++) {

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
    fitter->SetParameter(2, "a",      0.01, 0.1, 0,0);
    fitter->SetParameter(3, "b",      0.01, 0.1, 0,0);
    fitter->SetParameter(4, "theta", 1.7, 0.1, 0,3.30);
    // fitter->FixParameter(0);
    // fitter->FixParameter(1);

    Double_t arglist[1] = {0};
    fitter->ExecuteCommand("MIGRAD", arglist, 0);
   
    //Getparameters
    float a = fitter->GetParameter(2)/(1-fitter->GetParameter(3)*fitter->GetParameter(3));

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
   fitter->SetParameter(2, "a",      0.01, 0.1, 0,0);
   fitter->SetParameter(3, "b",      0.01, 0.1, 0,0);
   fitter->SetParameter(4, "theta", 1.7, 0.1, 0,3.30);
   // fitter->FixParameter(0);
   // fitter->FixParameter(1);

   Double_t arglist[1] = {0};
   fitter->ExecuteCommand("MIGRAD", arglist, 0);
   
   //Getparameters
   float a = fitter->GetParameter(2)/(1-fitter->GetParameter(3)*fitter->GetParameter(3));

   //Draw the circle on top of the points
   ell95 = new TEllipse(fitter->GetParameter(0),
			fitter->GetParameter(1),
			fitter->GetParameter(2),
			fitter->GetParameter(3), 
			0, 360,
			fitter->GetParameter(4) );
   ell95->SetLineColor(kRed);
   ell95->SetLineWidth(3);
   ell95->SetLineStyle(7);
   ell95->SetFillStyle(0);
   
 
  }
  cout<<" grou4 "<<endl;
  contour->GetXaxis()->SetLimits(-0.15,0.15);
  contour->GetYaxis()->SetLimits(-0.15,0.15);
  contour->GetXaxis()->SetTitle("f4_{Z}");
  contour->GetYaxis()->SetTitle("f5_{Z}");
  contour->Draw("ap");
  contour->GetYaxis()->SetRangeUser(-0.15,0.15);
  ell95->Draw("same");
  ell68->Draw("same");

}


void
ZZUseTree::DrawATGCYields() {


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
