#include "ZUseTree.hh"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooPlot.h>

using namespace RooFit;
using namespace std;

ClassImp(ZUseTree)


ZUseTree::ZUseTree():
UseTree()
{

  LoadDBWeight();
  LoadResponseDB();

  float ZIsoCuts[2][4][3]= { { {0.05, 0.06, 0.03}, {0.09, 0.07, 0.1},
			       {0.12, 0.09, 0.1}, {0.15, 2.0, 0.12} },
			     { {0.025, 0.025, 0.02}, {0.04, 0.05, 0.025},
			       {0.05, 0.06, 0.03}, {0.08, 0.06, 0.05} } };
  
  float ZIDCuts[2][4][4]= { { {0.01,0.004,0.03,0.025}, {0.01,0.004,0.06,0.04},
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
	  _IsoCuts[j][i][k] = ZIsoCuts[j][i][k];
	_IDCuts[j][i][k] = ZIDCuts[j][i][k];
      }
    }
  }

  
  histoManager.ConfigureContamination(0.05);

}

void
ZUseTree::SpecialZVariables(bool Response) {
  response = Response;
}

void
ZUseTree::FillZTree() {

  //Add Varaibles to be computed

  //Prepare suffix for outfiles
  string suffix = "METs";
  if(switchRMS && FillProf)
    suffix = "RMS_METs";

  if(Draw3on1)
    for(int i=0;i<6;i++)
      Suffix[i] = suffix;

  string metT[6]={"pf","tc","calo","caloT1","caloT2","pfT1"};
  string lep[2]={"l1","l2"};

  cout<<" Now book histograms "<<endl;
  
  //Vertex
  histoManager.AddVariable("Vertex",10,0,10,"number of vertices","Vertex");

  //Eta Phi Pt
    for(int ii=0;ii<2;ii++) {
      histoManager.AddVariable("Eta"+lep[ii],100,-3,3,"#eta "+lep[ii],"Eta"+lep[ii]);
    histoManager.AddVariable("Phi"+lep[ii],128,0,3.2,"#phi "+lep[ii],"Phi"+lep[ii]);
    histoManager.AddVariable("Pt"+lep[ii],200,0,100,lep[ii]+" p_{T} [GeV]","Pt"+lep[ii]);
    histoManager.AddVariable("PfPt"+lep[ii],200,0,100,lep[ii]+"pf p_{T} [GeV]","PfPt"+lep[ii]);
    histoManager.AddVariable("PfPtCor"+lep[ii],200,0,100,lep[ii]+"pf cor p_{T} [GeV]","PfPtCor"+lep[ii]);
    histoManager.AddVariable("EtSC"+lep[ii],200,0,100,lep[ii]+" SC E_{T} [GeV]","EtSC"+lep[ii]);
    histoManager.AddVariable("EtSCCor"+lep[ii],200,0,100,lep[ii]+"SC cor E_{T} [GeV]","EtSCCor"+lep[ii]);
    
    //ID
    histoManager.AddVariable("sigieie"+lep[ii],100,0,0.1,"#sigmai#etai#eta ("+lep[ii]+")","sigieie"+lep[ii]);
    histoManager.AddVariable("Deta"+lep[ii],100,-0.05,0.05,"#Delta#eta ("+lep[ii]+")","Deta"+lep[ii]);
    histoManager.AddVariable("Dphi"+lep[ii],100,-0.2,0.2,"#Delta#Phi ("+lep[ii]+")","Dphi"+lep[ii]);
    histoManager.AddVariable("HoE"+lep[ii],100,0,0.1,"H/E ("+lep[ii]+")","HoE"+lep[ii]);
    
    //iso
    histoManager.AddVariable("TrackIso"+lep[ii],500,0,1,"TrackIso ("+lep[ii]+")","TrackIso"+lep[ii]);
    histoManager.AddVariable("EcalIso"+lep[ii],500,0,1,"EcalIso ("+lep[ii]+")","EcalIso"+lep[ii]);
    histoManager.AddVariable("HcalIso"+lep[ii],500,0,1,"HcalIso ("+lep[ii]+")","HcalIso"+lep[ii]);
    
  }
    
  //Underlying event variables
  histoManager.AddVariable("dZv",200,0,100,"#Delta z [cm]","Zee_dZ_vertices");
  histoManager.AddVariable("genZPt",400,0,200,"q_{T} Z generated [GeV]","ZGen_Pt");
  histoManager.AddVariable("jetPt",400,0,200,"q_{T} jet [GeV]","Zee_Pt_Jet");
  histoManager.AddVariable("photonPt",400,0,200,"p_{T} photon [GeV]","Zee_Pt_Photon");
  histoManager.AddVariable("photonMult",10,0,10,"N photons","Zee_Photon_Mult");
  histoManager.AddVariable("jetMult",10,0,10,"N jets","Zee_Jet_Mult");

  histoManager.AddVariable("dPhiJetZ",100,0,3.14,"#Delta#phi(jet,Z) [rad]","Zee_dPhi_Jet");
  histoManager.AddVariable("dPhiPhotonZ",100,0,3.14,"#Delta#phi(jet,Z) [rad]","Zee_dPhi_Photon");

  histoManager.AddVariable("NPhot1",20,0,20,"nPhot1","Zee_NPot1_");
  histoManager.AddVariable("NPhot2",20,0,20,"nPhot1","Zee_NPot1_");
 
  //Underlying variables and quality variables
  histoManager.AddVariable("ZEta",100,-3,3,"#eta Z","Zee_Eta");
  histoManager.AddVariable("ZPhi",128,-4,4,"Z phi ","Zee_phi");
  histoManager.AddVariable("ZPt",400,0,200,"q_{T} Z [GeV]","Zee_PT");
  histoManager.AddVariable("ZPtUnc",400,0,200,"q_{T} Z [GeV]","Zee_PT");
  histoManager.AddVariable("ZPtCor",400,0,200,"cor p_{T} Z [GeV]","Zee_PT");
  histoManager.AddVariable("ZMass",400,0,200,"M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassUnc",400,0,200,"M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassCor",400,0,200,"cor M_{ee} [GeV]","Zee_MassCor");
  histoManager.AddVariable("ZPtControl",400,0,200,"p_{T} Z [GeV]","Zee_PT");
  histoManager.AddVariable("ZPtCorControl",400,0,200,"cor p_{T} Z [GeV]","Zee_PT");
  histoManager.AddVariable("ZMassControl",400,0,200,"M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassCorControl",400,0,200,"cor M_{ee} [GeV]","Zee_MassCor");
 
  //ZMass special Laser =================================
   histoManager.AddVariable("ZMassEBEEmInner",400,0,200,"EE- Inner M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassEBEEmOuter",400,0,200,"EE- Outer M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassEBEEpInner",400,0,200,"EE+ Inner M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassEBEEpOuter",400,0,200,"EE+ Outer M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassEBEEp",400,0,200,"EE+ M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassEBEEm",400,0,200,"EE- M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassEBEEmRight",400,0,200,"EE- Right M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassEBEEmLeft",400,0,200,"EE- Left M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassEBEEpRight",400,0,200,"EE+ Right M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassEBEEpLeft",400,0,200,"EE+ Left M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassEBEBp",400,0,200,"BB+ M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassEBEBm",400,0,200,"BB- M_{ee} [GeV]","Zee_Mass");

  histoManager.AddVariable("ZMassSC",400,0,200,"cor M_{ee} [GeV]","Zee_MassSC");
  histoManager.AddVariable("ZMassSCUnCor",400,0,200,"uncor M_{ee} [GeV]","Zee_MassSCUnCor");

  //=====================================================

  histoManager.AddVariable("ZPt_V1",400,0,200,"p_{T} Z [GeV]","Zee_PT");
  histoManager.AddVariable("ZPt_V2",400,0,200,"p_{T} Z [GeV]","Zee_PT");
  
  //SumET and MET decomposition
    histoManager.AddVariable("caloMETEB",400,0,200,"calo#slash{E}_{T} EB [GeV]","Zee_caloMET_EB");
  histoManager.AddVariable("caloMETEE",400,0,200,"calo#slash{E}_{T} EE [GeV]","Zee_caloMET_EE");
  histoManager.AddVariable("caloMETHB",400,0,200,"calo#slash{E}_{T} HB [GeV]","Zee_caloMET_HB");
  histoManager.AddVariable("caloMETHE",400,0,200,"calo#slash{E}_{T} HE [GeV]","Zee_caloMET_HE");
  histoManager.AddVariable("caloMETHF",400,0,200,"calo#slash{E}_{T} HF [GeV]","Zee_caloMET_HF");

  histoManager.AddVariable("caloMETPhiEB",360,0,360,"calo#slash{E}_{T} #Phi EB [°]","Zee_caloMETPhi_EB");
  histoManager.AddVariable("caloMETPhiEE",360,0,360,"calo#slash{E}_{T} #Phi EE [°]","Zee_caloMETPhi_EE");
  histoManager.AddVariable("caloMETPhiHB",360,0,360,"calo#slash{E}_{T} #Phi HB [°]","Zee_caloMETPhi_HB");
  histoManager.AddVariable("caloMETPhiHE",360,0,360,"calo#slash{E}_{T} #Phi HE [°]","Zee_caloMETPhi_HE");
  histoManager.AddVariable("caloMETPhiHF",360,0,360,"calo#slash{E}_{T} #Phi HF [°]","Zee_caloMETPhi_HF");

  histoManager.AddVariable("caloSumETEB",500,0,500,"calo #Sigma E_{T} EB [GeV]","Zee_caloSumET_EB");
  histoManager.AddVariable("caloSumETEE",500,0,500,"calo #Sigma E_{T} EE [GeV]","Zee_caloSumET_EE");
  histoManager.AddVariable("caloSumETHB",500,0,500,"calo #Sigma E_{T} HB [GeV]","Zee_caloSumET_HB");
  histoManager.AddVariable("caloSumETHE",500,0,500,"calo #Sigma E_{T} HE [GeV]","Zee_caloSumET_HE");
  histoManager.AddVariable("caloSumETHF",500,0,500,"calo #Sigma E_{T} HF [GeV]","Zee_caloSumET_HF");

  histoManager.AddVariable("caloSumETEBControl",500,0,500,"calo #Sigma E_{T} EB [GeV]","Zee_caloSumET_EB");
  histoManager.AddVariable("caloSumETEEControl",500,0,500,"calo #Sigma E_{T} EE [GeV]","Zee_caloSumET_EE");
  histoManager.AddVariable("caloSumETHBControl",500,0,500,"calo #Sigma E_{T} HB [GeV]","Zee_caloSumET_HB");
  histoManager.AddVariable("caloSumETHEControl",500,0,500,"calo #Sigma E_{T} HE [GeV]","Zee_caloSumET_HE");
  histoManager.AddVariable("caloSumETHFControl",500,0,500,"calo #Sigma E_{T} HF [GeV]","Zee_caloSumET_HF");

  histoManager.AddVariable("pfMETEB",400,0,200,"pf#slash{E}_{T} EB [GeV]","Zee_pfMET_EB");
  histoManager.AddVariable("pfMETEE",400,0,200,"pf#slash{E}_{T} EE [GeV]","Zee_pfMET_EE");
  histoManager.AddVariable("pfMETHB",400,0,200,"pf#slash{E}_{T} HB [GeV]","Zee_pfMET_HB");
  histoManager.AddVariable("pfMETHE",400,0,200,"pf#slash{E}_{T} HE [GeV]","Zee_pfMET_HE");
  histoManager.AddVariable("pfMETHF",400,0,200,"pf#slash{E}_{T} HF [GeV]","Zee_pfMET_HF");

  histoManager.AddVariable("pfMETPhiEB",360,0,360,"pf#slash{E}_{T} #Phi EB [#circ]","Zee_pfMETPhi_EB");
  histoManager.AddVariable("pfMETPhiEE",360,0,360,"pf#slash{E}_{T} #Phi EE [#circ]","Zee_pfMETPhi_EE");
  histoManager.AddVariable("pfMETPhiEEp",360,0,360,"pf#slash{E}_{T} #Phi EE+ [#circ]","Zee_pfMETPhi_EE");
  histoManager.AddVariable("pfMETPhiEEm",360,0,360,"pf#slash{E}_{T} #Phi EE- [#circ]","Zee_pfMETPhi_EE");
  histoManager.AddVariable("pfMETPhiHB",360,0,360,"pf#slash{E}_{T} #Phi HB [#circ]","Zee_pfMETPhi_HB");
  histoManager.AddVariable("pfMETPhiHE",360,0,360,"pf#slash{E}_{T} #Phi HE [#circ]","Zee_pfMETPhi_HE");
  histoManager.AddVariable("pfMETPhiHF",360,0,360,"pf#slash{E}_{T} #Phi HF [#circ]","Zee_pfMETPhi_HF");

  histoManager.AddVariable("pfSumETEB",500,0,500,"pf #Sigma E_{T} EB [GeV]","Zee_pfSumET_EB");
  histoManager.AddVariable("pfSumETEE",500,0,500,"pf #Sigma E_{T} EE [GeV]","Zee_pfSumET_EE");
  histoManager.AddVariable("pfSumETHB",500,0,500,"pf #Sigma E_{T} HB [GeV]","Zee_pfSumET_HB");
  histoManager.AddVariable("pfSumETHE",500,0,500,"pf #Sigma E_{T} HE [GeV]","Zee_pfSumET_HE");
  histoManager.AddVariable("pfSumETHF",500,0,500,"pf #Sigma E_{T} HF [GeV]","Zee_pfSumET_HF");

  histoManager.AddVariable("pfSumETEBControl",500,0,500,"pf #Sigma E_{T} EB [GeV]","Zee_pfSumET_EB");
  histoManager.AddVariable("pfSumETEEControl",500,0,500,"pf #Sigma E_{T} EE [GeV]","Zee_pfSumET_EE");
  histoManager.AddVariable("pfSumETHBControl",500,0,500,"pf #Sigma E_{T} HB [GeV]","Zee_pfSumET_HB");
  histoManager.AddVariable("pfSumETHEControl",500,0,500,"pf #Sigma E_{T} HE [GeV]","Zee_pfSumET_HE");
  histoManager.AddVariable("pfSumETHFControl",500,0,500,"pf #Sigma E_{T} HF [GeV]","Zee_pfSumET_HF");

  histoManager.Add2DVariable("pfMETETPhiEB","pf#slash{E}_{T} #Phi EB [#circ]","#slash{E}_{T} [GeV]","Zee_pfMETEtPhi_EB",10,0,360,12,0,60);
  histoManager.Add2DVariable("pfMETETPhiEE","pf#slash{E}_{T} #Phi EE [#circ]","#slash{E}_{T} [GeV]","Zee_pfMETEtPhi_EE",10,0,360,12,0,60);
  histoManager.Add2DVariable("pfMETETPhiHB","pf#slash{E}_{T} #Phi HB [#circ]","#slash{E}_{T} [GeV]","Zee_pfMETEtPhi_HB",10,0,360,12,0,60);
  histoManager.Add2DVariable("pfMETETPhiHE","pf#slash{E}_{T} #Phi EE [#circ]","#slash{E}_{T} [GeV]","Zee_pfMETEtPhi_HE",10,0,360,12,0,60);
  histoManager.Add2DVariable("pfMETETPhiHF","pf#slash{E}_{T} #Phi HF [#circ]","#slash{E}_{T} [GeV]","Zee_pfMETEtPhi_HF",10,0,360,12,0,60);

  histoManager.Add2DVariable("caloMETETPhiEB","calo#slash{E}_{T} #Phi EB [#circ]","#slash{E}_{T} [GeV]","Zee_caloMETEtPhi_EB",10,0,360,12,0,60);
  histoManager.Add2DVariable("caloMETETPhiEE","calo#slash{E}_{T} #Phi EE [#circ]","#slash{E}_{T} [GeV]","Zee_caloMETEtPhi_EE",10,0,360,12,0,60);
  histoManager.Add2DVariable("caloMETETPhiHB","calo#slash{E}_{T} #Phi HB [#circ]","#slash{E}_{T} [GeV]","Zee_caloMETEtPhi_HB",10,0,360,12,0,60);
  histoManager.Add2DVariable("caloMETETPhiHE","calo#slash{E}_{T} #Phi EE [#circ]","#slash{E}_{T} [GeV]","Zee_caloMETEtPhi_HE",10,0,360,12,0,60);
  histoManager.Add2DVariable("caloMETETPhiHF","calo#slash{E}_{T} #Phi HF [#circ]","#slash{E}_{T} [GeV]","Zee_caloMETEtPhi_HF",10,0,360,12,0,60);

  histoManager.AddVariable("pfMETEcal",400,0,200,"pf#slash{E}_{T} Ecal [GeV]","Zee_pfMET_Ecal");
  histoManager.AddVariable("pfMETHcal",400,0,200,"pf#slash{E}_{T} Hcal [GeV]","Zee_pfMET_Hcal");
  histoManager.AddVariable("pfMETAll",400,0,200,"pf#slash{E}_{T} [GeV]","Zee_pfMET_All");
  histoManager.AddVariable("caloMETAll",400,0,200,"calo#slash{E}_{T} [GeV]","Zee_caloMET_All");
  histoManager.AddVariable("caloMETEcal",400,0,200,"calo#slash{E}_{T} Ecal [GeV]","Zee_caloMET_Ecal");
  histoManager.AddVariable("caloMETHcal",400,0,200,"calo#slash{E}_{T} Hcal [GeV]","Zee_caloMET_Hcal");

  histoManager.AddVariable("pfMETPhiEcal",360,0,360,"pf#slash{E}_{T} #Phi Ecal [#circ]","Zee_pfMETPhi_Ecal");
  histoManager.AddVariable("pfMETPhiHcal",360,0,360,"pf#slash{E}_{T} #Phi Hcal [#circ]","Zee_pfMETPhi_Hcal");
  histoManager.AddVariable("caloMETPhiEcal",360,0,360,"calo#slash{E}_{T} #Phi Ecal [#circ]","Zee_caloMETPhi_Ecal");
  histoManager.AddVariable("caloMETPhiHcal",360,0,360,"calo#slash{E}_{T} #Phi Hcal [#circ]","Zee_caloMETPhi_Hcal");
  
  histoManager.AddVariable("pfSumETHF_V1",500,0,500,"pf #Sigma E_{T} HF [GeV]","Zee_pfSumET_HF");
  histoManager.AddVariable("pfSumETHF_V2",500,0,500,"pf #Sigma E_{T} HF [GeV]","Zee_pfSumET_HF");
  
  histoManager.AddVariable("diffPfTcMETX",400,-5,5,"(pf-tc)/pf MET X [GeV]","Zee_difftcpfX");
  histoManager.AddVariable("diffPfTcMETY",400,-5,5,"(pf-tc)/pf MET Y [GeV]","Zee_difftcpfY");
  histoManager.AddVariable("diffPfTcMETSFL",400,-5,5,"(pf-tc)/pf MET SframeL [GeV]","Zee_difftcpfSFL");
  histoManager.AddVariable("diffPfTcMETSFT",400,-5,5,"(pf-tc)/pf MET SframeT [GeV]","Zee_difftcpfSFT");

  //PF frac
  histoManager.AddVariable("pfEmNFrac",100,0,1,"pf EM N frac","Zee_pfEMNFrac");
  histoManager.AddVariable("pfEmCFrac",100,0,1,"pf EM C frac","Zee_pfEMCFrac");
  histoManager.AddVariable("pfMuCFrac",100,0,1,"pf Mu C frac","Zee_pfMuCFrac");
  histoManager.AddVariable("pfHadCFrac",100,0,1,"pf Had C frac","Zee_pfHadCFrac");
  histoManager.AddVariable("pfHadNFrac",100,0,1,"pf Had N frac","Zee_pfHadNFrac");
  
 
  for(int kk=0;kk<6;kk++) {

    //METs 
    histoManager.AddVariable(metT[kk]+"MET",400,0,200,metT[kk]+" #slash{E}_{T} [GeV]","Zee_MET_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METPhi",360,0,360,metT[kk]+" #slash{E}_{T} #Phi ","Zee_METPhi_"+Suffix[kk]);
    //METProj
    histoManager.AddVariable(metT[kk]+"METPara",400,-100,100,metT[kk]+" #slash{E}_{T,||} [GeV]","Zee_MET_Para_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METPerp",400,-100,100,metT[kk]+" #slash{E}_{T,#perp}   [GeV]","Zee_MET_Perp_"+Suffix[kk]);
  
    histoManager.AddVariable(metT[kk]+"METUpper",400,0,200,metT[kk]+" #slash{E}_{T} upper [GeV]","Zee_METUpper_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METLower",400,0,200,metT[kk]+" #slash{E}_{T} lower [GeV]","Zee_METLower_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METSFramePara",400,-100,100,metT[kk]+" #slash{E}_{T,||} SF [GeV]","Zee_MET_SFPara_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METSFramePerp",400,-100,100,metT[kk]+" #slash{E}_{T,#perp}  SF [GeV]","Zee_MET_SFPerp_"+Suffix[kk]);

    histoManager.AddVariable(metT[kk]+"METUpper_V1",400,0,200,metT[kk]+" #slash{E}_{T} upper [GeV]","Zee_METUpper_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METLower_V1",400,0,200,metT[kk]+" #slash{E}_{T} lower [GeV]","Zee_METLower_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METSFramePara_V1",400,-100,100,metT[kk]+" #slash{E}_{T,||} SF [GeV]","Zee_MET_SFPara_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METSFramePerp_V1",400,-100,100,metT[kk]+" #slash{E}_{T,#perp}  SF [GeV]","Zee_MET_SFPerp_"+Suffix[kk]);

    histoManager.AddVariable(metT[kk]+"METUpper_V2",400,0,200,metT[kk]+" #slash{E}_{T} upper [GeV]","Zee_METUpper_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METLower_V2",400,0,200,metT[kk]+" #slash{E}_{T} lower [GeV]","Zee_METLower_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METSFramePara_V2",400,-100,100,metT[kk]+" #slash{E}_{T,||} SF [GeV]","Zee_MET_SFPara_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METSFramePerp_V2",400,-100,100,metT[kk]+" #slash{E}_{T,#perp}  SF [GeV]","Zee_MET_SFPerp_"+Suffix[kk]);

    histoManager.AddVariable(metT[kk]+"dPhiMETZ",180,0,180,"#slash{E}_{T} #Delta#phi ("+metT[kk]+") [#circ ]","Zee_dPhi_MET_"+Suffix[kk]);
    //  histoManager.AddVariable(metT[kk]+"dPhiMETZLepCor",100,0,3.14,"#Delta#phi ("+metT[kk]+") [rad]","Zee_dPhi_Recoil_"+Suffix[kk]);

    histoManager.AddVariable(metT[kk]+"METPhiL1",360,-180,180,"dPhi MET L1","Zee_dphiMETL1_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METPhiL2",360,-180,180,"dPhi MET L2","Zee_dphiMETL2_"+Suffix[kk]);

    //Recoil
    histoManager.AddVariable(metT[kk]+"Recoil",400,0,200,"u_{T} ("+metT[kk]+") [GeV]","Zee_Recoil_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilCor",400,0,200,"u_{T} ("+metT[kk]+") [GeV]","Zee_Recoil_"+Suffix[kk]);
  
    histoManager.AddVariable(metT[kk]+"RecoilPhi",360,0,360,metT[kk]+" u_{T} #Phi ","Zee_RecoilPhi_"+Suffix[kk]);

    //RecoilProj
    histoManager.AddVariable(metT[kk]+"RecoilPara",800,-200,200,"u_{||} ("+metT[kk]+") [GeV]","Zee_Recoil_Para_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilPerp",800,-200,200,"u_{#perp}   ("+metT[kk]+") [GeV]","Zee_Recoil_Perp_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilParaCor",400,-200,200,"u_{||} ("+metT[kk]+") [GeV]","Zee_Recoil_Para_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilPerpCor",400,-200,200,"u_{#perp}   ("+metT[kk]+") [GeV]","Zee_Recoil_Perp_"+Suffix[kk]);
    if(kk==0) {
      histoManager.AddVariable("ParaBis",800,-200,200,"u_{||} ("+metT[kk]+") [GeV]","Zee_Recoil_Para_"+Suffix[kk]);
      histoManager.AddVariable("PerpBis",800,-200,200,"u_{#perp}   ("+metT[kk]+") [GeV]","Zee_Recoil_Perp_"+Suffix[kk]);
    }

    //SumET
    histoManager.AddVariable(metT[kk]+"SumEt",500,0,500,"#Sigma E_{T} ("+metT[kk]+") [GeV]","Zee_SumEt_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"SumEtControl",500,0,500,"#Sigma E_{T} ("+metT[kk]+") [GeV]","Zee_SumEt_"+Suffix[kk]);

    //MET X Y 
    histoManager.AddVariable(metT[kk]+"METX",400,-100,100,metT[kk]+" #slash{E}_{T,X} [GeV]","Zee_METX_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METXControl",400,-100,100,metT[kk]+" #slash{E}_{T,X} [GeV]","Zee_METX_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METX_V1",400,-100,100,metT[kk]+" #slash{E}_{T,X} [GeV]","Zee_METX_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METX_V2",400,-100,100,metT[kk]+" #slash{E}_{T,X} [GeV]","Zee_METX_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METY",400,-100,100,metT[kk]+" #slash{E}_{T,Y} [GeV]","Zee_METY_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METYControl",400,-100,100,metT[kk]+" #slash{E}_{T,Y} [GeV]","Zee_METY_"+Suffix[kk]);
  
    histoManager.AddVariable(metT[kk]+"METtt",400,-100,100,metT[kk]+" #slash{E}_{T,t} [GeV]","Zee_METX_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METtl",400,-100,100,metT[kk]+" #slash{E}_{T,l} [GeV]","Zee_METX_"+Suffix[kk]);

    histoManager.AddVariable(metT[kk]+"METXCor",400,-100,100,metT[kk]+" #slash{E}_{T,X} [GeV]","Zee_METX_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METXCorControl",400,-100,100,metT[kk]+" #slash{E}_{T,X} [GeV]","Zee_METX_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METXCor_V1",400,-100,100,metT[kk]+" #slash{E}_{T,X} [GeV]","Zee_METX_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METXCor_V2",400,-100,100,metT[kk]+" #slash{E}_{T,X} [GeV]","Zee_METX_"+Suffix[kk]);

    //Angle
    histoManager.AddVariable(metT[kk]+"dPhiRecoil",100,0,3.14,"#Delta#phi ("+metT[kk]+") [rad]","Zee_dPhi_Recoil_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"dPhiRecoilCor",100,0,3.14,"#Delta#phi ("+metT[kk]+") [rad]","Zee_dPhi_Recoil_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"dPhiRecoilControl",100,0,3.14,"#Delta#phi ("+metT[kk]+") [rad]","Zee_dPhi_Recoil_"+Suffix[kk]);

    //Corrected MET variable
    histoManager.AddVariable(metT[kk]+"diffMet",400,-1,1,metT[kk]+" diff #slash{E}_{T} [GeV]","Zee_METdiff_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METLepCor",400,0,200,metT[kk]+"cor #slash{E}_{T} [GeV]","Zee_METLepCor_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METParaLepCor",400,-200,200,metT[kk]+"cor #slash{E}_{T,||} [GeV]","Zee_MET_ParaLepCor_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METPerpLepCor",400,-200,200,metT[kk]+"cor #slash{E}_{T,#perp}   [GeV]","Zee_MET_PerpLepCor_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilLepCor",400,0,200,"cor u_{T} ("+metT[kk]+") [GeV]","Zee_RecoilLepCor_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilParaLepCor",400,-200,200,"cor u_{||} ("+metT[kk]+") [GeV]","Zee_Recoil_ParaLepCor_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilPerpLepCor",400,-200,200,"cor u_{#perp}   ("+metT[kk]+") [GeV]","Zee_Recoil_PerpLepCor_"+Suffix[kk]);


    float tB[6]={200,0,200,200,0,200};
    tB[0]=100; tB[1]=0; tB[2]=200; tB[3]=100; tB[4]= 0; tB[5]=3.14; 
    //  histoManager.Add2DVariable(metT[kk]+"dPhiRecoil2D","u_{T} ("+metT[kk]+") [GeV]","#Delta#phi ("+metT[kk]+") [rad] ","Zee_dPhi_vs_Recoil_"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    //  histoManager.Add2DVariable(metT[kk]+"dPhiRecoil2DCor","u_{T} ("+metT[kk]+") [GeV]","#Delta#phi ("+metT[kk]+") [rad] ","Zee_dPhi_vs_Recoil_"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
 

    //Ut/Pt
    histoManager.AddVariable(metT[kk]+"RecoilPt",200,0,2,"u_{T}/pt_{lepton} ("+metT[kk]+")","Zee_Recoil_pt_"+Suffix[kk]);
    histoManager.Add2DVariable(metT[kk]+"dPhiRecoilPt2D","u/q_{T} (pf) [GeV]","#Delta#phi(u_{T,"+metT[kk]+"},lepton) [rad] ","Zee_dPhi_vs_RoQ_"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
 
    //MET vs angle
    tB[0]=100; tB[1]=0; tB[2]=200; tB[3]=180; tB[4]= -180; tB[5]=180;
    /* histoManager.Add2DVariable(metT[kk]+"dPhiL22D"," #slash{E}_{T} [GeV]"," angle L2","Zee_dPhiL22D_"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    histoManager.Add2DVariable(metT[kk]+"dPhiL12D"," #slash{E}_{T} [GeV]"," angle L1","Zee_dPhiL12D_"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    histoManager.Add2DVariable(metT[kk]+"dPhiZ2D"," #slash{E}_{T} [GeV]"," angle Z","Zee_dPhiZ2D_"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);*/

    //Ut vs Qt
      tB[0]=100; tB[1]=0; tB[2]=200; tB[3]=200; tB[4]= -300; tB[5]=100; 
      // histoManager.Add2DVariable(metT[kk]+"UtVsQt","q_{T} [GeV]","u_{T} [GeV]","Zee_UtvsQt"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
      //  histoManager.Add2DVariable(metT[kk]+"UtVsQtCor","q_{T} [GeV]","u_{T} [GeV]","Zee_UtvsQt"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
  
    //Upara vs Qt
    histoManager.Add2DVariable(metT[kk]+"UparaVsQt","q_{T} [GeV]","u_{||} [GeV]","Zee_UparaVsQt"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    histoManager.Add2DVariable(metT[kk]+"UparaVsQtCor","q_{T} [GeV]","u_{||} [GeV]","Zee_UparapVsQt"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    histoManager.Add2DVariable(metT[kk]+"UparaVsQtControl","q_{T} [GeV]","u_{||} [GeV]","Zee_UparaVsQt"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    histoManager.Add2DVariable(metT[kk]+"UparaVsQtCorControl","q_{T} [GeV]","u_{||} [GeV]","Zee_UparapVsQt"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
      
    //Uperp vs Qt
    histoManager.Add2DVariable(metT[kk]+"UperpVsQt","q_{T} [GeV]","u_{#perp} [GeV]","Zee_UperpVsQt"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    histoManager.Add2DVariable(metT[kk]+"UperpVsQtCor","q_{T} [GeV]","u_{#perp} [GeV]","Zee_UperpVsQt"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    histoManager.Add2DVariable(metT[kk]+"UperpVsQtControl","q_{T} [GeV]","u_{#perp} [GeV]","Zee_UperpVsQt"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    histoManager.Add2DVariable(metT[kk]+"UperpVsQtCorControl","q_{T} [GeV]","u_{#perp} [GeV]","Zee_UperpVsQt"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
  
    double varBin[33] = {0,5,10,20,30,40,50,60,70,80,90,100,110,120,140,160,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000}; 
    histoManager.AddProfVariable(metT[kk]+"Response",32,varBin,"q_{T} [GeV]","Zee_response"+Suffix[kk]);
    
    double varBin2[34] = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,55,60,65,70,80,100,150,200};
    histoManager.AddProfVariable(metT[kk]+"CorPara",33,varBin2,"q_{T} [GeV]","Zee_CorPara"+Suffix[kk]);
    histoManager.AddProfVariable(metT[kk]+"CorPerp",33,varBin2,"q_{T} [GeV]","Zee_CorPerp"+Suffix[kk]);

    for(int iv=1;iv<5;iv++) {
      ostringstream os;
      os << iv;
      
      histoManager.AddVariable(metT[kk]+"RecoilPara_V"+os.str(),800,-200,200,"u_{||} ("+metT[kk]+") [GeV]","Zee_Recoil_Para_"+Suffix[kk]);
      histoManager.AddVariable(metT[kk]+"RecoilPerp_V"+os.str(),800,-200,200,"u_{#perp}   ("+metT[kk]+") [GeV]","Zee_Recoil_Perp_"+Suffix[kk]);
      histoManager.AddVariable(metT[kk]+"RecoilParaCor_V"+os.str(),400,-200,200,"u_{||} ("+metT[kk]+") [GeV]","Zee_Recoil_Para_"+Suffix[kk]);
      histoManager.AddVariable(metT[kk]+"RecoilPerpCor_V"+os.str(),400,-200,200,"u_{#perp}   ("+metT[kk]+") [GeV]","Zee_Recoil_Perp_"+Suffix[kk]);

      histoManager.AddVariable(metT[kk]+"RecoilX"+os.str(),800,-200,200,"u_{X} ("+metT[kk]+") [GeV]","Zee_Recoil_Para_"+Suffix[kk]);
      histoManager.AddVariable(metT[kk]+"RecoilY"+os.str(),800,-200,200,"u_{Y}   ("+metT[kk]+") [GeV]","Zee_Recoil_Perp_"+Suffix[kk]);

      histoManager.Add2DVariable(metT[kk]+"UperpVsQt_V"+os.str(),"q_{T} [GeV]","u_{#perp} [GeV]","Zee_UperpVsQt_V"+os.str()+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
      histoManager.Add2DVariable(metT[kk]+"UperpVsQtCor_V"+os.str(),"q_{T} [GeV]","u_{#perp} [GeV]","Zee_UperpVsQt_V"+os.str()+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
      histoManager.Add2DVariable(metT[kk]+"UparaVsQt_V"+os.str(),"q_{T} [GeV]","u_{||} [GeV]","Zee_UparaVsQt_V"+os.str()+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
      histoManager.Add2DVariable(metT[kk]+"UparaVsQtCor_V"+os.str(),"q_{T} [GeV]","u_{||} [GeV]","Zee_UparaVsQt_V"+os.str()+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    }

    //METProj
    histoManager.AddVariable(metT[kk]+"METPara_V1",400,-200,200,metT[kk]+" #slash{E}_{T,||} [GeV]","Zee_MET_Para_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METPara_V2",400,-200,200,metT[kk]+" #slash{E}_{T,||} [GeV]","Zee_MET_Para_"+Suffix[kk]);

    histoManager.AddVariable(metT[kk]+"METPerp_V1",400,-200,200,metT[kk]+" #slash{E}_{T,#perp} [GeV]","Zee_MET_Para_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METPerp_V2",400,-200,200,metT[kk]+" #slash{E}_{T,#perp} [GeV]","Zee_MET_Para_"+Suffix[kk]);

    //MET vs SumET
    /*  tB[0]=50; tB[1]=0; tB[2]=100; tB[3]=100; tB[4]= -100; tB[5]=100;
	histoManager.Add2DVariable(metT[kk]+"METXvsSumET","#Sigma E_{T} [GeV]","#slash{E}_{T} X [GeV]","METvsSumET"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
	histoManager.Add2DVariable(metT[kk]+"METYvsSumET","#Sigma E_{T} [GeV]","#slash{E}_{T} Y [GeV]","METvsSumET"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
	tB[3]=50; tB[4]= 0; tB[5]=100;
	histoManager.Add2DVariable(metT[kk]+"METvsSumET","#Sigma E_{T} [GeV]","#slash{E}_{T} [GeV]","METvsSumET"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);*/

    //METX vs METY
    tB[0]=200; tB[1]=-100; tB[2]=100; tB[3]=200; tB[4]= -100; tB[5]=100;
    /*   histoManager.Add2DVariable(metT[kk]+"METXvsMETY","#slash{E}_{T} X [GeV]","#slash{E}_{T} Y [GeV]","Zee_METXY_",tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    histoManager.Add2DVariable(metT[kk]+"METZFrame","#slash{E}_{T} Para [GeV]","#slash{E}_{T} Perp [GeV]","Zee_METZframe_",tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    histoManager.Add2DVariable(metT[kk]+"METTruthFrame","#slash{E}_{T} tl [GeV]","#slash{E}_{T} tt [GeV]","Zee_METtruth_",tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
    
    tB[0]=100; tB[1]=0; tB[2]=100; tB[3]=60; tB[4]= 0; tB[5]=360;
    histoManager.Add2DVariable(metT[kk]+"METDirNorm","#slash{E}_{T} [GeV]","Phi MET [#circ ]","Zee_phivsnorm_"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);*/

    //***************************************** Control samples
  
    //METs
    histoManager.AddVariable(metT[kk]+"MET_V1",400,0,200,metT[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"MET_V2",400,0,200,metT[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
  
    // histoManager.AddVariable(metT[kk]+"SumETRecoil_V1",500,0,500,"#sum E_{T} ("+metT[kk]+") [GeV]","We_SumET_"+Suffix[kk]);
    //  histoManager.AddVariable(metT[kk]+"SumETRecoil_V2",500,0,500,"#sum E_{T} ("+metT[kk]+") [GeV]","We_SumET_"+Suffix[kk]);

    histoManager.AddVariable(metT[kk]+"SumEt_V1",500,0,500,"#sum E_{T} ("+metT[kk]+") [GeV]","We_SumET_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"SumEt_V2",500,0,500,"#sum E_{T} ("+metT[kk]+") [GeV]","We_SumET_"+Suffix[kk]);
  
    histoManager.AddVariable(metT[kk]+"METControl",400,0,200,metT[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METPhiControl",360,0,360,metT[kk]+" #slash{E}_{T} #Phi ","Zee_METPhi_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METParaControl",400,-200,200,metT[kk]+" #slash{E}_{T,||} [GeV]","We_MET_Para_"+Suffix[kk]);  
    histoManager.AddVariable(metT[kk]+"METPerpControl",400,-200,200,metT[kk]+" #slash{E}_{T,#perp}   [GeV]","We_MET_Perp_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"dPhiRecoilCorControl",100,0,3.14,"#Delta#phi(u_{T,"+metT[kk]+"},lepton) [rad]","We_dPhi_Recoil_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilCorControl",400,0,200,"u_{T} ("+metT[kk]+") [GeV]","We_Recoil_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilParaCorControl",400,-200,200,"u_{||} ("+metT[kk]+") [GeV]","We_Recoil_Para_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilPerpCorControl",400,-200,200,"u_{#perp}   ("+metT[kk]+") [GeV]","We_Recoil_Perp_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilControl",400,0,200,"u_{T} ("+metT[kk]+") [GeV]","We_Recoil_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilParaControl",800,-200,200,"u_{||} ("+metT[kk]+") [GeV]","We_Recoil_Para_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilPerpControl",800,-200,200,"u_{#perp}   ("+metT[kk]+") [GeV]","We_Recoil_Perp_"+Suffix[kk]);
  
    histoManager.AddVariable(metT[kk]+"diffMetControl",400,-1,1,metT[kk]+" diff #slash{E}_{T} [GeV]","Zee_METdiff_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METLepCorControl",400,0,200,metT[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METParaLepCorControl",400,-200,200,metT[kk]+" #slash{E}_{T,||} [GeV]","We_MET_Para_"+Suffix[kk]);  
    histoManager.AddVariable(metT[kk]+"METPerpLepCorControl",400,-200,200,metT[kk]+" #slash{E}_{T,#perp}   [GeV]","We_MET_Perp_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilLepCorControl",400,0,200,"u_{T} ("+metT[kk]+") [GeV]","We_Recoil_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilParaLepCorControl",400,-200,200,"u_{||} ("+metT[kk]+") [GeV]","We_Recoil_Para_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilPerpLepCorControl",400,-200,200,"u_{#perp}   ("+metT[kk]+") [GeV]","We_Recoil_Perp_"+Suffix[kk]);


    //Uncertainties ********************************************

    histoManager.AddVariable(metT[kk]+"METUnc",400,0,200,metT[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METParaUnc",400,-200,200,metT[kk]+" #slash{E}_{T,||} [GeV]","We_MET_Para_"+Suffix[kk]);  
    histoManager.AddVariable(metT[kk]+"METPerpUnc",400,-200,200,metT[kk]+" #slash{E}_{T,#perp}   [GeV]","We_MET_Perp_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"AcopUnc",100,0,3.14,"#zeta ("+metT[kk]+") [rad]","We_Acop_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"dPhiRecoilUnc",100,0,3.14,"#Delta#phi(u_{T,"+metT[kk]+"},lepton) [rad]","We_dPhi_Recoil_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilUnc",400,0,200,"u_{T} ("+metT[kk]+") [GeV]","We_Recoil_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilParaUnc",400,-200,200,"u_{||} ("+metT[kk]+") [GeV]","We_Recoil_Para_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilPerpUnc",400,-200,200,"u_{#perp}   ("+metT[kk]+") [GeV]","We_Recoil_Perp_"+Suffix[kk]);


    //Response spec as W
    // histoManager.Add2DVariable(metT[kk]+"RecoilvsqT","Gen q_{T} [GeV]","<|u_{||}|> [GeV]","",21,-42,42,21,-42,42);
    //  histoManager.Add2DVariable(metT[kk]+"AbsRecoilvsqT","Gen q_{T} [GeV]","<|u_{||}|> [GeV]","",21,-42,42,21,-42,42);

  }
 
  //End Declaration, now preparation
  histoManager.PrepareHistos(name);
  histoManager.Prepare2DHistos(name);
  histoManager.PrepareProfiles(name);
  cout<<" End declaration, now analyse "<<endl;

  vector<vector<vector<float> > > Wghts;
 
  //if( MultiWeight || Remove2V || ShapeWeight)
  Wghts = GetFitWeight();
  
  if( ShapeWeight)
    NVert=1;

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

  //sélection
  for(int i=0;i<nt+1;i++) {
    
    //Déclaration des variables
    float IsoVar[2][3];
    float IDVar[2][6];
    float Mets[6][2];
    float MetProj[6][2];
    float Recoils[6][2];
    float RecoilProj[6][2];
    float Lepton[2][3];
    float PFLepton[2][3];
    float EtSC[2];

    float dPhiRecoilZ[6];
    float dPhiMETZ[6];
    float SumET[6];


    bool ConvRej1;
    bool InHit1;
    bool ExpInHit1;
    bool ConvRej2;
    bool InHit2;
    bool ExpInHit2;

    int NPhot1;
    int NPhot2;

    int Run;
    int Event;
    int AbsEvent;
    char sample[21];
    char fileName[21];

    float GenRecoil[2][2];
    float SpecRecoilProj[6][2];
    float CTRecoil[6][2];
    
    float Corrections[6];
    
    //Corrected MET variables
    float MetCor[6][2];
    float MetCorProj[6][2];
    float RecoilCor[6][2];
    float RecoilCorProj[6][2];

    float GenZPt;
    float ZPt;
    float ZMass;
    float DZvertices;

    float ZMassSC;
    float ZMassSCUnCor;

    //Underlying event
    float dPhiPhotonZ;
    float dPhiJetZ;
    int JetMultiplicity;
    int PhotonMultiplicity;
    float LeadingJet;
    float LeadingPhoton;
    
    int NVertex;
    
    //Decomposed MET SumET
    float caloMETEB[2];
    float caloMETEE[2];
    float caloMETHB[2];
    float caloMETHE[2];
    float caloMETHF[2];
    
    float caloSumETEB;
    float caloSumETEE;
    float caloSumETHB;
    float caloSumETHE;
    float caloSumETHF;

    float pfMETEB[2];
    float pfMETEE[2];
    float pfMETEEx[2]={0,0};
    float pfMETHB[2];
    float pfMETHE[2];
    float pfMETHF[2];
    
    float pfSumETEB;
    float pfSumETEE;
    float pfSumETHB;
    float pfSumETHE;
    float pfSumETHF;
    
    float pfEcalMET[2];
    float pfHcalMET[2];
    float caloEcalMET[2];
    float caloHcalMET[2];

    float cMuFrac;
    float nHadFrac;
    float cHadFrac;
    float nEmFrac;
    float cEmFrac;


    float Weight=0;
    float WeightSav=0;

    bool Sel[2]={true,true};
    int epart[2]={-1,-1};
    TLorentzVector lepton[2];
    int METXY[2][2]={{0,0},{0,0}};
    int ZXY[2][2]={{0,0},{0,0}};

    tChains[i]->SetBranchAddress("NVertex",&NVertex);

    tChains[i]->SetBranchAddress("IsoVar",IsoVar);
    tChains[i]->SetBranchAddress("IDVar",IDVar);
    tChains[i]->SetBranchAddress("Lepton",Lepton);
    tChains[i]->SetBranchAddress("PFLepton",PFLepton);
    tChains[i]->SetBranchAddress("EtSC",EtSC);
   
    tChains[i]->SetBranchAddress("Mets",Mets);
    tChains[i]->SetBranchAddress("MetProj",MetProj);
 
    tChains[i]->SetBranchAddress("SumET",SumET);

    //  tChains[i]->SetBranchAddress("Corrections",Corrections);
    
    tChains[i]->SetBranchAddress("BasicRecoils",Recoils);
    tChains[i]->SetBranchAddress("BasicRecoilProj",RecoilProj);
    tChains[i]->SetBranchAddress("RecoilProj",SpecRecoilProj);

    tChains[i]->SetBranchAddress("BasicdPhiRecoilZ",dPhiRecoilZ);
    tChains[i]->SetBranchAddress("dPhiMETZ",dPhiMETZ);
  
    tChains[i]->SetBranchAddress("InHit1",&InHit1);
    tChains[i]->SetBranchAddress("ExpInHit1",&ExpInHit1);
    tChains[i]->SetBranchAddress("ConvRej1",&ConvRej1);
    tChains[i]->SetBranchAddress("InHit2",&InHit2);
    tChains[i]->SetBranchAddress("ExpInHit2",&ExpInHit2);
    tChains[i]->SetBranchAddress("ConvRej2",&ConvRej2);

   

   
    tChains[i]->SetBranchAddress("Event",&Event);
    tChains[i]->SetBranchAddress("AbsEvent",&AbsEvent);
    tChains[i]->SetBranchAddress("Run",&Run);
    //  tChains[i]->SetBranchAddress("sample",&sample);
    tChains[i]->SetBranchAddress("fileName",&fileName);

    tChains[i]->SetBranchAddress("ZPt",&ZPt); //ZMass
    tChains[i]->SetBranchAddress("ZMass",&ZMass); 
    tChains[i]->SetBranchAddress("DZvertices",&DZvertices);

   
    //MET Corrected variables
    if(i==nt) {
      tChains[i]->SetBranchAddress("MetCor",MetCor);
      tChains[i]->SetBranchAddress("MetCorProj",MetCorProj);
      tChains[i]->SetBranchAddress("RecoilCor",RecoilCor);
      tChains[i]->SetBranchAddress("RecoilCorProj",RecoilCorProj);
      tChains[i]->SetBranchAddress("NPhot1",&NPhot1);
      tChains[i]->SetBranchAddress("NPhot2",&NPhot2);

      tChains[i]->SetBranchAddress("ZMassSC",&ZMassSC); 
      tChains[i]->SetBranchAddress("ZMassSCUnCor",&ZMassSCUnCor); 

    }

    //Decomposed MET/SumET
    tChains[i]->SetBranchAddress("CaloMETEB",caloMETEB);
    tChains[i]->SetBranchAddress("CaloMETEE",caloMETEE);
    tChains[i]->SetBranchAddress("CaloMETHB",caloMETHB);
    tChains[i]->SetBranchAddress("CaloMETHE",caloMETHE);
    tChains[i]->SetBranchAddress("CaloMETHF",caloMETHF);
    
    tChains[i]->SetBranchAddress("CaloSumETEB",&caloSumETEB);
    tChains[i]->SetBranchAddress("CaloSumETEE",&caloSumETEE);
    tChains[i]->SetBranchAddress("CaloSumETHB",&caloSumETHB);
    tChains[i]->SetBranchAddress("CaloSumETHE",&caloSumETHE);
    tChains[i]->SetBranchAddress("CaloSumETHF",&caloSumETHF);
    
    tChains[i]->SetBranchAddress("PfMETEB",pfMETEB);
    tChains[i]->SetBranchAddress("PfMETEE",pfMETEE);
    tChains[i]->SetBranchAddress("PfMETEEx",pfMETEEx);
    tChains[i]->SetBranchAddress("PfMETHB",pfMETHB);
    tChains[i]->SetBranchAddress("PfMETHE",pfMETHE);
    tChains[i]->SetBranchAddress("PfMETHF",pfMETHF);
    
    tChains[i]->SetBranchAddress("PfSumETEB",&pfSumETEB);
    tChains[i]->SetBranchAddress("PfSumETEE",&pfSumETEE);
    tChains[i]->SetBranchAddress("PfSumETHB",&pfSumETHB);
    tChains[i]->SetBranchAddress("PfSumETHE",&pfSumETHE);
    tChains[i]->SetBranchAddress("PfSumETHF",&pfSumETHF);

    tChains[i]->SetBranchAddress("PfEcalMET",pfEcalMET);
    tChains[i]->SetBranchAddress("PfHcalMET",pfHcalMET);
    tChains[i]->SetBranchAddress("CaloEcalMET",caloEcalMET);
    tChains[i]->SetBranchAddress("CaloHcalMET",caloHcalMET);

    tChains[i]->SetBranchAddress("cMuFrac",&cMuFrac);
    tChains[i]->SetBranchAddress("nHadFrac",&nHadFrac);
    tChains[i]->SetBranchAddress("cHadFrac",&cHadFrac);
    tChains[i]->SetBranchAddress("nEmFrac",&nEmFrac);
    tChains[i]->SetBranchAddress("cEmFrac",&cEmFrac);

    //Underlying event
    tChains[i]->SetBranchAddress("LeadingJet",&LeadingJet);
    tChains[i]->SetBranchAddress("LeadingPhoton",&LeadingPhoton);
    tChains[i]->SetBranchAddress("PhotonMultiplicity",&PhotonMultiplicity);
    tChains[i]->SetBranchAddress("JetMultiplicity",&JetMultiplicity);
    tChains[i]->SetBranchAddress("dPhiJetZ",&dPhiJetZ);
    tChains[i]->SetBranchAddress("dPhiPhotonZ",&dPhiPhotonZ);
	
    if(i==nt-1) {
      tChains[i]->SetBranchAddress("GenRecoil",GenRecoil);
      tChains[i]->SetBranchAddress("CTRecoil",CTRecoil);
      tChains[i]->SetBranchAddress("GenZPt",&GenZPt);
    }
    else {
      for(int jj=0;jj<5;jj++) {
	CTRecoil[jj][0]=0; CTRecoil[jj][1]=0;
      }
    }
    
    int EP=0;
    cout<<" beginning tree  " <<name[i]<<"   "<<tChains[i]->GetEntries()<<"   "<<Events.size()<<endl;
    int ent = tChains[i]->GetEntries();

    //Boucle et sélection, remplissage histos
    for(int ie=0;ie<ent;ie++) {
   
      if(ie%1000==0)
	fprintf(stdout, "Done : %d / %d\r", ie, ent );


      if(i!=nt)
	Weight = GetWeight(i, ie);
      else
	Weight = 1;
      
      /* if(NVert==1 && i!=nt)
	Weight *= Frac1V;
      else if(NVert==2 && i!=nt)
      Weight *= (1-Frac1V);*/

      //FIXME 
      if(i<nt-2) {
	if(NVert==1 )
	  Weight *= Frac1V;
	if(NVert==2 )
	  Weight *= (1-Frac1V);
      }

      if(name[i].substr(0,4)=="Z #r")
      	Weight *= SearchWeight(ZPt);

      WeightSav = Weight;

      tChains[i]->GetEntry(ie);      
   
      if(name[i].substr(0,4)=="data")
	{
	  bool doubleCount=false;
	 
	  for(size_t ik=0;ik<Events.size();ik++)
	    {
	      if( ((Events[ik]).first.first) == Run && ((Events[ik]).first.second) == Event)
	      	{doubleCount=true; break;}
	    }

	  if(doubleCount || (EventFilter && Run>EventNum ) )
	    { continue; }
	 
	}

   
      //Filling Histos
      
      //Passing selection
      for(int id=0;id<2;id++) {

	if(fabs(Lepton[id][1])>1.479)
	  {EP=1; }
	else
	  {EP=0; }
	
	if(name[i]!="data" || EP==0 || noDeltaCor) {
	  IDVar[id][4] = 0;
	  IDVar[id][5] =0;
	}
	
	Sel[id]=true;

	//first pt cut
	if(EtSC[id] < PTcut)  Sel[id]=false;
	
	if(convRej && !ExpInHit1 && !ExpInHit2)   Sel[id]=false;
	
	//Met Cut
	
	//And Selections
	//  cout<<IDVar[4]<<endl;
	if(IDVar[id][0] > _IDcuts[EP][0] && Cbool(skipCut, "sigieie"==observable.substr(0,7)) )  Sel[id]=false;
	if( fabs(IDVar[id][1]) > _IDcuts[EP][1] && Cbool(skipCut, "Deta"==observable.substr(0,4)) )  Sel[id]=false;
	if( fabs(IDVar[id][2]) > _IDcuts[EP][2] && Cbool(skipCut, "Dphi"==observable.substr(0,4)) )  Sel[id]=false;
	if(IDVar[id][3] > _IDcuts[EP][3] && Cbool(skipCut, "HoE"==observable.substr(0,3)) )  Sel[id]=false;
	
	if(IsoVar[id][0]/Lepton[id][0] > _Isocuts[EP][0] && Cbool(skipCut, "TrackIso"==observable.substr(0,8)) )  Sel[id]=false;
	if(IsoVar[id][1]/Lepton[id][0] > _Isocuts[EP][1] && Cbool(skipCut, "EcalIso"==observable.substr(0,7)) )  Sel[id]=false;
	if(IsoVar[id][2]/Lepton[id][0] > _Isocuts[EP][2] && Cbool(skipCut, "HcalIso"==observable.substr(0,7)) )  Sel[id]=false;
	
	if(MultiWeight) {
	  if(i<=1 && i!=nt) {
	    if(NVert==0)
	      Weight = (Wghts[0][0][1]+Wghts[0][1][1])*GetWeight(i, ie);
	    else if(NVert==1)
	      Weight = Wghts[0][0][1]*GetWeight(i, ie);
	    else if(NVert==2)
	      Weight = Wghts[0][1][1]*GetWeight(i, ie);
	  }
	  else if(i>1 && i!=nt) {
	    if(NVert==0)
	      Weight = (Wghts[0][0][0]+Wghts[0][1][0])*GetWeight(i, ie);
	    else if(NVert==1)
	      Weight = Wghts[0][0][0]*GetWeight(i, ie);
	    else if(NVert==2)
	      Weight = Wghts[0][1][0]*GetWeight(i, ie);
	  }
	  else if(i==nt)
	    Weight = 1;
	}
	
	if(EcalP!=EP && EcalP!=2 && EcalP!=3) Sel[id]=false;

	if(EcalP==3) 
	  epart[id]=EP;
	
	if(NVert==2 && NVertex<2 && ( i>=nt-2/*name[i].substr(0,4)=="data"|| name[i].substr(0,4)=="Z #r"*/) )
	  Sel[id]=false;
	else if(NVert==1 && NVertex>1 && ( i>=nt-2/*name[i].substr(0,4)=="data" || name[i].substr(0,4)=="Z #r"*/) )
	  Sel[id]=false;
	//	cout<<NVertex<<"    "<<name[i].substr(0,4)<<endl;
	if(!Sel[id]) break;
	
	//ID
	histoManager.fill("sigieie"+lep[id],i, IDVar[id][0],Weight);
	histoManager.fill("Deta"+lep[id],i, fabs(IDVar[id][1]-IDVar[id][4]),Weight);
	histoManager.fill("Dphi"+lep[id],i, IDVar[id][2]-IDVar[id][5],Weight);
	histoManager.fill("HoE"+lep[id],i, IDVar[id][3],Weight);
	
	//Iso
	histoManager.fill("TrackIso"+lep[id],i,IsoVar[id][0]/Lepton[id][0],Weight);
	histoManager.fill("EcalIso"+lep[id],i,IsoVar[id][1]/Lepton[id][0],Weight);
	histoManager.fill("HcalIso"+lep[id],i,IsoVar[id][2]/Lepton[id][0],Weight);
	

      }

      if(!Sel[0] || !Sel[1])
	continue;
      
      if(EcalP==3 && epart[0]==epart[1] ) continue;
      //if(ZPt>30) continue;
      //  if( ( (epart[1]==0) && (epart[0]==1 && (Lepton[0][1]<-2.1 || Lepton[0][1]>-1.479) ) || 
      //	  ( (epart[0]==0) && (epart[1]==1 && (Lepton[1][1]<-2.1 || Lepton[0][1]>-1.479) ) ) continue;

      // if(ZPt < 20 ) continue;
      /*  if( name[i].substr(0,4)=="data")
	  cout<<Run<<"    "<<Event<<"   "
	  <<Lepton[0][0]<<"   "<<Lepton[0][1]<<"   "<<Lepton[0][2]<<"   "
	  <<Lepton[1][0]<<"   "<<Lepton[1][1]<<"   "<<Lepton[1][2]<<endl;*/

      histoManager.fill("jetMult",i,JetMultiplicity,Weight);
      //if(JetMultiplicity == 0 ) continue;

      TVector2 mettmp(0,0);
      mettmp.SetMagPhi(Mets[0][0],Mets[0][1]);
      TVector2 mettmp2(0,0);
      mettmp2.SetMagPhi(Mets[1][0],Mets[1][1]);

      if( ( fabs(mettmp.Px() )  < METcut && fabs(mettmp.Py())  < METcut  && !invCut && ZMass > MTCut && ZMass < 120) || 
	  (Mets[0][0] < METcut && invCut && ZMass > MTCut && ZMass < 140 ) ) { 
	float cor;
	for(int id=0;id<2;id++) {
	  	  histoManager.fill("Pt"+lep[id],i,Lepton[id][0],Weight);
	  if(PFLepton[id][0]!=-1)
	       histoManager.fill("PfPt"+lep[id],i,PFLepton[id][0],Weight);
	     histoManager.fill("EtSC"+lep[id],i,EtSC[id],Weight);
	  if(fabs(Lepton[id][1])<1.5)
	    cor=1.008;
	  else
	    cor=1.027;
	  if(i==nt) {
	    histoManager.fill("EtSCCor"+lep[id],i,EtSC[id]*cor,Weight);
	    histoManager.fill("PfPtCor"+lep[id],i,PFLepton[id][0]*cor,Weight);
	  }
	  else {
	    histoManager.fill("EtSCCor"+lep[id],i,EtSC[id],Weight);
	    histoManager.fill("PfPtCor"+lep[id],i,PFLepton[id][0]*cor,Weight);
	  }
	  histoManager.fill("Eta"+lep[id],i,Lepton[id][1],Weight);
	  histoManager.fill("Phi"+lep[id],i,Lepton[id][2],Weight);
	  histoManager.fill("EcalIso"+lep[id],i,IsoVar[id][1]/Lepton[id][0],Weight);
	  histoManager.fill("TrackIso"+lep[id],i,IsoVar[id][0]/Lepton[id][0],Weight);
	  histoManager.fill("HcalIso"+lep[id],i,IsoVar[id][2]/Lepton[id][0],Weight);
	  

	  lepton[id].SetPtEtaPhiM(Lepton[id][0],Lepton[id][1],Lepton[id][2],0.000000511);
	}
	
	if(i==nt) {
	  histoManager.fill("NPhot1",i,NPhot1,Weight);
	  histoManager.fill("NPhot2",i,NPhot2,Weight);
	}
	else {
	  histoManager.fill("NPhot1",i,-100,Weight);
	  histoManager.fill("NPhot2",i,-100,Weight);
	}

	//Rebuild the Z
	TLorentzVector Zvect = lepton[0] + lepton[1];
	histoManager.fill("ZPhi",i,Zvect.Phi(),Weight);
	histoManager.fill("ZEta",i,Zvect.Eta(),Weight);
	if(Zvect.Px() < 0 ) {
	  if(Zvect.Py() < 0 )
	    ZXY[0][0]++;
	  else if(Zvect.Py() > 0)
	    ZXY[0][1]++;
	}
	else if(Zvect.Px() > 0) {
	  if(Zvect.Py() < 0 )
	    ZXY[1][0]++;
	  else if(Zvect.Py() > 0 ) 
	    ZXY[1][1]++;
	}

	TVector2 TruthFrame(0,0);
	TruthFrame = lepton[0].Vect().XYvector() - lepton[1].Vect().XYvector();
	TVector2 tl = TruthFrame.Unit();
	
  
	TVector2 tt = tl.Rotate(3.1415/2);
	tt *= cos(tt.DeltaPhi(lepton[1].Vect().XYvector()))/fabs( cos(tt.DeltaPhi(lepton[1].Vect().XYvector())) ); 
	if(  fabs( lepton[0].Vect().XYvector().DeltaPhi(lepton[1].Vect().XYvector()) ) < 3.1415/2. ) {
	  tt = Zvect.Vect().XYvector().Unit();
	  tl = tt.Rotate( -acos(-1)/2);
	  tl *= - cos(tl.DeltaPhi(lepton[1].Vect().XYvector()))/fabs( cos(tl.DeltaPhi(lepton[1].Vect().XYvector())) ); 
	}



	//ZPt
	histoManager.fill("ZPt",i,ZPt,Weight);
	histoManager.fill("ZMass",i,ZMass,Weight);
	

	//Separation of mass ======================
	
	if( Lepton[0][1]>=0 && Lepton[1][1]>=0) //BB+
	  histoManager.fill("ZMassEBEBp",i,ZMass,Weight);

	if( Lepton[0][1]<0 && Lepton[1][1]<0) //BB-
	  histoManager.fill("ZMassEBEBm",i,ZMass,Weight);

	if( ( (epart[1]==0) && (epart[0]==1 && (Lepton[0][1]>-2.1 && Lepton[0][1]<-1.479) ) ) || 
	    ( (epart[0]==0) && (epart[1]==1 && (Lepton[1][1]>-2.1 && Lepton[1][1]<-1.479) ) ) ) //EE - inner
	  histoManager.fill("ZMassEBEEmInner",i,ZMass,Weight);
	
	if( ( (epart[1]==0) && (epart[0]==1 && (Lepton[0][1]<-2.1) ) ) || 
	    ( (epart[0]==0) && (epart[1]==1 && (Lepton[1][1]<-2.1) ) ) ) //EE - outer
	  histoManager.fill("ZMassEBEEmOuter",i,ZMass,Weight);
	      
	if( ( (epart[1]==0) && (epart[0]==1 && (Lepton[0][1]<2.1 && Lepton[0][1]>1.479) ) ) || 
	    ( (epart[0]==0) && (epart[1]==1 && (Lepton[1][1]<2.1 && Lepton[1][1]>1.479) ) ) ) //EE + inner
	  histoManager.fill("ZMassEBEEpInner",i,ZMass,Weight);
	
	if( ( (epart[1]==0) && (epart[0]==1 && (Lepton[0][1]>2.1) ) ) || 
	    ( (epart[0]==0) && (epart[1]==1 && (Lepton[1][1]>2.1) ) ) ) //EE + outer
	  histoManager.fill("ZMassEBEEpOuter",i,ZMass,Weight);
	
	if( ( (epart[1]==0) && (epart[0]==1 && (Lepton[0][1]>1.479) ) ) || 
	    ( (epart[0]==0) && (epart[1]==1 && (Lepton[1][1]>1.479) ) ) ) // EE +
	  histoManager.fill("ZMassEBEEp",i,ZMass,Weight);
	
	if( ( (epart[1]==0) && (epart[0]==1 && (Lepton[0][1]<-1.479) ) ) || 
	    ( (epart[0]==0) && (epart[1]==1 && (Lepton[1][1]<-1.479) ) ) ) // EE -
	  histoManager.fill("ZMassEBEEm",i,ZMass,Weight);

	if( ( (epart[1]==0) && (epart[0]==1 && (Lepton[0][1]<-1.479 && (Lepton[0][2]<3.1415/2 || Lepton[0][2]>3.1415*3/2) ) ) ) || 
	    ( (epart[0]==0) && (epart[1]==1 && (Lepton[1][1]<-1.479 && (Lepton[1][2]<3.1415/2 || Lepton[1][2]>3.1415*3/2)) ) ) ) // EE - Right
	  histoManager.fill("ZMassEBEEmRight",i,ZMass,Weight);

	if( ( (epart[1]==0) && (epart[0]==1 && (Lepton[0][1]<-1.479 && (Lepton[0][2]>3.1415/2 && Lepton[0][2]<3.1415*3/2) ) ) ) || 
	    ( (epart[0]==0) && (epart[1]==1 && (Lepton[1][1]<-1.479 && (Lepton[1][2]>3.1415/2 && Lepton[1][2]<3.1415*3/2)) ) ) ) // EE - Left
	  histoManager.fill("ZMassEBEEmLeft",i,ZMass,Weight);

	if( ( (epart[1]==0) && (epart[0]==1 && (Lepton[0][1]>1.479 && (Lepton[0][2]<3.1415/2 || Lepton[0][2]>3.1415*3/2) ) ) ) || 
	    ( (epart[0]==0) && (epart[1]==1 && (Lepton[1][1]>1.479 && (Lepton[1][2]<3.1415/2 || Lepton[1][2]>3.1415*3/2)) ) ) ) // EE + Right
	  histoManager.fill("ZMassEBEEpRight",i,ZMass,Weight);

	if( ( (epart[1]==0) && (epart[0]==1 && (Lepton[0][1]>1.479 && (Lepton[0][2]>3.1415/2 && Lepton[0][2]<3.1415*3/2) ) ) ) || 
	    ( (epart[0]==0) && (epart[1]==1 && (Lepton[1][1]>1.479 && (Lepton[1][2]>3.1415/2 && Lepton[1][2]<3.1415*3/2)) ) ) ) // EE + Left
	  histoManager.fill("ZMassEBEEpLeft",i,ZMass,Weight);

	
	//=========================================

	//Corrected mass ===========================
	float Et[2]={Lepton[0][0],Lepton[1][0]};
	for(int jj=0;jj<2;jj++) {
	  if(fabs(Lepton[jj][1])>1.5)
	    Et[jj]*=1.027;
	  else
	    Et[jj]*=1.008;
	}
	if(i==nt-1)
	  histoManager.fill("ZMassCor",i,ComputeMinv(EtSC[0],EtSC[1],Lepton[0][1],Lepton[1][1],Lepton[0][2],Lepton[1][2]) ,Weight);
	else if(i==nt)
	  histoManager.fill("ZMassCor",i,ZMassSC,Weight);

	if(i==nt) {
	  histoManager.fill("ZMassSC",i,ZMassSC,Weight);
	  histoManager.fill("ZMassSCUnCor",i,ZMassSCUnCor,Weight);
	}
	else {
	  histoManager.fill("ZMassSC",i,ZMass,Weight);
	  histoManager.fill("ZMassSCUnCor",i,ZMass,Weight);
	}

	//=========================================
	
	//Vertex
	histoManager.fill("Vertex",i,NVertex,Weight);

	//Another frame
	TVector2 noiseAxis(1,0);
	noiseAxis = noiseAxis.Unit();
	noiseAxis = noiseAxis.Rotate( -34.86*acos(-1)/180. );
	TVector2 noiseAxisT = noiseAxis.Rotate(acos(-1)/2. );
	//Difference between mets
	histoManager.fill("diffPfTcMETX",i,(mettmp.Px()-mettmp2.Px())/Mets[0][0],Weight);
	histoManager.fill("diffPfTcMETY",i,(mettmp.Py()-mettmp2.Py())/Mets[0][0],Weight);
	histoManager.fill("diffPfTcMETSFL",i,(mettmp*noiseAxis - mettmp2*noiseAxis)/Mets[0][0],Weight);
	histoManager.fill("diffPfTcMETSFT",i,(mettmp*noiseAxisT - mettmp2*noiseAxisT)/Mets[0][0],Weight);
	

	TVector2 pft1mettmp(0,0);
	TVector2 lept1(0,0), lept2(0,0);

	pft1mettmp.SetMagPhi(Mets[5][0],Mets[5][1]);
	lept1.SetMagPhi(1,Lepton[0][2]);
	lept2.SetMagPhi(1,Lepton[1][2]);
	
	TVector2 rp = (lept1 + lept2).Unit();

	histoManager.fill("ParaBis",i, rp*pft1mettmp, Weight);
	histoManager.fill("PerpBis",i, (rp.Rotate(TMath::Pi()))*pft1mettmp, Weight);

	//Pffracs
	histoManager.fill("pfEmNFrac",i,nEmFrac,Weight);
	histoManager.fill("pfEmCFrac",i,cEmFrac,Weight);
	histoManager.fill("pfHadNFrac",i,nHadFrac,Weight);
	histoManager.fill("pfHadCFrac",i,cHadFrac,Weight);
	histoManager.fill("pfMuCFrac",i,cMuFrac,Weight);
	

	for(int kk=0;kk<6;kk++) {

	  Corrections[kk] = SearchESWeight(ZPt, kk);

	  if(MultiWeight) {
	    if(i<=1 && i!=nt) {
	      if(NVert==0)
		Weight = (Wghts[kk][0][1]+Wghts[kk][1][1])*GetWeight(i, ie);
	      else if(NVert==1)
		Weight = Wghts[kk][0][1]*GetWeight(i, ie);
	      else if(NVert==2)
		Weight = Wghts[kk][1][1]*GetWeight(i, ie);
	    }
	    else if(i>1 && i!=nt) {
	 
	      if(NVert==0)
		Weight =( Wghts[kk][0][0]+Wghts[kk][1][0])*GetWeight(i, ie);
	      else if(NVert==1)
		Weight = (Wghts[kk][0][0])*GetWeight(i, ie);
	      else if(NVert==2)
		Weight = (Wghts[kk][1][0])*GetWeight(i, ie);
	    }
	    else
	      Weight = 1;
	  }
	  //	cout<<ShapeWeight<<"   "<<i<<"  "<<nt<<"   "<<NVert<<endl;
	
	  if(ShapeWeight && i!=nt ) {
	    Weight = WeightSav;
	    /*
	      Weight = Weight*GetMCWeightFromDataShape(kk, SumET[kk],
	      metT[kk]+"SumEt_V1",
	      metT[kk]+"SumEt_V1", 10. );*/
	    Weight = Weight*GetMCWeightFromDataShape(5, ZPt,
						     "ZPt_V1",
						     "ZPt_V1", 5., i );
	
	    /*  Weight = Weight*GetMCWeightFromDataShape(0, caloSumETEE+caloSumETHF,
		"pfSumETHF_V1",
		"pfSumETHF_V1", 10.);*/
	
	  }
	  else if(i==nt)
	    Weight=1;
	  
	  //MET
	  //	if(kk==4)
	  //	  cout<<Mets[kk][0]<<endl;
	  histoManager.fill(metT[kk]+"MET",i,Mets[kk][0],Weight);
	  histoManager.fill(metT[kk]+"METPhi",i,Mets[kk][1]*180/3.1415,Weight);
	
	  //MET Proj
	  histoManager.fill(metT[kk]+"METPara",i,-MetProj[kk][0],Weight);
	  histoManager.fill(metT[kk]+"METPerp",i,-MetProj[kk][1],Weight);
	
	  // histoManager.fill2D(metT[kk]+"METDirNorm",i,Mets[kk][0],Mets[kk][1]*180/3.1415,Weight);

	  //	histoManager.fill2D(metT[kk]+"METZFrame",i,MetProj[kk][0],MetProj[kk][1],Weight);

	  histoManager.fill(metT[kk]+"dPhiMETZ",i,fabs(dPhiMETZ[kk]),Weight);


	  histoManager.fill(metT[kk]+"SumEt",i,SumET[kk]-Lepton[0][0]-Lepton[1][0],Weight);
	
	  mettmp.SetMagPhi(Mets[kk][0],Mets[kk][1]);

	  histoManager.fill(metT[kk]+"METX",i,mettmp.Px(),Weight);
	  histoManager.fill(metT[kk]+"METY",i,mettmp.Py(),Weight);

	  histoManager.fill(metT[kk]+"METPhiL1",i,mettmp.DeltaPhi(lepton[0].Vect().XYvector() ) * 180/3.14,Weight);
	  histoManager.fill(metT[kk]+"METPhiL2",i,mettmp.DeltaPhi(lepton[1].Vect().XYvector() ) * 180/3.14,Weight);

	  //	histoManager.fill2D(metT[kk]+"dPhiL12D",i,Mets[kk][0],mettmp.DeltaPhi(lepton[0].Vect().XYvector() ) * 180/3.14,Weight);
	  //	histoManager.fill2D(metT[kk]+"dPhiL22D",i,Mets[kk][0],mettmp.DeltaPhi(lepton[1].Vect().XYvector() ) * 180/3.14,Weight);
	  //	histoManager.fill2D(metT[kk]+"dPhiZ2D",i,Mets[kk][0],fabs(dPhiMETZ[kk]),Weight);

	  histoManager.fill(metT[kk]+"METtt",i, tt*mettmp,Weight);
	  histoManager.fill(metT[kk]+"METtl",i, tl*mettmp,Weight);

	  // histoManager.fill2D(metT[kk]+"METTruthFrame",i,tl*mettmp,tt*mettmp,Weight);

	  TVector2 mettmpCor(0,0);
	  mettmpCor.SetMagPhi(MetCor[kk][0],MetCor[kk][1]);

	  //Another frame!
	  histoManager.fill(metT[kk]+"METSFramePara",i,noiseAxis*mettmp,Weight);
	  histoManager.fill(metT[kk]+"METSFramePerp",i,noiseAxisT*mettmp,Weight);

	  if(mettmp.Py() > mettmp.Px()*3) //upper
	    histoManager.fill(metT[kk]+"METUpper",i,Mets[kk][0],Weight);
	  else //lower
	    histoManager.fill(metT[kk]+"METLower",i,Mets[kk][0],Weight);
	  /*	histoManager.fill2D(metT[kk]+"METXvsSumET",i,pfSumETHF,mettmp.Px(),NumberEntries[i]);
		histoManager.fill2D(metT[kk]+"METYvsSumET",i,pfSumETHF,mettmp.Py(),NumberEntries[i]);
		histoManager.fill2D(metT[kk]+"METvsSumET",i,pfSumETHF,Mets[kk][0],NumberEntries[i]);
	  */
	  //	histoManager.fill2D(metT[kk]+"METXvsMETY",i,noiseAxis*mettmp,noiseAxisT*mettmp,Weight);
	  if(kk==1) {
	    if(mettmp.Px() < 0 && mettmp.Px() > -15) {
	      if(mettmp.Py() < 0 && mettmp.Py() > -15)
		METXY[0][0]++;
	      else if(mettmp.Py() > 0 && mettmp.Py() < 15)
		METXY[0][1]++;
	    }
	    else if(mettmp.Px() > 0 && mettmp.Px() < 15) {
	      if(mettmp.Py() < 0 && mettmp.Py() > -15)
		METXY[1][0]++;
	      else if(mettmp.Py() > 0 && mettmp.Py() < 15) 
		METXY[1][1]++;
	    }
	  }	


	  //Corrected MET variables
	  /* if(i==nt) {
	    histoManager.fill(metT[kk]+"METLepCor",i,MetCor[kk][0],Weight);
	    histoManager.fill(metT[kk]+"METParaLepCor",i,MetCorProj[kk][0],Weight);
	    histoManager.fill(metT[kk]+"METPerpLepCor",i,MetCorProj[kk][1],Weight);
	    histoManager.fill(metT[kk]+"RecoilLepCor",i,RecoilCor[kk][0],Weight);
	    histoManager.fill(metT[kk]+"RecoilParaLepCor",i,RecoilCorProj[kk][0],Weight);
	    histoManager.fill(metT[kk]+"RecoilPerpLepCor",i,RecoilCorProj[kk][1],Weight);
	    histoManager.fill(metT[kk]+"diffMet",i,(MetCor[kk][0]-Mets[kk][0])/Mets[kk][0],Weight);

	    histoManager.fill(metT[kk]+"METXCor",i,mettmpCor.Px(),Weight);
	  }
	  else {
	    histoManager.fill(metT[kk]+"METLepCor",i,Mets[kk][0],Weight);
	    histoManager.fill(metT[kk]+"METParaLepCor",i,MetProj[kk][0],Weight);
	    histoManager.fill(metT[kk]+"METPerpLepCor",i,MetProj[kk][1],Weight);
	    histoManager.fill(metT[kk]+"RecoilLepCor",i,Recoils[kk][0],Weight);
	    histoManager.fill(metT[kk]+"RecoilParaLepCor",i,RecoilProj[kk][0],Weight);
	    histoManager.fill(metT[kk]+"RecoilPerpLepCor",i,RecoilProj[kk][1],Weight);
	    histoManager.fill(metT[kk]+"diffMet",i,-1000,Weight);

	    histoManager.fill(metT[kk]+"METXCor",i,mettmp.Px(),Weight);
	  }*/

	  //****** Cutted variables *********


	
	  //Underlying event
	  if(kk==0) {

	    if(LeadingJet> 15)
	      histoManager.fill("jetPt",i,LeadingJet,Weight);
	    if(LeadingPhoton> 15)
	      histoManager.fill("photonPt",i,LeadingPhoton,Weight);
	    histoManager.fill("photonMult",i,PhotonMultiplicity,Weight);
	 
	    histoManager.fill("dZv",i,DZvertices,Weight); 
	  
	    if(LeadingJet> 15)
	      histoManager.fill("dPhiJetZ",i,fabs(dPhiJetZ),Weight);
	    if(LeadingPhoton> 5)
	      histoManager.fill("dPhiPhotonZ",i,fabs(dPhiPhotonZ),Weight);
	  
	    histoManager.fill("pfMETPhiEB",i,pfMETEB[1]*180/3.1415,Weight);
	    histoManager.fill("pfMETPhiEE",i,pfMETEE[1]*180/3.1415,Weight);
	    histoManager.fill("pfMETPhiEEp",i,pfMETEEx[0]*180/3.1415,Weight);
	    histoManager.fill("pfMETPhiEEm",i,pfMETEEx[1]*180/3.1415,Weight);
	    histoManager.fill("pfMETPhiHB",i,pfMETHB[1]*180/3.1415,Weight);
	    histoManager.fill("pfMETPhiHE",i,pfMETHE[1]*180/3.1415,Weight);
	    histoManager.fill("pfMETPhiHF",i,pfMETHF[1]*180/3.1415,Weight);


	    /*
	    histoManager.fill("caloMETEB",i,caloMETEB[0],Weight);
	    histoManager.fill("caloMETEE",i,caloMETEE[0],Weight);
	    histoManager.fill("caloMETHB",i,caloMETHB[0],Weight);
	    histoManager.fill("caloMETHE",i,caloMETHE[0],Weight);
	    histoManager.fill("caloMETHF",i,caloMETHF[0],Weight);

	    histoManager.fill("caloMETPhiEB",i,caloMETEB[1]*180/3.1415,Weight);
	    histoManager.fill("caloMETPhiEE",i,caloMETEE[1]*180/3.1415,Weight);
	    histoManager.fill("caloMETPhiHB",i,caloMETHB[1]*180/3.1415,Weight);
	    histoManager.fill("caloMETPhiHE",i,caloMETHE[1]*180/3.1415,Weight);
	    histoManager.fill("caloMETPhiHF",i,caloMETHF[1]*180/3.1415,Weight);

	    histoManager.fill("caloSumETEB",i,caloSumETEB,Weight);
	    histoManager.fill("caloSumETEE",i,caloSumETEE,Weight);
	    histoManager.fill("caloSumETHB",i,caloSumETHB,Weight);
	    histoManager.fill("caloSumETHE",i,caloSumETHE,Weight);
	    histoManager.fill("caloSumETHF",i,caloSumETHF,Weight);

	    histoManager.fill("pfMETEB",i,pfMETEB[0],Weight);
	    histoManager.fill("pfMETEE",i,pfMETEE[0],Weight);
	    histoManager.fill("pfMETHB",i,pfMETHB[0],Weight);
	    histoManager.fill("pfMETHE",i,pfMETHE[0],Weight);
	    histoManager.fill("pfMETHF",i,pfMETHF[0],Weight);

	   
	    histoManager.fill("pfSumETEB",i,pfSumETEB,Weight);
	    histoManager.fill("pfSumETEE",i,pfSumETEE,Weight);
	    histoManager.fill("pfSumETHB",i,pfSumETHB,Weight);
	    histoManager.fill("pfSumETHE",i,pfSumETHE,Weight);
	    histoManager.fill("pfSumETHF",i,pfSumETHF,Weight);
	    *//*
	    histoManager.fill("pfMETEcal",i,pfEcalMET[0],Weight);
	    histoManager.fill("pfMETHcal",i,pfHcalMET[0],Weight);
	    TVector2 pfMET1(0,0); pfMET1.SetMagPhi(pfEcalMET[0],pfEcalMET[1]);
	    TVector2 pfMET2(0,0); pfMET2.SetMagPhi(pfHcalMET[0],pfHcalMET[1]);
	    TVector2 pfMETHF2(0,0); pfMETHF2.SetMagPhi(pfMETHF[0],pfMETHF[1]);
	    TVector2 pfMET = pfMET1 + pfMET2 - pfMETHF2;
	    histoManager.fill("pfMETAll",i,pfMET.Mod(),Weight);
	    histoManager.fill("pfMETPhiEcal",i,pfEcalMET[1]*180/3.1415,Weight);
	    histoManager.fill("pfMETPhiHcal",i,pfHcalMET[1]*180/3.1415,Weight);
	  
	    histoManager.fill("caloMETEcal",i,caloEcalMET[0],Weight);
	    histoManager.fill("caloMETHcal",i,caloHcalMET[0],Weight);
	    TVector2 caloMET1(0,0); caloMET1.SetMagPhi(caloEcalMET[0],caloEcalMET[1]);
	    TVector2 caloMET2(0,0); caloMET2.SetMagPhi(caloHcalMET[0],caloHcalMET[1]);
	    TVector2 caloMETHF2(0,0); caloMETHF2.SetMagPhi(caloMETHF[0],caloMETHF[1]);
	    TVector2 caloMET = caloMET1 + caloMET2 - caloMETHF2;
	    histoManager.fill("caloMETAll",i,caloMET.Mod(),Weight);
	    histoManager.fill("caloMETPhiEcal",i,caloEcalMET[1]*180/3.1415,Weight);
	    histoManager.fill("caloMETPhiHcal",i,caloHcalMET[1]*180/3.1415,Weight);*/
	    /*
	      histoManager.fill2D("pfMETETPhiEB",i,pfMETEB[1]*180/3.1415,pfMETEB[0],Weight);
	      histoManager.fill2D("pfMETETPhiEE",i,pfMETEE[1]*180/3.1415,pfMETEE[0],Weight);
	      histoManager.fill2D("pfMETETPhiHB",i,pfMETHB[1]*180/3.1415,pfMETHB[0],Weight);
	      histoManager.fill2D("pfMETETPhiHE",i,pfMETHE[1]*180/3.1415,pfMETHE[0],Weight);
	      histoManager.fill2D("pfMETETPhiHF",i,pfMETHF[1]*180/3.1415,pfMETHF[0],Weight);
	      histoManager.fill2D("caloMETETPhiEB",i,caloMETEB[1]*180/3.1415,caloMETEB[0],Weight);
	      histoManager.fill2D("caloMETETPhiEE",i,caloMETEE[1]*180/3.1415,caloMETEE[0],Weight);
	      histoManager.fill2D("caloMETETPhiHB",i,caloMETHB[1]*180/3.1415,caloMETHB[0],Weight);
	      histoManager.fill2D("caloMETETPhiHE",i,caloMETHE[1]*180/3.1415,caloMETHE[0],Weight);
	      histoManager.fill2D("caloMETETPhiHF",i,caloMETHF[1]*180/3.1415,caloMETHF[0],Weight);
	    */


	  }
	
	  //Uncertainties
	  if( name[i].substr(0,4)=="Z #r") {
	
	    float dbw = SearchWeight(ZPt);
	
	    //if( (Mets[kk][0] > METcut && !invCut && MT[kk] > MTCut) || 
	    //	(Mets[kk][0] < METcut && invCut && MT[kk] > MTCut ) ) { 
		
	    if(kk==0) {
	      histoManager.fill("ZPtUnc",i,ZPt,Weight*(1-dbw) );
	      histoManager.fill("ZMassUnc",i,ZMass,Weight*(1-dbw) );
	    }

		histoManager.fill(metT[kk]+"METUnc",i,Mets[kk][0],Weight*(1-dbw) );
		histoManager.fill(metT[kk]+"METParaUnc",i,MetProj[kk][0],Weight*(1-dbw) );
		histoManager.fill(metT[kk]+"METPerpUnc",i,MetProj[kk][1],Weight*(1-dbw) );
		
		//histoManager.fill(metT[kk]+"dPhiRecoilUnc",i,(dPhiRecoilLep[kk]*3.14/180),Weight*(1-dbw) );
		histoManager.fill(metT[kk]+"RecoilUnc",i,Recoils[kk][0],Weight*(1-dbw) );
		histoManager.fill(metT[kk]+"RecoilParaUnc",i,RecoilProj[kk][0],Weight*(1-dbw) );
		histoManager.fill(metT[kk]+"RecoilPerpUnc",i,RecoilProj[kk][1],Weight*(1-dbw) );
		
		//	}
	  }

	
	  //Recoil
	  histoManager.fill(metT[kk]+"Recoil",i,Recoils[kk][0],Weight);
	  histoManager.fill(metT[kk]+"RecoilPhi",i,Recoils[kk][1]*180/3.1415,Weight);

	  float response = SearchESWeight(ZPt, kk);
	  histoManager.fill(metT[kk]+"RecoilCor",i,Recoils[kk][0]*response,Weight);
	  //Recoil Proj
	  histoManager.fill(metT[kk]+"RecoilPara",i,RecoilProj[kk][0],Weight);
	  histoManager.fill(metT[kk]+"RecoilPerp",i,RecoilProj[kk][1],Weight);
  	  histoManager.fill(metT[kk]+"RecoilParaCor",i,RecoilProj[kk][0]*response,Weight);
	  histoManager.fill(metT[kk]+"RecoilPerpCor",i,RecoilProj[kk][1]*response,Weight);

	  //Angles
	  histoManager.fill(metT[kk]+"dPhiRecoil",i,fabs(dPhiRecoilZ[kk]*3.14/180),Weight);
	  histoManager.fill(metT[kk]+"dPhiRecoilCor",i,fabs(dPhiRecoilZ[kk]*3.14/180),Weight);
		  /*	  histoManager.fill2D(metT[kk]+"dPhiRecoil2D",i,Recoils[kk][0],fabs(dPhiRecoilZ[kk]*3.14/180),Weight);
		  histoManager.fill2D(metT[kk]+"dPhiRecoil2DCor",i,Recoils[kk][0]*Corrections[kk],fabs(dPhiRecoilZ[kk]*3.14/180),Weight);*/
	
	  //Response and other plots
	  //	  histoManager.fill2D(metT[kk]+"UtVsQt",i,ZPt,Recoils[kk][0],Weight);
	  //	  histoManager.fill2D(metT[kk]+"UtVsQtCor",i,ZPt,Recoils[kk][0]*Corrections[kk],Weight);
	  //Veto pour ttbar
	  if(Run == 142191 && AbsEvent == 25487594 ) continue;

	  if(ZPt>1) {
	    histoManager.fill2D(metT[kk]+"UparaVsQt",i,ZPt,RecoilProj[kk][0],Weight);
	    histoManager.fill2D(metT[kk]+"UparaVsQtCor",i,ZPt,RecoilProj[kk][0]*response,Weight);
	    histoManager.fill2D(metT[kk]+"UperpVsQt",i,ZPt,RecoilProj[kk][1],Weight);
	    histoManager.fill2D(metT[kk]+"UperpVsQtCor",i,ZPt,RecoilProj[kk][1]*response,Weight);

	    histoManager.fillProf(metT[kk]+"Response",i,ZPt,-RecoilProj[kk][0]/ZPt,Weight);

	    histoManager.fillProf(metT[kk]+"CorPara",i,ZPt,-RecoilProj[kk][0],Weight);
	    histoManager.fillProf(metT[kk]+"CorPerp",i,ZPt,-RecoilProj[kk][1],Weight);

	    ostringstream os;
	    int nv = NVertex;
	    if(nv>4)
	      nv=4;
	    if(nv==0)
	      nv=1;

	    os << nv;
	    
	    histoManager.fill(metT[kk]+"RecoilPara_V"+os.str(),i,RecoilProj[kk][0],Weight);
	    histoManager.fill(metT[kk]+"RecoilPerp_V"+os.str(),i,RecoilProj[kk][1],Weight);
	    histoManager.fill(metT[kk]+"RecoilParaCor_V"+os.str(),i,RecoilProj[kk][0]*Corrections[kk],Weight);
	    histoManager.fill(metT[kk]+"RecoilPerpCor_V"+os.str(),i,RecoilProj[kk][1]*Corrections[kk],Weight);
	    

	    histoManager.fill2D(metT[kk]+"UparaVsQt_V"+os.str(),i,ZPt,RecoilProj[kk][0],Weight);
	    histoManager.fill2D(metT[kk]+"UparaVsQtCor_V"+os.str(),i,ZPt,RecoilProj[kk][0]*Corrections[kk],Weight);
	    histoManager.fill2D(metT[kk]+"UperpVsQt_V"+os.str(),i,ZPt,RecoilProj[kk][1],Weight);
	    histoManager.fill2D(metT[kk]+"UperpVsQtCor_V"+os.str(),i,ZPt,RecoilProj[kk][1]*Corrections[kk],Weight);

	    TVector2 rectmp(0,0);
	    rectmp.SetMagPhi(Recoils[kk][0],Recoils[kk][1]);
	    // cout<<Recoils[kk][0]<<"   "<<Recoils[kk][1]<<"    "<<rectmp.Px()<<"   "<<rectmp.Py()<<endl;
	    histoManager.fill(metT[kk]+"RecoilX"+os.str(),i,rectmp.Px(),Weight);
	    histoManager.fill(metT[kk]+"RecoilY"+os.str(),i,rectmp.Py(),Weight);

	  }  

	  // if(name[i].substr(0,4)=="data" && ZPt>=80 && ZPt<=90 && kk>3)
	  //  cout<<ZPt<<"   "<<RecoilProj[kk][1]*response<<"     "<<Event<<"    "<<Run<<"   "<<kk<<endl;

	  //Spec response as W
	  //histoManager.fill2D(metT[kk]+"RecoilvsqT",i,ZPt,-RecoilProj[kk][0],Weight);
	  // histoManager.fill2D(metT[kk]+"AbsRecoilvsqT",i,ZPt,fabs(RecoilProj[kk][0]),Weight);
	  
	  if(name[i].substr(0,4)=="data" && kk==5)
	    {
	   
	      //   if(RecoilProj[kk][1]*response<-70 && RecoilProj[kk][1]*response>-100  && ZPt>=80 && ZPt<=90 )
	      {
		  ostringstream os;
		  os << AbsEvent;
		
		  std::pair<int,int> tmp(Run,AbsEvent);
		  string t2(fileName); //string t1(sample);
		  std::pair<string,string> tmp2( os.str(), t2 );
		  std::pair< std::pair<int,int> , std::pair<string,string> > tmp3(tmp,tmp2);
		  Events.push_back(tmp3);
		}
	    }
	
	} //end loop MET

	
      } //End if condition >         
      
      if(name[i].substr(0,1)=="Z")
	{ S++; 
	}
      else if(name[i]!="data")
	{B++;/* GenRecoil[0]=0; GenRecoil[1]=0;*/ }
      
      NumberEntries[i]++;
	       
    }//End events
    if(name[i].substr(0,4)=="data") {
      cout<<endl;
      cout<<" -- "<<METXY[0][0]<<endl;
      cout<<" -+ "<<METXY[0][1]<<endl;
      cout<<" +- "<<METXY[1][0]<<endl;
      cout<<" ++ "<<METXY[1][1]<<endl;
      cout<<" And the Z "<<endl;
      cout<<" -- "<<ZXY[0][0]<<endl;
      cout<<" -+ "<<ZXY[0][1]<<endl;
      cout<<" +- "<<ZXY[1][0]<<endl;
      cout<<" ++ "<<ZXY[1][1]<<endl;

    }

     } //End datasets
 
  cout<<" End Filling , S="<<S<<" ;  B= "<<B<<endl;
  cout<<" detail "<<endl;
  for(int i=0;i<nt+1;i++) {
    cout<<" --->  "<<name[i]<<"  "<< NumberEntries[i]<<endl;
  }

  //Cleaning 2D Histos
  /* {
     histoManager.Cleaning2DHisto(NEntries );
     }*/

  Histos = histoManager.GetHistos();
  Histos2D = histoManager.GetHistos2D();
  Profiles = histoManager.GetProfiles();

  for(size_t i=0;i<Histos.size();i++)
    Weighted.push_back(false);

  cout<<"Histos "<<Histos.size()<<"   "<<Histos2D.size()<<endl;

  //  Delete();

}



void
ZUseTree::PrepareDatasets() {

  int nQCD=0;
  int nZ=0;

  //We
  // reposi=;
  if(!vetoE) { 
    //datasets.push_back("ttbar");
 
    /*   datasets.push_back("GamJet_PT1530");
    datasets.push_back("GamJet_PT3050");
    datasets.push_back("GamJet_PT5080");
    datasets.push_back("GamJet_PT80120");
    if(QCDType=="EM") {
      datasets.push_back("QCD_EM2030");
      datasets.push_back("QCD_EM3080");
      datasets.push_back("QCD_EM80170");
      datasets.push_back("QCD_bc2030");
      datasets.push_back("QCD_bc3080");
      datasets.push_back("QCD_bc80170");
      nQCD = 5;
    }
    else if(QCDType=="QCD"){
      cout<<" *** Using inclusive QCD"<<endl;
      datasets.push_back("QCD");
      nQCD = 1;
    }
    else {cout<<" bad QCD sample !"<<endl; return;}
  
    datasets.push_back("Z_2t_PU");
    datasets.push_back("W_en_PU");
    
    nZ=1;

    if(WType=="Pythia")
      { 
	datasets.push_back("Z_2e_PU"); //Z_2e_TuneZ2
      }
    else if(WType=="Powheg") {
      datasets.push_back("Z_2e_PU");
    }
    else if(WType=="Spec") {
      datasets.push_back("Z_2e_PU_v2");
    }
  
    //Combining
    if(QCDType=="EM") {
      int dttmp[20] ={0,1,1,2,2,2,2,2,3,3,4}; 
      for(int i=0;i<(nQCD+5);i++) {
	FillAddWeight(datasets[i]);
	dt.push_back(dttmp[i]);
      }
    }
    else if(QCDType=="QCD"){
      int dttmp[20] ={0,0,0,0,0,1,2,2,2,3,3,3,3,3,3};
      for(int i=0;i<(nQCD+14);i++) {
	FillAddWeight(datasets[i]);
	dt.push_back(dttmp[i]);
      }
    }
  
    nt=5;

    colors.push_back(kRed+1);
    colors.push_back(kMagenta+4);
    colors.push_back(kViolet-5);
    colors.push_back(kOrange+7); //+1
 
 
    colors.push_back(kOrange-2);
    colors.push_back(kBlack);//put to kBlack for PAS, otherwise kOrange+3
  
    string name_tmp[6]={"ttbar","#gamma+jet","QCD","EWK","Z #rightarrow ee","data"}; // "t#bar{t}",
    for(int i=0;i<nt+1;i++)
      name.push_back(name_tmp[i]);
*/
    colors.push_back(kBlack);
    
    Sorder.push_back("Z #rightarrow ee");
    Sorder.push_back("EWK");
    Sorder.push_back("QCD");
    Sorder.push_back("#gamma+jet");
    Sorder.push_back("ttbar");

    for(int i=0;i<datasets.size();i++)
      FillAddWeight(datasets[i]);

    name.push_back("data");
  }
  else {
    if(WType=="Pythia")
      {
	datasets.push_back("Z_2e_TuneZ2");

	int dttmp[7] ={0,1,0,0,0,0,1};
	for(int i=0;i<1;i++) {
	  FillAddWeight(datasets[i]);
	  dt.push_back(dttmp[i]);
	}
      }

    else if(WType=="Powheg") {
      datasets.push_back("Z_2e_PU");

      FillAddWeight(datasets[0]);
      dt.push_back(0);
    
    }   
    nt=1;

    colors.push_back(kOrange-2);
    //colors.push_back(kGreen-2);
    //  colors.push_back(kBlue-2);
    colors.push_back(kBlack);//put to kBlack for PAS, otherwise kOrange+3

    string name_tmp[3]={"Z #rightarrow ee","data"}; //
    for(int i=0;i<nt+1;i++)
      name.push_back(name_tmp[i]);

    Sorder.push_back("Z #rightarrow ee");

  }

  
}




TH1F* 
ZUseTree::GetVtxCorrection(string obs,int nds, int bin) {

  //Protection
  if(NVert!=1) {
    cout<<" Warning, no requirement of only 1 Vtx event !!!!!"<<endl;
  }
  
  float contamination = 0.05; //6.9 %
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
ZUseTree::GetFitWeight( ) {

  string metT[6]={"pf","tc","calo","caloT1","caloT2","pfT1"};
 
  // histoManager.PrepareHistos(name);

  vector<int> NumberEntries(nt+1,0);
  string Nvert;
 
  bool  fillBoth= false;

  cout<<" ==================================================== "<<endl;
  cout<<" ================= En cours de weight =============== "<<endl;

  //sÃ©lection
  for(int i=0;i<nt+1;i++) {

    
    //Déclaration des variables
    float IsoVar[2][3];
    float IDVar[2][6];
    float Mets[6][2];
    float MetProj[6][2];
    float Recoils[6][2];
    float RecoilProj[6][2];
    float Lepton[2][3];
    float EtSC[2];

    float dPhiRecoilZ[6];
    float dPhiMETZ[6];

    float SumET[6];

    bool ConvRej1;
    bool InHit1;
    bool ExpInHit1;
    bool ConvRej2;
    bool InHit2;
    bool ExpInHit2;

    int Run;
    int Event;
    int AbsEvent;
    //char sample[21];
    char fileName[21];

    float GenRecoil[2][2];
    float SpecRecoilProj[6][2];
    float CTRecoil[6][2];
    
    float Corrections[6];
    
    //Corrected MET variables
    float MetCor[6][2];
    float MetCorProj[6][2];
    float RecoilCor[6][2];
    float RecoilCorProj[6][2];

    float GenZPt;
    float ZPt;
    float ZMass;
    float DZvertices;
    //Underlying event
    float dPhiPhotonZ;
    float dPhiJetZ;
    int JetMultiplicity;
    int PhotonMultiplicity;
    float LeadingJet;
    float LeadingPhoton;
    
    //Decomposed MET SumET
    float caloSumETEB;
    float caloSumETEE;
    float caloSumETHB;
    float caloSumETHE;
    float caloSumETHF;
    float pfSumETEB;
    float pfSumETEE;
    float pfSumETHB;
    float pfSumETHE;
    float pfSumETHF;

    int NVertex;

    float Weight=0;

    bool Sel[2]={true,true};
    int epart[2]={-1,-1};

    tChains[i]->SetBranchAddress("NVertex",&NVertex);

    tChains[i]->SetBranchAddress("IsoVar",IsoVar);
    tChains[i]->SetBranchAddress("IDVar",IDVar);
    tChains[i]->SetBranchAddress("Lepton",Lepton);
    tChains[i]->SetBranchAddress("EtSC",EtSC);
 
    tChains[i]->SetBranchAddress("SumET",SumET);

    tChains[i]->SetBranchAddress("Mets",Mets);
    tChains[i]->SetBranchAddress("MetProj",MetProj);
 
    //    tChains[i]->SetBranchAddress("Corrections",Corrections);
    
    tChains[i]->SetBranchAddress("BasicRecoils",Recoils);
    tChains[i]->SetBranchAddress("BasicRecoilProj",RecoilProj);
    tChains[i]->SetBranchAddress("RecoilProj",SpecRecoilProj);

    tChains[i]->SetBranchAddress("BasicdPhiRecoilZ",dPhiRecoilZ);
    tChains[i]->SetBranchAddress("dPhiMETZ",dPhiMETZ);
  
    tChains[i]->SetBranchAddress("InHit1",&InHit1);
    tChains[i]->SetBranchAddress("ExpInHit1",&ExpInHit1);
    tChains[i]->SetBranchAddress("ConvRej1",&ConvRej1);
    tChains[i]->SetBranchAddress("InHit2",&InHit2);
    tChains[i]->SetBranchAddress("ExpInHit2",&ExpInHit2);
    tChains[i]->SetBranchAddress("ConvRej2",&ConvRej2);

    tChains[i]->SetBranchAddress("AbsEvent",&AbsEvent);
    tChains[i]->SetBranchAddress("Run",&Run);
    //    tChains[i]->SetBranchAddress("sample",&sample);
    tChains[i]->SetBranchAddress("fileName",&fileName);
    tChains[i]->SetBranchAddress("Event",&Event);

    tChains[i]->SetBranchAddress("ZPt",&ZPt); //ZMass
    tChains[i]->SetBranchAddress("ZMass",&ZMass); 
    tChains[i]->SetBranchAddress("DZvertices",&DZvertices);

    //MET Corrected variables
    tChains[i]->SetBranchAddress("MetCor",MetCor);
    tChains[i]->SetBranchAddress("MetCorProj",MetCorProj);
    tChains[i]->SetBranchAddress("RecoilCor",RecoilCor);
    tChains[i]->SetBranchAddress("RecoilCorProj",RecoilCorProj);

    //Underlying event
    tChains[i]->SetBranchAddress("LeadingJet",&LeadingJet);
    tChains[i]->SetBranchAddress("LeadingPhoton",&LeadingPhoton);
    tChains[i]->SetBranchAddress("PhotonMultiplicity",&PhotonMultiplicity);
    tChains[i]->SetBranchAddress("JetMultiplicity",&JetMultiplicity);
    tChains[i]->SetBranchAddress("dPhiJetZ",&dPhiJetZ);
    tChains[i]->SetBranchAddress("dPhiPhotonZ",&dPhiPhotonZ);
	
    //Decomposed MET
    tChains[i]->SetBranchAddress("CaloSumETEB",&caloSumETEB);
    tChains[i]->SetBranchAddress("CaloSumETEE",&caloSumETEE);
    tChains[i]->SetBranchAddress("CaloSumETHB",&caloSumETHB);
    tChains[i]->SetBranchAddress("CaloSumETHE",&caloSumETHE);
    tChains[i]->SetBranchAddress("CaloSumETHF",&caloSumETHF); 
    tChains[i]->SetBranchAddress("PfSumETEB",&pfSumETEB);
    tChains[i]->SetBranchAddress("PfSumETEE",&pfSumETEE);
    tChains[i]->SetBranchAddress("PfSumETHB",&pfSumETHB);
    tChains[i]->SetBranchAddress("PfSumETHE",&pfSumETHE);
    tChains[i]->SetBranchAddress("PfSumETHF",&pfSumETHF);

    if(i==nt-1) {
      tChains[i]->SetBranchAddress("GenRecoil",GenRecoil);
      tChains[i]->SetBranchAddress("CTRecoil",CTRecoil);
      tChains[i]->SetBranchAddress("GenZPt",&GenZPt);
    }
    else {
      for(int jj=0;jj<6;jj++) {
	CTRecoil[jj][0]=0; CTRecoil[jj][1]=0;
      }
    }
    
    int EP=0;
    cout<<" beginning tree  " <<name[i]<<"   "<<tChains[i]->GetEntries()<<endl;
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
      
      if(name[i].substr(0,4)=="data")
	{
	  fillBoth= false;
	  bool doubleCount=false;
	  for(size_t ik=0;ik<Events.size();ik++)
	    {
	      if( ((Events[ik]).first.first) == Run && ((Events[ik]).first.second) == Event)
		{doubleCount=true; break;}
	    }
	  
	  if(doubleCount || (EventFilter && Run<EventNum ) )
	    { continue; }
	}

      //Filling Histos
      
      //Passing selection
      for(int id=0;id<2;id++) {

	if(fabs(Lepton[id][1])>1.479)
	  {EP=1; }
	else
	  {EP=0; }
	
	if(name[i]!="data" || EP==0 || noDeltaCor) {
	  IDVar[id][4] = 0;
	  IDVar[id][5] =0;
	}
	
	Sel[id]=true;

	//first pt cut
	if(EtSC[id] < PTcut)  Sel[id]=false;
	
	if(convRej && !ExpInHit1 && !ExpInHit2)   Sel[id]=false;
	
	//Met Cut
	
	//And Selections
	//  cout<<IDVar[4]<<endl;
	if(IDVar[id][0] > _IDcuts[EP][0] )  Sel[id]=false;
	if( fabs(IDVar[id][1]-IDVar[id][4]) > _IDcuts[EP][1] )  Sel[id]=false;
	if( fabs(IDVar[id][2]-IDVar[id][5]) > _IDcuts[EP][2] )  Sel[id]=false;
	if(IDVar[id][3] > _IDcuts[EP][3] )  Sel[id]=false;
	
	if(IsoVar[id][0]/Lepton[id][0] > _Isocuts[EP][0] )  Sel[id]=false;
	if(IsoVar[id][1]/Lepton[id][0] > _Isocuts[EP][1] )  Sel[id]=false;
	if(IsoVar[id][2]/Lepton[id][0] > _Isocuts[EP][2] )  Sel[id]=false;
	
	/*	if(MultiWeight) {
		if(i<=1 && i!=nt) {
		if(NVert==0)
		Weight = (Weights[0][0][1]+Weights[0][1][1])*GetWeight(i, ie);
		else if(NVert==1)
		Weight = Weights[0][0][1]*GetWeight(i, ie);
		else if(NVert==2)
		Weight = Weights[0][1][1]*GetWeight(i, ie);
		}
		else if(i>1 && i!=nt) {
		if(NVert==0)
		Weight = (Weights[0][0][0]+Weights[0][1][0])*GetWeight(i, ie);
		else if(NVert==1)
		Weight = Weights[0][0][0]*GetWeight(i, ie);
		else if(NVert==2)
		Weight = Weights[0][1][0]*GetWeight(i, ie);
		}
		else if(i==nt)
		Weight = 1;
		}*/
	//	Weight = 1; //FIXME

	if(EcalP!=EP && EcalP!=2 && EcalP!=3) Sel[id]=false;
	
	if(EcalP==3)  
	  epart[id]=EP;
	  

	/*	if(NVert==2 && NVertex<2 && name[i].substr(0,4)=="data")
		Sel[id]=false;
		else if(NVert==1 && NVertex>1 && name[i].substr(0,4)=="data")
		Sel[id]=false;*/
	


	if(!Sel[id]) break;

      }
      
      if(!Sel[0] || !Sel[1])
	continue;
     
      if(EcalP==3 && epart[0]==epart[1] ) continue;
      
      //  if(JetMultiplicity == 0 ) continue;

      if(NVertex==1 && i>=nt-2/*name[i].substr(0,4)=="data"*/)
	Nvert="1";
      else if(NVertex>1 && i>=nt-2/*name[i].substr(0,4)=="data"*/)
	Nvert="2";
      else if(i<nt-2/*name[i].substr(0,4)!="data"*/)
	fillBoth=true;
      
      if(i>=nt-2) {
	fillBoth=false;
      }
      
      for(int kk=0;kk<6;kk++) {

	Corrections[kk] = SearchESWeight(ZPt, kk);
	  
	if( (Mets[kk][0] > METcut && !invCut && ZMass > MTCut && ZMass < 120) || 
	    (Mets[kk][0] < METcut && invCut && ZMass > MTCut ) ) {
	
	  if(!fillBoth) {
		    
	    
	    
	    histoManager.fill(metT[kk]+"MET_V"+Nvert,i,Mets[kk][0],Weight);
	    TVector2 mettmp(0,0);
	    mettmp.SetMagPhi(Mets[kk][0],Mets[kk][1]);
	    histoManager.fill(metT[kk]+"METX_V"+Nvert,i,Mets[kk][0],Weight);

	    //Another frame!
	    TVector2 noiseAxis(0.6,-0.2);
	    noiseAxis = noiseAxis.Unit();
	    TVector2 noiseAxisT = noiseAxis.Rotate(acos(-1)/2. );
	    histoManager.fill(metT[kk]+"METSFramePara_V"+Nvert,i,noiseAxis*mettmp,Weight);
	    histoManager.fill(metT[kk]+"METSFramePerp_V"+Nvert,i,noiseAxisT*mettmp,Weight);

	    if(mettmp.Py() > mettmp.Px()*3) //upper
	      histoManager.fill(metT[kk]+"METUpper_V"+Nvert,i,Mets[kk][0],Weight);
	    else //lower
	      histoManager.fill(metT[kk]+"METLower_V"+Nvert,i,Mets[kk][0],Weight);

	    TVector2 mettmpCor(0,0);
	    mettmpCor.SetMagPhi(MetCor[kk][0],MetCor[kk][1]);
	    histoManager.fill(metT[kk]+"METXCor_V"+Nvert,i,mettmpCor.Px(),Weight);
	    
	    histoManager.fill(metT[kk]+"SumEt_V"+Nvert,i,SumET[kk]-Lepton[0][0]-Lepton[1][0],Weight);

	    if(kk==0) {
	      histoManager.fill("ZPt_V"+Nvert,i,ZPt,Weight);
	      histoManager.fill("pfSumETHF_V"+Nvert,i,caloSumETEE+caloSumETHF,Weight);
	    }
	  
	    if(NVertex>1 && i>=nt-2/*name[i].substr(0,4)=="data"*/) {
	   
	      histoManager.fill(metT[kk]+"METControl",i,Mets[kk][0],Weight);
	      histoManager.fill(metT[kk]+"METPhiControl",i,Mets[kk][1]*180/3.1415,Weight);
	      histoManager.fill(metT[kk]+"METParaControl",i,MetProj[kk][0],Weight);
	      histoManager.fill(metT[kk]+"METPerpControl",i,MetProj[kk][1],Weight);
	    
	      histoManager.fill(metT[kk]+"diffMetControl",i,(MetCor[kk][0]-Mets[kk][0])/Mets[kk][0],Weight);
	      histoManager.fill(metT[kk]+"METLepCorControl",i,MetCor[kk][0],Weight);
	      histoManager.fill(metT[kk]+"METParaLepCorControl",i,-MetCorProj[kk][0],Weight);
	      histoManager.fill(metT[kk]+"METPerpLepCorControl",i,-MetCorProj[kk][1],Weight);
	      histoManager.fill(metT[kk]+"RecoilLepCorControl",i,RecoilCor[kk][0],Weight);
	      histoManager.fill(metT[kk]+"RecoilParaLepCorControl",i,RecoilCorProj[kk][0],Weight);
	      histoManager.fill(metT[kk]+"RecoilPerpLepCorControl",i,RecoilCorProj[kk][1],Weight);
	      histoManager.fill(metT[kk]+"dPhiRecoilControl",i,fabs(dPhiRecoilZ[kk]*3.14/180),Weight);
	      histoManager.fill(metT[kk]+"dPhiRecoilCorControl",i,fabs(dPhiRecoilZ[kk]*3.14/180),Weight);

	      float response = SearchESWeight(ZPt, kk);

	      histoManager.fill(metT[kk]+"RecoilCorControl",i,Recoils[kk][0]*response,Weight);
	      histoManager.fill(metT[kk]+"RecoilParaCorControl",i,RecoilProj[kk][0]*response,Weight);
	      histoManager.fill(metT[kk]+"RecoilPerpCorControl",i,RecoilProj[kk][1]*response,Weight);
	      histoManager.fill(metT[kk]+"RecoilControl",i,Recoils[kk][0],Weight);
	      histoManager.fill(metT[kk]+"RecoilParaControl",i,RecoilProj[kk][0],Weight);
	      histoManager.fill(metT[kk]+"RecoilPerpControl",i,RecoilProj[kk][1],Weight);

	      histoManager.fill(metT[kk]+"SumEtControl",i,SumET[kk]-Lepton[0][0]-Lepton[1][0],Weight);
	    

	      histoManager.fill(metT[kk]+"METXControl",i,mettmp.Px(),Weight);
	      histoManager.fill(metT[kk]+"METYControl",i,mettmp.Py(),Weight);

	      histoManager.fill(metT[kk]+"METXCorControl",i,mettmp.Px(),Weight);

	      if(ZPt>1) {
		histoManager.fill2D(metT[kk]+"UparaVsQtControl",i,ZPt,RecoilProj[kk][0],Weight);
		histoManager.fill2D(metT[kk]+"UparaVsQtCorControl",i,ZPt,RecoilProj[kk][0]*Corrections[kk],Weight);
		histoManager.fill2D(metT[kk]+"UperpVsQtControl",i,ZPt,RecoilProj[kk][1],Weight);
		histoManager.fill2D(metT[kk]+"UperpVsQtCorControl",i,ZPt,RecoilProj[kk][1]*Corrections[kk],Weight);
	      }

	      if(kk==0) {
		histoManager.fill("caloSumETEBControl",i,caloSumETEB,Weight);
		histoManager.fill("caloSumETEEControl",i,caloSumETEE,Weight);
		histoManager.fill("caloSumETHBControl",i,caloSumETHB,Weight);
		histoManager.fill("caloSumETHEControl",i,caloSumETHE,Weight);
		histoManager.fill("caloSumETHFControl",i,caloSumETHF,Weight);
		histoManager.fill("pfSumETEBControl",i,pfSumETEB,Weight);
		histoManager.fill("pfSumETEEControl",i,pfSumETEE,Weight);
		histoManager.fill("pfSumETHBControl",i,pfSumETHB,Weight);
		histoManager.fill("pfSumETHEControl",i,pfSumETHE,Weight);
		histoManager.fill("pfSumETHFControl",i,pfSumETEE,Weight);
	      }

	      //ZPt
	      histoManager.fill("ZPtControl",i,ZPt,Weight);
	      histoManager.fill("ZMassControl",i,ZMass,Weight);

	      //Corrected mass ===========================
	      float Et[2]={Lepton[0][0],Lepton[1][0]};
	      for(int jj=0;jj<2;jj++) {
		if(fabs(Lepton[jj][1])>1.5)
		  Et[jj]*=1.027;
		else
		  Et[jj]*=1.008;
	      }
	      if(i==nt)
		histoManager.fill("ZMassCorControl",i,ComputeMinv(Et[0],Et[1],Lepton[0][1],Lepton[1][1],Lepton[0][2],Lepton[1][2]) ,Weight);
	      else
		histoManager.fill("ZMassCorControl",i,ZMass,Weight);
	      //=========================================
	  
	    }
	  }
	  else {
	    
	    if(kk==0) {
	      histoManager.fill("ZPt_V1",i,ZPt,Weight);
	      //    histoManager.fill("pfSumETHF_V1",i,caloSumETEE+caloSumETHF,Weight);
	    }

	    TVector2 mettmp(0,0);
	    mettmp.SetMagPhi(Mets[kk][0],Mets[kk][1]);
	    TVector2 noiseAxis(0.6,-0.2);
	    noiseAxis = noiseAxis.Unit();
	    TVector2 noiseAxisT = noiseAxis.Rotate(acos(-1)/2. );
	    histoManager.fill(metT[kk]+"METSFramePara_V1",i,noiseAxis*mettmp,Weight);
	    histoManager.fill(metT[kk]+"METSFramePerp_V1",i,noiseAxisT*mettmp,Weight);

	    if(mettmp.Py() > mettmp.Px()*3) //upper
	      histoManager.fill(metT[kk]+"METUpper_V1",i,Mets[kk][0],Weight);
	    else //lower
	      histoManager.fill(metT[kk]+"METLower_V1",i,Mets[kk][0],Weight);


	    histoManager.fill(metT[kk]+"SumEt_V1",i,SumET[kk]-Lepton[0][0]-Lepton[1][0],Weight);
	    histoManager.fill(metT[kk]+"SumEt_V2",i,SumET[kk]-Lepton[0][0]-Lepton[1][0],Weight);
	    
	    //MET
	    histoManager.fill(metT[kk]+"MET_V1",i,Mets[kk][0],Weight);
	    histoManager.fill(metT[kk]+"MET_V2",i,Mets[kk][0],Weight);
	    //MET Proj
	    histoManager.fill(metT[kk]+"METPara_V1",i,MetProj[kk][0],Weight);
	    histoManager.fill(metT[kk]+"METPara_V2",i,MetProj[kk][0],Weight);

	    histoManager.fill(metT[kk]+"METPerp_V1",i,MetProj[kk][1],Weight);

	  }
	}
      }

      if(name[i].substr(0,4)=="data")
	{
	  std::pair<int,int> tmp(Run,Event);
	  string t2(fileName);  //,t1(sample);
	  std::pair<string,string> tmp2("data",t2);
	  std::pair< std::pair<int,int> , std::pair<string,string> > tmp3(tmp,tmp2);
	  Events.push_back(tmp3);
	}
      
    }//End events

  }// End datasets
 
    
  Events.clear();

  Histos = histoManager.GetHistos();
    
  //Compute fraction of events containing PU 
  int S1 =  histoManager.FindNVar("pfMET_V1");
  int S2 =  histoManager.FindNVar("pfMET_V2");

  int n1= histoManager.access(S1,nt);
  int n2= histoManager.access(S2,nt); cout<<" s "<<S1<<"  "<<S2<<"    "<<n1<<"   "<<n2<<endl;
  TH1F* h1 = (TH1F*)Histos[ n1 ]->Clone();
  TH1F* h2 = (TH1F*)Histos[ n2 ]->Clone();
  
  float integ1 = h1->Integral(0, h1->GetNbinsX() +1 );
  float integ2 = h2->Integral(0, h2->GetNbinsX() +1 );

  cout<<" Compute fraction of events containing PU "<<endl;
  cout<<" Total "<<integ1+integ2<<endl;
  cout<<" fraction 1 vtx "<<integ1<<" ----> "<<integ1/(integ1+integ2)<<endl;
  cout<<" fraction 2 vtx "<<integ2<<" ----> "<<integ2/(integ1+integ2)<<endl;
  Frac1V = integ1/(integ1+integ2);
  //======================================

  vector<vector<vector<float> > > Wghts(6,vector<vector<float> >(2,vector<float>(2,1)));
  /* if(!Remove2V) {
     for(int i=0;i<2;i++) {
     for(int kk=0;kk<5;kk++) {
     ostringstream os;
     os << i+1;
     TString Name = (TString)metT[kk]+"MET_V"+ os.str();
     Weights[kk][i]= GetAlpha(Name);
     }
     }
     }*/  
    
  return Wghts;
    
    
}
 

float ZUseTree::ComputeMinv(float Et1, float Et2, float eta1, float eta2, float phi1, float phi2) {

  float Minv=0;

  float E1=ConversionPtE(Et1,eta1);
  float E2=ConversionPtE(Et2,eta2);

  float l1[3]={E1,eta1,phi1};
  float l2[3]={E2,eta2,phi2};
  
  float x1,y1,z1,x2,y2,z2;

  Conversion_REP_carte(l1,x1,y1,z1);
  Conversion_REP_carte(l2,x2,y2,z2);
  TVector3 L1(x1,y1,z1);
  TVector3 L2(x2,y2,z2);
  
  TVector2 p2_(0,0);
  float pz_(0);
  float E_(0);
  
  p2_ = L1.XYvector() + L2.XYvector();
  E_  = E1+E2;
  pz_ = z1+z2;
  TVector3 p3_( p2_.X(), p2_.Y(), pz_ );
  float m2_= pow(E_,2)-p3_.Mag2();
  if( m2_<0 ) m2_=0;
  Minv = sqrt(m2_);

  return Minv;

}




float ZUseTree::ConversionEtaTheta(float eta) //Convertion eta en theta
{
  float theta;
  theta = 2*atan(exp(-eta));
  return theta;
}

float ZUseTree::ConversionThetaEta(float theta) //Conversion theta en eta
{
  float eta;
  eta = -log(tan((theta)/2));
  return eta;
}



float ZUseTree::ConversionEpt(float energy,float eta ) //Convertion Energie en pt
{

  float pt =  energy*sin(ConversionEtaTheta(eta));
  return pt;
}

float ZUseTree::ConversionPtE(float pt, float eta) //Conversion pt en energie
{
  float energy = pt/sin(ConversionEtaTheta(eta));
  return energy;

}

void ZUseTree::Conversion_REP_carte(float coord[3], float& x, float& y, float &z) //Conversion Rho/eta/phi en coordonnées cartésiennes
{
  float Theta = ConversionEtaTheta(coord[1]);
  x = coord[0]*sin(Theta)*cos(coord[2]);
  y = coord[0]*sin(Theta)*sin(coord[2]);
  z = coord[0]*cos(Theta);
  return;
}



/*
  void ZUseTree::FillZResponse(string obs, string observable2 ) {

  // string observable="tcAbsRecoilvsGen";
  // string observable2="pfAbsSpecRecoilvsGen";

  int Obs = histoManager.FindNVar2D(obs);
      
  if(Obs==-1) {
  cout<<" Be careful, no such observable "<<obs<<endl;
  abort();
  }

  XTitle = (histoManager.FindLeg2D( obs )).first;
  YTitle = (histoManager.FindLeg2D( obs )).second;
  TH1F* emptyHisto= new TH1F(("bidon"+obs).c_str(),"bidon",10000,-10000,10000);
  emptyHisto->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
  emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);
  emptyHisto->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
  emptyHisto->GetYaxis()->SetNdivisions(Ydiv[0],Ydiv[1],Ydiv[2]);
  emptyHisto->GetXaxis()->SetTitle(XTitle.c_str());
  emptyHisto->GetYaxis()->SetTitle(YTitle.c_str());
  emptyHisto->GetYaxis()->SetTitleOffset(1.1);
  emptyHisto->SetFillColor(0);
  emptyHisto->SetLineColor(0);
         
  int NW = histoManager.access2D(Obs,nt);
  TGraph* WGraph= Histos2D[ NW ];
      
  double x,y;
  for(int i=0;i<WGraph->GetN();i++)
  {
  WGraph->GetPoint(i,x,y);
  //	  cout<<" x "<<x<<" y "<<y<<endl;
  }

  vector<float> temp = histoManager.GetTemplate2D(obs);
  
  TGraphAsymmErrors* Wprof =  histoManager.GraphReduction(WGraph,temp[0],temp[1],temp[2],"ZR","");
  TGraphAsymmErrors* Wprof2 = histoManager.GraphReduction(WGraph,temp[0],temp[1],temp[2],"ZRE","s");
  TPolyLine* polline = PolyLineFromGraph( (TGraph*)Wprof2 );
  polline->SetFillColor(kYellow-9);
      
  Wprof->SetLineWidth(2);
  Wprof->SetLineColor(2);
  Wprof->SetMarkerColor(2);
  Wprof->GetYaxis()->SetRangeUser(-42,42);
  emptyHisto->GetYaxis()->SetRangeUser(-42,50);
  emptyHisto->GetXaxis()->SetRangeUser(-42,42);
  emptyHisto->GetXaxis()->SetTitleOffset(1);
  Wprof->GetXaxis()->SetLimits(-42,42);
  emptyHisto->GetXaxis()->SetTitle("Generated Z p_{T,||} [GeV]");
  Wprof->GetYaxis()->SetTitle(YTitle.c_str() );
  TLine* line = new TLine(-42,-42,42,42);
  TLine* line4 = new TLine(-42,42,0,0);
  TLine* line2 = new TLine(-42,0,42,0);
  TLine* line3 = new TLine(0,-50,0,50);

  //Prepare the pads
  TPad* pad_[2];
  size_t n_ = 2;
  Double_t e_ = 7; 
  Double_t t_ = 30;
  Double_t b_ = 80;
  Double_t g_ = 85;
  Double_t d_ = 30;
  Double_t h_ = 480;
  Double_t w_ = 480;
  Double_t hh_ = b_ + h_ + e_ ;
  Double_t W_ = g_ + w_ + d_;
  Double_t H_ = t_ + h_*n_ + 2*(n_-1)*e_ + b_ ;
  c2=new TCanvas("c2","Test",300,300,W_,H_);
  Double_t xlow_= 0.;
  Double_t ylow_=0.;
  Double_t xup_=1.;
  Double_t yup_=0.;
  for(size_t i=0;i<n_;i++) {

  TString padname_("pad_");
  padname_+=i;
  yup_ = ylow_ + hh_ /H_;
  pad_[i] = new TPad( padname_, padname_, 
  xlow_, ylow_, xup_, yup_,
  kWhite,0,0);
  ylow_ += (h_ + 2*e_)/H_;
  pad_[i]->SetLeftMargin(  g_/W_ );
  pad_[i]->SetRightMargin( d_/W_ );
  pad_[i]->SetTopMargin(  e_/hh_ );
  pad_[i]->SetBottomMargin( b_/hh_ );
  pad_[i]->SetFillColor(0);
  pad_[i]->SetTickx(1);
  pad_[i]->SetTicky(1);
  pad_[i]->SetFrameFillStyle(0);
  pad_[i]->SetFrameLineWidth(Wline);
  pad_[i]->SetFrameBorderMode(0);
  }

  pad_[1]->Draw();
  pad_[1]->cd();
	  
  emptyHisto->GetXaxis()->SetLabelSize(0);
  emptyHisto->GetXaxis()->SetTitleSize(0);
  emptyHisto->Draw();
  polline->Draw("Fsame");
  Wprof->Draw("PE");
  line->SetLineColor(kGreen+2);
  line->SetLineWidth(2);
  line4->SetLineColor(kGreen+2);
  line4->SetLineWidth(2);
  line->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  // line4->Draw("same");
    
      
  //*******************************************
  c2->cd();
  pad_[0]->Draw();
  pad_[0]->cd();
  int Obs2 = histoManager.FindNVar2D(observable2);
      
  if(Obs2==-1) {
  cout<<" Be careful, no such observable "<<observable2<<endl;
  abort();
  }

  XTitle = (histoManager.FindLeg2D( observable2 )).first;
      
  int NW2 = histoManager.access2D(Obs2,nt);
  TGraph* WGraph2= Histos2D[ NW2 ];
      
  vector<float> temp2 = histoManager.GetTemplate2D(observable2);
  TGraphAsymmErrors* WprofGen2 =  histoManager.GraphReduction(WGraph2,temp2[0],temp2[1],temp2[2],"WGen3","");
  TGraphAsymmErrors* WprofGen3 = histoManager.GraphReduction(WGraph2,temp2[0],temp2[1],temp2[2],"WGenE2","s");
  TPolyLine* polline2 = PolyLineFromGraph( (TGraph*)WprofGen3 );
  polline2->SetFillColor(kYellow-9);

  TH1F* emptyHisto2 = (TH1F*)emptyHisto->Clone();
  emptyHisto2->GetXaxis()->SetLabelSize(0.05);
  emptyHisto2->GetXaxis()->SetTitleSize(0.05);
  emptyHisto2->GetYaxis()->SetTitleSize(0);
  emptyHisto2->Draw();
  polline2->Draw("Fsame");
  WprofGen2->SetLineWidth(2);
  WprofGen2->SetLineColor(4);
  WprofGen2->SetMarkerStyle(29);
  WprofGen2->SetMarkerColor(4);
  WprofGen2->Draw("PE");

  TLegend* legend=new TLegend(0.6,0.6,0.8,0.8);
  line->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  //    line4->Draw("same");
  legend->SetFillColor(0);
  legend->SetLineColor(0);
  legend->SetTextSize(0.04);

  legend->AddEntry(Wprof,"EGamma electron","lp"); //"EGamma electron"
  legend->AddEntry(WprofGen2,"CaloTowers","lp"); //"pf electron"
  legend->AddEntry(polline2,"CL 68%","f");
  legend->Draw("same");
     
      
  cout<< " Comparaison****** "<<endl;
  double x1,y1,x2,y2;
  for(int i=0;i<WprofGen2->GetN();i++) {
  Wprof->GetPoint(i,x1,y1);
  WprofGen2->GetPoint(i,x2,y2);

  cout<<" x "<<x1<<"<>"<<x2<<"  y "<<y1<<"<>"<<y2<<endl;
      
  }

  }
*/



void
ZUseTree::LoadDBWeight() {

  ifstream in("/home/mmarionn/Documents/CMS/CMSDatabase/Z_Resbos.txt", ios::in );  

  vector<vector<float> > tmpw;//(18,vector<float>(3,0));

  if(in) {
    cout<<" Loading Database "<<endl;
    while(!in.eof()) {
      vector<float> tmpv(3,0);
      in >> tmpv[0] >> tmpv[1] >> tmpv[2];
      tmpw.push_back(tmpv);
    }
  }
  else 
    cout<<" No DB loaded !! "<<endl;
  
  DBWeights = tmpw;

}

void ZUseTree::LoadResponseDB() {

  ifstream file("/home/mmarionn/Documents/CMS/CMSDatabase/ResponseZMC_v2.db", ios::in);
  ifstream fileD("/home/mmarionn/Documents/CMS/CMSDatabase/ResponseZData.db", ios::in);

  vector<vector<float> > vtmp2;
  vector<float> vtmp(3,0);
  for(int i=0;i<6;i++) {
    MCResponseCor.push_back(vtmp2);
    DataResponseCor.push_back(vtmp2);
  }

  if(file) {

    int met, m, M;
    double cor;

    while(!file.eof()) {
      file >> met >> m >> M >> cor;
      vtmp[0] = m;
      vtmp[1] = M;
      vtmp[2] = cor;
      MCResponseCor[met].push_back(vtmp);
    }
  }
  else cout<<" missing db!!! "<<endl;


  if(fileD) {

    int met, m, M;
    double cor;

    while(!fileD.eof()) {
      fileD >> met >> m >> M >> cor;
      vtmp[0] = m;
      vtmp[1] = M;
      vtmp[2] = cor;
      DataResponseCor[met].push_back(vtmp);
    }
  }
  else cout<<" missing db!!! "<<endl;
}


float 
ZUseTree::SearchWeight(float Zpt) {
  
  float w=1;
  
  for(int i=0;i<DBWeights.size();i++) {
    if(Zpt >= DBWeights[i][0] && Zpt < DBWeights[i][1])
      { w = DBWeights[i][2]; break;}
  }

  return w;
}

float 
ZUseTree::SearchESWeight(float Zpt, int met) {
  
  float w=1;
  
  for(int i=0;i<MCResponseCor[met].size();i++) {
    if(Zpt >= MCResponseCor[met][i][0] && Zpt < MCResponseCor[met][i][1])
      { w = MCResponseCor[met][i][2]; break;}
  }
  return w;
}


void
ZUseTree::PlotResolutionPerVtx(string met, string type,bool isMC) {


  string obs = met+"U"+type+"VsQt_V";

  
  c2=new TCanvas("c2","Test",300,300,600,600);

  TH1F* emptyHisto= new TH1F(("bidon"+obs).c_str(),"bidon",1000,0,1000);
  emptyHisto->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
  emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);
  emptyHisto->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
  emptyHisto->GetYaxis()->SetNdivisions(Ydiv[0],Ydiv[1],Ydiv[2]);
  emptyHisto->GetXaxis()->SetTitle(XTitle.c_str());
  emptyHisto->GetYaxis()->SetTitle(YTitle.c_str());
  emptyHisto->SetFillColor(0);
  emptyHisto->SetLineColor(0);
  emptyHisto->Fill(0.,0.);
  emptyHisto->Draw();

  int Nbin=26;
  double BinVar[27] = {0,10,20,30,40,50,60,80,100,150,200,250,300,350,
		       400,450,500,550,600,650,700,750,800,850,900,950,1000};

  int mType[4] = {29,20,21,22};//,23};
  //int mTypeMC[4] = {30,24,25,26};//,32};
  int mcolor[4]={kBlue+1,kGreen+1,kOrange-2,/*kOrange+7,*/kRed+1};


  TLegend* legend=new TLegend(0.6,0.6,0.8,0.8);
  legend->SetFillColor(0);
  legend->SetLineColor(0);
  legend->SetTextSize(0.04);
  
  for(int iv=1;iv<5;iv++) {

    ostringstream os;
    os<< iv;

    int Obs =  histoManager.FindNVar2D(obs+ os.str() );
    int Ndata= histoManager.access2D(Obs,nt);
    if(isMC)
      Ndata =  histoManager.access2D(Obs,nt-1);

    // int NSignal= histoManager.access2D(Obs,nt-1);

    TH2F* HData = (TH2F*)Histos2D[ Ndata ]->Clone();
    TH1D* RMSData = histoManager.ConvertHistoToRMSPlotVarBin(HData,Nbin,BinVar,obs+os.str(),Remove2V,false);

    // TH2F* HSignal = (TH2F*)Histos2D[ NSignal ]->Clone();
    // TH1D* RMSSignal = histoManager.ConvertHistoToRMSPlotVarBin(HSignal,Nbin,BinVar,obs+os.str(),Remove2V);

    RMSData->SetLineWidth(LineWidth);
    RMSData->SetMarkerSize(MarkerSize);
    RMSData->SetMarkerStyle(mType[iv-1]);
    RMSData->SetMarkerColor(mcolor[iv-1]);
    RMSData->SetLineColor(mcolor[iv-1]);

    /* RMSSignal->SetLineWidth(LineWidth);
    RMSSignal->SetMarkerSize(MarkerSize);
    RMSSignal->SetMarkerStyle(mTypeMC[iv-1]);
    RMSSignal->SetMarkerColor(mcolor[iv-1]);
    RMSSignal->SetLineColor(mcolor[iv-1]);
    RMSSignal->SetLineStyle(2);*/

    RMSData->Draw("same");
    //RMSSignal->Draw("same");
    legend->AddEntry(RMSData, (met+"RMS "+type+"  "+os.str()+" Vtx").c_str(), "lp");
    // legend->AddEntry(RMSSignal, ("MC: "+met+"RMS "+type+"  "+os.str()+" Vtx").c_str(), "lp");
  }

  legend->Draw("same");
  cmsPrel(Lumi);
  
}

void
ZUseTree::PlotGlobalResolutionPerVtx( string type,bool isMC) {

  
  string mets[4]={"caloT2","tc","pf","pfT1"};
  int colors[4]={kOrange+7, 896,38, kBlue+1};
  int mType[4] = {29,20,21,22};

  c2=new TCanvas("c2","Test",300,300,600,600);

 
  TLegend* legend=new TLegend(0.6,0.6,0.8,0.8);
  legend->SetFillColor(0);
  legend->SetLineColor(0);
  legend->SetTextSize(0.04);


  TH1F* emptyHisto=new TH1F("emp","",4,0,5);
  emptyHisto->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
  emptyHisto->GetYaxis()->SetTitle(("RMS(u_{"+type+"} )").c_str());
  emptyHisto->GetYaxis()->SetTitleOffset(1.1);
  emptyHisto->GetXaxis()->SetTitle("Number of vertices");
  emptyHisto->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
  emptyHisto->GetYaxis()->SetNdivisions(Ydiv[0],Ydiv[1],Ydiv[2]);
  emptyHisto->Draw();

  for(int im=0;im<4;im++) {

    string obs = mets[im]+"Recoil"+type+"_V";
    TGraphErrors* RMSPerVtx= new TGraphErrors(4);
    RMSPerVtx->SetName(obs.c_str() );

   

    for(int iv=1;iv<5;iv++) {

      ostringstream os;
      os<< iv;

      int Obs =  histoManager.FindNVar(obs+ os.str() );
      int Ndata= histoManager.access(Obs,nt);
      if(isMC)
	Ndata =  histoManager.access(Obs,nt-1);

      TH1F* HData = (TH1F*)Histos[ Ndata ]->Clone();

      RMSPerVtx->SetPoint(iv-1,iv,HData->GetRMS());
      RMSPerVtx->SetPointError(iv-1,0.5,HData->GetRMSError());
      //   cout<<iv<<"   "<<obs<<"  "<<Ndata<<"  "<<HData->GetRMS()<<"   "<<HData->GetRMSError()<<endl;
    }
    
    RMSPerVtx->SetMarkerStyle(mType[im]);
    RMSPerVtx->SetLineStyle(im);
    RMSPerVtx->SetMarkerColor(colors[im]);
    RMSPerVtx->SetMarkerSize(1.5);
    RMSPerVtx->SetLineColor(colors[im]);
    
    if(im==0)
      RMSPerVtx->Draw("P");
    else
      RMSPerVtx->Draw("P");
    
    legend->AddEntry(RMSPerVtx, (mets[im]+"MET").c_str(), "lp");

  }
  legend->SetShadowColor(0);
  legend->Draw();

 cmsPrel(Lumi);

}




void
ZUseTree::PlotGlobalResponse( bool isMC) {

  
  string mets[6]={"calo","caloT1","caloT2","tc","pf","pfT1"};
  int colors[6]={kRed+1,kOrange+7,kOrange, 896,38, kBlue+1};
  int mType[6] = {33,30,29,20,21,22};

  c2=new TCanvas("c2","Test",300,300,600,600);

 
  TLegend* legend=new TLegend(0.6,0.6,0.8,0.8);
  legend->SetFillColor(0);
  legend->SetLineColor(0);
  legend->SetTextSize(0.04);


  TH1F* emptyHisto=new TH1F("emp","",1000,0,1000);
  emptyHisto->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
  emptyHisto->GetYaxis()->SetTitle( "<u_{||}>/q_{T}" );
  emptyHisto->GetYaxis()->SetTitleOffset(1.1);
  emptyHisto->GetXaxis()->SetTitle("q_{T} [GeV]");
  emptyHisto->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
  emptyHisto->GetYaxis()->SetNdivisions(Ydiv[0],Ydiv[1],Ydiv[2]);
  emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);
      
  emptyHisto->Draw();

  for(int im=0;im<6;im++) {

    string obs = mets[im]+"Response";
    
    int Obs =  histoManager.FindNVarP( obs );
    int Ndata= histoManager.accessProf(Obs,nt);
    if(isMC)
      Ndata =  histoManager.accessProf(Obs,nt-1);

    TProfile* dataProf = (TProfile*)Profiles[Ndata]->Clone();
    
    dataProf->SetLineWidth(LineWidth);
    dataProf->SetMarkerSize(MarkerSize);
    dataProf->SetMarkerStyle(mType[im]);
    dataProf->SetMarkerColor(colors[im]);
    dataProf->SetLineColor(colors[im]);
    dataProf->SetLineStyle(im);
    dataProf->Draw("same");
    
    legend->AddEntry(dataProf, (mets[im]+"MET").c_str(), "lp");


    for(int ib=0;ib<dataProf->GetNbinsX();ib++) {
      if(dataProf->GetBinContent(ib+1) != 0)
	cout<<mets[im]<<"   "<<	dataProf->GetBinLowEdge(ib+1)<<"   "
	    << dataProf->GetBinLowEdge(ib+1) + dataProf->GetBinWidth(ib+1)<<"    "<< 1./dataProf->GetBinContent(ib+1)<<endl;

    }
    
  }
  legend->SetShadowColor(0);
  legend->Draw();

  TLine* line= new TLine(RangeX[0],1,RangeX[1],1);
  line->SetLineStyle(2);
  line->SetLineColor(kGreen+1);
  line->Draw("same");

  cmsPrel(Lumi);

}


void ZUseTree::PlotResolutionUsingRooFit(int Nvtx, string type,string proj,bool isMC) {


  TTree* ZMC=new TTree("ZMC","MC");
  TTree* ZData=new TTree("ZData","data");

  float uparaMC,uperpMC,uparaD, uperpD, qTMC, qTD;

  ZMC->Branch("Upara",&uparaMC,"uparaMC/F");
  ZMC->Branch("Uperp",&uperpMC,"uperpMC/F");
  ZMC->Branch("qT",&qTMC,"qTMC/F");

  ZData->Branch("Upara",&uparaD,"uparaD/F");
  ZData->Branch("Uperp",&uperpD,"uperpD/F");
  ZData->Branch("qT",&qTD,"qTD/F");

  float RecoilProjMC[6][2],  RecoilProjData[6][2];
  int NVertexMC, NVertexData;
  float ZPtMC, ZPtData;

  int Nmet=0;
  if(type=="tc")
    Nmet=1;
  if(type=="calo")
    Nmet=2;
  if(type=="caloT1")
    Nmet=3;
  if(type=="caloT2")
    Nmet=4;
  if(type=="pfT1")
    Nmet=5;

  tChains[nt-1]->SetBranchAddress("BasicRecoilProj",RecoilProjMC);
  tChains[nt-1]->SetBranchAddress("NVertex",&NVertexMC);
  tChains[nt-1]->SetBranchAddress("ZPt",&ZPtMC);

  tChains[nt]->SetBranchAddress("BasicRecoilProj",RecoilProjData);
  tChains[nt]->SetBranchAddress("NVertex",&NVertexData);
  tChains[nt]->SetBranchAddress("ZPt",&ZPtData);

  cout<<" Begin MC tree"<<endl;

  for(int imc=0;imc< tChains[nt-1]->GetEntries();imc++) {
    
    if(NVertexMC!=Nvtx) continue;

    uparaMC = RecoilProjMC[Nmet][0];
    uperpMC = RecoilProjMC[Nmet][1];
    qTMC = ZPtMC;
    
    ZMC->Fill();

  }
  cout<<" End MC tree, begin data tree" <<endl;
  for(int id=0;id< tChains[nt]->GetEntries();id++) {
   
    if(NVertexData!=Nvtx) continue;

    uparaD = RecoilProjData[Nmet][0];
    uperpD = RecoilProjData[Nmet][1];
    qTD = ZPtData;
    
    ZData->Fill();
  }
  cout<<" End data tree "<<endl;

  double BinVar[21] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
  TCanvas* c34= new TCanvas("RooFitCanvas","RooFitCanvas",1200,1200);
  c34->Divide(5,4);
  TCanvas* c35= new TCanvas("RooFitCanvasData","RooFitCanvasData",1200,1200);
  c35->Divide(5,4);

  TH1D* ParaResoMC=new TH1D("pararesoMC","parareso",20,0,200);
  TH1D* ParaResoD=new TH1D("pararesoD","parareso",20,0,200);

  for(int iqt=0;iqt<20;iqt++) {

    RooRealVar qTcat("qT","qT",BinVar[iqt],BinVar[iqt+1]);
    
    RooRealVar upara("upara","upara",-400,100);
      
    RooDataSet MCpara("MCdata","MCData",ZMC, RooArgSet(qTcat,upara) );
    RooDataSet datapara("MCdata","MCData",ZData, RooArgSet(qTcat,upara) );

    RooRealVar sig_mc("sigmc","sigmc",10,0,70);
    RooRealVar mean_mc("sigmc","sigmc",iqt*(-10),-500,100);
    RooGaussian gausMC("GausMC","GausMC",upara,mean_mc,sig_mc);

    RooRealVar sig_dat("sigdat","sigdat",10,0,70);
    RooRealVar mean_dat("sigdat","sigdat",iqt*(-10),-500,100);
    RooGaussian gausDAT("GausDAT","GausDAT",upara,mean_dat,sig_dat);
      
    c34->cd(iqt+1);
    RooFitResult* result = gausMC.fitTo(MCpara,RooFit::SumW2Error(kFALSE),RooFit::PrintLevel(-1) );
    RooPlot* frameNN = upara.frame() ;
    MCpara.plotOn(frameNN) ;
    gausMC.plotOn(frameNN) ;
    gausMC.paramOn(frameNN);
    frameNN->Draw();

    c35->cd(iqt+1);
    RooFitResult* result2 = gausDAT.fitTo(datapara,RooFit::SumW2Error(kFALSE),RooFit::PrintLevel(-1) );
    RooPlot* frameNN2 = upara.frame() ;
    datapara.plotOn(frameNN2) ;
    gausDAT.plotOn(frameNN2) ;
    gausDAT.paramOn(frameNN2);
    frameNN2->Draw();

    ParaResoMC->SetBinContent(iqt+1,sig_mc.getVal() );
    ParaResoMC->SetBinError(iqt+1,sig_mc.getError() );
    ParaResoD->SetBinContent(iqt+1,sig_dat.getVal() );
    ParaResoD->SetBinError(iqt+1,sig_dat.getError() );
  }

  TCanvas* c23=new TCanvas("c23","resoviaroofit",600,600); 
  ParaResoMC->Draw();
  ParaResoD->Draw("same");

}





void
ZUseTree::SaveResolutionPerVtx(string met, string type,bool isMC, bool isCor=false) {


  string obs= met+"U"+type+"VsQt_V";

  if(isCor)
    obs= met+"U"+type+"VsQtCor_V";
  
  int Nbin=13;
  double BinVar[14] = {0,5,10,20,30,40,50,60,80,100,120,140,160,200};

  int mType[4] = {29,20,21,22};
  int mcolor[4]={kBlue+1,kGreen+1,kOrange-2,kRed+1};

  string name = obs;
  if(isMC)
    name += "_MC";

  TFile out( (name+".root").c_str(),"RECREATE");

  for(int iv=1;iv<5;iv++) {

    ostringstream os;
    os<< iv;

    int Obs =  histoManager.FindNVar2D(obs+ os.str() );
    int Ndata= histoManager.access2D(Obs,nt);
    if(isMC)
      Ndata =  histoManager.access2D(Obs,nt-1);
    TH2F* HData = (TH2F*)Histos2D[ Ndata ]->Clone();
    TH1D* RMSData = histoManager.ConvertHistoToRMSPlotVarBin(HData,Nbin,BinVar,obs+os.str(),Remove2V,false);

    RMSData->SetLineWidth(LineWidth);
    RMSData->SetMarkerSize(MarkerSize);
    RMSData->SetMarkerStyle(mType[iv-1]);
    RMSData->SetMarkerColor(mcolor[iv-1]);
    RMSData->SetLineColor(mcolor[iv-1]);

    RMSData->Write();
  }

  out.Close();

}




void ZUseTree::GetRecoilCorrections(string comp, string met) {

  c2=new TCanvas("c2","RecoilCor",300,300,600,600);

 TH1F* emptyHisto=new TH1F("emp","",1000,0,1000);
  emptyHisto->GetYaxis()->SetRangeUser(RangeY[0],RangeY[1]);
  emptyHisto->GetYaxis()->SetTitle( "<u_{||}>/q_{T}" );
  emptyHisto->GetYaxis()->SetTitleOffset(1.1);
  emptyHisto->GetXaxis()->SetTitle("q_{T} [GeV]");
  emptyHisto->GetXaxis()->SetNdivisions(Xdiv[0],Xdiv[1],Xdiv[2]);
  emptyHisto->GetYaxis()->SetNdivisions(Ydiv[0],Ydiv[1],Ydiv[2]);
  emptyHisto->GetXaxis()->SetRangeUser(RangeX[0],RangeX[1]);
  emptyHisto->Draw();

  string obs = met+"Cor"+comp;
    
  int Obs =  histoManager.FindNVarP( obs );
  int Ndata= histoManager.accessProf(Obs,nt);
  int NMC= histoManager.accessProf(Obs,nt-1);

  TProfile* dataProf = (TProfile*)Profiles[Ndata]->Clone();
  TProfile* mcProf = (TProfile*)Profiles[NMC]->Clone();

  vector<TProfile*> cloneP; //Noise profiles
  for(int i=0;i<nt-1;i++) {
    int Nhisto = histoManager.accessProf(Obs,nt-1-i);
    cloneP.push_back((TProfile*)Profiles[Nhisto]->Clone() );
  }

  dataProf->Draw("same");
  
  mcProf->SetLineColor(kBlue);
  mcProf->SetMarkerColor(kBlue);
  mcProf->Draw("same");

  for(int ib=0;ib<dataProf->GetNbinsX();ib++) {

    cout<<0<<"   "<<dataProf->GetBinLowEdge(ib+1)<<"   "
	<<dataProf->GetBinLowEdge(ib+1) + dataProf->GetBinWidth(ib+1)<<"   "
      //	<<dataProf->GetBinContent(ib+1)<<"   "
      //	<<mcProf->GetBinContent(ib+1)<<"   "
      	<<dataProf->GetBinContent(ib+1)/mcProf->GetBinContent(ib+1)<<endl;
  }



}





void
ZUseTree::PrintResponseSummary() {

  
  string mets[6]={"calo","caloT1","caloT2","tc","pf","pfT1"};
 
  int tmp = histoManager.FindNVarP("caloResponse");
  int tmpn= histoManager.accessProf(tmp,nt);
  TProfile* tmpp = (TProfile*)Profiles[tmpn]->Clone();
  

  for(int ib=0;ib<tmpp->GetNbinsX();ib++) {

    for(int im=0;im<6;im++) {

      string obs = mets[im]+"Response";
    
      int Obs =  histoManager.FindNVarP( obs );
      int Ndata= histoManager.accessProf(Obs,nt);
   
      int NMC =  histoManager.accessProf(Obs,nt-1);

      TProfile* dataProf = (TProfile*)Profiles[Ndata]->Clone();
      TProfile* MCProf = (TProfile*)Profiles[NMC]->Clone();
    
   
      if(dataProf->GetBinContent(ib+1) != 0) {
	
	if(im==0)
	  cout<<fixed<< setprecision(0)<<" $[ "<<	dataProf->GetBinLowEdge(ib+1)<<" - "
	      << dataProf->GetBinLowEdge(ib+1) + dataProf->GetBinWidth(ib+1)<<" ]$ & ";

	cout<<fixed<< setprecision(2)<< 1./dataProf->GetBinContent(ib+1)<<" $\\pm$ "<< setprecision(2)<<1/(dataProf->GetBinContent(ib+1)) *dataProf->GetBinError(ib+1)/dataProf->GetBinContent(ib+1)<<" & ";
	float error = MCProf->GetBinContent(ib+1)/dataProf->GetBinContent(ib+1) * sqrt( pow(1/(dataProf->GetBinContent(ib+1)) *dataProf->GetBinError(ib+1)/dataProf->GetBinContent(ib+1),2) + pow(1/(MCProf->GetBinContent(ib+1)) *MCProf->GetBinError(ib+1)/MCProf->GetBinContent(ib+1),2 ) );
	cout<<fixed<< setprecision(2)<< MCProf->GetBinContent(ib+1)/dataProf->GetBinContent(ib+1)<<" $\\pm$ "<< setprecision(2)<<error;

	if(im!=5)
	  cout<<" & ";
	else
	  cout<<" \\\\ "<<endl;

      }
      delete dataProf;
      delete MCProf;
    }
  }

  delete tmpp;

}

float ZUseTree::Rapidity(float p, float mass, float pt, float eta) {

  //compute pz
  float pz = sqrt(p*p-pt*pt);

  //compute E
  float E = sqrt(p*p + mass*mass);

  //and now rapidity
  float y = 0.5 * log( (E+pz) / (E-pz) );

  return y;

}
