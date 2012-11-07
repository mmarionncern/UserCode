#include "WUseTree.hh"

#include "TVector2.h"

using namespace std;


ClassImp(WUseTree)


WUseTree::WUseTree():
UseTree()
{

  LoadDBWeight();
  LoadDBResponse();
  LoadRecoilCorrections();

 float WIsoCuts[2][4][3]= { { {0.05, 0.06, 0.03}, {0.09, 0.07, 0.1},
			       {0.12, 0.09, 0.1}, {0.15, 2.0, 0.12} },
			     { {0.025, 0.025, 0.02}, {0.04, 0.05, 0.025},
			       {0.05, 0.06, 0.03}, {0.08, 0.06, 0.05} } };
  
  float WIDCuts[2][4][4]= { { {0.01,0.004,0.03,0.025}, {0.01,0.004,0.06,0.04},
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
	  _IsoCuts[j][i][k] = WIsoCuts[j][i][k];
	_IDCuts[j][i][k] = WIDCuts[j][i][k];
      }
    }
  }

  getRatio =false;

  Sorder.push_back("W #rightarrow e#nu");
  Sorder.push_back("EWK");
  Sorder.push_back("QCD");
  Sorder.push_back("#gamma+jet");

  //for(int i=0;i<80;i++)
  //  VtxCorSyst.push_back(0);

}

void
WUseTree::FillWTree() {

  //Add Varaibles to be computed

  //Prepare suffix for outfiles
  string suffix = "METs";
  if(switchRMS && FillProf)
    suffix = "RMS_METs";

  if(Draw3on1)
    for(int i=0;i<5;i++)
      Suffix[i] = suffix;

  string metT[6]={"pf","tc","calo","caloT1","caloT2","pfT1"};
  string metName[6]={"PF","TC","Calo","Calo Type I","Calo Type II","PF Type I"};
  //Vertex
  histoManager.AddVariable("Vertex",10,0,10,"number of vertices","Vertex");

  //Eta Phi Pt
  histoManager.AddVariable("Eta",100,-3,3,"#eta","Eta");
  histoManager.AddVariable("Phi",128,0,3.2,"#phi","Phi");
  histoManager.AddVariable("Pt",200,0,100,"p_{T}","Pt");
  histoManager.AddVariable("EtSC",200,0,100,"p_{T}","Pt");

  histoManager.AddVariable("PtControl",200,0,100,"p_{T}","Pt");
  histoManager.AddVariable("PtUnc",200,0,100,"p_{T}","Pt");

  histoManager.AddVariable("ExpInHit",2,0,2,"expInHit","expinHit");

  //ID
  histoManager.AddVariable("sigieie",200,0,0.1,"#sigma_{i#etai#eta}","sigieie");
  histoManager.AddVariable("Deta",100,-0.05,0.05,"#Delta#eta","Deta");
  histoManager.AddVariable("Dphi",100,-0.2,0.2,"#Delta#Phi","Dphi");
  histoManager.AddVariable("HoE",100,0,0.1,"H/E","HoE");
  
 //iso
  histoManager.AddVariable("TrackIso",500,0,1,"TrackIso","TrackIso");
  histoManager.AddVariable("EcalIso",500,0,1,"EcalIso","EcalIso");
  histoManager.AddVariable("HcalIso",500,0,1,"HcalIso","HcalIso");
    
  //Underlying event variables
  histoManager.AddVariable("dZv",400,-20,20,"#Delta z [cm]","We_dZ_vertices");
  histoManager.AddVariable("zv2",400,-100,100,"# z vertex2 [cm]","We_Zv2_vertices");
  histoManager.AddVariable("genWPt",400,0,200,"p_{T} W generated [GeV]","WGen_Pt");
  histoManager.AddVariable("jetPt",400,0,200,"p_{T} jet [GeV]","We_Pt_Jet");
  histoManager.AddVariable("photonPt",400,0,200,"p_{T} photon [GeV]","We_Pt_Photon");
  histoManager.AddVariable("photonMult",10,0,10,"N photons","We_Photon_Mult");
  histoManager.AddVariable("jetMult",10,0,10,"N jets","We_Jet_Mult");

  for(int kk=0;kk<6;kk++) {

  //METs
  histoManager.AddVariable(metT[kk]+"MET",400,0,200,metName[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
  //METProj
  histoManager.AddVariable(metT[kk]+"METPara",400,-200,200,metName[kk]+" #slash{E}_{T,||} [GeV]","We_MET_Para_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"METPerp",400,-200,200,metName[kk]+" #slash{E}_{T,#perp}   [GeV]","We_MET_Perp_"+Suffix[kk]);
  
  //Acoplanarity pfAcop
  histoManager.AddVariable(metT[kk]+"Acop",100,0,3.14,"#zeta ("+metName[kk]+") [rad]","We_Acop_"+Suffix[kk]);
  
  //Recoil
  histoManager.AddVariable(metT[kk]+"Recoil",400,0,200,metName[kk]+" u_{T} [GeV]","We_Recoil_"+Suffix[kk]);
  // histoManager.AddVariable(metT[kk]+"RecoilCor",400,0,200,metName[kk]+" u_{T} [GeV]","We_Recoil_Cor_"+Suffix[kk]);
  
  //RecoilProj
  histoManager.AddVariable(metT[kk]+"RecoilPara",400,-200,200,metName[kk]+" u_{||} [GeV]","We_Recoil_Para_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"RecoilPerp",400,-200,200,metName[kk]+" u_{#perp}  [GeV]","We_Recoil_Perp_"+Suffix[kk]);
 histoManager.AddVariable(metT[kk]+"RecoilParaCor",400,-200,200,metName[kk]+" u_{||} [GeV]","We_Recoil_Para_Cor_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"RecoilPerpCor",400,-200,200,metName[kk]+" u_{#perp}  [GeV]","We_Recoil_Perp_Cor_"+Suffix[kk]);
  
  histoManager.Add2DVariable(metT[kk]+"RecoilParaVtx","number of vertices ","RMS( u_{||} ) [GeV]","We_Recoil_Para_Vtx_"+Suffix[kk],10,0,10,200,-100,100);
  histoManager.Add2DVariable(metT[kk]+"RecoilPerpVtx","number of vertices ","RMS( u_{#perp}  ) [GeV]","We_Recoil_Perp_Vtx_"+Suffix[kk],10,0,10,200,-100,100);
  histoManager.Add2DVariable(metT[kk]+"RecoilCorParaVtx","number of vertices ","RMS( u_{||} ) [GeV]","We_Recoil_Para_Vtx_"+Suffix[kk],10,0,10,200,-100,100);
  histoManager.Add2DVariable(metT[kk]+"RecoilCorPerpVtx","number of vertices ","RMS( u_{#perp}  ) [GeV]","We_Recoil_Perp_Vtx_"+Suffix[kk],10,0,10,200,-100,100);

  //SumET
  histoManager.AddVariable(metT[kk]+"SumET",500,0,500,"#Sigma E_{T} ("+metT[kk]+") [GeV]","We_SumET_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"SumETRecoil",500,0,500,"#Sigma E_{T} ("+metT[kk]+") [GeV]","We_SumETRecoil_"+Suffix[kk]);
  histoManager.Add2DVariable(metT[kk]+"SumETRecoil2D","u_{T} ("+metT[kk]+") [GeV]","#Sigma E_{T} ("+metT[kk]+") [GeV]","We_SumETRecoil2D_"+Suffix[0],11,0,110,25,0,500);

  //MT
  histoManager.AddVariable(metT[kk]+"MT",400,0,200,"M_{T} ("+metT[kk]+") [GeV]","We_MT_"+Suffix[kk]);
 
  //MT vs Acop
  histoManager.Add2DVariable(metT[kk]+"AcopMT","#zeta ("+metT[kk]+") [rad]","M_T [GeV]","We_MT_vs_Acop_"+Suffix[kk],5,0,3.2,40,0,200);
  
  //Angle
  histoManager.AddVariable(metT[kk]+"dPhiRecoil",42,0,3.15,metName[kk]+" #Delta#phi(e,u_{T}) [rad]","We_dPhi_Recoil_"+Suffix[kk]);
  
  //Underlying variables and quality variables
  histoManager.AddVariable(metT[kk]+"WPt",400,0,200,"p_{T} W ("+metT[kk]+") [GeV]","We_PT_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"dPhiJetW",100,0,3.14,"#Delta#phi(jet,W ("+metT[kk]+") ) [rad]","We_dPhi_Jet_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"dPhiPhotonW",100,0,3.14,"#Delta#phi(jet,W ("+metT[kk]+") ) [rad]","We_dPhi_Photon_"+Suffix[kk]);

  //Corrected variables ======================
    histoManager.AddVariable(metT[kk]+"METCor",400,0,200,metName[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METParaCor",400,-200,200,metName[kk]+" #slash{E}_{T,||} [GeV]","We_MET_Para_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METPerpCor",400,-200,200,metName[kk]+" #slash{E}_{T,#perp}   [GeV]","We_MET_Perp_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilCor",400,0,200,metName[kk]+" u_{T} [GeV]","We_Recoil_Cor_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"MTCor",400,0,200,"M_{T} ("+metT[kk]+") [GeV]","We_MT_"+Suffix[kk]);

    histoManager.AddVariable(metT[kk]+"METCorControl",400,0,200,metName[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METParaCorControl",400,-200,200,metName[kk]+" #slash{E}_{T,||} [GeV]","We_MET_Para_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"METPerpCorControl",400,-200,200,metName[kk]+" #slash{E}_{T,#perp}   [GeV]","We_MET_Perp_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"RecoilCorControl",400,0,200,metName[kk]+" u_{T} [GeV]","We_Recoil_Cor_"+Suffix[kk]);
    histoManager.AddVariable(metT[kk]+"MTCorControl",400,0,200,"M_{T} ("+metT[kk]+") [GeV]","We_MT_"+Suffix[kk]);
  //==========================================


  float tB[6];
  if(FillProf)
    { tB[0]=100; tB[1]=0; tB[2]=200; }
  else
    {  tB[0]=40; tB[1]=0; tB[2]=200; tB[3]=50; tB[4]= 0; tB[5]=3.14; }
  histoManager.Add2DVariable(metT[kk]+"dPhiRecoil2D","u_{T} ("+metT[kk]+") [GeV]","#Delta#phi ("+metT[kk]+") [rad] ","We_dPhi_vs_Recoil_"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);

  //Ut/Pt
  histoManager.AddVariable(metT[kk]+"RecoilPt",200,0,2,"u_{T}/pt_{lepton} ("+metT[kk]+")","We_Recoil_pt_"+Suffix[kk]);
    
  histoManager.Add2DVariable(metT[kk]+"dPhiRecoilPt2D","u/q_{T} (pf) [GeV]","#Delta#phi ("+metT[kk]+") [rad] ","We_dPhi_vs_RoQ_"+Suffix[kk],tB[0],tB[1],tB[2],tB[3],tB[4],tB[5]);
 
  //Special Variables for WGen study ********
  histoManager.Add2DVariable(metT[kk]+"RecoilvsGen","Gen q_{T} [GeV]","u_{||} [GeV]","",21,-42,42,21,-42,42);
   histoManager.Add2DVariable(metT[kk]+"SpecRecoilvsGen","Gen q_{T} [GeV]","u_{||} [GeV]","",21,-42,42,21,-42,42);
   histoManager.Add2DVariable(metT[kk]+"AbsRecoilvsGen","Gen q_{T} [GeV]","<|u_{||}|> [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable(metT[kk]+"AbsSpecRecoilvsGen","Gen q_{T} [GeV]","<|u_{||}|> [GeV]","",21,-42,42,21,-42,42);
  
  histoManager.Add2DVariable(metT[kk]+"RecoilEta","#eta","u_{T} [GeV]","",25,-3,3,22,0,110);

  //  Histomanager.Add2DVariable(met[kk]+"GenVsAnglep",,);
  //  Histomanager.Add2DVariable(met[kk]+"GenVsAnglep",,);

  //***************************************** Control samples
  
  //METs
  histoManager.AddVariable(metT[kk]+"MET_V1",400,0,200,metT[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"MET_V2",400,0,200,metT[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"Recoil_V1",400,0,200,metT[kk]+" u_{T} [GeV]","We_Recoil_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"Recoil_V2",400,0,200,metT[kk]+" u_{T} [GeV]","We_Recoil_"+Suffix[kk]);
  
  histoManager.AddVariable(metT[kk]+"SumETRecoil_V1",500,0,500,"#sum E_{T} ("+metT[kk]+") [GeV]","We_SumET_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"SumETRecoil_V2",500,0,500,"#sum E_{T} ("+metT[kk]+") [GeV]","We_SumET_"+Suffix[kk]);

  histoManager.AddVariable(metT[kk]+"METControl",400,0,200,metT[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"MET_V1Control",400,0,200,metT[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"METParaControl",400,-200,200,metT[kk]+" #slash{E}_{T,||} [GeV]","We_MET_Para_"+Suffix[kk]);  
  histoManager.AddVariable(metT[kk]+"METPerpControl",400,-200,200,metT[kk]+" #slash{E}_{T,#perp}   [GeV]","We_MET_Perp_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"AcopControl",100,0,3.14,"#zeta ("+metT[kk]+") [rad]","We_Acop_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"dPhiRecoilControl",100,0,3.14,"#Delta#phi ("+metT[kk]+") [rad]","We_dPhi_Recoil_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"RecoilControl",400,0,200,"u_{T} ("+metT[kk]+") [GeV]","We_Recoil_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"RecoilParaControl",400,-200,200,"u_{||} ("+metT[kk]+") [GeV]","We_Recoil_Para_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"RecoilPerpControl",400,-200,200,"u_{#perp}   ("+metT[kk]+") [GeV]","We_Recoil_Perp_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"MTControl",400,0,200,"M_{T} ("+metT[kk]+") [GeV]","We_MT_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"SumETRecoil_V1Control",500,0,500,"#sum E_{T} ("+metT[kk]+") [GeV]","We_SumET_"+Suffix[kk]);

  //Uncertainties ********************************************

  histoManager.AddVariable(metT[kk]+"METUnc",400,0,200,metT[kk]+" #slash{E}_{T} [GeV]","We_MET_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"METParaUnc",400,-200,200,metT[kk]+" #slash{E}_{T,||} [GeV]","We_MET_Para_"+Suffix[kk]);  
  histoManager.AddVariable(metT[kk]+"METPerpUnc",400,-200,200,metT[kk]+" #slash{E}_{T,#perp}   [GeV]","We_MET_Perp_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"AcopUnc",100,0,3.14,"#zeta ("+metT[kk]+") [rad]","We_Acop_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"dPhiRecoilUnc",100,0,3.14,"#Delta#phi ("+metT[kk]+") [rad]","We_dPhi_Recoil_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"RecoilUnc",400,0,200,"u_{T} ("+metT[kk]+") [GeV]","We_Recoil_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"RecoilParaUnc",400,-200,200,"u_{||} ("+metT[kk]+") [GeV]","We_Recoil_Para_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"RecoilPerpUnc",400,-200,200,"u_{#perp}   ("+metT[kk]+") [GeV]","We_Recoil_Perp_"+Suffix[kk]);
  histoManager.AddVariable(metT[kk]+"MTUnc",400,0,200,"M_{T} ("+metT[kk]+") [GeV]","We_MT_"+Suffix[kk]);


  }

  //Special variables ***********************************
  
  histoManager.AddVariable("GenParaP",400,-200,200,"u_{||} [GeV]","We_GenPara_plus");
  histoManager.AddVariable("GenParaM",400,-200,200,"u_{||} [GeV]","We_GenPara_minus");
  histoManager.AddVariable("GenPerpP",400,-200,200,"u_{#perp}  [GeV]","We_GenPerp_plus");
  histoManager.AddVariable("GenPerpM",400,-200,200,"u_{#perp}  [GeV]","We_GenPerp_minus");

  histoManager.Add2DVariable("GenAngleP","Gen q_{T} [GeV]","dphi plus","",40,0,200,90,-3.14,3.14);
  histoManager.Add2DVariable("GenAngleM","Gen q_{T} [GeV]","dphi minus","",40,0,200,90,-3.14,3.14);


  histoManager.Add2DVariable("pfRecoilvsGenR","Gen q_{T} [GeV]","u_{||} [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("tcRecoilvsGenR","Gen q_{T} [GeV]","u_{||} [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("pfSpecRecoilvsGenR","Gen q_{T} [GeV]","u_{||} [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("tcSpecRecoilvsGenR","Gen q_{T} [GeV]","u_{||} [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("pfAbsRecoilvsGenR","Gen q_{T} [GeV]","<|u_{||}|> [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("tcAbsRecoilvsGenR","Gen q_{T} [GeV]","<|u_{||}|> [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("pfAbsSpecRecoilvsGenR","Gen q_{T} [GeV]","<|u_{||}|> [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("tcAbsSpecRecoilvsGenR","Gen q_{T} [GeV]","<|u_{||}|> [GeV]","",21,-42,42,21,-42,42);

  histoManager.Add2DVariable("pfCTRecoilvsGenR","Gen q_{T} [GeV]","u_{||} [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("tcCTRecoilvsGenR","Gen q_{T} [GeV]","u_{||} [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("pfCTRecoilvsGen","Gen q_{T} [GeV]","u_{||} [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("tcCTRecoilvsGen","Gen q_{T} [GeV]","u_{||} [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("pfAbsCTRecoilvsGenR","Gen q_{T} [GeV]","<|u_{||}|> [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("tcAbsCTRecoilvsGenR","Gen q_{T} [GeV]","<|u_{||}|> [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("pfAbsCTRecoilvsGen","Gen q_{T} [GeV]","<|u_{||}|> [GeV]","",21,-42,42,21,-42,42);
  histoManager.Add2DVariable("tcAbsCTRecoilvsGen","Gen q_{T} [GeV]","<|u_{||}|> [GeV]","",21,-42,42,21,-42,42);

  //************************************************************

 //End Declaration, now preparation
  histoManager.PrepareHistos(name);
  histoManager.Prepare2DHistos(name);
 

  vector<vector<vector<float> > > Wghts;
  if( MultiWeight || Remove2V || ShapeWeight)
     Wghts = GetFitWeight();
  
  string Epart="EE";
  if(EcalP==0)
    Epart = "EB";
  if(EcalP==2)
    Epart="Combined";

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
    float IsoVar[3];
    float IDVar[6];
    float Mets[6][2];
    float MetProj[6][2];
    float Recoils[6][2];
    float RecoilProj[6][2];
    float MT[6];
    float SumET[6];
    float Lepton[3];
    float EtSC;

    float dPhiRecoilLep[6];
    float dPhiMETLep[6];

    bool VetoE;
    bool ConvRej;
    bool InHit;
    bool ExpInHit;

    int Run;
    int Event;

    float GenCharge;
    float GenRecoil[2][2];
    float SpecRecoilProj[6][2];
    float CTRecoil[6][2];
    
    float GenWPt[2];
    float GenEPt[2];
    float WPt[6];
    float DZvertices;
    float zvertex2;
    //Underlying event
    float dPhiPhotonW[6];
    float dPhiJetW[6];
    int JetMultiplicity;
    int PhotonMultiplicity;
    float LeadingJet;
    float LeadingPhoton;
    
    int NVertex;

    float Weight=1;

    tChains[i]->SetBranchAddress("NVertex",&NVertex);

    tChains[i]->SetBranchAddress("IsoVar",IsoVar);
    tChains[i]->SetBranchAddress("IDVar",IDVar);
    tChains[i]->SetBranchAddress("Lepton",Lepton);
    tChains[i]->SetBranchAddress("EtSC",&EtSC);
 
    tChains[i]->SetBranchAddress("Mets",Mets);
    tChains[i]->SetBranchAddress("MetProj",MetProj);
 
    tChains[i]->SetBranchAddress("BasicRecoils",Recoils);
    tChains[i]->SetBranchAddress("BasicRecoilProj",RecoilProj);
    tChains[i]->SetBranchAddress("RecoilProj",SpecRecoilProj);

    tChains[i]->SetBranchAddress("MT",MT);

    tChains[i]->SetBranchAddress("SumET",SumET);

    tChains[i]->SetBranchAddress("BasicdPhiRecoilLep",dPhiRecoilLep);
    tChains[i]->SetBranchAddress("dPhiMETLep",dPhiMETLep);
  
    tChains[i]->SetBranchAddress("VetoE",&VetoE);
    tChains[i]->SetBranchAddress("InHit",&InHit);
    tChains[i]->SetBranchAddress("ExpInHit",&ExpInHit);
    tChains[i]->SetBranchAddress("ConvRej",&ConvRej);

    tChains[i]->SetBranchAddress("Run",&Run);
    tChains[i]->SetBranchAddress("AbsEvent",&Event);

    tChains[i]->SetBranchAddress("WPt",WPt);
    tChains[i]->SetBranchAddress("DZvertices",&DZvertices);
    tChains[i]->SetBranchAddress("zv2",&zvertex2);

    //Underlying event
    tChains[i]->SetBranchAddress("LeadingJet",&LeadingJet);
    tChains[i]->SetBranchAddress("LeadingPhoton",&LeadingPhoton);
    tChains[i]->SetBranchAddress("PhotonMultiplicity",&PhotonMultiplicity);
    tChains[i]->SetBranchAddress("JetMultiplicity",&JetMultiplicity);
    tChains[i]->SetBranchAddress("dPhiJetW",dPhiJetW);
    tChains[i]->SetBranchAddress("dPhiPhotonW",dPhiPhotonW);
	
    if(i==nt-1) {
      tChains[i]->SetBranchAddress("GenCharge",&GenCharge);
      tChains[i]->SetBranchAddress("GenRecoil",GenRecoil);
      tChains[i]->SetBranchAddress("CTRecoil",CTRecoil);
      tChains[i]->SetBranchAddress("GenWPt",GenWPt);
      tChains[i]->SetBranchAddress("GenEPt",GenEPt);
    }
    else {
      for(int jj=0;jj<5;jj++) {
	CTRecoil[jj][0]=0; CTRecoil[jj][1]=0;
      }
      GenCharge=-100;
      GenRecoil[0][0]=-1000; GenRecoil[0][1]=-1000;
      GenRecoil[1][0]=-1000; GenRecoil[1][1]=-1000;
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

      /* if(i<nt-2) {
	if(NVert==1 )
	  Weight *= Frac1V;
	if(NVert==2 )
	Weight *= (1-Frac1V);
	}*/
 /*  if(NVert==1 && i!=nt && i!=nt-1 )
     Weight *= Frac1V;*/

      //Reweight PT
      if( name[i].substr(0,4)=="W #r")
        Weight *=  SearchWeight(GenWPt[0]);

      tChains[i]->GetEntry(ie);      

      //Throw multijet events //FIXME
      //if(JetMultiplicity!=1) continue;

      if(name[i].substr(0,4)=="data")
	{
	  bool doubleCount=false;
	  for(size_t ik=0;ik<Events.size();ik++)
	    {
	      if( ((Events[ik]).first.first) == Run && ((Events[ik]).first.second) == Event)
		{doubleCount=true; break;}
	    }
	  if(doubleCount || (EventFilter && Event<=EventNum ) )
	    continue;
	}

   
      //Filling Histos
      if(fabs(Lepton[1])>1.479)
	{EP=1; }
      else
	{EP=0; }

      if(name[i]!="data" || EP==0 || noDeltaCor) {
	IDVar[4] = 0;
	IDVar[5] =0;
      }

      //first pt cut
      if(EtSC < PTcut) continue;
  
      if(convRej && !ExpInHit)  continue;
      if(vetoE && !VetoE) continue;
  
      //Met Cut
       
      //And Selections
      //  cout<<IDVar[4]<<endl;
      if(IDVar[0] > _IDcuts[EP][0] && Cbool(skipCut, "sigieie"==Nm1Var) ) continue;
      if( fabs(IDVar[1]/*-IDVar[4]*/) > _IDcuts[EP][1] && Cbool(skipCut, "Deta"==Nm1Var) ) continue;
      if( fabs(IDVar[2]/*-IDVar[5]*/) > _IDcuts[EP][2] && Cbool(skipCut, "Dphi"==Nm1Var) ) continue;
      if(IDVar[3] > _IDcuts[EP][3] && Cbool(skipCut, "HoE"==Nm1Var) ) continue;
   
      if(IsoVar[0]/Lepton[0] > _Isocuts[EP][0] && Cbool(skipCut, "TrackIso"==Nm1Var) ) continue;
      if(IsoVar[1]/Lepton[0] > _Isocuts[EP][1] && Cbool(skipCut, "EcalIso"==Nm1Var) ) continue;
      if(IsoVar[2]/Lepton[0] > _Isocuts[EP][2] && Cbool(skipCut, "HcalIso"==Nm1Var) ) continue;
    
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

    
     
      if(EcalP!=EP && EcalP!=2) continue;
      
      if(NVert==2 && NVertex<2/* && ( i>=nt-2 /* name[i].substr(0,4)=="data"  || name[i].substr(0,4)=="W #r"*/ )
	continue;
      else if(NVert==1 && NVertex>1 /*&& ( i>=nt-2 /*name[i].substr(0,4)=="data"  || name[i].substr(0,4)=="W #r")*/ )
	continue;

      //Vertex
      histoManager.fill("Vertex",i,NVertex,Weight);

      histoManager.fill("Phi",i,Lepton[2],Weight);
     
      histoManager.fill("EtSC",i,EtSC,Weight);

      //ID
      histoManager.fill("sigieie",i, IDVar[0],Weight);
      histoManager.fill("Deta",i, fabs(IDVar[1]-IDVar[4]),Weight);
       histoManager.fill("Dphi",i, IDVar[2]-IDVar[5],Weight);
       histoManager.fill("HoE",i,IDVar[3] ,Weight);
    
      //Iso
       histoManager.fill("TrackIso",i,IsoVar[0]/Lepton[0],Weight);
      histoManager.fill("EcalIso",i,IsoVar[1]/Lepton[0],Weight);
      histoManager.fill("HcalIso",i,IsoVar[2]/Lepton[0],Weight);
    
      for(int kk=0;kk<6;kk++) {

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

	if(ShapeWeight && i!=nt/* && NVert>1*/) {
	  Weight = Weight*GetMCWeightFromDataShape(kk, SumET[kk]-Lepton[0],
						   metT[kk]+"SumETRecoil_V1",
						   metT[kk]+"SumETRecoil_V1", 10. );
	  //  cout<<Weight<<"     ";
	  /*  Weight = Weight*GetMCWeightFromDataShape(kk, Recoils[kk][0],
						   metT[kk]+"Recoil_V1",
						   metT[kk]+"Recoil_V1", 4. );*/
	  //  cout<<Weight<<"    "<<endl;
	}



      //MET
      	histoManager.fill(metT[kk]+"MET",i,Mets[kk][0],Weight);
       //FIXME
	//      continue;
      
      //MET Proj
      histoManager.fill(metT[kk]+"METPara",i,MetProj[kk][0],Weight);
      histoManager.fill(metT[kk]+"METPerp",i,MetProj[kk][1],Weight);
     
      //WPt
      histoManager.fill(metT[kk]+"WPt",i,WPt[kk],Weight);
  

      //Corrected variables ======================
      /*  if(i==nt-1 && kk==0) {
	TVector2 rectmp(0,0);
	rectmp.SetMagPhi(Recoils[kk][0],Recoils[kk][1]);
	TVector2 leptmp(0,0);
	leptmp.SetMagPhi(Lepton[0],Lepton[2]);
	TVector2 GenWpttmp(0,0);
	GenWpttmp.SetMagPhi(GenWPt[0],GenWPt[1]);

	//Unit wPt vector
	TVector2 WPtUnit = GenWpttmp.Unit();
      
	//Corrections
	TVector2 correc(0,0);

	float upara = WPtUnit*rectmp;
	TVector2 tt = WPtUnit.Rotate(acos(-1)/2);
	float uperp = tt*rectmp;

	float paraoff = GetParaRecoilCor(kk,NVertex,GenWPt[0]);
	float perpoff = GetPerpRecoilCor(kk,NVertex,GenWPt[0]);
	float parares = GetParaRecoilCorRes(kk,NVertex,GenWPt[0]);
	float perpres = GetPerpRecoilCorRes(kk,NVertex,GenWPt[0]);

	float uparac = (upara*paraoff)*parares;
	float uperpc = (uperp+perpoff)*perpres;
	float mag = sqrt(uparac*uparac+uperpc*uperpc);

	TVector2 partmp(0,0); 
	partmp.SetMagPhi(uparac, WPtUnit.Phi() + acos(-1) );
	TVector2 pertmp(0,0); 
	pertmp.SetMagPhi(fabs(uperpc), (WPtUnit.Phi() + acos(-1)/2.*uperpc/fabs(uperpc) ) );
	
	correc = partmp + pertmp;

	TVector2 metcor(0,0);
	metcor -= correc + leptmp;

	float MTcor = sqrt(2*Lepton[0]*metcor.Mod()* (1-cos(metcor.DeltaPhi(leptmp))) );
	//+++++++++++

	histoManager.fill(metT[kk]+"METCor",i,metcor.Mod(),Weight);
	histoManager.fill(metT[kk]+"RecoilCor",i, mag,Weight);
	histoManager.fill(metT[kk]+"MTCor",i,MTcor,Weight);
      } 
      else {*/
	histoManager.fill(metT[kk]+"METCor",i,Mets[kk][0],Weight);
	histoManager.fill(metT[kk]+"RecoilCor",i,Recoils[kk][0] ,Weight);
	histoManager.fill(metT[kk]+"MTCor",i, MT[kk],Weight);
	//  }

      //==========================================




      //****** Cutted variables *********

      if( (Mets[kk][0] > METcut && !invCut && MT[kk] > MTCut) || 
	  (Mets[kk][0] < METcut && invCut && MT[kk] > MTCut ) ) { 

	//Underlying event
	if(kk==5) {

	  histoManager.fill("Pt",i,Lepton[0],Weight);
	  histoManager.fill("Eta",i,Lepton[1],Weight);
	  
	  histoManager.fill("ExpInHit",i,ExpInHit,Weight);

	  if(LeadingJet> 15)
	    histoManager.fill("jetPt",i,LeadingJet,Weight);
	  if(LeadingPhoton> 15)
	    histoManager.fill("photonPt",i,LeadingPhoton,Weight);
	  histoManager.fill("photonMult",i,PhotonMultiplicity,Weight);
	  histoManager.fill("jetMult",i,JetMultiplicity,Weight);
	  // if(fabs(DZvertices)>2)
	    histoManager.fill("dZv",i,DZvertices,Weight); 

	  histoManager.fill("zv2",i,zvertex2,Weight); 
	}      
	if(LeadingJet> 15)
	  histoManager.fill(metT[kk]+"dPhiJetW",i,fabs(dPhiJetW[kk]),Weight);
	if(LeadingPhoton> 5)
	  histoManager.fill(metT[kk]+"dPhiPhotonW",i,fabs(dPhiPhotonW[kk]),Weight);

	//SumET
	histoManager.fill(metT[kk]+"SumET",i,SumET[kk],Weight);
	histoManager.fill(metT[kk]+"SumETRecoil",i,SumET[kk]-Lepton[0],Weight);
	histoManager.fill2D(metT[kk]+"SumETRecoil2D",i,Recoils[kk][0],SumET[kk]-Lepton[0],Weight);

	//	if(kk==0 && Recoils[kk][0]<10 && Recoils[kk][0]>=5)
	//	  cout<<" Weight "<<Weight <<"   :  "<<SumET[kk]-Lepton[0]<<endl;

	//Acoplanarity
	histoManager.fill(metT[kk]+"Acop",i,3.14-fabs(dPhiMETLep[kk]*3.14/180),Weight);

	float dbrcor = SearchResponse(WPt[kk],kk);

	  //Recoil
	  histoManager.fill(metT[kk]+"Recoil",i,Recoils[kk][0],Weight);
	  histoManager.fill(metT[kk]+"RecoilCor",i,Recoils[kk][0]*dbrcor,Weight);
	  //Recoil Proj
	  histoManager.fill(metT[kk]+"RecoilPara",i,RecoilProj[kk][0],Weight);
	  histoManager.fill(metT[kk]+"RecoilPerp",i,RecoilProj[kk][1],Weight);
	  histoManager.fill(metT[kk]+"RecoilParaCor",i,RecoilProj[kk][0]*dbrcor,Weight);
	  histoManager.fill(metT[kk]+"RecoilPerpCor",i,RecoilProj[kk][1]*dbrcor,Weight);

	  histoManager.fill2D(metT[kk]+"RecoilParaVtx",i,NVertex,RecoilProj[kk][0],Weight);
	  histoManager.fill2D(metT[kk]+"RecoilPerpVtx",i,NVertex,RecoilProj[kk][1],Weight);
	 
	
	  
	  histoManager.fill2D(metT[kk]+"RecoilCorParaVtx",i,NVertex,RecoilProj[kk][0]*dbrcor,Weight);
	  histoManager.fill2D(metT[kk]+"RecoilCorPerpVtx",i,NVertex,RecoilProj[kk][1]*dbrcor,Weight);
	  histoManager.fill2D(metT[kk]+"RecoilEta",i,Lepton[1],Recoils[kk][0],Weight);
	   
	  //MT
	  histoManager.fill(metT[kk]+"MT",i,MT[kk],Weight);
	  histoManager.fill2D(metT[kk]+"AcopMT",i,3.14-fabs(dPhiMETLep[kk]*3.14/180),MT[kk],Weight);

	  //Angles
	  histoManager.fill(metT[kk]+"dPhiRecoil",i,fabs(dPhiRecoilLep[kk]*3.14/180),Weight);
	  histoManager.fill2D(metT[kk]+"dPhiRecoil2D",i,Recoils[kk][0],fabs(dPhiRecoilLep[kk]*3.14/180),Weight);
	

	  //Uncertainties
	  if( name[i].substr(0,4)=="W #r") {
	    
	    float dbw = SearchWeight(GenWPt[0]);
	   
	    histoManager.fill(metT[kk]+"METUnc",i,Mets[kk][0],Weight*(1-dbw) ); 
	    histoManager.fill(metT[kk]+"METParaUnc",i,MetProj[kk][0],Weight*(1-dbw) );
	    histoManager.fill(metT[kk]+"METPerpUnc",i,MetProj[kk][1],Weight*(1-dbw) );
	     
	    if( (Mets[kk][0] > METcut && !invCut && MT[kk] > MTCut) || 
		(Mets[kk][0] < METcut && invCut && MT[kk] > MTCut ) ) { 
	     
	      if(kk==0)
		histoManager.fill("PtUnc",i,Lepton[0],Weight*(1-dbw) );
	     
	      histoManager.fill(metT[kk]+"AcopUnc",i,3.14-(dPhiMETLep[kk]*3.14/180),Weight*(1-dbw) );
	      histoManager.fill(metT[kk]+"dPhiRecoilUnc",i,(dPhiRecoilLep[kk]*3.14/180),Weight*(1-dbw) );
	      histoManager.fill(metT[kk]+"RecoilUnc",i,Recoils[kk][0],Weight*(1-dbw) );
	      histoManager.fill(metT[kk]+"RecoilParaUnc",i,RecoilProj[kk][0],Weight*(1-dbw) );
	      histoManager.fill(metT[kk]+"RecoilPerpUnc",i,RecoilProj[kk][1],Weight*(1-dbw) );
	      histoManager.fill(metT[kk]+"MTUnc",i,MT[kk],Weight*(1-dbw) );
	     
	    }
	  }
	  
    } //End if condition >
  } //end loop MET

      //Ut/pt
      histoManager.fill("pfRecoilPt",i,Recoils[0][0]/Lepton[0],Weight);
      histoManager.fill("tcRecoilPt",i,Recoils[1][0]/Lepton[0],Weight);
      histoManager.fill("caloRecoilPt",i,Recoils[2][0]/Lepton[0],Weight);
      histoManager.fill("caloT1RecoilPt",i,Recoils[3][0]/Lepton[0],Weight);
      histoManager.fill("caloT2RecoilPt",i,Recoils[4][0]/Lepton[0],Weight);

      histoManager.fill2D("pfdPhiRecoilPt2D",i,Recoils[0][0]/Lepton[0],fabs(dPhiRecoilLep[0]),Weight);
      histoManager.fill2D("tcdPhiRecoilPt2D",i,Recoils[1][0]/Lepton[0],fabs(dPhiRecoilLep[1]),Weight);
      histoManager.fill2D("calodPhiRecoilPt2D",i,Recoils[2][0]/Lepton[0],fabs(dPhiRecoilLep[2]),Weight);
      histoManager.fill2D("caloT1dPhiRecoilPt2D",i,Recoils[3][0]/Lepton[0],fabs(dPhiRecoilLep[3]),Weight);
      histoManager.fill2D("caloT2dPhiRecoilPt2D",i,Recoils[4][0]/Lepton[0],fabs(dPhiRecoilLep[4]),Weight);

      //Ut/pt

      if(name[i].substr(0,4)=="data")
	{
	  std::pair<int,int> tmp(Run,Event);
	  std::pair<string,string> tmp2("data","data");
	  std::pair< std::pair<int,int> , std::pair<string,string> > tmp3(tmp,tmp2);
	  Events.push_back(tmp3);
	}
      if(name[i].substr(0,1)=="W")
	{ S++; 
	}
      else if(name[i]!="data")
	{B++;/* GenRecoil[0]=0; GenRecoil[1]=0;*/ }
      


      //Special variables

      float angle = GenWPt[1]-GenEPt[1];
      if(angle> acos(1.)*2 )
	angle = angle-2*acos(1.);
      else if(angle < acos(1.)*2 )
	angle = 2*acos(1.)-angle;
      
      histoManager.fill("genWPt",i,GenWPt[0],Weight);
    
	if(GenCharge==1) {
	  histoManager.fill("GenParaP",i,GenRecoil[0][0],Weight);
	  histoManager.fill("GenPerpP",i,GenRecoil[0][1],Weight);
	  histoManager.fill2D("GenAngleP",i,GenWPt[0],angle,Weight);
	} else {
	  histoManager.fill("GenParaM",i,GenRecoil[0][0],Weight);
	  histoManager.fill("GenPerpM",i,GenRecoil[0][1],Weight);
	  histoManager.fill2D("GenAngleP",i,GenWPt[0],angle,Weight);
	}

      histoManager.fill2D("pfRecoilvsGen",i,GenRecoil[0][0],RecoilProj[0][0],Weight);
      histoManager.fill2D("tcRecoilvsGen",i,GenRecoil[0][0],RecoilProj[1][0],Weight);
      histoManager.fill2D("caloRecoilvsGen",i,GenRecoil[0][0],RecoilProj[2][0],Weight);
      histoManager.fill2D("caloT1RecoilvsGen",i,GenRecoil[0][0],RecoilProj[3][0],Weight);
      histoManager.fill2D("caloT2RecoilvsGen",i,GenRecoil[0][0],RecoilProj[4][0],Weight);

      histoManager.fill2D("pfAbsRecoilvsGen",i,GenRecoil[0][0],fabs(RecoilProj[0][0]),Weight);
      histoManager.fill2D("tcAbsRecoilvsGen",i,GenRecoil[0][0],fabs(RecoilProj[1][0]),Weight);
      histoManager.fill2D("caloAbsRecoilvsGen",i,GenRecoil[0][0],fabs(RecoilProj[2][0]),Weight);
      histoManager.fill2D("caloT1AbsRecoilvsGen",i,GenRecoil[0][0],fabs(RecoilProj[3][0]),Weight);
      histoManager.fill2D("caloT2AbsRecoilvsGen",i,GenRecoil[0][0],fabs(RecoilProj[4][0]),Weight);

      histoManager.fill2D("pfRecoilvsGenR",i,GenRecoil[1][0],RecoilProj[0][0],Weight);
      histoManager.fill2D("tcRecoilvsGenR",i,GenRecoil[1][0],RecoilProj[1][0],Weight);
      histoManager.fill2D("pfAbsRecoilvsGenR",i,GenRecoil[1][0],fabs(RecoilProj[0][0]),Weight);
      histoManager.fill2D("tcAbsRecoilvsGenR",i,GenRecoil[1][0],fabs(RecoilProj[1][0]),Weight);

      histoManager.fill2D("pfCTRecoilvsGen",i,GenRecoil[0][0],CTRecoil[0][0],Weight);
      histoManager.fill2D("tcCTRecoilvsGen",i,GenRecoil[0][0],CTRecoil[1][0],Weight);
      histoManager.fill2D("pfAbsCTRecoilvsGen",i,GenRecoil[0][0],fabs(CTRecoil[0][0]),Weight);
      histoManager.fill2D("tcAbsCTRecoilvsGen",i,GenRecoil[0][0],fabs(CTRecoil[1][0]),Weight);
      histoManager.fill2D("pfCTRecoilvsGenR",i,GenRecoil[1][0],CTRecoil[0][0],Weight);
      histoManager.fill2D("tcCTRecoilvsGenR",i,GenRecoil[1][0],CTRecoil[1][0],Weight);
      histoManager.fill2D("pfAbsCTRecoilvsGenR",i,GenRecoil[1][0],fabs(CTRecoil[0][0]),Weight);
      histoManager.fill2D("tcAbsCTRecoilvsGenR",i,GenRecoil[1][0],fabs(CTRecoil[1][0]),Weight);
      
      /*     if( SpecRecoilProj[0][0] != -100000 )
	{
	  histoManager.fill2D("pfSpecRecoilvsGen",i,GenRecoil[0][0],SpecRecoilProj[0][0],Weight);
	  histoManager.fill2D("tcSpecRecoilvsGen",i,GenRecoil[0][0],SpecRecoilProj[1][0],Weight);
	  histoManager.fill2D("caloSpecRecoilvsGen",i,GenRecoil[0][0],SpecRecoilProj[2][0],Weight);
	  histoManager.fill2D("caloT1SpecRecoilvsGen",i,GenRecoil[0][0],SpecRecoilProj[3][0],Weight);
	  histoManager.fill2D("caloT2SpecRecoilvsGen",i,GenRecoil[0][0],SpecRecoilProj[4][0],Weight);

	  histoManager.fill2D("pfAbsSpecRecoilvsGen",i,GenRecoil[0][0],fabs(SpecRecoilProj[0][0]),Weight);
	  histoManager.fill2D("tcAbsSpecRecoilvsGen",i,GenRecoil[0][0],fabs(SpecRecoilProj[1][0]),Weight);
	  histoManager.fill2D("caloAbsSpecRecoilvsGen",i,GenRecoil[0][0],fabs(SpecRecoilProj[2][0]),Weight);
	  histoManager.fill2D("caloAbsT1SpecRecoilvsGen",i,GenRecoil[0][0],fabs(SpecRecoilProj[3][0]),Weight);
	  histoManager.fill2D("caloAbsT2SpecRecoilvsGen",i,GenRecoil[0][0],fabs(SpecRecoilProj[4][0]),Weight);

	  histoManager.fill2D("pfSpecRecoilvsGenR",i,GenRecoil[1][0],SpecRecoilProj[0][0],Weight);
	  histoManager.fill2D("tcSpecRecoilvsGenR",i,GenRecoil[1][0],SpecRecoilProj[1][0],Weight);
	  histoManager.fill2D("pfAbsSpecRecoilvsGenR",i,GenRecoil[1][0],fabs(SpecRecoilProj[0][0]),Weight);
	  histoManager.fill2D("tcAbsSpecRecoilvsGenR",i,GenRecoil[1][0],fabs(SpecRecoilProj[1][0]),Weight);

	}*/
     
      NumberEntries[i]++;
	       
    }//End events
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
 
}



TH1F* WUseTree::GetVtxCorrection(string obs,int nds, int bin) {

  //Protection
  if(NVert!=1) {
    cout<<" Warning, no requirement of only 1 Vtx event !!!!!"<<endl;
    abort();
  }
  
  float contamination = 0.05;
  float Systcontamination = 0.02;

  int Obs = histoManager.FindNVar(obs);
  int Ndata= histoManager.access(Obs,nds);
  TH1F* dataTMP = (TH1F*)Histos[ Ndata ]->Clone();

  float Int = dataTMP->Integral(0, dataTMP->GetNbinsX()+1);
  cout<<" vertex correction "<<obs<<endl;
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
  
  cout<<obs<<" Contamination  "<<contamination<<" ;  "<<Nevts<<"  ;  "<<Int<<"   ;  "<<Nevts/Int<<"   "<<endl;
  getRatio =true;
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



void
WUseTree::PrepareDatasets() {

  int nQCD=0;
  int nW=0;

  //We
  //  reposi=;
  /*
  // datasets.push_back("GamJet_PT20_TuneZ2");
  datasets.push_back("GamJet_PT1530_TuneZ2");
  datasets.push_back("GamJet_PT3050_TuneZ2");
  datasets.push_back("GamJet_PT5080_TuneZ2");
  datasets.push_back("GamJet_PT80120_TuneZ2");
  if(QCDType=="EM") {
    datasets.push_back("QCD_EM2030");
    datasets.push_back("QCD_EM3080_TuneZ2_v2");
    datasets.push_back("QCD_EM80170_TuneZ2_v2");
    datasets.push_back("QCD_bc2030");
    datasets.push_back("QCD_bc3080");
    datasets.push_back("QCD_bc80170");
    nQCD = 6;
  }
  else if(QCDType=="QCD"){
    cout<<" *** Using inclusive QCD"<<endl;
    datasets.push_back("QCD");
    nQCD = 1;
  }
  else {cout<<" bad QCD sample !"<<endl; return;}
  datasets.push_back("W_tn_TuneZ2");
  datasets.push_back("Z_2e_PU_TuneZ2");
  
  if(WType=="Pythia")
    { datasets.push_back("W_en_TuneZ2");nW=1;}
  else if(WType=="Pythia_PU")
    { datasets.push_back("W_en_PU");nW=1;}
  else if(WType=="Powheg") {
    datasets.push_back("W_en_minus_powheg_TuneZ2");
    datasets.push_back("W_en_plus_powheg_TuneZ2");
    nW=2;
  }   
  else if(WType=="Powheg_PU") {
    datasets.push_back("W_en_minus_powheg_PU_TuneZ2");
    datasets.push_back("W_en_plus_powheg_PU_TuneZ2");
    nW=2;
  }   
  
  //Combining
  if(QCDType=="EM") {
    int dttmp[12] ={0,0,0,0,1,1,1,1,1,1,2,2}; 
    for(int i=0;i<(nQCD+2+4);i++) {
      FillAddWeight(datasets[i]);
      dt.push_back(dttmp[i]);
    }
  }
  else if(QCDType=="QCD"){
    int dttmp[9] ={0,0,0,0,0,1,2,2};
    for(int i=0;i<(nQCD+7);i++) {
      FillAddWeight(datasets[i]);
      dt.push_back(dttmp[i]);
    }
  }
  
  //Now the Ws
  for(int i=0;i<nW;i++) {
    FillAddWeight(datasets[i+nQCD+2+4]);
    dt.push_back(3);
  }
  nt=4;
  
  colors.push_back(kMagenta+4);
  colors.push_back(kViolet-5);
  colors.push_back(kOrange+7); //+1
  colors.push_back(kOrange-2);*/

    
  colors.push_back(kBlack);//put to kBlack for PAS, otherwise kOrange+3
  
  // string name_tmp[5]={"#gamma+jet","QCD","EWK","W #rightarrow e#nu","data"}; //
  // for(int i=0;i<nt+1;i++)
  // name.push_back(data);
  for(int i=0;i<datasets.size();i++)
    FillAddWeight(datasets[i]);

  name.push_back("data");
  
}



void
WUseTree::LoadDBWeight() {

  ifstream in("/home/mmarionn/Documents/CMS/CMSDatabase/W_Resbos.txt", ios::in );  

  vector<vector<float> > tmpw(50,vector<float>(3,0));

  if(in) {
    cout<<" Loading Database "<<endl;
    for(int i=0;i<50;i++) {
      in >> tmpw[i][0] >> tmpw[i][1] >> tmpw[i][2];
      cout<<"   "<< tmpw[i][0]<<"   "<<tmpw[i][1]<<"   "<<tmpw[i][2]<<endl;
    }
  }
  else 
    cout<<" No DB loaded !! "<<endl;
  
  DBWeights = tmpw;

}

float 
WUseTree::SearchWeight(float Wpt) {
  
  float w=1;
  
  for(int i=0;i<50;i++) {
    if(Wpt >= DBWeights[i][0] && Wpt < DBWeights[i][1])
      { w = DBWeights[i][2]; break;}
  }
  
  return w;

}



float 
WUseTree::SearchResponse(float Wpt,int met) {
  
  float w =1 ;
  for(int i=0;i<MCResponseCor[met].size();i++) {
    if(Wpt >= MCResponseCor[met][i][0] && Wpt < MCResponseCor[met][i][1])
      { w = MCResponseCor[met][i][2]; break;}
  }
  /*  for(int unsigned i=0;i<ResponseDB[met].size()-1;i++) {
    if(Wpt >= ResponseDB[met][i].first && Wpt < ResponseDB[met][i+1].first)
      { w = ResponseDB[met][i].second; break; }
      }*/
  
  return w;

}




vector<vector<vector<float> > >  
WUseTree::GetFitWeight() {

  string metT[6]={"pf","tc","calo","caloT1","caloT2","pfT1"};
 
 // histoManager.PrepareHistos(name);

 vector<int> NumberEntries(nt+1,0);
 string Nvert;
 
 bool  fillBoth= false;

 cout<<" ==================================================== "<<endl;
 cout<<" ================= En cours de weight =============== "<<endl;

  //sÃ©lection
 for(int i=0;i<nt+1;i++) {

    //DÃ©claration des variables
    float IsoVar[3];
    float IDVar[6];
    float Mets[6][2];
    float MetProj[6][2];
    float Recoils[6][2];
    float RecoilProj[6][2];
    float MT[6];
    float SumET[6];
    float Lepton[3];
    float EtSC;

    float WPt[6];

    float dPhiRecoilLep[6];
    float dPhiMETLep[6];

    bool VetoE;
    bool ConvRej;
    bool InHit;
    bool ExpInHit;

    int NVertex;

    float Weight=0;

    float GenWPt[2];

    int Run, Event;

    int JetMultiplicity;
    tChains[i]->SetBranchAddress("JetMultiplicity",&JetMultiplicity);
    tChains[i]->SetBranchAddress("NVertex",&NVertex);

    tChains[i]->SetBranchAddress("IsoVar",IsoVar);
    tChains[i]->SetBranchAddress("IDVar",IDVar);
    tChains[i]->SetBranchAddress("Lepton",Lepton);
    tChains[i]->SetBranchAddress("EtSC",&EtSC);
 
    tChains[i]->SetBranchAddress("WPt",WPt);

    tChains[i]->SetBranchAddress("Mets",Mets);
    tChains[i]->SetBranchAddress("MetProj",MetProj);
    tChains[i]->SetBranchAddress("MT",MT);

    tChains[i]->SetBranchAddress("SumET",SumET);

    tChains[i]->SetBranchAddress("VetoE",&VetoE);
    tChains[i]->SetBranchAddress("InHit",&InHit);
    tChains[i]->SetBranchAddress("ExpInHit",&ExpInHit);
    tChains[i]->SetBranchAddress("ConvRej",&ConvRej);
    
    tChains[i]->SetBranchAddress("Run",&Run);
    tChains[i]->SetBranchAddress("AbsEvent",&Event);

    tChains[i]->SetBranchAddress("BasicRecoils",Recoils);
    tChains[i]->SetBranchAddress("BasicRecoilProj",RecoilProj);
    tChains[i]->SetBranchAddress("BasicdPhiRecoilLep",dPhiRecoilLep);
    tChains[i]->SetBranchAddress("dPhiMETLep",dPhiMETLep);

    tChains[i]->SetBranchAddress("GenWPt",GenWPt);
  
    int EP=0;
    cout<<" beginning tree with vertices  " <<name[i]<<"   "<<tChains[i]->GetEntries()<<endl;
    int ent = tChains[i]->GetEntries();
    //Boucle et sÃ©lection, remplissage histos
    for(int ie=0;ie<ent;ie++) {
   
      if(ie%1000==0)
	fprintf(stdout, "Done : %d / %d\r", ie, ent );

      if(i!=nt)
	Weight = GetWeight(i, ie);
      else
	Weight = 1;
      
      //Reweight PT
      if( name[i].substr(0,4)=="W #r")
        Weight *=  SearchWeight(GenWPt[0]);

      //FIXME...
      if(NVert==1 && i!=nt && i!=nt-1 )
	Weight *= 0.3;

      tChains[i]->GetEntry(ie);      
      
      //Throw multijet events //FIXME
      //if(JetMultiplicity!=1) continue;

      
      if(name[i].substr(0,4)=="data")
	{
	  fillBoth= false;
	  bool doubleCount=false;
	  for(size_t ik=0;ik<Events.size();ik++)
	    {
	      if( ((Events[ik]).first.first) == Run && ((Events[ik]).first.second) == Event)
		{doubleCount=true; break;}
	    }
	  
	  if(doubleCount || (EventFilter && Event<=EventNum ) )
	    continue;
	}

      //Filling Histos
      if(fabs(Lepton[1])>1.479)
	{EP=1; }
      else
	{EP=0; }
      
      if(name[i]!="data" || EP==0 || noDeltaCor) {
	IDVar[4] = 0;
	IDVar[5] =0;
      }

      //first pt cut
      if(EtSC < PTcut) continue;
  
      if(convRej && !ExpInHit)  continue;
      if(vetoE && !VetoE) continue;
    
      //Met Cut
        
      //And Selections
      if(IDVar[0] > _IDcuts[EP][0] && Cbool(skipCut, "sigieie"==observable) ) continue;
      if( fabs(IDVar[1]/*-IDVar[4]*/) > _IDcuts[EP][1] && Cbool(skipCut, "Deta"==observable) ) continue;
      if( fabs(IDVar[2]/*-IDVar[5]*/) > _IDcuts[EP][2] && Cbool(skipCut, "Dphi"==observable) ) continue;
      if(IDVar[3] > _IDcuts[EP][3] && Cbool(skipCut, "HoE"==observable) ) continue;
   
      if(IsoVar[0]/Lepton[0] > _Isocuts[EP][0] && Cbool(skipCut, "TrackIso"==observable) ) continue;
      if(IsoVar[1]/Lepton[0] > _Isocuts[EP][1] && Cbool(skipCut, "EcalIso"==observable) ) continue;
      if(IsoVar[2]/Lepton[0] > _Isocuts[EP][2] && Cbool(skipCut, "HcalIso"==observable) ) continue;
    
      if(EcalP!=EP && EcalP!=2) continue;
     
      // if(name[i].substr(0,4)=="W #r")
      //	fillBoth= false;

      /*  if(NVertex==1 && ( name[i].substr(0,4)=="data" || name[i].substr(0,4)=="W #r"))
	 Nvert="1";
      else if(NVertex>1 && ( name[i].substr(0,4)=="data" || name[i].substr(0,4)=="W #r" ) )
	Nvert="2";
      else if( ( name[i].substr(0,4)!="data" && name[i].substr(0,4)!="W #r") )
      fillBoth=true;*/

      if(NVertex==1 /* && i>=nt-2/*name[i].substr(0,4)=="data"*/)
	Nvert="1";
      else if(NVertex>1/* && i>=nt-2/*name[i].substr(0,4)=="data"*/)
	Nvert="2";
      /*  else if(i<nt-2/*name[i].substr(0,4)!="data")
	fillBoth=true;*/
      
	  //  if(i>=nt-2) {
	fillBoth=false;
	//   }
      
      for(int kk=0;kk<6;kk++) {

	//Reweight PT
	if( name[i].substr(0,4)=="W #r")
	  Weight *=  SearchWeight(WPt[kk]);

	if(!fillBoth) {
	  //MET
	    histoManager.fill(metT[kk]+"MET_V"+Nvert,i,Mets[kk][0],Weight);
	    
	   
	  //MET Proj
	  //  histoManager.fill(metT[kk]+"METPara_V"+Nvert,i,MetProj[kk][0],Weight);
	  
	  if( (Mets[kk][0] > METcut && !invCut && MT[kk] > MTCut) || 
	      (Mets[kk][0] < METcut && invCut && MT[kk] > MTCut ) ) { 
	  
	    histoManager.fill(metT[kk]+"Recoil_V"+Nvert,i,Recoils[kk][0],Weight);
	    histoManager.fill(metT[kk]+"SumETRecoil_V"+Nvert,i,SumET[kk]-Lepton[0],Weight);
	  }
	  
	  if(NVertex>1/* && i>=nt-2*//* && name[i].substr(0,4)=="data"*/) {
	   
	    histoManager.fill(metT[kk]+"METControl",i,Mets[kk][0],Weight);
	    histoManager.fill(metT[kk]+"MET_V1Control",i,Mets[kk][0],Weight);
	    histoManager.fill(metT[kk]+"METParaControl",i,MetProj[kk][0],Weight);
	    histoManager.fill(metT[kk]+"METPerpControl",i,MetProj[kk][1],Weight);
	

	    if( (Mets[kk][0] > METcut && !invCut && MT[kk] > MTCut) || 
		(Mets[kk][0] < METcut && invCut && MT[kk] > MTCut ) ) { 
	      
	      if(kk==0)
		histoManager.fill("PtControl",i,Lepton[0],Weight);
	      histoManager.fill(metT[kk]+"AcopControl",i,3.14-fabs(dPhiMETLep[kk]*3.14/180),Weight);
	      histoManager.fill(metT[kk]+"dPhiRecoilControl",i,fabs(dPhiRecoilLep[kk]*3.14/180),Weight);
	      histoManager.fill(metT[kk]+"RecoilControl",i,Recoils[kk][0],Weight);

	      float dbrcor = SearchResponse(WPt[kk],kk);
	      histoManager.fill(metT[kk]+"RecoilCorControl",i,Recoils[kk][0]*dbrcor,Weight);

	      histoManager.fill(metT[kk]+"RecoilParaControl",i,RecoilProj[kk][0],Weight);
	      histoManager.fill(metT[kk]+"RecoilPerpControl",i,RecoilProj[kk][1],Weight);
	      histoManager.fill(metT[kk]+"MTControl",i,MT[kk],Weight);
	      histoManager.fill(metT[kk]+"SumETRecoil_V1Control",i,SumET[kk]-Lepton[0],Weight);
	    }
	    
	    //Corrected variables ======================
	    /*  if( kk==0) {
	      TVector2 rectmp(0,0);
	      rectmp.SetMagPhi(Recoils[kk][0],Recoils[kk][1]);
	      TVector2 leptmp(0,0);
	      leptmp.SetMagPhi(Lepton[0],Lepton[2]);
	      TVector2 GenWpttmp(0,0);
	      GenWpttmp.SetMagPhi(GenWPt[0],GenWPt[1]);
	      
	      //Unit wPt vector
	      TVector2 WPtUnit = GenWpttmp.Unit();
	      
	      //Corrections
	      TVector2 correc(0,0);
	      
	      float upara = WPtUnit*rectmp;
	      TVector2 tt = WPtUnit.Rotate(acos(-1)/2);
	      float uperp = tt*rectmp;
	      
	      float paraoff = GetParaRecoilCor(kk,NVertex,GenWPt[0]);
	      float perpoff = GetPerpRecoilCor(kk,NVertex,GenWPt[0]);
	      float parares = GetParaRecoilCorRes(kk,NVertex,GenWPt[0]);
	      float perpres = GetPerpRecoilCorRes(kk,NVertex,GenWPt[0]);
	      
	      float uparac = (upara*paraoff)*parares;
	      float uperpc = (uperp+perpoff)*perpres;
	      float mag = sqrt(uparac*uparac+uperpc*uperpc);
	      
	      TVector2 partmp(0,0); 
	      partmp.SetMagPhi(uparac, WPtUnit.Phi() + acos(-1) );
	      TVector2 pertmp(0,0); 
	      pertmp.SetMagPhi(fabs(uperpc), (WPtUnit.Phi() + acos(-1)/2.*uperpc/fabs(uperpc) ) );
	      
	      correc = partmp + pertmp;
	      
	      TVector2 metcor(0,0);
	      metcor -= correc + leptmp;
	      
	      float MTcor = sqrt(2*Lepton[0]*metcor.Mod()* (1-cos(metcor.DeltaPhi(leptmp))) );
	      //+++++++++++
	      
	      histoManager.fill(metT[kk]+"METCorControl",i,metcor.Mod(),Weight);
	      histoManager.fill(metT[kk]+"RecoilCorControl",i, mag,Weight);
	      histoManager.fill(metT[kk]+"MTCorControl",i,MTcor,Weight);
	    } 
	    else {*/
	      histoManager.fill(metT[kk]+"METCorControl",i,Mets[kk][0],Weight);
	      histoManager.fill(metT[kk]+"RecoilCorControl",i,Recoils[kk][0] ,Weight);
	      histoManager.fill(metT[kk]+"MTCorControl",i, MT[kk],Weight);
	   
	      //	    }

      //==========================================


	  }
	}
	else {
	  //MET
	  histoManager.fill(metT[kk]+"MET_V1",i,Mets[kk][0],Weight);
	  histoManager.fill(metT[kk]+"MET_V2",i,Mets[kk][0],Weight);
	
	

	  if( (Mets[kk][0] > METcut && !invCut && MT[kk] > MTCut) || 
	      (Mets[kk][0] < METcut && invCut && MT[kk] > MTCut ) ) { 
	  
	    Nvert="1";
	    histoManager.fill(metT[kk]+"SumETRecoil_V"+Nvert,i,SumET[kk]-Lepton[0],Weight);
	    histoManager.fill(metT[kk]+"Recoil_V1",i,Recoils[kk][0],Weight);
	    histoManager.fill(metT[kk]+"Recoil_V2",i,Recoils[kk][0],Weight);
	  }

	  //MET Proj
	  // histoManager.fill(metT[kk]+"METPara_V1"+Nvert,i,MetProj[kk][0],Weight);
	  //histoManager.fill(metT[kk]+"METPara_V2"+Nvert,i,MetProj[kk][0],Weight);
	}
      }
    }//End events
  }// End datasets
  

 Histos = histoManager.GetHistos();

 Events.clear();
    
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
  if(!Remove2V) {
    for(int i=0;i<2;i++) {
      for(int kk=0;kk<6;kk++) {
	ostringstream os;
	os << i+1;
	TString NameO = (TString)metT[kk]+"MET_V"+ os.str();
	Wghts[kk][i]= GetAlpha(NameO);
      }
    }
  }  
  
  return Wghts;
}




void
WUseTree::LoadDBResponse() {
  /*
  TFile* dbfile = new TFile("../Response/GammaJet2Aug.root","READ");

  vector<vector<pff> > tmpdbr;
  vector< pff > tmppair;

  vector<string> ndb;
  ndb.push_back("pfMET");
  ndb.push_back("tcMET");
  ndb.push_back("CaloMET");
  ndb.push_back("CaloMETtypeI");
  ndb.push_back("CaloMETtypeII");
 
  for(size_t im=0;im<ndb.size();im++)
    {
      TH1F* tmp = (TH1F*)dbfile->Get(ndb[im].c_str() );
      tmpdbr.push_back(tmppair);
  
      for(int ib=0;ib<tmp->GetNbinsX();ib++)
	{
	  pff tp(tmp->GetBinLowEdge(ib), 1./tmp->GetBinContent(ib) );
	  if(tmp->GetBinContent(ib) == 0)
	    tp.second = 1;
	  
	  tmpdbr[im].push_back(tp);
	}
    }

  ResponseDB = tmpdbr;
  */

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
WUseTree::GetMCWeightFromDataShape(int met, float var, string shref, string shvar, double bin) {

  if(shape[met]==NULL) {
    cout<<" Utilisation W MCWeight shape "<<endl;
    cout<<" MET "<<met<<"  encore nul "<<"   "<<shref<<"   "<<shvar<<"    "<<endl;
    int Sref =  histoManager.FindNVar(shref);
    int Svar =  histoManager.FindNVar(shvar);
 
    int nhref= histoManager.access(Sref,nt);
 
    TH1F* href = (TH1F*)Histos[ nhref ]->Clone();

    cout<<" numbers "<<" ref  "<<nhref<<endl;

    vector<TH1F*> hVars;
    for(int i=0;i<nt;i++) {
      int nhvar= histoManager.access(Svar,i);
      hVars.push_back( (TH1F*)Histos[nhvar]->Clone() );
    }

    float NEWK=0,NQCD=0,NData=0;

    NEWK = hVars[3]->Integral(0, hVars[3]->GetNbinsX()+1) + hVars[2]->Integral(0, hVars[2]->GetNbinsX()+1);
    NQCD = hVars[0]->Integral(0, hVars[0]->GetNbinsX()+1) + hVars[1]->Integral(0, hVars[1]->GetNbinsX()+1);
    NData = href->Integral(0, href->GetNbinsX()+1);
    
    if(FitR && !MultiWeight && !NoData) {
      //Get MT distribution and determine a************************
      TString obsMT;
      if(shref.substr(0,2)=="pf")
	obsMT="pf";
      else if(shref.substr(0,2)=="tc")
	obsMT="tc";
      else if(shref.substr(0,6)=="caloT1")
	obsMT="caloT1";
      else if(shref.substr(0,6)=="caloT2")
	obsMT="caloT2";
      else
	obsMT="calo";
      cout<<" NEWK "<<NEWK<<"  ;NQCD  "<<NQCD<<"  ;NData  "<<NData<<endl;
      vector<float> a = GetAlpha(obsMT+"MET_V1");
      //***********************************************************
      float NormEWK2 = a[0];
      float NormQCD2 = a[1];

      float R = NData/(NormEWK2*NEWK+NormQCD2*NQCD);
      for(int i=0;i<nt;i++) {
	if(hVars[i]->Integral()!=0) {

	  if(i>1)
	    hVars[i]->Scale( NormEWK2*R );
	  else 
	    hVars[i]->Scale( NormQCD2*R );
   
	  float integ = hVars[i]->Integral(0, hVars[i]->GetNbinsX()+1 ); 
	  cout<<" ---> After reweight "<<name[i] << " -> "<<integ<<endl;
	}
      }
    }

      vector<float> templ = histoManager.GetTemplate(Svar);
      TH1F* hvar =new TH1F("tmphistosforweight","Observable",(int)templ[0],templ[1],templ[2]);
      for(int i=0;i<nt;i++) { 
	if(i==0) {
	  hvar = (TH1F*)hVars[i]->Clone();
	}
	else {
	  hvar->Add(hVars[i],1);
	}
      }

      href->Rebin( /*Bin*pas*/bin/href->GetBinWidth(1) );
      hvar->Rebin( /*Bin*pas*/bin/hvar->GetBinWidth(1) );
 
      float integref = href->Integral(0, href->GetNbinsX() +1 );
      float integvar = hvar->Integral(0, href->GetNbinsX() +1 );
      cout<<" integral  "<<integref<<"   "<<integvar<<endl;
      href->Scale( 1./integref );
      hvar->Scale( 1./integvar );

      shape[met] = HistoManager::ShapeWeight(href,hvar);
      delete href;
      delete hvar;
  }

  float w = shape[met]->GetBinContent( shape[met]->FindBin(var) );
  /*if(met==1)
    cout<<shape[met]->FindBin(var)<<"   "<<weight<<endl;*/
  int ii= shape[met]->FindBin(var);
  while(w==0 && ii!=0 )
    { w =  shape[met]->GetBinContent( ii ); ii--; }
  
  return w;
}


float WUseTree::GetPerpRecoilCorRes(int met, int Nvtx, float GenWqT) {

  
  if(met!=0) return 1; //FIXME, need all corrections

  float datval = sqrt( pow(0.56*sqrt(GenWqT + 17) + 2 ,2) + (Nvtx-1)*pow( 2.37 ,2 ) );
  float mcval = sqrt( pow(0.33*sqrt(GenWqT + 17) + 3 ,2) + (Nvtx-1)*pow( 3.12 ,2 ) );

  return datval/mcval;

}


float WUseTree::GetParaRecoilCorRes(int met, int Nvtx, float GenWqT) {

  if(met!=0) return 1; //FIXME, need all corrections
  
  float datval = sqrt( pow(0.8*sqrt(GenWqT + 1) + 3 ,2) + (Nvtx-1)*pow( 2.99 ,2 ) );
  float mcval = sqrt( pow(0.82*sqrt(GenWqT + 1) + 2.5 ,2) + (Nvtx-1)*pow( 3.06 ,2 ) );
  
  return datval/mcval;
    
}

float WUseTree::GetPerpRecoilCor(int met, int Nvtx, float GenWqT) {

  return -0.16;

}


float WUseTree::GetParaRecoilCor(int met, int Nvtx, float GenWqT) {

  if(met!=0) return 1;
  else 
    return  SearchRecoilCor(GenWqT);
  

}


void WUseTree::LoadRecoilCorrections() {

  ifstream in("/home/mmarionn/Documents/CMS/CMSDataBase/ZMCResponse_v2.db", ios::in );  

  vector<float> tmpw(3,0);

  if(in) {
    cout<<" Loading Database "<<endl;
    
    float m,M,cor;
    int met;
    
    while(!in.eof()) {
      in >> met >> m >> M >> cor;
      tmpw[0] = m;
      tmpw[1] = M;
      tmpw[2] = cor;
      RecoilCors.push_back(tmpw);
    }
  }
  else 
    cout<<" No response DB loaded !! "<<endl;
  
}


float WUseTree::SearchRecoilCor(float GenWpt) {

  
  float w=1;
  
  for(int i=0;i<RecoilCors.size();i++) {
    if(GenWpt >= RecoilCors[i][0] && GenWpt < RecoilCors[i][1])
      { w = RecoilCors[i][2]; break;}
  }
  
  return w;
  
}
