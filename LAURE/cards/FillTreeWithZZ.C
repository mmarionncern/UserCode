{

  if(Recompute)
    ZZUseTree Analyse;
  
  //paramètres généraux ********************* paramètres généraux
  string repository="ZZAnalysis";
  if(Recompute) {
    Analyse.SetTreeName("ZZ_2e2n/");
  }


  vector<string> data;
  
  /* data.push_back("EG_Nov4_132440_141000");
  data.push_back("EG_Nov4_141001_143000");
  data.push_back("EG_Nov4_143001_144114");*/
  /*data.push_back("Electron_Nov4_146240_147000");
  data.push_back("Electron_Nov4_147001_148000");
  data.push_back("Electron_Nov4_148001_149000");
  data.push_back("Electron_Nov4_149001_149442");*/
  
  data.push_back("EG_Dec22_132440_141000");
  data.push_back("EG_Dec22_141001_143000");
  data.push_back("EG_Dec22_143001_144114");
  data.push_back("Electron_Dec22_146240_147000");
  data.push_back("Electron_Dec22_147001_148000");
  data.push_back("Electron_Dec22_148001_149000");
  data.push_back("Electron_Dec22_149001_149442");
  
  data.push_back("Mu_Dec22_132440_141000");
  data.push_back("Mu_Dec22_141001_143000");
  data.push_back("Mu_Dec22_143001_144114");
  data.push_back("Mu_Dec22_146240_147000");
  data.push_back("Mu_Dec22_147001_148000");
  data.push_back("Mu_Dec22_148001_149000");
  data.push_back("Mu_Dec22_149001_149442");
  
  
  bool MCOnly = false;

  bool EventFilter=false;
  int EventNumber=139466;

  string QCDType="";
  string ZType="Powheg";

  //Variable à étudier ********************** Variable à étudier
  string observable="ZCateg";
  
  bool Draw3on1 =false;
  vector<string> SevObs;
  SevObs.push_back("Ptl1");
  SevObs.push_back("Ptl2");
  SevObs.push_back("Ptl1");

  bool Fill2DHistos=false;
  bool FillProfile=false;
  bool response=false;

  //Binning & titre ************************* Binning & titre
  string YTitle="nombre d'événements ";
  int Binning=4;
  int AddBinBkg=1; //BinB = binning*AddBin
  double RangeY[2]={0.001,2};
  double RangeX[2]={-1,500};
  int Xdiv[3]={6,5,0};
  int Ydiv[3]={5,5,0}; //Nlabel /  sous-Div /ssdiv
  bool SuperColor=false;
  bool OverFlowBin=true;
  bool UnderFlowBin=false;
  bool ShowDMCRatio=false;
  float MarkerSize=0.8;
  float LineWidth=2;

  bool BasicProf=true;
  bool switchRMS=false;
  string errorOpt="s";

  bool savePlot = false;
  //Paramètres de l'analyse ****************** Paramètres de l'analyse
  bool N1Plot=false;
  string Nm1Var = "rejection #gamma";

  int isEndcaps=2; //1=EE, 0=EB , 2=combined

  bool ConvRejection=true;
  bool NoBackground=false;
  int WPIso = 80;
  int WPID = 80;

  bool ForceCut=false;

  //                     sieie  deta  dphi  hoe
  float IdCuts[2][4] ={{ -1,   -1,   -1,   -1 }, //EB
		       { -1,   1000,  -1,   -1 }}; //EE
  //                       trk  ecal  hcal
  float IsoCuts[2][3] ={{ -1,    -1,   -1  },   //EB
			{ -1,    -1,   -1  } }; //EE
 
  int NVertex=0;
  float PTcut=20;
  string METType="pfMET";
  float METcut=4.2; //4.2 // 4.5
  float MassCut=60;
  bool inversionMassCut = false;
  bool inversionMetCut = false;
  bool noDeltaCor=false;
  string Norm="";


  float lumi=36; //pb-1
 
  //paramètres additionnels
  bool UsePdgId = true; //remove for WW cross check
  int pdgID = 11;
  int MaxJet = 1; //+1
  bool fixJet = false;
  bool inverseJet = false;
  int MaxPhoton = 100; //+1
  int MaxLepton = 1; //+1
  float NormMETCutHigh = 2; //2
  float NormMETCutLow = -0.25; //0.25
  float dPhiJMCut = 20; //20
  float NNoutMETCut = -0.9; //0.85
  float NNoutCut = 0.15; //0.15
  bool bTag=true;

  bool useAccForEff=true;
  //Lumis pb-1 & KFactors ************************************
  map<string,float> Lumis;

  Lumis[ "GamJet_PT1530" ] = 5.97 ;
  Lumis[ "GamJet_PT3050" ] = 61.41 ;
  Lumis[ "GamJet_PT5080" ] = 376.69 ;
  Lumis[ "GamJet_PT80120" ] = 2343.95 ;
  Lumis[ "GamJet_PT120170" ] = 11947.74 ;
  Lumis[ "QCD_EM2030" ] = 15.04 ;
  Lumis[ "QCD_EM3080" ] = 18.55 ;
  Lumis[ "QCD_EM80170" ] = 57.87 ;
  Lumis[ "QCD_bc2030" ] = 16.98 ;
  Lumis[ "QCD_bc3080" ] = 14.58;
  Lumis[ "QCD_bc80170" ] = 111.47;
  Lumis[ "ttbar" ] = 9634.02;          //
  Lumis[ "T_schan" ] = 499966.66;      //
  Lumis[ "T_tchan" ] = 23050.48;       //
  Lumis[ "TW" ] = 46871.30;            //
  Lumis[ "VGam" ] = 6371.06;           //
  Lumis[ "WJets" ] = 615.05;           //
  Lumis[ "Z_2t" ] = 860.04;            //
  Lumis[ "Z_2m" ] = 1238.49;           //
  Lumis[ "Z_2e" ] = 1238.53;           //
  Lumis[ "ZJets" ] = 2265.1;
  Lumis[ "ZZ_2l2n" ] = 736455.81 ;     //
  Lumis[ "WZ_3ln" ] = 323529.41;       //
  Lumis[ "WW_2l2n" ] = 37931.04;       //
  Lumis[ "gg_WW_2l2n" ] = 715175 ;     //

  /*NOV4  Lumis[ "GamJet_PT1530" ] = 6.03 ;
  Lumis[ "GamJet_PT3050" ] = 62.15;
  Lumis[ "GamJet_PT5080" ] = 376.7;
  Lumis[ "GamJet_PT80120" ] = 2349.2;
  Lumis[ "QCD_EM2030" ] = 13.55;
  Lumis[ "QCD_EM3080" ] = 17.45;
  Lumis[ "QCD_EM80170" ] = 60.2;
  Lumis[ "QCD_bc2030" ] = 20.8;
  Lumis[ "QCD_bc3080" ] = 8.1;
  Lumis[ "QCD_bc80170" ] = 33.5;
  Lumis[ "WW_2l2n" ] = 37500;
  Lumis[ "WZ_3ln" ] = 321600;
  Lumis[ "ZZ_2l2n" ] = 634218;
  Lumis[ "ZZ_2l2n_PU" ] = 736456;
  Lumis[ "WGam" ] = 7110 ;
  Lumis[ "ZGam_2l" ] = 22600;
  Lumis[ "WJets_PU" ] = 615.6;
  Lumis[ "Z_2t_PU" ] = 1683;
  Lumis[ "Z_2e_PU" ] = 1604.30;
  Lumis[ "Z_2e_PU_Powheg" ] = 1998990/1666.4/0.95;
  Lumis[ "Z_2e_16PU" ] = 1703450/1300;
  Lumis[ "Z_2m" ] = 1604.30;
  Lumis[ "ZJets" ] = 641;
  Lumis[ "ttbar" ] = 10600;*/

  map<string,float> KFactors;
  //Nov4
  // KFactors[ "Z_2t_PU" ] = 1.3*.95;
  // KFactors[ "Z_2e_PU" ] = 1.3*.95;
  // KFactors[ "Z_2e_16PU" ] = 1.3*.95/1.2;
  //KFactors[ "Z_2e_PU_Powheg" ] = 1.3*.95;
  // KFactors[ "Z_2m" ] = 1.3;
  // KFactors[ "ZJets" ] = 1.28*.95;
  // KFactors[ "WJets_PU" ] = 1.3*.95;
  // KFactors[ "WW" ] = 1.53;
  // KFactors[ "WW_2l2n" ] = 1.53;
  // KFactors[ "WZ" ] = 1.71;
  // KFactors[ "WZ_3ln" ] = 1.71;
  // KFactors[ "ZZ" ] = 1.37;
  // KFactors[ "ZZ_4l" ] = 1.37;
  // KFactors[ "ZZ_2l2n" ] = 1.37;
  // KFactors[ "ttbar" ] = 1.75;
  // KFactors[ "T_tchan" ] = 1.75;
  // KFactors[ "T_schan" ] = 1.75;
  // KFactors[ "TW" ] = 1.75;


  KFactors[ "Z_2t" ] = 1.032*0.956; // FEWZ->r mll_20->60 : 0.6098 ->931/984=0.95 (NNLOkF=1.03)
  KFactors[ "Z_2e" ] = 1.032*0.956;
  KFactors[ "Z_2m" ] = 1.032*0.956;
  KFactors[ "ZJets" ] = 1.31;
  KFactors[ "ZZ_2l2n" ] = 1; //1.37
  KFactors[ "WW_2l2n" ] = 0.97; //1.48 //2.9 /   4.3 pb tot  -> 1.48
  KFactors[ "WZ_3ln" ] = 1; //1.73
  KFactors[ "VGam" ] = 1.11; //0.95 
  KFactors[ "WJets" ] = 1.27;
  KFactors[ "ttbar" ] = 173/121.; //1.3
  KFactors[ "T_tchan" ] = 64.6/21; //1.75
  KFactors[ "T_schan" ] = 4.6/0.99; //1.41
  KFactors[ "TW" ] = 10.6/10.56; //1.75           //all top : 153.55 / 173  -> 1.13
 
//MonteCarlo Samples ************************** MC samples
  if(Recompute) {
    Analyse.AddMCSample( "ZZ_2l2n"        ,                        "ZZ",   kViolet+7 );
    Analyse.AddMCSample( "WZ_3ln"         ,                        "WZ",   kViolet-2 );
    Analyse.AddMCSample( "WW_2l2n"        , "WW #rightarrow (l#nu)(l#nu)",    kTeal-7  );
    Analyse.AddMCSample( "ttbar"          ,                  "t#bar{t}",      kRed+1 );
    Analyse.AddMCSample( "T_schan"        ,                          "",           0 );
    Analyse.AddMCSample( "T_tchan"        ,                          "",           0 );
    Analyse.AddMCSample( "TW"             ,                          "",           0 );
    /*Analyse.AddMCSample( "GamJet_PT1530"  ,                       "QCD",  kMagenta+4 );
    Analyse.AddMCSample( "GamJet_PT3050"  ,                          "",           0 );
    Analyse.AddMCSample( "GamJet_PT5080"  ,                          "",           0 );
    Analyse.AddMCSample( "GamJet_PT80120" ,                          "",           0 );
    Analyse.AddMCSample( "QCD_bc2030"     ,                          "",           0 );
    Analyse.AddMCSample( "QCD_bc3080"     ,                          "",           0 );
    Analyse.AddMCSample( "QCD_bc80170"    ,                          "",           0 );
    Analyse.AddMCSample( "QCD_EM2030"     ,                          "",           0 );
    Analyse.AddMCSample( "QCD_EM3080"     ,                          "",           0 );
    Analyse.AddMCSample( "QCD_EM80170"    ,                          "",           0 );*/
    Analyse.AddMCSample( "VGam"           ,                   "V#gamma",   kOrange+7 );
    Analyse.AddMCSample( "WJets"          ,                    "W+jets",     kBlue+1 ); //W+jets
    Analyse.AddMCSample( "ZJets"          ,    "Z #rightarrow ll + X",   kOrange-2 );
    /* Analyse.AddMCSample( "Z_2t"           ,      "Z #rightarrow ll + X",   kOrange-2 );
    if(UsePdgId) {
      if(pdgID==11)
	Analyse.AddMCSample( "Z_2e"       ,                          "",   kOrange-2 );
      else
	Analyse.AddMCSample( "Z_2m"       ,                          "",  kOrange-2 );
    }
    else {
      Analyse.AddMCSample( "Z_2e"         ,                          "",   kOrange-7 );
      Analyse.AddMCSample( "Z_2m"         ,                          "",           0 );
    }
    */
    
  }

  //Lines *********************************** Lines
  vector<int> Lines;
  
  //*********************************************************************
  //Execution macro ******************************************************
  
  Analyse.ConfigureData(EventFilter,EventNumber,MCOnly);
  Analyse.ConfigureLumi(Lumis, KFactors,lumi);
  Analyse.ConfigureEndAnalysis(ForceCut, WPIso, WPID, IdCuts,IsoCuts,PTcut,METcut,inversionMassCut,METType,MassCut,NVertex,Norm);
  Analyse.ConfigurePlots(Binning,AddBinBkg,OverFlowBin, UnderFlowBin,RangeY,RangeX, YTitle,Xdiv,Ydiv,MarkerSize,LineWidth,Lines,SevObs,Draw3on1,SuperColor,switchRMS,errorOpt,FillProfile,BasicProf,ShowDMCRatio);
  Analyse.SavePlot(savePlot);
  if(Recompute) {
    Analyse.AddVariables(UsePdgId, pdgID, MaxJet, MaxPhoton, MaxLepton, dPhiJMCut, NormMETCutHigh, NormMETCutLow, NNoutMETCut,fixJet,inverseJet,bTag,inversionMetCut, NNoutCut, useAccForEff);

  Analyse.EndAnalysis(data,observable,isEndcaps,ConvRejection,NoBackground,
		      N1Plot,Nm1Var,noDeltaCor,QCDType,ZType, repository);
  Recompute=false;
  Analyse.FillZZTree();
  }
 
  Analyse.PlotDistribution(observable,Fill2DHistos,FillProfile);

}
