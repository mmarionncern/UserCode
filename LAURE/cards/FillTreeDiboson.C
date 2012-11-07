{

  if(Recompute)
    DibosonUseTree Analyse;
  
  //paramètres généraux ********************* paramètres généraux
  string repository="Diboson";

  vector<string> MCAddTrees;
  vector<string> LepAddTrees;
  if(Recompute) {
    Analyse.SetTreeName("ZZ_2e2n/ZZtuple");
    Analyse.SetTreeName("ZZ_2m2n/ZZtuple");

    MCAddTrees.push_back("ZZ_2e2n/VV_mcTruth");
    MCAddTrees.push_back("ZZ_2m2n/VV_mcTruth");
    
    LepAddTrees.push_back("ZZ_2e2n/ZZLeptons");
    //LepAddTrees.push_back("ZZ_2m2n/ZZLeptons");
  }



  vector<string> data;
  
  data.push_back("DoubleEle_160404_162000");
  data.push_back("DoubleEle_162001_163000");
  data.push_back("DoubleEle_163001_163869");
  data.push_back("DoubleEle_165071_165400");
  data.push_back("DoubleEle_165401_165700");
  data.push_back("DoubleEle_165543_165900");
  data.push_back("DoubleEle_165901_166200");
  data.push_back("DoubleEle_166201_166502");
  data.push_back("DoubleEle_166503_166600");
  data.push_back("DoubleEle_166601_166700");
  data.push_back("DoubleEle_166701_166861");
  data.push_back("DoubleEle_166862_167000");
  data.push_back("DoubleEle_167001_167151");
  data.push_back("DoubleEle_167152_167450");
  data.push_back("DoubleEle_167451_167784");
  data.push_back("DoubleEle_167785_167913");

  data.push_back("DoubleMu_160404_162000");
  data.push_back("DoubleMu_162001_163000");
  data.push_back("DoubleMu_163001_163869");
  data.push_back("DoubleMu_165071_165400");
  data.push_back("DoubleMu_165401_165700");
  data.push_back("DoubleMu_165543_165900");
  data.push_back("DoubleMu_165901_166200");
  data.push_back("DoubleMu_166201_166502");
  data.push_back("DoubleMu_166503_166600");
  data.push_back("DoubleMu_166601_166700");
  data.push_back("DoubleMu_166701_166861");
  data.push_back("DoubleMu_166862_167000");
  data.push_back("DoubleMu_167001_167151");
  data.push_back("DoubleMu_167152_167450");
  data.push_back("DoubleMu_167451_167784");
  data.push_back("DoubleMu_167785_167913");
  
  
  bool MCOnly = false;

  bool EventFilter=false;
  int EventNumber=139466;

  string QCDType="";
  string ZType="Powheg";

  //Variable à étudier ********************** Variable à étudier
  string observable="ZMassPresel";
  
  bool Draw3on1 =false;
  vector<string> SevObs;
  SevObs.push_back("");
  SevObs.push_back("J1DR2");
  SevObs.push_back("J2DR2");

  bool Fill2DHistos=false;
  bool FillProfile=false;
  bool response=false;

  //Binning & titre ************************* Binning & titre
  string YTitle="number of events ";
  int Binning=1;
  int AddBinBkg=1; //BinB = binning*AddBin
  double RangeY[2]={0.03,500000};
  double RangeX[2]={60,120};
  int Xdiv[3]={6,5,0};
  int Ydiv[3]={5,5,0}; //Nlabel /  sous-Div /ssdiv
  bool SuperColor=false;
  bool OverFlowBin=false;
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
  bool InvCut=false;
  string Nm1Var = "met";

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
  float METcut=50;
  float MassCut=60;
  bool noDeltaCor=false;
  string Norm="";

  float lumi=1070; //pb-1
 
  //paramètres additionnels
  bool UsePdgId = false; //remove for WW cross check
  int pdgID = 13;
  int MaxJet = 1; 
  string typeJet = "<";
  int MaxLepton = 1; 
  float balMax = 1.8;
  float balMin = 0.4;
  float dPhiJMCut = 20;
  float dPhiZMCut = 60;
  bool bTag=false; //true= b jet in the event
  float advMassMin = 80;
  float advMassMax = 100;

  bool useAccForEff=false;
  //Lumis pb-1 & KFactors ************************************
  map<string,float> Lumis;


  Lumis[ "ttbar" ] = 9621.55;          //
  Lumis[ "T_schan" ] = 499966.66;      //
  Lumis[ "T_tchan" ] = 23050.48;       //
  Lumis[ "TW" ] = 46346.30;            //
  Lumis[ "VGam" ] = 6371.06;           //
  Lumis[ "WJets" ] = 11542792/24640.;           //
  Lumis[ "Z_2t" ] = 1236.29;            //
  Lumis[ "Z_2m" ] = 1185.97;           //
  Lumis[ "Z_2e" ] = 1234.37;           //
  Lumis[ "ZJets" ] = 970.6;
  Lumis[ "ZZ_2l2n" ] = 736455.81 ;     //
  Lumis[ "WZ_3ln" ] = 323529.41;       //
  Lumis[ "WW_2l2n" ] = 37931.04;       //
  Lumis[ "gg_WW_2l2n" ] = 715175 ;     //
  Lumis[ "ZZ_2m_Pythia" ] = 1070 ;     //

  map<string,float> KFactors;
  KFactors[ "Z_2t" ] = 1.032*0.956; // FEWZ->r mll_20->60 : 0.6098 ->931/984=0.95 (NNLOkF=1.03)
  KFactors[ "Z_2e" ] = 1.032*0.956;
  KFactors[ "Z_2m" ] = 1.032*0.956;
  KFactors[ "ZJets" ] = 1.31;
  KFactors[ "ZZ_2l2n" ] = 1.12; //1.37
  KFactors[ "WW_2l2n" ] = 1; //1.48 //2.9 /   4.3 pb tot  -> 1.48
  KFactors[ "WZ_3ln" ] = 1; //1.73
  KFactors[ "VGam" ] = 1.11; //0.95 
  KFactors[ "WJets" ] = 1.27;
  KFactors[ "ttbar" ] = 173/121.; //1.3
  KFactors[ "T_tchan" ] = 64.6/21; //1.75
  KFactors[ "T_schan" ] = 4.6/0.99; //1.41
  KFactors[ "TW" ] = 10.6/10.56; //1.75           //all top : 153.55 / 173  -> 1.13
 
  //MonteCarlo Samples ************************** MC samples
  if(Recompute) {
    Analyse.AddMCSample( "ZZ_2l2n"        , "ZZ #rightarrow (l#nu)(l#nu)",   kOrange-2 );
    Analyse.AddMCSample( "WZ_3ln"         , "WZ #rightarrow (ll)(l#nu)",   kOrange+7 ); //WZ #rightarrow (ll)(l#nu)
    Analyse.AddMCSample( "WW_2l2n"        , "WW #rightarrow (l#nu)(l#nu)", kRed+1    ); //WW #rightarrow (l#nu)(l#nu)
    Analyse.AddMCSample( "T_schan"        ,                  "single-t",   kViolet+1 );
    Analyse.AddMCSample( "T_tchan"        ,                          "",           0 );
    Analyse.AddMCSample( "TW"             ,                          "",           0 );
    Analyse.AddMCSample( "ttbar"          ,                  "t#bar{t}",  kMagenta+3 );
    // /*Analyse.AddMCSample( "GamJet_PT1530"  ,                       "QCD",  kMagenta+4 );
    // Analyse.AddMCSample( "GamJet_PT3050"  ,                          "",           0 );
    // Analyse.AddMCSample( "GamJet_PT5080"  ,                          "",           0 );
    // Analyse.AddMCSample( "GamJet_PT80120" ,                          "",           0 );
    // Analyse.AddMCSample( "QCD_bc2030"     ,                          "",           0 );
    // Analyse.AddMCSample( "QCD_bc3080"     ,                          "",           0 );
    // Analyse.AddMCSample( "QCD_bc80170"    ,                          "",           0 );
    // Analyse.AddMCSample( "QCD_EM2030"     ,                          "",           0 );
    // Analyse.AddMCSample( "QCD_EM3080"     ,                          "",           0 );
    // Analyse.AddMCSample( "QCD_EM80170"    ,                          "",           0 );*/
    // Analyse.AddMCSample( "WJets"          ,                    "W+jets",     kBlue+1 ); //W+jets
    Analyse.AddMCSample( "VGam"           ,                   "V#gamma",   kGray+2);
    // //Analyse.AddMCSample( "ZJets"          ,    "Z #rightarrow ll + X",   kGray+2 );
     Analyse.AddMCSample( "Z_2t"           ,"Z #rightarrow #tau#tau",   kGray+1 );
     if(UsePdgId) {
       if(pdgID==11)
     	Analyse.AddMCSample( "Z_2e"       ,"Z #rightarrow ee",   kGray+1 );
       else
     	Analyse.AddMCSample( "Z_2m"       ,"Z #rightarrow #mu#mu",  kGray+1 );
     }
     else {
       Analyse.AddMCSample( "Z_2e"         ,"Z #rightarrow ee,#mu#mu",    kGray );
       Analyse.AddMCSample( "Z_2m"         ,                          "",          0 );
      }
    
    
  }

  //Lines *********************************** Lines
  vector<int> Lines;
  
  //*********************************************************************
  //Execution macro ******************************************************
  
  Analyse.ConfigureData(EventFilter,EventNumber,MCOnly);
  Analyse.ConfigureLumi(Lumis, KFactors,lumi);
  Analyse.ConfigureEndAnalysis(ForceCut, WPIso, WPID, IdCuts,IsoCuts,PTcut,METcut,
			       "",MassCut,NVertex,Norm);
  
  Analyse.ConfigurePlots(Binning,AddBinBkg,OverFlowBin, UnderFlowBin,RangeY,RangeX,YTitle,
			 Xdiv,Ydiv,MarkerSize,LineWidth,Lines,SevObs,Draw3on1,SuperColor,
			 switchRMS,errorOpt,FillProfile,BasicProf,ShowDMCRatio);
  Analyse.SavePlot(savePlot);

  if(Recompute) {
    Analyse.AddVariables(UsePdgId, pdgID, MaxJet, MaxLepton, dPhiJMCut, 
			 dPhiZMCut, balMin, balMax, bTag, 
			 advMassMin, advMassMax, useAccForEff, typeJet );

    Analyse.EndAnalysis(data,observable,isEndcaps,ConvRejection,NoBackground,InvCut,
			N1Plot,Nm1Var,noDeltaCor,QCDType,ZType, repository);

    Analyse.SetAdditionnalTrees( "MCTruth", MCAddTrees );
    Analyse.SetAdditionnalTrees( "Leptons", LepAddTrees );

  Recompute=false;
  Analyse.FillZZTree();
  }
 
  Analyse.PlotDistribution(observable,Fill2DHistos,FillProfile);

}
