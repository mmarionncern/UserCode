{

  if(Recompute)
    PFCor Analyse;

  //paramètres généraux ********************* paramètres généraux
  // string analyse="OZee";
  vector<string> data;
  //  map<string, vector<int> > MC;


  data.push_back("Electron_Dec22_146240_147000");
  data.push_back("Electron_Dec22_147001_148000");
  data.push_back("Electron_Dec22_148001_149000");
  data.push_back("Electron_Dec22_149001_149442");
  /*  data.push_back("Electron_Dec22_147001_148000");
  data.push_back("Electron_Dec22_148001_149000");
  data.push_back("Electron_Dec22_149001_149442");*/

  bool MCOnly = false;

  bool EventFilter=false;
  int EventNumber=147205;

  string QCDType="EM";
  string ZType="Powheg";

  //Variable à étudier ********************** Variable à étudier
  string observable="ZMass";
  bool N1Plot=false;
 
  bool Draw3on1 =true;
  vector<string> SevObs;
  SevObs.push_back("RawZMass");
  SevObs.push_back("ZMassCor");
  SevObs.push_back("ZMass");

  bool Fill2DHistos=false;
  bool FillProfile=false;
  bool response=true;

  //Binning & titre ************************* Binning & titre
  string YTitle="number of events / 1 GeV ";
  int Binning=2;
  int AddBinBkg=2; //BinB = binning*AddBin
  double RangeY[2]={0.1,10000};
  double RangeX[2]={70,120};
  int Xdiv[3]={7,5,0};
  int Ydiv[3]={7,5,0}; //Nlabel /  sous-Div /ssdiv 
  bool SuperColor=false;
  bool OverFlowBin=false;
  bool UnderFlowBin=false;
  float MarkerSize=0.8;
  float LineWidth=1;

  bool BasicProf=true;
  bool switchRMS=true;
  string errorOpt="";

  //Paramètres de l'analyse ****************** Paramètres de l'analyse

  int isEndcaps=0; //1=EE, 0=EB , 2=combined, 3= EBEE for Z

  bool ConvRejection=false;
  bool NoBackground=false;
  int WPIso = 95;
  int WPID = 95;

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
  float METcut=1000;
  float MTcut=70;
  bool inversionMETCut = true;
  bool noDeltaCor=false;
  string Norm="Norm";


  float lumi=36; //pb-1 
  //Lumis pb-1 & KFactors ************************************
  map<string,float> Lumis;
  Lumis[ "GamJet_PT20_TuneZ2" ] = 2366.4;
  Lumis[ "QCD_EM2030_TuneZ2" ] = 21.61;
  Lumis[ "QCD_EM3080_TuneZ2" ] = 11.43;
  Lumis[ "QCD_EM80170_TuneZ2" ] = 50.2;
  Lumis[ "QCD_bc2030_TuneZ2" ] = 20.3;
  Lumis[ "QCD_bc3080_TuneZ2" ] = 14.4;
  Lumis[ "QCD_bc80170_TuneZ2" ] = 111;
  Lumis[ "Z_2t_TuneZ2" ] = 1579;
  Lumis[ "Z_2e_powheg_TuneZ2" ] = 9;
  Lumis[ "W_en_TuneZ2" ] = 377;
  Lumis[ "Z_2e_TuneZ2" ] = 1633;
  // Lumis[ "Z_2e_PU_TuneZ2" ] = 1600; //old one
  Lumis[ "Z_2e_PU_TuneZ2" ] = 1492;
   // Lumis[ "ttbar" ] = 6700; //old one
  
  Lumis[ "GamJet_PT1530" ] = 6.03 ;
  Lumis[ "GamJet_PT3050" ] = 62.15;
  Lumis[ "GamJet_PT5080" ] = 376.7;
  Lumis[ "GamJet_PT80120" ] = 2349.2;
  Lumis[ "QCD_EM2030" ] = 13.55;
  Lumis[ "QCD_EM3080" ] = 17.45;
  Lumis[ "QCD_EM80170" ] = 60.2;
  Lumis[ "QCD_bc2030" ] = 20.8;
  Lumis[ "QCD_bc3080" ] = 8.1;
  Lumis[ "QCD_bc80170" ] = 33.5;
  Lumis[ "Z_2t_PU" ] = 1547;
  Lumis[ "W_en_PU" ] = 325;
  Lumis[ "Z_2e_PU" ] = 1604.30;
  Lumis[ "ttbar" ] = 10600;

  map<string,float> KFactors;
  KFactors[ "Z_2t_PU" ] = 1.3;
  KFactors[ "Z_2e_PU" ] = 1.3;
  KFactors[ "W_en_PU" ] = 1.3;
  KFactors[ "ttbar" ] = 2;
  //KFactors[ "W_en_PU" ] = 1.3;

  //MonteCarlo Samples ************************** MC samples
  if(Recompute) {
    /*  Analyse.AddMCSample( "ttbar" ,          0,         "t#bar{t}" ,      kRed+1 );
    Analyse.AddMCSample( "GamJet_PT1530" ,  1,        "#gamma+jet",  kMagenta+4 );
    Analyse.AddMCSample( "GamJet_PT3050" ,  1,                  "",           0 );
    Analyse.AddMCSample( "GamJet_PT5080" ,  1,                  "",           0 );
    Analyse.AddMCSample( "GamJet_PT80120" , 1,                  "",           0 );
    Analyse.AddMCSample( "QCD_bc2030" ,     2,               "QCD",   kViolet-5 );
    Analyse.AddMCSample( "QCD_bc3080" ,     2,                  "",          0  );
    Analyse.AddMCSample( "QCD_bc80170" ,    2,                  "",          0  );
    Analyse.AddMCSample( "QCD_EM2030" ,     2,                  "",          0  );
    Analyse.AddMCSample( "QCD_EM3080" ,     2,                  "",          0  );
    Analyse.AddMCSample( "QCD_EM80170" ,    2,                  "",          0  );
    Analyse.AddMCSample( "Z_2t_PU" ,        3,               "EWK",   kOrange+7 );
    Analyse.AddMCSample( "W_en_PU" ,        3,                  "",           0 );
    Analyse.AddMCSample( "Z_2e_PU" ,        4,  "Z #rightarrow ee",   kOrange-2 );*/
  }
  //Lines *********************************** Lines
  vector<int> Lines;
 
  //*********************************************************************
  //Execution macro ******************************************************
  

  Analyse.ConfigureData(EventFilter,EventNumber,MCOnly);
  Analyse.ConfigureLumi(Lumis, KFactors, lumi);
  Analyse.ConfigureEndAnalysis(ForceCut, WPIso, WPID, IdCuts,IsoCuts,PTcut,METcut,inversionMETCut,METType,MTcut,NVertex,Norm);
  Analyse.ConfigurePlots(Binning,AddBinBkg,OverFlowBin, UnderFlowBin,RangeY,RangeX, YTitle,Xdiv,Ydiv,MarkerSize,LineWidth,Lines,SevObs,Draw3on1,SuperColor,switchRMS,errorOpt,FillProfile,BasicProf);
  if(Recompute) {
    Analyse.EndAnalysis(data,observable,isEndcaps,ConvRejection,NoBackground,
			N1Plot, noDeltaCor,QCDType,ZType);
    Recompute=false;
    
    Analyse.FillTree();
  }
  Analyse.PlotDistribution(observable,Fill2DHistos,FillProfile);

  
}
