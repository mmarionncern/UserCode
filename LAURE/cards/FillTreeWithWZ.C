{

  if(Recompute)
    WZUseTree Analyse;

  //paramètres généraux ********************* paramètres généraux
  string analyse="OZee";
  vector<string> data;
  
  data.push_back("EG_Nov4_132440_141000");
  data.push_back("EG_Nov4_141001_143000");
  data.push_back("EG_Nov4_143001_144114");
  data.push_back("Electron_Nov4_146240_147000");
  data.push_back("Electron_Nov4_147001_148000");
  data.push_back("Electron_Nov4_148001_149000");
  data.push_back("Electron_Nov4_149001_149442");


  bool MCOnly = false;

  bool EventFilter=false;
  int EventNumber=139466;

  string QCDType="";
  string ZType="Powheg";

  //Variable à étudier ********************** Variable à étudier
  string observable="LeptonMult";
  bool N1Plot=false;
 
  bool Draw3on1 =true;
  vector<string> SevObs;
  SevObs.push_back("LElecEta");
  SevObs.push_back("dPhiMETL1");
  SevObs.push_back("dPhiMETL2");

  bool Fill2DHistos=false;
  bool FillProfile=false;
  bool response=false;

  //Binning & titre ************************* Binning & titre
  string YTitle="nevts  ";
  int Binning=2;
  double RangeY[2]={0.00001,50};
  double RangeX[2]={-180,300};
  int Xdiv[3]={6,5,0};
  int Ydiv[3]={5,5,0}; //Nlabel /  sous-Div /ssdiv
  bool SuperColor=false;
  float MarkerSize=0.8;
  float LineWidth=2;

  bool BasicProf=true;
  bool switchRMS=true;
  string errorOpt="s";

  //Paramètres de l'analyse ****************** Paramètres de l'analyse

  int isEndcaps=2; //1=EE, 0=EB , 2=combined

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
  float METcut=00;
  float MTcut=60;
  bool inversionMETCut = false;
  bool noDeltaCor=false;
  string Norm="";


  float lumi=36; //pb-1
 
  //Lumis pb-1 & KFactors ************************************
  map<string,float> Lumis;
  Lumis[ "WW" ] = 4392;
  Lumis[ "WW_2l" ] = 75460;
  Lumis[ "WZ" ] = 11250;
  Lumis[ "WZ_3l" ] = 318700;
  Lumis[ "ZZ" ] = 33806;
  Lumis[ "ZZ_2l2n" ] = 570200;
  Lumis[ "ZZ_4l" ] = 1700000;
  Lumis[ "WGam" ] = 4614;
  Lumis[ "ZGam" ] = 13760;
  Lumis[ "ZJets" ] = 452*3;
  Lumis[ "Z_2t_powheg" ] = 100/0.9;
  Lumis[ "Z_2t" ] = 1683;
  Lumis[ "Z_2e_powheg" ] = 1058*3;
  Lumis[ "Z_2e" ] = 1956*3;
  Lumis[ "Z_2m_powheg" ] = 1084;
  Lumis[ "Z_2m" ] = 10000;
  Lumis[ "ttbar" ] = 6700;

  map<string,float> KFactors;
  KFactors[ "Z_2t" ] = 1.3;
  KFactors[ "Z_2e" ] = 1.3;
  KFactors[ "Z_2m" ] = 1.3;
  KFactors[ "WW" ] = 1.53;
  KFactors[ "WW_2l" ] = 1.53;
  KFactors[ "WZ" ] = 1.71;
  KFactors[ "WZ_3l" ] = 1.71;
  KFactors[ "ZZ" ] = 1.37;
  KFactors[ "ZZ_4l" ] = 1.37;
  KFactors[ "ZZ_2l2n" ] = 1.37;
  KFactors[ "ttbar" ] = 1.75;
 
  //Lines *********************************** Lines
  vector<int> Lines;
  /*Lines.push_back(70);
    Lines.push_back(80);
    Lines.push_back(90);
    Lines.push_back(95);*/




  // string XTitle="pfMET";

  //*********************************************************************
  //Execution macro ******************************************************
 
  Analyse.ConfigureData(EventFilter,EventNumber,MCOnly);
  Analyse.ConfigureLumi(Lumis, KFactors,lumi);
  Analyse.ConfigureEndAnalysis(ForceCut, WPIso, WPID, IdCuts,IsoCuts,PTcut,METcut,inversionMETCut,METType,MTcut,NVertex,Norm);
  Analyse.ConfigurePlots(Binning,RangeY,RangeX, YTitle,Xdiv,Ydiv,MarkerSize,LineWidth,Lines,SevObs,Draw3on1,SuperColor,switchRMS,errorOpt,FillProfile,BasicProf);
  
  if(Recompute) {
  Analyse.EndAnalysis(data,observable,isEndcaps,ConvRejection,NoBackground,
		      N1Plot,noDeltaCor,QCDType,ZType);
  Recompute=false;
  Analyse.FillWZTree();
   }
  
  Analyse.PlotDistribution(observable,Fill2DHistos,FillProfile);

}
