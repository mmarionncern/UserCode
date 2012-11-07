 {

   if(Recompute)
     THEANALYSIS Analyse;
  

   //paramètres généraux ********************* paramètres généraux
   string repository="THEPATH"; //MET/53Xv1 //root://eoscms//eos/cms/store/user/mmarionn/Leptoquarks/LQ2012/

   vector<string> data;
  
   Analyse.SetTreeName("rootTupleTree/tree");
 
   data.push_back("TauPlusX_A_190456_193621");
   data.push_back("TauPlusX_A_r_190782_190949");
   data.push_back("TauPlusX_B_190456_196531");
   data.push_back("TauPlusX_C_1_197770_198913");
   data.push_back("TauPlusX_C_2_198934_201678");
   data.push_back("TauPlusX_C_2_201679_202016");
   data.push_back("TauPlusX_C_2_202017_202305");
   data.push_back("TauPlusX_C_2_202306_203002");
   
   //skimming procedure
   bool skim=true;

   bool MCOnly = true;

   bool EventFilter=false; //should be called run filter and run number...
   int EventNumber=193621;

   string QCDType=""; //??? obsolete
   string ZType="Powheg"; //??? obsolete

   //Variable à étudier ********************** Variable à étudier
   string observable="mass";
  
   bool Draw3on1 =false;
   vector<string> SevObs;
   SevObs.push_back("mass");
   SevObs.push_back("Tmass");
   SevObs.push_back("qT");

   bool Fill2DHistos=false;
   bool FillProfile=false;
   bool response=false;

   //Binning & titre ************************* Binning & titre
   string YTitle="number of events";
   int Binning=4;
   int AddBinBkg=2; //BinB = binning*AddBin
   double RangeY[2]={1,1000000};
   double RangeX[2]={0.,200.};
   bool logYscale =true;
   int Xdiv[3]={6,5,0};
   int Ydiv[3]={5,5,0}; //Nlabel /  sous-Div /ssdiv
   bool SuperColor=false;
   bool OverFlowBin=false;
   bool UnderFlowBin=false;
   bool ShowDMCRatio=true;
   bool AddSystematics=true;
   float MarkerSize=0.9;
   float LineWidth=2;

   bool BasicProf=true;
   bool switchRMS=true;
   string errorOpt="s";

   bool savePlot = false;

   //Paramètres de l'analyse ****************** Paramètres de l'analyse
   bool N1Plot=false;
   bool InvCut=false;
   string Nm1Var = "HLT";

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


   // Normalization type

   string Norm="";
   float lumi=10259; //pb-1  //5049 //5295
 
   //paramètres additionnels

   //Lumis( or XSections ) pb-1 & KFactors ************************************
   bool useXS=true;

   map<string,float> LumisXS;

   //via XSect
   LumisXS [ "ZJets_2l" ] = 3503.71;
   LumisXS [ "WJets_ln" ] = 36257.2;
   LumisXS [ "TTbarJets" ] = 234;
   LumisXS [ "T_tWchan" ] = 11.1;
   LumisXS [ "T_tchan" ] = 30.7;
   LumisXS [ "T_schan" ] = 1.76;
   LumisXS [ "Tbar_tWchan" ] = 11.1;
   LumisXS [ "Tbar_tchan" ] = 56.4;
   LumisXS [ "Tbar_schan" ] = 3.79;
   LumisXS [ "ZZ_4l" ] = 0.0843;
   LumisXS [ "WZ_3ln" ] = 1.05;
   LumisXS [ "WWJets_2l2n" ] = 5.998;
   LumisXS [ "WW" ] = 55.47;
   LumisXS [ "ZZ_2l2n" ] = 0.33;
   LumisXS [ "ZZ" ] = 8.27;
   LumisXS [ "WZ" ] = 33.59;
   
   //via lumi
   // LumisXS [ "ZJets_2l_1" ] = 8459.9; //lumi
   // LumisXS [ "ZJets_2l_2" ] = 8459.9;
   // LumisXS [ "TTbarJets_1" ] = 25364.06;
   // LumisXS [ "WZ_3ln" ] = ;
   // LumisXS [ "WW_2ln" ] = ;
   

  map<string,float> KFactors;
  KFactors [ "ZJets_2l" ] = 1.;

  KFactors [ "TTbarJets" ] = 1.00;
  KFactors [ "T_tWchan" ] = 1.;
  KFactors [ "T_tchan" ] = 1.;
  KFactors [ "Tbar_tWchan" ] = 1.;
  KFactors [ "Tbar_tchan" ] = 1.;
  KFactors [ "ZZ_4l" ] = 1.;
  KFactors [ "WZ_3ln" ] = 1.;
  KFactors [ "WWJets_2l2n" ] = 1.;
  KFactors [ "ZZ_2l2n" ] = 1.;
  KFactors [ "ZZ" ] = 1.;
  KFactors [ "WZ" ] = 1.;

 

 
  //MonteCarlo Samples ************************** MC samples
  if(Recompute) {
   //  Analyse.AddMCSample( "T_tWchan"               ,             "single-t", kMagenta+1 );
//     Analyse.AddMCSample( "Tbar_tWchan"            ,                     "",          0 );
    Analyse.AddMCSample( "MAINDS"                ,                     "test",          0 );
 //    Analyse.AddMCSample( "Tbar_tchan"             ,                     "",          0 );
 //    Analyse.AddMCSample( "T_schan"                ,                     "",          0 );
//     Analyse.AddMCSample( "Tbar_schan"             ,                     "",          0 );
//     Analyse.AddMCSample( "TTbarJets"              ,             "t#bar{t}",  kViolet+1 );
//     Analyse.AddMCSample( "WZ"                     ,                   "VV",     kRed+1 );
//     Analyse.AddMCSample( "WW"                     ,                     "",          0 );
//     Analyse.AddMCSample( "ZZ"                     ,                     "",          0 );
//     Analyse.AddMCSample( "WJets_ln"               ,               "W+jets",          0 );
//     Analyse.AddMCSample( "ZJets_2l"               ,     "Z #rightarrow ll",  kOrange-3 );
    
  }

  //Lines *********************************** Lines
  vector<int> Lines;
  
  //*********************************************************************
  //Execution macro ******************************************************
  
  Analyse.ConfigureData(EventFilter,EventNumber,MCOnly);
  Analyse.ConfigureLumi(LumisXS, KFactors,lumi,useXS);
  Analyse.ConfigureEndAnalysis(ForceCut, WPIso, WPID, IdCuts,IsoCuts,PTcut,METcut,
			       "",MassCut,NVertex,Norm);
  
  Analyse.ConfigurePlots(Binning,AddBinBkg,OverFlowBin, UnderFlowBin,RangeY,RangeX,YTitle,
			 Xdiv,Ydiv,MarkerSize,LineWidth,Lines,SevObs,Draw3on1,SuperColor,
			 switchRMS,errorOpt,FillProfile,BasicProf,ShowDMCRatio,logYscale,AddSystematics);
  Analyse.SavePlot(savePlot);

  if(Recompute) {
    Analyse.AddVariables(skim);
  
    Analyse.EndAnalysis(data,observable,isEndcaps,ConvRejection,NoBackground,InvCut,
			N1Plot,Nm1Var,noDeltaCor,QCDType,ZType, repository);
    
    //Bloody René Brun
    bool * rtmp= const_cast<bool*> pr;
    *rtmp = false;
    Analyse.FillMETTree();
  }
 
  //  if(skimming)
  // gProcessLine->(".q");

  Analyse.PlotDistribution(observable,Fill2DHistos,FillProfile);

 }
