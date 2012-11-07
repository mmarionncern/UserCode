 {

   if(Recompute)
     LeptoQuarkUseTree Analyse;
  

   //paramètres généraux ********************* paramètres généraux
   string repository="LeptoQuark";

   vector<string> data;
  
   Analyse.SetTreeName("tree");
   data.push_back("data");
   // data.push_back("RootNtuple-DATA-TauX-PromptRecoV4_20120214_233210"); //861.6
   // data.push_back("RootNtuple-DATA-TauX-May10_20120202_115736"); //163.9 
   // data.push_back("RootNtuple-DATA-TauX-2011B_20120215_080824"); //2486
   // data.push_back("RootNtuple-DATA-TauX-05Aug2011_20120202_115952"); //336.6
   // data.push_back("RootNtuple-DATA-TauX-03Oct2011_20120220_110223"); //649.3

   //skimming procedure
   bool skim=true;

   bool MCOnly = true;

   bool EventFilter=false;
   int EventNumber=139466;

   string QCDType="";
   string ZType="Powheg";

   //Variable à étudier ********************** Variable à étudier
   string observable="ST_TauIso";
  
   bool Draw3on1 =false;
   vector<string> SevObs;
   SevObs.push_back("TLMass_Presel");
   SevObs.push_back("TLCharge_PreselElec");
   SevObs.push_back("TLCharge_PreselMu");

   bool Fill2DHistos=false;
   bool FillProfile=false;
   bool response=false;

   //Binning & titre ************************* Binning & titre
   string YTitle="number of events";
   int Binning=5;
   int AddBinBkg=1; //BinB = binning*AddBin
   double RangeY[2]={0.8,50};
   double RangeX[2]={-180,1000};
   int Xdiv[3]={6,5,0};
   int Ydiv[3]={5,5,0}; //Nlabel /  sous-Div /ssdiv
   bool SuperColor=false;
   bool OverFlowBin=false;
   bool UnderFlowBin=false;
   bool ShowDMCRatio=false;
   float MarkerSize=1.1;
   float LineWidth=2;

   bool BasicProf=true;
   bool switchRMS=false;
   string errorOpt="s";

   bool savePlot = false;
   //Paramètres de l'analyse ****************** Paramètres de l'analyse
   bool N1Plot=true;
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
   string Norm="";

   float lumi=4497.4; //pb-1
 
   //paramètres additionnels
   // bool UsePdgId = true; //remove for WW cross check
   // int pdgID = 13;
   float ptlep = 30;
   float pttau = 50;
   int chr = -1;
   int nbtag=2;
  

   bool useAccForEff=false;
   //Lumis( or XSections ) pb-1 & KFactors ************************************
   bool useXS=false;

   map<string,float> LumisXS;
   // LumisXS [ "RootNtuple-MC-DY_20120203_192111" ] = 2321;
   // LumisXS [ "RootNtuple-MC-QCD_20120203_192035" ] = 296.9e6;
   // LumisXS [ "RootNtuple-MC-QCD-EL_20120203_191710" ] = 0;
   // LumisXS [ "RootNtuple-MC-TT_20120220_213535" ] = 164.4;
   // LumisXS [ "RootNtuple-MC-WJETS_20120220_214711" ] = 24640;
   // LumisXS [ "WW" ] = 47;
   // LumisXS [ "WZ" ] = 18.2;
   // LumisXS [ "ZZ" ] = 7.41;
   // LumisXS [ "Mass350" ] = 0.5;
   // LumisXS [ "Mass450" ] = 0.5;

  // LumisXS [ "skim_TTbar" ] = 0.5;
  
  //
 
  

  map<string,float> KFactors;
  // KFactors [ "RootNtuple-MC-DY_20120203_192111" ] = 1.;
  // KFactors [ "RootNtuple-MC-QCD_20120203_192035" ] = 1.;
  // KFactors [ "RootNtuple-MC-QCD-EL_20120203_191710" ] = 1.;
  // KFactors [ "RootNtuple-MC-TT_20120220_213535" ] = 1.;
  // KFactors [ "RootNtuple-MC-WJETS_20120220_214711" ] = 1.;
  // KFactors [ "WW" ] = 1.;
  // KFactors [ "WZ" ] = 1.;
  // KFactors [ "ZZ" ] = 1.;
  // KFactors [ "Mass350" ] = 1.;
  // KFactors [ "Mass450" ] = 1.;



 
  //MonteCarlo Samples ************************** MC samples
  if(Recompute) {
      Analyse.AddMCSample( "ZZ"                                  ,       "VV",    kRed+1 );
      Analyse.AddMCSample( "WW"                                  ,         "",         0 );
      Analyse.AddMCSample( "WZ"                                  ,         "",         0 );
    //  // Analyse.AddMCSample( "RootNtuple-MC-TT_20120220_213535"    , "t#bar{t}", kViolet+1 );
    //  // Analyse.AddMCSample( "RootNtuple-MC-WJETS_20120220_214711" ,   "W+jets", kOrange+7 ); 
    //  // Analyse.AddMCSample( "RootNtuple-MC-DY_20120203_192111"    ,       "DY", kOrange-2 );
    // // Analyse.AddMCSample( "RootNtuple-MC-QCD_20120203_192035"   ,      "QCD",   kBlue+1 ); 
    // // Analyse.AddMCSample( "RootNtuple-MC-QCD-EL_20120203_191710"   ,         "",         0 ); 
    // Analyse.AddMCSample( "skim_TTbar"                                  ,       "VV",    kRed+1 );
    // Analyse.AddMCSample( "skim_WJets"                                  ,       "cdc",    kBlue+1 );
    // Analyse.AddMCSample( "Mass350"                             ,   "LQ 350",   kGray+1 );
    // Analyse.AddMCSample( "Mass450"                             ,   "LQ 450",   kGray+3 );

    // Analyse.AddMCSample( "LQ_Mass550"     ,   "LQ 550",   kBlue+1 );
    // Analyse.AddMCSample( "skim_ZZ"        ,       "VV",    kRed+1 );
    // Analyse.AddMCSample( "skim_WW"        ,         "",         0 );
    // Analyse.AddMCSample( "skim_WZ"        ,         "",         0 );
    // Analyse.AddMCSample( "skim_TTbar"     , "t#bar{t}", kViolet+1 );
    // Analyse.AddMCSample( "skim_WJets"     ,   "W+jets", kOrange+7 );
    // Analyse.AddMCSample( "skim_DYJets"    ,       "DY", kOrange-2 ); 
    //    skim_QCD_Pt_20to30_BCtoE
    // skim_QCD_Pt_20to30_EMEnriched
    // skim_QCD_Pt_30to80_BCtoE
    // skim_QCD_Pt_30to80_EMEnriched
    // skim_QCD_Pt_80to170_BCtoE
  
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
			 switchRMS,errorOpt,FillProfile,BasicProf,ShowDMCRatio);
  Analyse.SavePlot(savePlot);

  if(Recompute) {
    Analyse.AddVariables(/*UsePdgId, pdgID*/ nbtag, ptlep, pttau, chr, skim );
  
    Analyse.EndAnalysis(data,observable,isEndcaps,ConvRejection,NoBackground,InvCut,
			N1Plot,Nm1Var,noDeltaCor,QCDType,ZType, repository);
    
    //Bloody René Brun
    bool * rtmp= const_cast<bool*> pr;
    *rtmp = false;
    Analyse.FillLQTree();
  }
 
  Analyse.PlotDistribution(observable,Fill2DHistos,FillProfile);

 }
