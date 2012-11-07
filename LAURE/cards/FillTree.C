{

 if(Recompute)
   WUseTree Analyse;

  //paramètres généraux ********************* paramètres généraux
  string analyse="OWe";
  vector<string> data;
  //data.push_back("FullLumi");
  /*   data.push_back("EG_Nov4_132440_141000");
  data.push_back("EG_Nov4_141001_143000");
  data.push_back("EG_Nov4_143001_144114");
  data.push_back("Electron_Nov4_146240_147000");
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
  
  bool MCOnly = false;

  bool EventFilter=false;
  int EventNumber=139466;

 string QCDType="EM";
 string WType="Pythia_PU";

 //Variable à étudier ********************** Variable à étudier
 string observable="Eta";
 
 bool Draw3on1 =false;
 vector<string> SevObs;
 SevObs.push_back("tcdPhiRecoil");
 SevObs.push_back("pfdPhiRecoil");
 SevObs.push_back("pfT1dPhiRecoil");

 bool Fill2DHistos=false;
 bool FillProfile=false;

 //Binning & titre ************************* Binning & titre
 //mean #Delta#Phi(U_{T}, lepton) [rad]
 string YTitle="number of events / 4 Gev";
 //#Delta(U_{T},lepton) [rad]
 int Binning=1;
 int AddBinBkg=2; //BinB = binning*AddBin
 double RangeY[2]={0.1,2000};
 double RangeX[2]={-156,48};
 int Xdiv[3]={5,6,0};
 int Ydiv[3]={5,6,0}; //Nlabel /  sous-Div /ssdiv
 bool SuperColor=false;
 bool OverFlowBin=true;
 bool UnderFlowBin=true;
 bool ShowDMCRatio=false;
 float MarkerSize=0.8;
 float LineWidth=1;

 bool BasicProf=false;
 bool switchRMS=false;
 string errorOpt="";

 //Paramètres de l'analyse ****************** Paramètres de l'analyse
 bool N1Plot=true;
 string Nm1Var = "sigieie";

 int isEndcaps=2; //1=EE, 0=EB , 2=combined

 bool ConvRejection=true;
 bool VetoE=true;
 int WPIso = 80;
 int WPID = 80;

 bool ForceCut=false;

//                     sieie  deta  dphi  hoe
 float IdCuts[2][4] ={{ -1,   -1,   -1,   -1 }, //EB
		      { -1,   1000,  -1,   -1 }}; //EE
//                       trk  ecal  hcal
float IsoCuts[2][3] ={{ -1,    1000,   -1  },   //EB
		      { -1,    1000,   -1  } }; //EE
 
 int NVertex=1;
 float PTcut=25;
 string METType="pfMET";
 float METcut=0;
 float MTcut=0;
 bool inversionMETCut = false;
 bool noDeltaCor=false;

 string Norm="Fit";

 float lumi=36; //pb-1
 
 //Lumis pb-1 & KFactors ************************************
  map<string,float> Lumis;
  Lumis[ "GamJet_PT20_TuneZ2" ] = 2366.4; //filter 0.0064
  Lumis[ "GamJet_PT1530_TuneZ2" ] = 0.88;
  Lumis[ "GamJet_PT3050_TuneZ2" ] = 11.97;
  Lumis[ "GamJet_PT5080_TuneZ2" ] = 66.67;
  Lumis[ "GamJet_PT80120_TuneZ2" ] = 190.8;
  /* Lumis[ "QCD_EM2030_TuneZ2" ] = 21.61;
  Lumis[ "QCD_EM3080_TuneZ2" ] = 11.43;
  Lumis[ "QCD_EM80170_TuneZ2" ] = 50.2;
  Lumis[ "QCD_bc2030_TuneZ2" ] = 20.3;
  Lumis[ "QCD_bc3080_TuneZ2" ] = 14.4;
  Lumis[ "QCD_bc80170_TuneZ2" ] = 111;
  // Lumis[ "Z_2t_TuneZ2" ] = 1579;
  Lumis[ "Z_2e_powheg_TuneZ2" ] = 9;
  Lumis[ "W_en_TuneZ2" ] = 377;*/
  //  Lumis[ "W_en_PU" ] = 816.2;
  Lumis[ "W_en_minus_powheg_PU_TuneZ2" ] = 358.7;
  Lumis[ "W_en_plus_powheg_PU_TuneZ2" ] = 525.7;
  Lumis[ "W_en_minus_powheg_TuneZ2" ] = 358.7;
  Lumis[ "W_en_plus_powheg_TuneZ2" ] = 525.7;
  Lumis[ "W_tn_TuneZ2" ] = 325.2;
  Lumis[ "Z_2e_TuneZ2" ] = 1633;
  Lumis[ "Z_2e_PU_TuneZ2" ] = 1600;
  // Lumis[ "ttbar" ] = 6700;
  
  Lumis[ "Z_2t_TuneZ2" ] = 77;
  Lumis[ "W_tn_TuneZ2" ] = 50.6;
  Lumis[ "Z_2e_PU_TuneZ2" ] = 77;
  Lumis[ "W_en_TuneZ2" ] = 32.5;
  //  Lumis[ "W_en_PU_TuneZ2_v2" ] = 276.3;
  Lumis[ "QCD_EM3080_TuneZ2_v2" ] = 3.9;
  Lumis[ "QCD_EM80170_TuneZ2_v2" ] = 43.95;
  Lumis[ "QCD_bc2030_TuneZ2" ] = 20.3;
  Lumis[ "QCD_bc3080_TuneZ2" ] = 14.4;
  Lumis[ "QCD_bc80170_TuneZ2" ] = 111;
 
  Lumis[ "ttbar" ] = 10600;
  Lumis[ "GamJet_PT1530" ] = 8.94;
  Lumis[ "GamJet_PT3050" ] = 621.5;
  Lumis[ "GamJet_PT5080" ] = 376.7;
  Lumis[ "GamJet_PT80120" ] = 2349.2;
  Lumis[ "QCD_EM2030" ] = 13.55;
  Lumis[ "QCD_EM3080" ] = 17.45;
  Lumis[ "QCD_EM80170" ] = 60.2;
  Lumis[ "QCD_bc2030" ] = 20.8;
  Lumis[ "QCD_bc3080" ] = 8.1;
  Lumis[ "QCD_bc80170" ] = 33.5;
  Lumis[ "Z_2e_PU" ] = 1604;
  Lumis[ "W_tn_PU" ] = 627.1;
  Lumis[ "W_en_PU" ] = 816.2;
  
  
  map<string,float> KFactors;
  KFactors[ "ttbar" ] = 2;
  KFactors[ "Z_2t" ] = 1.3;
  KFactors[ "Z_2e" ] = 1.3;
  KFactors[ "W_tn" ] = 1.3;
  KFactors[ "W_en" ] = 1.3;
  KFactors[ "W_en_PU" ] = 1.3;
  KFactors[ "W_en_PU_TuneZ2_v2" ] = 1.3;
  KFactors[ "W_tn_PU_TuneZ2" ] = 1.3;
  KFactors[ "Z_2e_PU_TuneZ2" ] = 1.3;
  KFactors[ "Z_t2_TuneZ2" ] = 1.3;
  
  
  //MonteCarlo Samples ************************** MC samples
  if(Recompute) {
    Analyse.AddMCSample( "ttbar" ,                     "t#bar{t}",      kRed+1 );
    Analyse.AddMCSample( "GamJet_PT1530" ,           "#gamma+jet",  kMagenta+4 );
    Analyse.AddMCSample( "GamJet_PT3050" ,                     "",           0 );
    Analyse.AddMCSample( "GamJet_PT5080" ,                     "",           0 );
    Analyse.AddMCSample( "GamJet_PT80120" ,                    "",           0 );
    Analyse.AddMCSample( "QCD_bc2030" ,                     "QCD",   kViolet-5 );
    Analyse.AddMCSample( "QCD_bc3080" ,                        "",          0  );
    Analyse.AddMCSample( "QCD_bc80170" ,                       "",          0  );
    Analyse.AddMCSample( "QCD_EM2030" ,                        "",          0  );
    Analyse.AddMCSample( "QCD_EM3080" ,                        "",          0  );
    Analyse.AddMCSample( "QCD_EM80170" ,                       "",          0  );
    Analyse.AddMCSample( "W_tn_PU" ,                        "EWK",   kOrange+7 );
    Analyse.AddMCSample( "Z_2e_PU" ,                           "",           0 );
    Analyse.AddMCSample( "W_en_PU" ,         "W #rightarrow e#nu",   kOrange-2 );
  }

  //Lines *********************************** Lines
  vector<int> Lines;



 //*********************************************************************
 //Execution macro ******************************************************
 

 
 Analyse.ConfigureData(EventFilter,EventNumber,MCOnly);
 Analyse.ConfigureLumi(Lumis, KFactors, lumi);
 Analyse.ConfigureEndAnalysis(ForceCut, WPIso, WPID, IdCuts,IsoCuts,PTcut,METcut,inversionMETCut,METType,MTcut,NVertex,Norm);
 Analyse.ConfigurePlots(Binning,AddBinBkg,OverFlowBin, UnderFlowBin,RangeY,RangeX, YTitle,Xdiv,Ydiv,MarkerSize,LineWidth,Lines,SevObs,Draw3on1,SuperColor,switchRMS,errorOpt,FillProfile,BasicProf,ShowDMCRatio);
 if(Recompute) {
   Analyse.EndAnalysis(data,observable,isEndcaps,ConvRejection,VetoE,
		       N1Plot, Nm1Var, noDeltaCor,QCDType,WType);
   Recompute=false;

   Analyse.FillWTree();
 }
 
 Analyse.PlotDistribution(observable,Fill2DHistos,FillProfile);
}
