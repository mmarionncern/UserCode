{

 if(Recompute)
   WUseTree Analyse;

  //paramètres généraux ********************* paramètres généraux
  string repository="Wen";
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
 string observable="pfT1MET";
 
 bool Draw3on1 =false;
 vector<string> SevObs;
 SevObs.push_back("tcdPhiRecoil");
 SevObs.push_back("pfdPhiRecoil");
 SevObs.push_back("pfT1dPhiRecoil");

 bool Fill2DHistos=false;
 bool FillProfile=false;

 //Binning & titre ************************* Binning & titre
 //mean #Delta#Phi(U_{T}, lepton) [rad]
 string YTitle="number of events / 2 GeV";
 //#Delta(U_{T},lepton) [rad]
 int Binning=4;
 int AddBinBkg=1; //BinB = binning*AddBin
 double RangeY[2]={0.001,6000};
 double RangeX[2]={0,120};
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
 bool N1Plot=false;
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

 string Norm="Remove2V";

 float lumi=36; //pb-1
 
 //Lumis pb-1 & KFactors ************************************
  map<string,float> Lumis;
  
  Lumis[ "GamJet_PT1530" ] = 5.97 ;
  Lumis[ "GamJet_PT3050" ] = 61.41 ;
  Lumis[ "GamJet_PT5080" ] = 376.69 ;
  Lumis[ "GamJet_PT80120" ] = 2343.95 ;
  Lumis[ "GamJet_PT120170" ] = 11947.74 ;
  Lumis[ "QCD_EM2030" ] = 15.04 ;///2./1.5;
  Lumis[ "QCD_EM3080" ] = 18.55 ;
  Lumis[ "QCD_EM3080_small" ] = 18.55;///2.;
  Lumis[ "QCD_EM80170" ] = 57.87 ;///2./1.5;
  Lumis[ "QCD_bc2030" ] = 16.98 ;///2./1.5;
  Lumis[ "QCD_bc3080" ] = 14.58 ;///2./1.5;
  Lumis[ "QCD_bc80170" ] = 111.47;///1.5;
  Lumis[ "ttbar" ] = 9634.02;
  Lumis[ "WJets" ] = 615.05;
  Lumis[ "W_en_Fall10" ] = 430000./6153.;
  Lumis[ "W_en_PU_TuneZ2" ] = 2117658./6153;
  Lumis[ "W_tn_Fall10" ] = 430000./7899.;
  Lumis[ "WJets_1" ] = 615.05/3.55;
  Lumis[ "WJets_2" ] = 615.05/3.;
  Lumis[ "WJets_3" ] = 615.05/3.;
  Lumis[ "Z_2t" ] = 860.04;
  Lumis[ "Z_2m" ] = 1238.49;
  Lumis[ "Z_2e" ] = 1238.53;

  map<string,float> KFactors;
  KFactors[ "ttbar" ] = 1.3;
  KFactors[ "Z_2t" ] = 1.03;
  KFactors[ "Z_2e" ] = 1.03;
  KFactors[ "W_tn" ] = 1.3;
  KFactors[ "W_en" ] = 1.3;
  KFactors[ "W_en_PU" ] = 1.3;
  KFactors[ "W_en_PU_TuneZ2_v2" ] = 1.3;
  KFactors[ "W_tn_PU_TuneZ2" ] = 1.3;
  KFactors[ "Z_2e_PU_TuneZ2" ] = 1.3;
  KFactors[ "Z_t2_TuneZ2" ] = 1.3;
  //KFactors[ "WJets_1" ] = 0.4;
  //KFactors[ "WJets_2" ] = 0.4;
  KFactors[ "W_en_Fall10" ] = 1.3;
  KFactors[ "W_en_PU_TuneZ2" ] = 1.3*0.98;
  KFactors[ "W_tn_Fall10" ] = 1.3;
  
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
    Analyse.AddMCSample( "QCD_EM3080_small" ,                  "",          0  );
    Analyse.AddMCSample( "QCD_EM80170" ,                       "",          0  );
    Analyse.AddMCSample( "Z_2t" ,                           "EWK",   kOrange+7 );
    Analyse.AddMCSample( "Z_2e" ,                              "",           0 );
    Analyse.AddMCSample( "W_tn_Fall10" ,                       "",           0 );
    Analyse.AddMCSample( "W_en_PU_TuneZ2",      "W #rightarrow e#nu",   kOrange-2 );
    //  Analyse.AddMCSample( "WJets_2",            "",   0 );
    //  Analyse.AddMCSample( "WJets_3",            "",   0 );
  }

  //Lines *********************************** Lines
  vector<int> Lines;



 //*********************************************************************
 //Execution macro ******************************************************
 

 
 Analyse.ConfigureData(EventFilter,EventNumber,MCOnly);
 Analyse.ConfigureLumi(Lumis, KFactors, lumi);
 Analyse.ConfigureEndAnalysis(ForceCut, WPIso, WPID, IdCuts,IsoCuts,PTcut,
			      METcut,METType,MTcut,NVertex,Norm);
 Analyse.ConfigurePlots(Binning,AddBinBkg,OverFlowBin, UnderFlowBin,RangeY,RangeX, YTitle,Xdiv,Ydiv,MarkerSize,LineWidth,Lines,SevObs,Draw3on1,SuperColor,switchRMS,errorOpt,FillProfile,BasicProf,ShowDMCRatio);
 if(Recompute) {
   Analyse.EndAnalysis(data,observable,isEndcaps,ConvRejection,VetoE, false,
		       N1Plot, Nm1Var, noDeltaCor,QCDType,WType, repository);
   Recompute=false;

   Analyse.FillWTree();
 }
 
 Analyse.PlotDistribution(observable,Fill2DHistos,FillProfile);
}
