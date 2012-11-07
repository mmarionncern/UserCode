#include "PFCorUseTree.hh"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"

using namespace RooFit;
using namespace std;

ClassImp(PFCorUseTree)


PFCorUseTree::PFCorUseTree():
UseTree()
{
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
  //  float t;
  
  histoManager.ConfigureContamination(0.05);

  treenames.clear();
  treenames.push_back( "PFCorVariables");

}

void
PFCorUseTree::FillTree() {

  //Add Varaibles to be computed

  //Prepare suffix for outfiles

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
  histoManager.AddVariable("RawEtPFSC",200,0,100,"Raw SC E_{T} [GeV]","RawEtSC");
  histoManager.AddVariable("EtPFSCCor",200,0,100,"SC cor E_{T} [GeV]","EtSCCor");
  histoManager.AddVariable("EtPFSC",200,0,100," SC E_{T} [GeV]","EtSC");
    
  //Underlying event variables
  histoManager.AddVariable("dZv",200,0,100,"#Delta z [cm]","Zee_dZ_vertices");

  //Underlying variables and quality variables
  histoManager.AddVariable("ZEta",100,-3,3,"#eta Z","Zee_Eta");
  histoManager.AddVariable("ZPhi",128,-4,4,"Z phi ","Zee_phi");
  histoManager.AddVariable("ZPt",400,0,200,"q_{T} Z [GeV]","Zee_PT");
  histoManager.AddVariable("ZPtCor",400,0,200,"cor p_{T} Z [GeV]","Zee_PT");
  histoManager.AddVariable("ZMass",400,0,200,"M_{ee} [GeV]","Zee_Mass");
  histoManager.AddVariable("ZMassCor",400,0,200,"cor M_{ee} [GeV]","Zee_MassCor");
  histoManager.AddVariable("RawZMass",400,0,200,"raw M_{ee} [GeV]","Zee_RawMass");
  
  histoManager.AddProfVariable("RatioECor",100,0,400,"raw vs cor","ratio raw/cor");
  histoManager.AddProfVariable("RatioE",100,0,400,"raw vs old","ratio raw/old");
 
  //End Declaration, now preparation
  histoManager.PrepareHistos(name);
  histoManager.Prepare2DHistos(name);
  histoManager.PrepareProfiles(name);
  cout<<" End declaration, now analyse "<<endl;

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

  //sélection
  for(int i=0;i<nt+1;i++) {
    
    float ZPt;
    float ZMass;
  
    //Leptons
    float PFLepton[2][3]; //ok
    float Lepton[2][3]; 
    float PfSCEt[2];
    float PfSCE[2];
    float PfSCESE[2];
    float PfSCRawE[2];
    float PfSCGeo[2][2];
    float PfSCShape[2][2];

    int charge[2];
    float eop[2];
    float IDVar[2][4];
    float IsoVar[2][3];
    int pdgId[2];

    //Selection
    bool inHit1;
    bool expInHit1;
    bool ConvRej1;
    bool inHit2;
    bool expInHit2;
    bool ConvRej2;
 
    //Event
    int Run;
    int Event;
    int AbsEvent;
    char sampleName[30];
    char fileName[128];
    int time;

    int Nvertex;
    float dZ;

    float Weight;

    bool Sel[2]={true,true};
    int epart[2]={-1,-1};
    TLorentzVector lepton[2];
  
    //Event
    tChains[i]->SetBranchAddress("Run",&Run);
    tChains[i]->SetBranchAddress("Event",&Event);
    tChains[i]->SetBranchAddress("AbsEvent",&AbsEvent);
    tChains[i]->SetBranchAddress("sample",&sampleName); 
    tChains[i]->SetBranchAddress("fileName",&fileName); 
    tChains[i]->SetBranchAddress("timestamp",&time);


    //Selection
    tChains[i]->SetBranchAddress("InHit1",&inHit1);
    tChains[i]->SetBranchAddress("ExpInHit1",&expInHit1);
    tChains[i]->SetBranchAddress("ConvRej1",&ConvRej1);

    tChains[i]->SetBranchAddress("InHit2",&inHit2);
    tChains[i]->SetBranchAddress("ExpInHit2",&expInHit2);
    tChains[i]->SetBranchAddress("ConvRej2",&ConvRej2);

    //Z
    tChains[i]->SetBranchAddress("ZPt",&ZPt);
    tChains[i]->SetBranchAddress("ZMass",&ZMass);
      
    //DeltaZ between vertices
    tChains[i]->SetBranchAddress("DZvertices",&dZ);
    //Lepton
    tChains[i]->SetBranchAddress("IDVar",IDVar);
    tChains[i]->SetBranchAddress("IsoVar",IsoVar); 
    tChains[i]->SetBranchAddress("Lepton",Lepton);
    tChains[i]->SetBranchAddress("PFLepton",PFLepton);
	
    tChains[i]->SetBranchAddress("PfSCRawE",PfSCRawE);
    tChains[i]->SetBranchAddress("PfSCE",PfSCE);
    tChains[i]->SetBranchAddress("PfSCESE",PfSCESE);
    tChains[i]->SetBranchAddress("PfSCEt",PfSCEt);
    tChains[i]->SetBranchAddress("PfSCGeo",PfSCGeo);
    tChains[i]->SetBranchAddress("PfSCShape",PfSCShape);

    tChains[i]->SetBranchAddress("pdgId",pdgId);
    tChains[i]->SetBranchAddress("Charge",charge);
    tChains[i]->SetBranchAddress("EoP",eop);

    tChains[i]->SetBranchAddress("NVertex",&Nvertex);

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

	if(fabs(PfSCGeo[id][0])>1.479)
	  {EP=1; }
	else
	  {EP=0; }

	Sel[id]=true;

	//first pt cut
	if(PfSCEt[id] < PTcut)  Sel[id]=false;
	
	//	if(convRej && !expInHit1 && !expInHit2)   Sel[id]=false;
	
	//Met Cut
	
	//And Selections

	if(IDVar[id][0] > _IDcuts[EP][0] && Cbool(skipCut, "sigieie"==observable.substr(0,7)) )  Sel[id]=false;
	if( fabs(IDVar[id][1]) > _IDcuts[EP][1] && Cbool(skipCut, "Deta"==observable.substr(0,4)) )  Sel[id]=false;
	if( fabs(IDVar[id][2]) > _IDcuts[EP][2] && Cbool(skipCut, "Dphi"==observable.substr(0,4)) )  Sel[id]=false;
	if(IDVar[id][3] > _IDcuts[EP][3] && Cbool(skipCut, "HoE"==observable.substr(0,3)) )  Sel[id]=false;
	
	if(IsoVar[id][0]/Lepton[id][0] > _Isocuts[EP][0] && Cbool(skipCut, "TrackIso"==observable.substr(0,8)) )  Sel[id]=false;
	if(IsoVar[id][1]/Lepton[id][0] > _Isocuts[EP][1] && Cbool(skipCut, "EcalIso"==observable.substr(0,7)) )  Sel[id]=false;
	if(IsoVar[id][2]/Lepton[id][0] > _Isocuts[EP][2] && Cbool(skipCut, "HcalIso"==observable.substr(0,7)) )  Sel[id]=false;
	

	if(EcalP!=EP && EcalP!=2 && EcalP!=3) Sel[id]=false;

	if(EcalP==3) 
	  epart[id]=EP;
	
	//ID
	histoManager.fill("sigieie"+lep[id],i, IDVar[id][0],Weight);
	histoManager.fill("Deta"+lep[id],i, fabs(IDVar[id][1]),Weight);
	histoManager.fill("Dphi"+lep[id],i, IDVar[id][2],Weight);
	histoManager.fill("HoE"+lep[id],i, IDVar[id][3],Weight);
	
	//Iso
	histoManager.fill("TrackIso"+lep[id],i,IsoVar[id][0]/Lepton[id][0],Weight);
	histoManager.fill("EcalIso"+lep[id],i,IsoVar[id][1]/Lepton[id][0],Weight);
	histoManager.fill("HcalIso"+lep[id],i,IsoVar[id][2]/Lepton[id][0],Weight);
	

      }

      if(!Sel[0] || !Sel[1])
	continue;
      
      if(EcalP==3 && epart[0]==epart[1] ) continue;
      


      if( (!invCut && ZMass > MTCut && ZMass < 140) || ( invCut && ZMass > MTCut && ZMass < 140 ) ) { 
	for(int id=0;id<2;id++) {
	  histoManager.fill("EtPFSC",i,PfSCEt[id],Weight);
	  histoManager.fill("RawEtPFSC",i,PfSCRawE[id],Weight);
	  
	  float brPFSC = computeBR( PfSCShape[id][0], PfSCShape[id][1] );
	  float corE = GetCalibCorEnergy(brPFSC, PfSCRawE[id] +  PfSCESE[id]  , PfSCGeo[id][0] );
	  float corEt = ConversionEpt(corE, PfSCGeo[id][0] );

	  histoManager.fill("EtPFSCCor",i,corEt,Weight);
	 
	  histoManager.fill("Eta"+lep[id],i,PfSCGeo[id][0],Weight);
	  histoManager.fill("Phi"+lep[id],i,PfSCGeo[id][1],Weight);
	  histoManager.fill("EcalIso"+lep[id],i,IsoVar[id][1]/Lepton[id][0],Weight);
	  histoManager.fill("TrackIso"+lep[id],i,IsoVar[id][0]/Lepton[id][0],Weight);
	  histoManager.fill("HcalIso"+lep[id],i,IsoVar[id][2]/Lepton[id][0],Weight);
	  //	  if(fabs(PfSCGeo[id][0])>1.5) 
	    /*  cout<<PfSCRawE[id]<<"    "<<GetBremCorEE(brPFSC,PfSCGeo[id][0])<<"  -->  "
		<<PfSCRawE[id]/GetBremCorEE(brPFSC,PfSCGeo[id][0])<<"  ---->::   "
		<<1/CalibCorEB(ConversionEpt(PfSCRawE[id]/GetBremCorEE(brPFSC,PfSCGeo[id][0]),PfSCGeo[id][0]),PfSCGeo[id][0])
		<<"  -->  "<<corE<<"    "<<PfSCE[id]<<endl;*/

	  histoManager.fillProf("RatioECor",i,PfSCRawE[id],corE/PfSCRawE[id],Weight);
	  histoManager.fillProf("RatioE",i,PfSCRawE[id],PfSCE[id]/PfSCRawE[id],Weight);

	  lepton[id].SetPtEtaPhiM(corEt,PfSCGeo[id][0],PfSCGeo[id][1],0.000000511);
	}

	//Rebuild the Z
	TLorentzVector Zvect = lepton[0] + lepton[1];
	histoManager.fill("ZPhi",i,Zvect.Phi(),Weight);
	histoManager.fill("ZEta",i,Zvect.Eta(),Weight);

	//	if( (fabs(PfSCGeo[0][0])<1.5 && fabs(PfSCGeo[0][0])>1.1) || 
	//	    (fabs(PfSCGeo[1][0])<1.5 && fabs(PfSCGeo[1][0])>1.1) ) {
	
	  //	if( computeBR( PfSCShape[0][0], PfSCShape[0][1] ) > 2.5 || 
	  //  computeBR( PfSCShape[1][0], PfSCShape[1][1] ) > 2.5 ) {
	  //ZPt
	  histoManager.fill("ZPt",i,ZPt,Weight);
	  histoManager.fill("ZMass",i,ZMass,Weight);
	
	  //Uncorrected mass ===============================
	  float rawEt1 = ConversionEpt(PfSCRawE[0],PfSCGeo[0][0]);
	  float rawEt2 = ConversionEpt(PfSCRawE[1],PfSCGeo[1][0]);
	  histoManager.fill("RawZMass",i,ComputeMinv(rawEt1,rawEt2,PfSCGeo[0][0],PfSCGeo[1][0],PfSCGeo[0][1],PfSCGeo[1][1]),Weight);
	  //================================================

	  //Corrected mass ===========================
	  histoManager.fill("ZMassCor",i,ComputeMinv(lepton[0].Pt(),lepton[1].Pt(),PfSCGeo[0][0],PfSCGeo[1][0],PfSCGeo[0][1],PfSCGeo[1][1]) ,Weight);
	  //=========================================
	
	  //	}

	//Vertex
	histoManager.fill("Vertex",i,Nvertex,Weight);

	  if(name[i].substr(0,4)=="data")
	    {
	   
	      if(ZPt>=80 && ZPt<=90 )
		{
		  ostringstream os;
		  os << AbsEvent;
		
		  std::pair<int,int> tmp(Run,Event);
		  string t2(fileName); //,t1(sample);
		  std::pair<string,string> tmp2( os.str(), t2 );
		  std::pair< std::pair<int,int> , std::pair<string,string> > tmp3(tmp,tmp2);
		  Events.push_back(tmp3);
		}
	    }

      } //End if condition >         
      
      NumberEntries[i]++;
	       
    }//End events
    
  } //End datasets
 
  cout<<" End Filling , S="<<S<<" ;  B= "<<B<<endl;
  cout<<" detail "<<endl;
  for(int i=0;i<nt+1;i++) {
    cout<<" --->  "<<name[i]<<"  "<< NumberEntries[i]<<endl;
  }

  Histos = histoManager.GetHistos();
  Histos2D = histoManager.GetHistos2D();
  Profiles = histoManager.GetProfiles();

  for(size_t i=0;i<Histos.size();i++)
    Weighted.push_back(false);

  cout<<"Histos "<<Histos.size()<<"   "<<Histos2D.size()<<endl;
  

}



void
PFCorUseTree::PrepareDatasets() {
  
  reposi="PFSCCor";
    
    colors.push_back(kBlack);
    
    for(int unsigned i=0;i<datasets.size();i++)
      FillAddWeight(datasets[i]);

    name.push_back("data");
    
}


float PFCorUseTree::ComputeMinv(float Et1, float Et2, float eta1, float eta2, float phi1, float phi2) {

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




float PFCorUseTree::ConversionEtaTheta(float eta) //Convertion eta en theta
{
  float theta;
  theta = 2*atan(exp(-eta));
  return theta;
}

float PFCorUseTree::ConversionThetaEta(float theta) //Conversion theta en eta
{
  float eta;
  eta = -log(tan((theta)/2));
  return eta;
}



float PFCorUseTree::ConversionEpt(float energy,float eta ) //Convertion Energie en pt
{

  float pt =  energy*sin(ConversionEtaTheta(eta));
  return pt;
}

float PFCorUseTree::ConversionPtE(float pt, float eta) //Conversion pt en energie
{
  float energy = pt/sin(ConversionEtaTheta(eta));
  return energy;

}

void PFCorUseTree::Conversion_REP_carte(float coord[3], float& x, float& y, float &z) //Conversion Rho/eta/phi en coordonnées cartésiennes
{
  float Theta = ConversionEtaTheta(coord[1]);
  x = coord[0]*sin(Theta)*cos(coord[2]);
  y = coord[0]*sin(Theta)*sin(coord[2]);
  z = coord[0]*cos(Theta);
  return;
}

float 
PFCorUseTree::GetCalibCorEnergy(float brLinear,float e, float eta) {

  if(fabs(eta) < 1.479 )
    {return  GetCalibCorEB(brLinear, e, eta);}
  else
    {return  GetCalibCorEE(brLinear, e, eta);}

}


float 
PFCorUseTree::GetCalibCorEB(float brLinear, float e, float eta) {

  float et = ConversionEpt(e,eta);
  float bremCorE = e/GetBremCorEB(brLinear, et);
  float bremCorEt = ConversionEpt(bremCorE,eta);
  float corE = bremCorE/CalibCorEB(bremCorEt, eta) ;

  return corE;
}



float 
PFCorUseTree::GetCalibCorEE(float brLinear, float e, float eta) {

  float bremCorE = e/GetBremCorEE(brLinear, eta);
  float bremCorEt = ConversionEpt(bremCorE,eta);
  float corE = bremCorE/CalibCorEE(bremCorEt, eta);
  // cout<<corE<<"    "<<bremCorE<<"    "<<e<<endl;
  return corE;
}



float 
PFCorUseTree::GetBremCorEB(float brLinear, float et) {

  if ( brLinear == 0 ) return 1;
  
  if ( brLinear < 0.6 ) brLinear = 1.1;
  if ( brLinear > 8 ) brLinear = 8.0;

    float p0 = -0.0255975;
    float p1 = 0.0576727;
    float p2 = 0.975442;
    float p3 = -0.000546394;
    float p4 = 1.26147;

    if(et<25) { //For Low PT electrons
       p0 = -0.02025;
       p1 = 0.04537;
       p2 = 0.9728;
       p3 = -0.0008962;
       p4 = 1.172;
    }


    float threshold = p4;
    float y = p0*threshold*threshold + p1*threshold + p2;
    float yprime = 2*p0*threshold + p1;
    float a = p3;
    float b = yprime - 2*a*threshold;
    float c = y - a*threshold*threshold - b*threshold;

    float fCorr = 1;
    if ( brLinear < threshold ) 
      fCorr = p0*brLinear*brLinear + p1*brLinear + p2;
    else 
      fCorr = a*brLinear*brLinear + b*brLinear + c;
    
    return fCorr;
}



float 
PFCorUseTree::GetBremCorEE(float brLinear, float eta) {

  if ( brLinear == 0 ) return 1;

 if ( brLinear > 6.5 ) brLinear = 6.5;
    if ( brLinear < 0.9 ) brLinear = 0.9;

    // ============= Fixed Matrix With Preshower SC
    float p0 = -0.0692932;  
    float p1 = 0.101776; 
    float p2 = 0.995338;  
    float p3 = -0.00236548;
    float p4 = 0.874998;  
    if(fabs(eta) >1.653) {
      p0 = -0.0750184; 
      p1 = 0.147000;  
      p2 = 0.923165;  
      p3 = 0.000474665; 
      p4 = 1.10782;  
    }

    float threshold = p4;
    float y = p0*threshold*threshold + p1*threshold + p2;
    float yprime = 2*p0*threshold + p1;
    float a = p3;
    float b = yprime - 2*a*threshold;
    float c = y - a*threshold*threshold - b*threshold;

    float fCorr = 1;
    if ( brLinear < threshold ) 
      fCorr = p0*brLinear*brLinear + p1*brLinear + p2;
    else 
      fCorr = a*brLinear*brLinear + b*brLinear + c;

    return fCorr;

}



float 
PFCorUseTree::CalibCorEB(float et, float eta) {

  float fCorr = 1;
  
  float c0 = 1.004;
  float c1 = -1.536;
  float c2 = 22.88;

  // fEtEta
  float c4 = 0.3555;
  float c5 = 0.6227;
  float c6 = 14.65;
  float c6b = 2051;

  // final fitting
  float c7 = 1.081;  // curve point in eta distribution
  float c8 = 7.6;     // sharpness of the curve
  float c3 = -0.00181;


  float p0 = c0 + c1/(et + c2) -1.467 /(et*et);
  float p1 = c4/(et + c5) + c6/(et*et + c6b);

  fCorr = p0 + p1*atan(c8*(c7 - fabs(eta))) + c3*fabs(eta);
  /*
  if(et<25) { //Low PT Electrons
    c0 = 0.9932;
    c1 = -0.5444;
    c2 = 0.5438;
    
    // fEtEta
    c4 = 0.7109;
    c5 = 7.645;
    c6 = 0.2904;
    
    
    // final fitting
    c7 = 1.081;  // curve point in eta distribution
    c8 = 7.6;     // sharpness of the curve
    c3 = -0.00181;
    
    
    p0 = c0 + c1/( et ) + c2/( et*et);
    p1 =  c4/(et +c5 ) + c6/( et*et);
    
    fCorr = p0 + p1*atan(c8*(c7 - fabs(eta))) + c3*fabs(eta);
  }*/

  return fCorr;
}

float 
PFCorUseTree::CalibCorEE(float et, float eta) {

  float fCorr = 1;

  double p0 = 1.153 - 16.5975/(5.668 + et);
  double p1 = -0.1772 + 16.22/(7.326 + et );
  double p2 = 0.0483 - 4.068/(9.406 + et );

  fCorr = p0 + p1*fabs(eta) + p2*fabs(eta)*fabs(eta);
  return fCorr;

}



TH1F*
PFCorUseTree::GetVtxCorrection(string obs, int nds, int bin) {
  return NULL; 
  
}

vector<vector<vector<float> > >
PFCorUseTree::GetFitWeight() {
  vector<vector<vector<float> > > tmp;
  return tmp;

}
