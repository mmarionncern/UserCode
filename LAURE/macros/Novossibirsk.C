#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <map>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMinuit.h"

using namespace std;
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>

#include <TH1.h>

#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooStringVar.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooAbsBinning.h>
#include <RooArgList.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooDataHist.h>
#include <RooRandom.h>
#include <RooSimultaneous.h>
#include <RooHist.h>
// models
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooGExpModel.h>
#include <RooVoigtian.h>
#include <RooCBShape.h>
#include <RooArgusBG.h>
#include <RooBreitWigner.h>
#include <RooFFTConvPdf.h>
#include <RooNumConvPdf.h>
#include <RooNovosibirsk.h>
#include <RooLandau.h>

#include <RooMsgService.h>

using namespace RooFit;


double novosibirsk(double* x, double* p)  {

  //p0 : norm
  //p1 : mean
  //p2 : width 
  //p3 : tail

  double qa, qb, qc, qx, qy;

  if(TMath::Abs(p[3]) < 1.e-7) 
    qc = 0.5*TMath::Power(((x[0]-p[1])/p[2]),2);
  else {
    qa = p[3]*sqrt(log(4.));
    qb = sinh(qa)/qa;
    qx = (x[0]-p[1])/p[2]*qb;
    qy = 1.+p[3]*qx;
  
    //---- Cutting curve from right side

    if( qy > 1.E-7) 
      qc = 0.5*(TMath::Power((log(qy)/p[3]),2) + p[3]*p[3]);
    else
      qc = 15.0;
  }

  //---- Normalize the result

  return p[0]*exp(-qc);
  
}


vector<vector<float> >  Novossibirsk(string name) {


  RooRealVar pmet( "pmet","pmet",0,100,"GeV" );

  TFile* mcFile = new TFile(name.c_str(),"READ");
  TH1F* Zee = (TH1F*)mcFile->Get("Zll");


  RooDataHist data_MC("data_MC","dataset with x",pmet, Zee ) ;
  

  RooRealVar* peak;

  //  if(Nvtx==0)
  peak= new RooRealVar( "peak","peak", 0/*,-2,2*/, "GeV" );
  /*  else
    peak= new RooRealVar( "peak","peak", 0,-5,5, "GeV" );*/
  RooRealVar width( "width", "width",1,0,10, "GeV" );
  RooRealVar tail( "tail", "tail",0.2,0,10, "GeV" );
  RooNovosibirsk Nov( "novo","novosi", pmet , *peak , width, tail );

  TCanvas* c3 =new TCanvas( ("c"+name).c_str(),("frameNV"+name).c_str(),600,600);

  RooFitResult* result2 = Nov.fitTo(data_MC,RooFit::SumW2Error(kFALSE) );

  RooPlot* frame = pmet.frame() ;
 
  data_MC.plotOn(frame) ;
  Nov.plotOn(frame) ;
  Nov.paramOn(frame);
  frame->Draw();
  frame->GetYaxis()->SetRangeUser(0.00001,1000);
  frame->GetXaxis()->SetRangeUser(0,60);

  RooCurve* func = (RooCurve*) frame->findObject("novo_Norm[pmet]");
  TH1F* nov = func->GetHistogram();
  // nov->Draw("same");
  //cout<<nov<<endl;
  
  double x,y;
  int ii;

  for(int i=0;i<func->GetN();i++) {

    func->GetPoint(i,x,y);
    //cout<<i<<"  "<<x<<"  "<<y<<endl;
    
    if(x>=4.2)
      { ii = i; break;}
  }

  cout<<" ----> Integral Novo "<<"   " <<func->Integral( ii , func->GetN() )<<" --> "<<func->Integral(ii, func->GetN() )/func->Integral(1, func->GetN() )<<endl;


  vector<float> Values;
  vector<float> Errors;
  Values.push_back( peak->getVal() );
  Values.push_back( width.getVal() );
  Values.push_back( tail.getVal() );
  Errors.push_back( peak->getError() );
  Errors.push_back( width.getError() );
  Errors.push_back( tail.getError() );

  vector<vector<float> > vals;
  vals.push_back(Values);
  vals.push_back(Errors);

  return vals;
}


void Argus() {


  RooRealVar pmet( "scpmet","scpmet",0,1000 );
 

  TFile* mcFile = new TFile("SigniCorProjMETEnd.root","READ");
  TH1F* ZZ = (TH1F*)mcFile->Get("dib");
  ZZ->Rebin(4);

  TH1F* ZZ2 = (TH1F*)ZZ->Clone();
  /*ZZ2->Reset("icem");
  
  for(int i=0;i<ZZ->GetNbinsX()+2;i++) {
    ZZ2->SetBinContent(i, ZZ->GetBinContent(ZZ->GetNbinsX() +1 -i) );
    ZZ2->SetBinError(i, ZZ->GetBinError(ZZ->GetNbinsX() +1 -i) );
  }*/

  RooDataHist data_MC("data_MC","dataset with x",pmet, ZZ2 ) ;
  
  /*  RooRealVar peak( "peak","peak", 2,0,10, "GeV" );
      RooRealVar width( "width", "width",1,0,4, "GeV" );
      RooRealVar tail( "tail", "tail",0.2,0,10, "GeV" );
  RooNovosibirsk argus( "novo","novosi", pmet , peak , width, tail );
  */

  /*  RooRealVar argmean("argmean","argus mean parameter",-20.0,-100.,-1.) ;
  RooRealVar argpar("argpar","argus shape parameter",-20.0,-100.,-1.) ;
  RooArgusBG argus("background","Argus PDF",pmet, argmean ,argpar) ;*/
  
  RooRealVar cb_bias_MC( "cb_b_MC", "bias",0.0, -2.0, 2.0 );
  RooRealVar cb_width_MC( "cb_w_MC","width", 2,0.0,3 ); 
  RooRealVar cb_alpha_MC( "cb_a_MC","alpha", 0.94,0.0,2.0 ); 
  RooRealVar cb_power_MC( "cb_n_MC","power", 5/*, 0.5, 5.0*/ ); 
  RooCBShape argus( "cb_pdf_MC", "A  Crystal Ball Lineshape", 
			pmet,cb_bias_MC,  cb_width_MC, 
			cb_alpha_MC, cb_power_MC );


  TCanvas* c3 =new TCanvas( "c","pot",600,600);

  RooFitResult* result2 = argus.fitTo(data_MC,RooFit::SumW2Error(kFALSE) );

  RooPlot* frame = pmet.frame() ;
 
  data_MC.plotOn(frame) ;
  argus.plotOn(frame) ;
  argus.paramOn(frame);
  frame->Draw();

  


}


void MakeGraphs(int n=5) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
  //  int n=5;

  TGraphErrors* mean=new TGraphErrors(n-1);
  TGraphErrors* width=new TGraphErrors(n-1);
  TGraphErrors* tail=new TGraphErrors(n-1);

  mean->GetXaxis()->SetLimits(0,n+1);
  width->GetXaxis()->SetLimits(0,n+1);
  tail->GetXaxis()->SetLimits(0,n+1);

  vector<vector<float> > vals(3,vector<float>(3,0));

  
  for(int i=1;i<n;i++) {
    
    ostringstream os;
    os<<i;

    string name= "rootfiles/ProjMET_V"+os.str()+".root";
    
    vals = Novossibirsk(name);
    
    mean->SetPoint( i-1, i , vals[0][0]  );
    mean->SetPointError( i-1, 0.1 , vals[1][0]  );

    width->SetPoint( i-1, i , vals[0][1]  );
    width->SetPointError( i-1, 0.1 , vals[1][1]  );
    
    tail->SetPoint( i-1, i , vals[0][2]  );
    tail->SetPointError( i-1, 0.1 , vals[1][2]  );

    
  }

  TCanvas* c= new TCanvas("c","c",1200,400);
  c->Divide(3);

  c->cd(1);
  mean->Draw("AP");
  mean->GetXaxis()->SetTitle("Nvtx");
  mean->GetYaxis()->SetTitle("Mean");
  TF1* poly=new TF1("poly","pol1",0,60);
  mean->Fit("poly");

  c->cd(2);
  width->Draw("AP");
  width->GetXaxis()->SetTitle("Nvtx");
  width->GetYaxis()->SetTitle("Width [GeV]");
  TF1* square=new TF1("square","sqrt([1]*[1]*(x-1) + [0]*[0])",60);
  square->SetParName(0,"#sigma_{HP}");
  square->SetParName(1,"#sigma_{PU}");
  width->Fit("square");


  c->cd(3);
  tail->Draw("AP");
  tail->GetXaxis()->SetTitle("Nvtx");
  tail->GetYaxis()->SetTitle("Tail [GeV]");
  //  TF1* expf=new TF1("expf","expo",0,60);
  TF1* expf=new TF1("expf","[0] + (x-1)*[1]",0,60);
  tail->Fit("expf");

  int N=0;

  for (int k=0;k<4;k++) {
    N=k;

    float tau = expf->GetParameter(0) + N*expf->GetParameter(1) ;
    float sigma = sqrt( pow(square->GetParameter(0),2) + N*pow(square->GetParameter(1),2 ) ) ;
    
    float A = sqrt(log(4));
    
    float sigp = (pow(tau,4)*sigma*sigma*A*A)/( pow(sinh(tau*A)*tau,2) );
    
    cout<<k<<"   "<<sigma<<"  sigma prime    "<<sqrt(sigp)<<endl;
  }

}





void MakeGraphsGaus(int n=5) {

  //  int n=5;

  TGraphErrors* mean=new TGraphErrors(n-1);
  TGraphErrors* width=new TGraphErrors(n-1);
 
  mean->GetXaxis()->SetLimits(0,n+1);
  width->GetXaxis()->SetLimits(0,n+1);
 
  vector<vector<float> > vals(3,vector<float>(3,0));
  TCanvas* cc= new TCanvas("cc","cc",400,400);
  int N=0;
  for(int i=1;i<n;i++) {
    
    ostringstream os;
    os<<i;
    //ProjMET/PowhegFall10/  ProjMET/goodFiles_Fall10_powheg/
    string name= "ProjMET/goodFiles_Fall10_powheg/ProjMET_V"+os.str()+".root";
    
    TFile* mcFile = new TFile(name.c_str(),"READ");
    TH1F* Zee = (TH1F*)mcFile->Get("data");
    
    Zee->GetXaxis()->SetLabelOffset(0.005);
    Zee->GetYaxis()->SetLabelOffset(0.005);
    
    Zee->GetXaxis()->SetLabelSize(0.05);
    Zee->GetXaxis()->SetTitleSize(0.05);
    
    Zee->GetYaxis()->SetLabelSize(0.05);
    Zee->GetYaxis()->SetTitleSize(0.05);
    
    Zee->GetXaxis()->SetNdivisions(6,5,0);
    Zee->GetYaxis()->SetNdivisions(6,5,0);
    
    Zee->GetXaxis()->SetTitleOffset(1);
    Zee->GetYaxis()->SetTitleOffset(1.2);
    
    Zee->GetYaxis()->SetRangeUser(10,400000);
    
    Zee->SetLineColor(1);

    Zee->GetYaxis()->SetTitle("nombre d'evenement");
    Zee->GetXaxis()->SetTitle("proj#slash{E}_{T} [GeV]");
    
    TF1* gaus=new TF1("f","gaus",0,50);
    gaus->SetLineColor(kBlue+1);
    gaus->SetLineWidth(2);
    gaus->SetLineStyle(2);
    if(i!=1)
      Zee->Fit("f","N0Q");
    else
      Zee->Fit("f","Q"); 

    // if(i!=3) {
      
      mean->SetPoint( i-1, i , gaus->GetParameter(1)  );
      mean->SetPointError( i-1, 0.1 , gaus->GetParError(1)   );
      
      width->SetPoint( i-1, i , gaus->GetParameter(2)  );
      width->SetPointError( i-1, 0.1 , gaus->GetParError(2) );
      
      cout<<" ----> Integral Gaus "<<i<<"   " <<gaus->Integral(31, 10000)<<endl;
      N++;
      //}

      /*  delete gaus;
    delete Zee;
    delete mcFile;*/
    
  }

  cc->SetLogy(1);

  gStyle->SetOptFit(0);
  TCanvas* c= new TCanvas("c","c",800,400);
  c->Divide(2);

  c->cd(1);
  width->Draw("AP");
  width->GetXaxis()->SetTitle("nombre de vertex");
  width->GetYaxis()->SetTitle("#sigma_{G} [GeV]");
  width->GetYaxis()->SetTitleOffset(1.13);
  width->GetXaxis()->SetNdivisions(7,5,0);
  TF1* square=new TF1("square","sqrt([1]*[1]*(x-1) + [0]*[0])",0,60);
  square->SetLineStyle(2);
  square->SetLineWidth(2);
  square->SetLineColor(kBlue+1);
  square->SetParameter(0,7.4);
  square->SetParameter(1,3);
  square->SetParName(0,"#sigma_{HP}");
  square->SetParName(1,"#sigma_{PU}");
  width->Fit("square","Q");

  

  /* TF1* square2=new TF1("square","sqrt([1]*[1]*(x-1) + [0]*[0])",0,60);
  square2->SetLineColor(kGreen+1);
  square2->SetParameter(0,7.36);
  square2->SetParameter(1,3.1);
  square2->Draw("same+");
*/
  double matrix[2][2];
  gMinuit->mnemat(&matrix[0][0],2);
  cout<<endl;
  cout<<" double SigM[2][2]={{"<<matrix[0][0]<<"  ,   "<<matrix[0][1]<<"},"<<endl;
  cout<<"{"<<matrix[1][0]<<"  ,   "<<matrix[1][1]<<"}};"<<endl;
  cout<<endl;
  cout<<"  double SigE[2] = {"<<square->GetParError(0)<<","<<square->GetParError(1)<<"};"<<endl;
  cout<<"  double SigV[2] = {"<<square->GetParameter(0)<<","<<square->GetParameter(1)<<"};"<<endl;
  cout<<endl;

  c->cd(2);
  mean->Draw("AP");
  mean->GetXaxis()->SetTitle("nombre de vertex");
  mean->GetYaxis()->SetTitle("x_{0} [GeV]");
  mean->GetYaxis()->SetTitleOffset(1.13);
  mean->GetXaxis()->SetNdivisions(7,5,0);
  //  TF1* expf=new TF1("expf","expo",0,60);
  TF1* expf=new TF1("expf","[0] + (x-1)*[1]",0,60);
  expf->SetLineStyle(2);
  expf->SetLineWidth(2);
  expf->SetLineColor(kBlue+1);
  expf->SetParameter(0,-4);
  expf->SetParameter(1,1);
  mean->Fit("expf","Q");

  gMinuit->mnemat(&matrix[0][0],2);
  cout<<endl;
  cout<<" double MeanM[2][2]={{"<<matrix[0][0]<<"  ,   "<<matrix[0][1]<<"},"<<endl;
  cout<<"{"<<matrix[1][0]<<"  ,   "<<matrix[1][1]<<"}};"<<endl;
  cout<<endl;
  cout<<"  double meanE[2] = {"<<expf->GetParError(0)<<","<<expf->GetParError(1)<<"};"<<endl;
  cout<<"  double MeanV[2] = {"<<expf->GetParameter(0)<<","<<expf->GetParameter(1)<<"};"<<endl;
  cout<<endl;

  // int N=0;

  // for (int k=0;k<4;k++) {
  //   N=k;

  //   float tau = expf->GetParameter(0) + N*expf->GetParameter(1) ;
  //   float sigma = sqrt( pow(square->GetParameter(0),2) + N*pow(square->GetParameter(1),2 ) ) ;
    
  //   float A = sqrt(log(4));
    
  //   float sigp = (pow(tau,4)*sigma*sigma*A*A)/( pow(sinh(tau*A)*tau,2) );
    
  //   cout<<k<<"   "<<sigma<<"  sigma prime    "<<sqrt(sigp)<<endl;
  // }

  //  cc->SaveAs("~/Documents/CMS/Redaction/These/Figures/ZZ/ParaProjMET.eps");
  // c->SaveAs("~/Documents/CMS/Redaction/These/Figures/ZZ/ParaVtx.eps");


}




void SingleGaus() {
  
  string name= "rootfiles/SigniCorProjMET.root";
    
  TFile* mcFile = new TFile(name.c_str(),"READ");
  TH1F* Zee = (TH1F*)mcFile->Get("Zll");
  TH1F* totaMC = (TH1F*)mcFile->Get("totalObs");
  TH1F* data = (TH1F*)mcFile->Get("data");
  
  Zee->GetXaxis()->SetLabelOffset(0.005);
  Zee->GetYaxis()->SetLabelOffset(0.005);

  Zee->GetXaxis()->SetLabelSize(0.05);
  Zee->GetXaxis()->SetTitleSize(0.05);

  Zee->GetYaxis()->SetLabelSize(0.05);
  Zee->GetYaxis()->SetTitleSize(0.05);

  Zee->GetXaxis()->SetNdivisions(6,5,0);
  Zee->GetYaxis()->SetNdivisions(6,5,0);

  Zee->GetXaxis()->SetTitleOffset(1);
  Zee->GetYaxis()->SetTitleOffset(1);

  Zee->GetYaxis()->SetRangeUser(0.001,400);
  Zee->SetLineColor(1);

  TCanvas* c = new TCanvas("c222","MC",600,600);

  TH1F* b=new TH1F("b","b",10,0,10);
 
  Zee->Draw("");
  Zee->GetXaxis()->SetRangeUser(0,6);
  Zee->GetYaxis()->SetRangeUser(0.001,1000);

  /* TF1* gaus=new TF1("f","gaus",0,50);
  gaus->SetLineColor(kRed+1);
  gaus->SetLineWidth(2);
  Zee->Fit("f","+");*/
 
  TF1* novo=new TF1("f2",novosibirsk,0,50,4);
  novo->SetParameter(0,600);
  novo->FixParameter(1,0);
  novo->SetParameter(2,0.4);
  novo->SetParameter(3,0.2);
  novo->SetLineColor(kBlue+1);
  novo->SetLineWidth(2);
  Zee->Fit("f2","+");

  // cout<<gaus->Integral(0, 1000.)<<" ----> Integral Gaus "<<"   " <<gaus->Integral(4.2, 1000.)<<" --> "<<gaus->Integral(4.2, 1000.)/gaus->Integral(0, 1000.)<<endl;

    
  c->RedrawAxis();

  TH1F* dataClean = (TH1F*)data->Clone();
  TH1F* otherProc = (TH1F*)totaMC->Clone();
  otherProc->Add(Zee, -1);

  dataClean->Add(otherProc, -1);
  
  for(int i=0;i<dataClean->GetNbinsX();i++) {
    if(dataClean->GetBinContent(i) < 0) {
      dataClean->SetBinContent(i, 0.);
      dataClean->SetBinError(i, 100000.);
    }

  }

  TF1* novoD=new TF1("nD",novosibirsk,0,50,4);
  novoD->SetParameter(0,600);
  novoD->FixParameter(1,0);
  novoD->SetParameter(2,0.4);
  novoD->SetParameter(3,0.2);
  novoD->SetLineColor(kGreen+1);
  novoD->SetLineWidth(2);

  
  TCanvas* c2 = new TCanvas("c2222","data",600,600);
  dataClean->Draw();
  dataClean->GetXaxis()->SetRangeUser(0,6);
  dataClean->GetYaxis()->SetRangeUser(0.001,1000);
  dataClean->Fit("nD","+");
  

  cout<<novo->Integral(0, 1000.)<<" ----> Integral Novo "<<"   " <<novo->Integral(4.2, 1000.)<<" --> "<<novo->Integral(4.2, 1000.)/novo->Integral(0, 1000.)<<endl;
  cout<<novoD->Integral(0, 1000.)<<" ----> Integral NovoData "<<"   " <<novoD->Integral(4.2, 1000.)<<" --> "<<novoD->Integral(4.2, 1000.)/novoD->Integral(0, 1000.)<<endl;

  TCanvas* c23 = new TCanvas("c22221","data",600,600);
  Zee->Draw();
  dataClean->SetLineColor(kRed+1);
  dataClean->Draw("same");
  novo->Draw("same");
  novoD->Draw("same");
  //  Novossibirsk(name);

}



void CompareElMu() {
  
  string name1= "CorProjMET.root";
  string name2= "CorProjMET_Mu.root";
    
  TFile* mcFile = new TFile(name1.c_str(),"READ");
  TH1F* Zee = (TH1F*)mcFile->Get("Zll");
  TH1F* totaMC = (TH1F*)mcFile->Get("totalObs");
  TH1F* data = (TH1F*)mcFile->Get("data");

  TFile* mcFile2 = new TFile(name2.c_str(),"READ");
  TH1F* Zmm = (TH1F*)mcFile2->Get("Zll");
  TH1F* totaMC2 = (TH1F*)mcFile2->Get("totalObs");
  TH1F* data2 = (TH1F*)mcFile2->Get("data");

  
  Zee->GetXaxis()->SetLabelOffset(0.005);
  Zee->GetYaxis()->SetLabelOffset(0.005);

  Zee->GetXaxis()->SetLabelSize(0.05);
  Zee->GetXaxis()->SetTitleSize(0.05);

  Zee->GetYaxis()->SetLabelSize(0.05);
  Zee->GetYaxis()->SetTitleSize(0.05);

  Zee->GetXaxis()->SetNdivisions(6,5,0);
  Zee->GetYaxis()->SetNdivisions(6,5,0);

  Zee->GetXaxis()->SetTitleOffset(1);
  Zee->GetYaxis()->SetTitleOffset(1);

  Zee->GetYaxis()->SetRangeUser(0.001,400);
  Zee->SetLineColor(1);


  Zmm->GetXaxis()->SetLabelOffset(0.005);
  Zmm->GetYaxis()->SetLabelOffset(0.005);

  Zmm->GetXaxis()->SetLabelSize(0.05);
  Zmm->GetXaxis()->SetTitleSize(0.05);

  Zmm->GetYaxis()->SetLabelSize(0.05);
  Zmm->GetYaxis()->SetTitleSize(0.05);

  Zmm->GetXaxis()->SetNdivisions(6,5,0);
  Zmm->GetYaxis()->SetNdivisions(6,5,0);

  Zmm->GetXaxis()->SetTitleOffset(1);
  Zmm->GetYaxis()->SetTitleOffset(1);

  Zmm->GetYaxis()->SetRangeUser(0.001,400);
  Zmm->SetLineColor(1);



  TCanvas* c = new TCanvas("c222","MC",600,600);

  TH1F* b=new TH1F("b","b",10,0,10);
  Zee->GetXaxis()->SetRangeUser(0,60);
  Zee->GetYaxis()->SetRangeUser(0.001,2);
  Zee->Scale( 1./Zee->Integral() );
  Zee->Draw(" L");

  Zmm->SetLineColor(kRed+1);
  Zmm->Scale( 1./Zmm->Integral() );
  Zmm->DrawNormalized("same L");
  Zmm->GetXaxis()->SetRangeUser(0,30);
  Zmm->GetYaxis()->SetRangeUser(0.001,500);


  TF1* novo=new TF1("f2",novosibirsk,0,50,4);
  novo->SetParameter(0,0.1);
  novo->FixParameter(1,0);
  novo->SetParameter(2,0.4);
  novo->SetParameter(3,0.2);
  novo->SetLineColor(kBlue+1);
  novo->SetLineWidth(2);
  Zee->Fit("f2","N0");





  TF1* novo2=new TF1("f3",novosibirsk,0,50,4);
  novo2->SetParameter(0,0.1);
  novo2->FixParameter(1,0);
  novo2->SetParameter(2,0.4);
  novo2->SetParameter(3,0.2);
  novo2->SetLineColor(kGreen+1);
  novo2->SetLineWidth(2);
  Zmm->Fit("f3","N0");

  novo->Draw("same");
  novo2->Draw("same");

  TLegend* leg= new TLegend(0.50,0.5,0.7,0.7);
  leg->AddEntry(Zee,"Zee","pl");
  leg->AddEntry(novo,"fit Zee","pl");
  leg->AddEntry(Zmm,"Zmm","pl");
  leg->AddEntry(novo2,"fit Zmm","pl");
  c->RedrawAxis();
  leg->Draw();
  
  TH1F* dataClean = (TH1F*)data->Clone();
  TH1F* otherProc = (TH1F*)totaMC->Clone();
  otherProc->Add(Zee, -1);

  dataClean->Add(otherProc, -1);


  TH1F* dataClean2 = (TH1F*)data2->Clone();
  TH1F* otherProc2 = (TH1F*)totaMC2->Clone();
  otherProc2->Add(Zmm, -1);

  dataClean2->Add(otherProc2, -1);
  
  for(int i=0;i<dataClean2->GetNbinsX();i++) {
    if(dataClean2->GetBinContent(i) < 0) {
      dataClean2->SetBinContent(i, 0.);
      dataClean2->SetBinError(i, 0.);
    }

    if(dataClean->GetBinContent(i) < 0) {
      dataClean->SetBinContent(i, 0.);
      dataClean->SetBinError(i, 0.);
    }

  }

  TCanvas* c2 = new TCanvas("c2222","data",600,600);

  dataClean->GetXaxis()->SetRangeUser(0,60);
  dataClean->GetYaxis()->SetRangeUser(0.001,2500);
  dataClean->DrawNormalized();


  dataClean2->SetMarkerColor(kRed+1);
  dataClean2->SetLineColor(kRed+1);
  dataClean2->DrawNormalized("same");




}
