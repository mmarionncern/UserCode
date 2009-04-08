

#ifndef mmUtil_cc
#define mmUtil_cc

#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <map>

#include <Math/VectorUtil.h>
//#include <TMath.h>
//#include <TVector3.h>
//#include <TLorentzVector.h>

using namespace std;

static const float myPi = acos(-1);
static const float mZ0 = 91.188;

int diveucli(int Event_cut_flag, int i) //Division euclidienne
{
  int reste = Event_cut_flag% (int)pow(2.,i);
  int quotient = (Event_cut_flag - reste)/ (int)pow(2.,i);
  return quotient;
}

double mmDeltaR(double eta1, double phi1, double eta2, double phi2, double pt1, double pt2, bool option) //Renvoie Dr , si option=true, renvoie Dr sur espace eta/phi/pt
{
  float dr,deta,dphi,dmass;

  deta = eta1 - eta2;

  if(phi1<0.)
    {
      if(phi2>0.)
	dphi=phi2-(phi1+2*myPi);
      else
	dphi=phi2-phi1;
    }		  
  else
    {
      if(phi2>0.)
	dphi=phi2-phi1;
      else
	dphi=phi2+2*myPi-phi1;  
    }
  if(option==1)
    {dmass= (pt1-pt2)/pt1;
      dr=sqrt( deta*deta + dphi*dphi + dmass*dmass );
    }
  else  dr=sqrt( deta*deta + dphi*dphi);
	
  return dr;
}

float ConversionEtaTheta(float eta) //Convertion eta en theta
{
  float theta;
  theta = 2*atan(exp(-eta));
  return theta;
}

float ConversionThetaEta(float theta) //Conversion theta en eta
{
  float eta;
  eta = -log(tan((theta)/2));
  return eta;
}



float ConversionEpt(float energy,float eta ) //Convertion Energie en pt
{

  float pt =  energy*sin(ConversionEtaTheta(eta));
  return pt;
}

float ConversionPtE(float pt, float eta) //Conversion pt en energie
{
  float energy = pt/sin(ConversionEtaTheta(eta));
  return energy;

}

void Conversion_REP_carte(float coord[3], float& x, float& y, float &z) //Conversion Rho/eta/phi en coordonnées cartésiennes
{
  float Theta = ConversionEtaTheta(coord[1]);
  x = coord[0]*sin(Theta)*cos(coord[2]);
  y = coord[0]*sin(Theta)*sin(coord[2]);
  z = coord[0]*cos(Theta);
  return;
}

float ConversionXY_phi(float carte[3]) //renvoie phi sur -pi/pi a partir de coordonnées cartesiennes
{
  float phi;
  if(carte[0]<0.)
    {
      if(carte[1]>0.) //cadran 2
	phi = atan(carte[1]/carte[0]) + myPi;
      else           //cadran 3
	phi = atan(carte[1]/carte[0]) - myPi ;
    }
  else
    phi = atan(carte[1]/carte[0]);
	
  return phi;
}

float Conversion_phiPP_phi2P(float phiPP) //renvoie phi sur 0/2pi à partir de phi sur -pi/pi
{
  if(phiPP<0.)
    {
      return 2*myPi-phiPP;
    }
  else
    return phiPP;

}


void Conversion_carte_REP(float carte[3], float& R, float& E, float& P)//Conversion coordonnées cartesienne dans rho/eta/phi
{
  float Theta;
  if(carte[2]>=0.)
    Theta = abs(atan( sqrt(carte[0]*carte[0] + carte[1]*carte[1])  /carte[2]));
  else
    Theta = myPi - abs(atan( sqrt(carte[0]*carte[0] + carte[1]*carte[1])  /carte[2]));
  if(carte[1]<0)
    {
      if(carte[0]<=0)
	P =myPi + atan(carte[0] / carte[1]);
      else
	P =myPi/2 +  atan(carte[0] / carte[1]);
    }
  else
    {
      if(carte[0]<0)
	P = 2*myPi + atan(carte[0] / carte[1]);
      else
	P = atan(carte[0] / carte[1]);
    }
  
  E = ConversionThetaEta(Theta);
  R = sqrt(carte[0]*carte[0] + carte[1]*carte[1] + carte[2]*carte[2]);
 
  return;
}


float Angle_2vecteurs(float Vect1[3], float Vect2[3],bool cartesien)//Renvoie l'angle formé par deux vecteurs, coordonnées de départ en REP ou Cartesien valide, mais juste à préciser
{
  if(!cartesien)
    {
      Conversion_REP_carte(Vect1,Vect1[0],Vect1[1],Vect1[2]);
      Conversion_REP_carte(Vect2,Vect2[0],Vect2[1],Vect2[2]);

    }
  double Norme1=0,Norme2=0,PS=0;
  float theta;
  Norme1 = sqrt(Vect1[0]*Vect1[0] + Vect1[1]*Vect1[1] +Vect1[2]*Vect1[2] );
  Norme2 = sqrt(Vect2[0]*Vect2[0] + Vect2[1]*Vect2[1] +Vect2[2]*Vect2[2] );
  PS = Vect1[0]*Vect2[0] + Vect1[1]*Vect2[1] +Vect1[2]*Vect2[2];
  
  theta=acos(PS/(Norme1*Norme2));
  
  return theta;
}


#endif
