#ifndef __VCLQ_
#define __VCLQ_

#include <iostream>
#include <vector>
#include <map>
#include <string>
//#include <EDataType.h>

#include "core/VarClass.hh"

#include <TROOT.h>
#include <TLorentzVector.h>
#include <TBranch.h>
#include <TBranchSTL.h>
#include <TClass.h>

#include <TObjArray.h>
#include <TTree.h>
#include <TBits.h>
#include <TChain.h>

using namespace std;

class VarClassLQ : public VarClass
{

public:

  VarClassLQ();
  ~VarClassLQ();
 
public:

  //4-Vectors
  TLorentzVector tau;
  TLorentzVector lep;
  TLorentzVector jet1;
  TLorentzVector jet2;
    
  TLorentzVector dilep;
  TLorentzVector dijet;
  TLorentzVector tjet1;
  TLorentzVector ljet1;
  TLorentzVector tjet2;
  TLorentzVector ljet2;

  TLorentzVector tljet1;
  TLorentzVector tljet2;
  TLorentzVector tjj;
  TLorentzVector ljj;

  TLorentzVector tljj;

  int jetLead;
  int jetTrail;
  int bJetLead;
  int bJetTrail;

  bool noBTag;

  int tauL;
  int lepL;

  float zpv;
  vector<size_t> v_electrons_vbtf95;
  vector<size_t> v_electrons;


  ClassDef(VarClassLQ,0)
};

#endif
