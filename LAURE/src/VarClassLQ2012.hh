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

class VarClassLQ2012 : public VarClass
{

public:

  VarClassLQ2012();
  ~VarClassLQ2012();
 
public:

  TLorentzVector el;
  TLorentzVector tau;
  

  ClassDef(VarClassLQ2012,0)
};

#endif
