#ifndef __VC_
#define __VC_

#include <iostream>
#include <vector>
#include <map>
#include <string>
//#include <EDataType.h>

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

typedef map<string,vector<int>* > mapVI;
typedef map<string,vector<unsigned int>* > mapVUI;
typedef map<string,vector<double>* > mapVD;
typedef map<string,vector<float>* > mapVF;
typedef map<string,vector<string>* > mapVS;
typedef map<string,vector<bool>* > mapVB;


typedef map<string,vector<int>* >::iterator itMapVI;
typedef map<string,vector<unsigned int>* >::iterator itMapVUI;
typedef map<string,vector<double>* >::iterator itMapVD;
typedef map<string,vector<float>* >::iterator itMapVF;
typedef map<string,vector<string>* >::iterator itMapVS;
typedef map<string,vector<bool>* >::iterator itMapVB;

typedef map<string,int > mapI;
typedef map<string,unsigned int > mapUI;
typedef map<string,double > mapD;
typedef map<string,float > mapF;
typedef map<string,bool > mapB;
typedef map<string,string > mapS;

typedef map<string,TBits* > mapTB;

typedef map<string, int >::iterator itMapI;
typedef map<string,unsigned int >::iterator itMapUI;
typedef map<string,double >::iterator itMapD;
typedef map<string,float >::iterator itMapF;
typedef map<string,bool >::iterator itMapB;
typedef map<string,string >::iterator itMapS;

typedef map<string,TBits* >::iterator itMapTB;

class VarClass
{

public:

  VarClass();
  ~VarClass();

  void ResetVariables();

  // void AssociateVariables();
  // void AssociateAddresses(TChain* tree);

  // void InitMaps();

  bool isUsefulVar(string name);
  void InitVar(string name);

  void FinalizeInit() {_init= false; };

  int getI(string name, int idx=0);
  unsigned int getUI(string name, int idx=0);
  bool getB(string name, int idx=0);
  double getD(string name, int idx=0);
  float getF(string name, int idx=0);
  string getS(string name, int idx=0);

  unsigned int getSize(string name);

  void BuildTree(TTree* tree, bool bypass);
  void RegisterBranch(TTree* tree, string name, string type, EDataType t);
 
private:

  //All variables
  mapVI varmVI;
  mapVUI varmVUI;
  mapVD varmVD;
  mapVF varmVF;
  mapVB varmVB;
  mapVS varmVS;
  

  mapI varmI;
  mapUI varmUI;
  mapS varmS;
  mapB varmB;
  mapD varmD;
  mapF varmF;

  mapTB varmTB;

  itMapVI itVI;
  itMapVUI itVUI;
  itMapVD itVD;
  itMapVF itVF;
  itMapVB itVB;
  itMapVS itVS;
  itMapI itI;
  itMapUI itUI;
  itMapS itS;
  itMapB itB;
  itMapD itD;
  itMapF itF;
  itMapTB itTB;

private:
  //list of variables
  vector<string> _varnames;

protected:
  bool _init;

  ClassDef(VarClass,0)
};

#endif
