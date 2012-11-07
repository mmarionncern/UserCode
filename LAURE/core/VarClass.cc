#include "VarClass.hh"



using namespace std;

ClassImp(VarClass)

VarClass::VarClass() {
  _init=true;
   
}



void VarClass::ResetVariables() {

  varmVI.clear();
  varmVUI.clear();
  varmVD.clear();
  varmVF.clear();
  varmVB.clear();
  varmVS.clear();
  
  varmI.clear();
  varmUI.clear();
  varmS.clear();
  varmB.clear();
  varmD.clear();
  varmF.clear();
  
  varmTB.clear();
  
  _varnames.clear();

  _init=true;

}


VarClass::~VarClass() {
 

}



bool
VarClass::isUsefulVar(string name) {
  
  for(size_t i=0;i<_varnames.size();i++) {
    if(name == _varnames[i] )
      {return true;}
  }
  
  return false;

}

void
VarClass::InitVar(string name) {
 
  for(size_t i=0;i<_varnames.size();i++) {
    if(name == _varnames[i] )
      {return;}
  }
  //  cout<<" initializtin "<<name<<endl;
  _varnames.push_back(name);
  return;
}

int
VarClass::getI(string name, int idx) {

  if(_init) {
    InitVar(name);
 
    return 1;
  }
  
  //cout<<"start reading "<<name<<endl;

  itI = varmI.find( name );
  if(itI == varmI.end() ) {
    itVI = varmVI.find( name );
    
    if(itVI == varmVI.end() ) {
      cout<<" error, no such variable "<<name<<endl;
    }
    else {
      //      cout<<" found var "<<name<<"   "<<((*itVI).second)->size()<<endl;
      return (*((*itVI).second))[idx];
    }
  }
  else {
    return (*itI).second;
  }
  return 0;
}

unsigned int
VarClass::getUI(string name, int idx) {

  if(_init) {
    InitVar(name);
    return 1;
  }

  // cout<<"start reading "<<name<<endl;
  
  itUI = varmUI.find( name );
  if(itUI == varmUI.end() ) {
    itVUI = varmVUI.find( name );
    
    if(itVUI == varmVUI.end() ) {
      cout<<" error, no such variable "<<name<<endl;
    }
    else {
      return (*((*itVUI).second))[idx];
    }
  }
  else {
    return (*itUI).second;
  }
  return 0;
}

bool
VarClass::getB(string name, int idx) {

  if(_init) {
    InitVar(name);
    return 1;
  }
  
  // cout<<"start reading "<<name<<endl;

  itB = varmB.find( name );
  if(itB == varmB.end() ) {
    itVB = varmVB.find( name );
    
    if(itVB == varmVB.end() ) {
      cout<<" error, no such variable "<<name<<endl;
    }
    else {
      return (*((*itVB).second))[idx];
    }
  }
  else {
    return (*itB).second;
  }
  return 0;
}

double
VarClass::getD(string name, int idx) {

  if(_init) {
    InitVar(name);
    return 1;
  }
  
  // cout<<"start reading "<<name<<endl;

  itD = varmD.find( name );
  if(itD == varmD.end() ) {
    itVD = varmVD.find( name );
    
    if(itVD == varmVD.end() ) {
      cout<<" error, no such variable "<<name<<endl;
    }
    else {
      return (*((*itVD).second))[idx];
    }
  }
  else {
    return (*itD).second;
  }
  return 0;
}

float
VarClass::getF(string name, int idx) {

  if(_init) {
    InitVar(name);
    return 1;
  }

  //  cout<<"start reading "<<name<<endl;
  
  itF = varmF.find( name );
  if(itF == varmF.end() ) {
    itVF = varmVF.find( name );
    
    if(itVF == varmVF.end() ) {
      cout<<" error, no such variable "<<name<<endl;
    }
    else {
      return (*((*itVF).second))[idx];
    }
  }
  else {
    return (*itF).second;
  }

  return 0;
}

string
VarClass::getS(string name, int idx) {

  if(_init) {
    InitVar(name);
    return "";
  }
  
  // cout<<"start reading "<<name<<endl;

  itS = varmS.find( name );
  if(itS == varmS.end() ) {
    itVS = varmVS.find( name );
    
    if(itVS == varmVS.end() ) {
      cout<<" error, no such variable "<<name<<endl;
    }
    else {
      return (*((*itVS).second))[idx];
    }
  }
  else {
    return (*itS).second;
  }

  return "";

}

unsigned int
VarClass::getSize(string name) {

 if(_init) {
    InitVar(name);
    return 1;
  }
  

  itVS = varmVS.find( name );
  if(itVS != varmVS.end() ) {
    return (*itVS).second->size();
  }
  itVD = varmVD.find( name );
  if(itVD != varmVD.end() ) {
    return (*itVD).second->size();
  }
  itVI = varmVI.find( name );
  if(itVI != varmVI.end() ) {
    return (*itVI).second->size();
  }
  itVF = varmVF.find( name );
  if(itVF != varmVF.end() ) {
    return (*itVF).second->size();
  }
  itVUI = varmVUI.find( name );
  if(itVUI != varmVUI.end() ) {
    return (*itVUI).second->size();
  }
  else{ 
    cout<<"Error for var "<<name<<endl;
    return  0;
  }
}


void 
VarClass::BuildTree(TTree* tree, bool bypass) {

  TObjArray* branchs =  tree->GetListOfBranches();
  string name;
  
  EDataType t;
  TClass* cc;
  string type;
  
  for(int ib=0;ib<branchs->GetEntries();ib++) {
    name = (string)( ((*branchs)[ib])->GetName());
    ((TBranchSTL*)((*branchs)[ib]))->GetExpectedType(cc,t);
    if(t==-1)
      type = (string)(cc->GetName());

    //by default, status disabled
    if(!bypass)
      tree->SetBranchStatus( name.c_str() , 0);
    
    if( isUsefulVar(name )) {
       //Status enabled
      tree->SetBranchStatus( name.c_str() , 1);
      // cout<<name<<"  "<<isUsefulVar(name )<<"   "<<type<<endl;
      //Register branch
      RegisterBranch(tree, name, type, t );
    }

  
  
  }

}


void 
VarClass::RegisterBranch(TTree* tree, string name, string type, EDataType t) {

  //vector or container first
  if(t==-1) {

    if(type=="vector<int>") { //vector<int>

      if( varmVI.find(name) !=varmVI.end() ) {
	cout<<" Warning, "<<name<<" already registered"<<endl;
	return;
      }
      varmVI[ name ] =NULL;
      tree->SetBranchAddress( name.c_str() , &(varmVI[ name ]) );

    }
    if(type=="vector<unsigned int>") { //vector<unsigned int>

      if( varmVUI.find(name) !=varmVUI.end() ) {
	cout<<" Warning, "<<name<<" already registered"<<endl;
	return;
      }
      varmVUI[ name ] =NULL;
      tree->SetBranchAddress( name.c_str() , &(varmVUI[ name ]) );

    }
   if(type=="vector<float>") { //vector<float>

     if( varmVF.find(name) !=varmVF.end() ) {
       cout<<" Warning, "<<name<<" already registered"<<endl;
       return;
     }
     varmVF[ name ] =NULL;
     tree->SetBranchAddress( name.c_str() , &(varmVF[ name ]) );

    }
   if(type=="vector<double>") { //vector<double
     
     if( varmVD.find(name) !=varmVD.end() ) {
       cout<<" Warning, "<<name<<" already registered"<<endl;
       return;
     }
     varmVD[ name ] =NULL;
     tree->SetBranchAddress( name.c_str() , &(varmVD[ name ]) );

    }
   if(type=="vector<bool>") { //vector<bool>
     
     if( varmVB.find(name) !=varmVB.end() ) {
       cout<<" Warning, "<<name<<" already registered"<<endl;
       return;
     }
     varmVB[ name ] =NULL;
     tree->SetBranchAddress( name.c_str() , &(varmVB[ name ]) );
     
    }
   if(type=="string") { //string
     
     if( varmS.find(name) !=varmS.end() ) {
       cout<<" Warning, "<<name<<" already registered"<<endl;
       return;
     }
     varmS[ name ] ="";
     tree->SetBranchAddress( name.c_str() , &(varmS[ name ]) );
     
   }
  if(type=="vector<string>") { //vector<string>
    
    if( varmVS.find(name) !=varmVS.end() ) {
       cout<<" Warning, "<<name<<" already registered"<<endl;
       return;
     }
    varmVS[ name ] =NULL;
    tree->SetBranchAddress( name.c_str() , &(varmVS[ name ]) );
     
   }
   if(type=="TBits") { //TBits
     
     if( varmTB.find(name) !=varmTB.end() ) {
       cout<<" Warning, "<<name<<" already registered"<<endl;
       return;
     }
     varmTB[ name ] =NULL;
     tree->SetBranchAddress( name.c_str() , &(varmTB[ name ]) );
   }

  }
  else if(t==3) { //int
    
    if( varmI.find(name) !=varmI.end() ) {
      cout<<" Warning, "<<name<<" already registered"<<endl;
      return;
    }
    varmI[ name ] =0;
    tree->SetBranchAddress( name.c_str() , &(varmI[ name ]) );

  }
 else if(t==13) { //unsigned int
    
    if( varmUI.find(name) !=varmUI.end() ) {
      cout<<" Warning, "<<name<<" already registered"<<endl;
      return;
    }
    varmUI[ name ] =0;
    tree->SetBranchAddress( name.c_str() , &(varmUI[ name ]) );

  }
  else if(t==5) { //float
   
    if( varmF.find(name) !=varmF.end() ) {
      cout<<" Warning, "<<name<<" already registered"<<endl;
      return;
    }
    varmF[ name ] =0.;
    tree->SetBranchAddress( name.c_str() , &(varmF[ name ]) );

  }
  else if(t==18) { //bool
   
    if( varmB.find(name) !=varmB.end() ) {
      cout<<" Warning, "<<name<<" already registered"<<endl;
      return;
    }
    varmB[ name ] =0;
    tree->SetBranchAddress( name.c_str() , &(varmB[ name ]) );

  }
  else if(t==8) { //double

    if( varmD.find(name) !=varmD.end() ) {
      cout<<" Warning, "<<name<<" already registered"<<endl;
      return;
    }
    varmD[ name ] =0.;
    tree->SetBranchAddress( name.c_str() , &(varmD[ name ]) );

  }

}
