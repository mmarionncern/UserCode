#ifndef __CombUtils
#define __CombUtils


#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <map>

using namespace std;

typedef map<string,int> mapSel;
typedef map<string,int>::iterator mapSelIter;

class CombUtils{

public:
  CombUtils();
  ~CombUtils();

  static bool exist(mapSel a, mapSel b);
  
  vector<mapSel> buildSelections(int Nleptons,
			 string typeLepton,
			 vector<string> sel);

private:
  void f(unsigned prof, int prof_max, 
    std::vector<string> perm,
    std::vector<bool> marked,
    vector<string> objs,
    vector<string>& out,
    bool mark);
  
  void perm(int nl, vector<string> objs,
	    vector<string>& out, bool mark );


};

CombUtils::CombUtils() {
}

CombUtils::~CombUtils() {
}


bool CombUtils::exist(mapSel a, mapSel b) {

  //bool e=false;

  mapSelIter itS1;
  mapSelIter itS2;
  for(itS1=a.begin();itS1!=a.end();itS1++) {

    vector<bool> v;
    itS2 = b.find( (*itS1).first );
    if(itS2 == b.end() ) return false; //single sel does not exist
    else {
      if( itS2->second != itS1->second)
	return false;
    }
  }
  return true;
}


void CombUtils::f(
    unsigned prof,  
    int prof_max, 
    std::vector<string> perm,
    std::vector<bool> marked,
    vector<string> objs,
    vector<string>& out,
    bool mark
){

  if (prof < (int unsigned)prof_max){
      for(unsigned i=0;i<objs.size();++i){
	if (marked[i] && mark) continue;

            std::vector<string> perm_suivant = perm;
            perm_suivant.push_back( objs[i] );
            std::vector<bool> marked_suivant = marked;
	    marked_suivant[i] = true;
            f(prof+1,prof_max,perm_suivant,marked_suivant, objs, out, mark);
	 
        }
    }else{ // final number reached

      out.push_back("");
      for(unsigned i=0;i<perm.size();++i)
	{
	  out[ out.size()-1 ] +=perm[i];
	}

    }
}

void CombUtils::perm(int nl, vector<string> objs, vector<string>& out, bool mark ){
 
  std::vector<string> perm0;
  std::vector<bool> marked0(10,false);
  f(0,nl,perm0,marked0,objs,out, mark);
}

vector<mapSel> CombUtils::buildSelections(int NPart, string typeLepton, vector<string> sel) {

  assert((int unsigned)NPart==sel.size());
  assert(NPart<10); //limiting the combinatory...
  
  
  vector<string> pdg;
  if(typeLepton=="l") {
    pdg.push_back("e");
    pdg.push_back("m");
    pdg.push_back("t");
  }
  if(typeLepton=="L") {
    pdg.push_back("e");
    pdg.push_back("m");
  }
  if(typeLepton=="e") {
    pdg.push_back("e");
  }
  if(typeLepton=="m") {
    pdg.push_back("m");
  }
  if(typeLepton=="t") {
    pdg.push_back("t");
  }
  if(typeLepton=="p") {
    pdg.push_back("p");
  }
  if(typeLepton=="j") {
    pdg.push_back("j");
  }
  if(typeLepton=="h") {
    pdg.push_back("h");
  }

  //int Nlepton=NPart;
  vector<string> combs;
  bool mark =false;
  //Nlepton / possibilities(l or L) / vector to fill out / double comb
  perm(NPart, pdg, combs, mark);

  mark=true;
  vector<string> selections;
 //Nlepton / possibilities(sels) / vector to fill out /  double comb
  perm(NPart, sel, selections, mark);

  
  vector<mapSel> sels;
  mapSelIter itS; 

  for(int unsigned i=0;i<combs.size();i++) {
    for(int unsigned j=0;j<selections.size();j++) {
      mapSel smap; 
      for(int k=0;k<NPart;k++) {
	string t=combs[i].substr(k,1)+selections[j].substr(k*3,3);
	
	itS = smap.find( t );
	if(itS == smap.end() )  
	  smap[ t ]=1;
	else
	  smap[ t ]+=1;
      }
      sels.push_back(smap);
    }
  }
  
  vector<mapSel> cSels;
  for(int unsigned i=0;i<sels.size();i++) {
    bool exists=false;
    for(int unsigned j=0;j<cSels.size();j++) {
      if(exist(cSels[j], sels[i]) )
	{ exists=true; break;}
    }
    if(!exists)
      cSels.push_back(sels[i]);
  }
  
//   for(int i=0;i<cSels.size();i++) {
    
//     cout<<" selection : ";
//     for(itS=cSels[i].begin();itS!=cSels[i].end();itS++)
//       cout<<(*itS).first<<" -> "<<(*itS).second<<"     ";
//     cout<<endl;

//   }

  return cSels;
}




#endif
