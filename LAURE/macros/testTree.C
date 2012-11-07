{

  TFile *_file0 = TFile::Open("../data/LeptoQuark/skim_TTbar.root");
  // ZZ_2m2n->cd();

  // cout<<ZZtuple->GetEntries()<<endl;  

  EDataType t;
  TClass* cc;

  TTree* tree=(TTree*)_file0->Get("tree");
  cout<<" tree "<<tree<<endl;
  TObjArray* br = tree->GetListOfBranches();
  cout<<br<<endl;
  int n = br->GetEntries();
  string name;
  //  cout<<<<endl;
  for(int ib=0;ib<n;ib++) {
    // name = (string)( ((*br)[ib])->GetName());
    ((TBranchSTL*)((*br)[ib]))->GetExpectedType(cc,t);
    cout<<((*br)[ib])->GetName()<<"   "<<t<<"   "<<cc<<endl;
    if(t==-1)
      cout<<cc->GetName()<<endl;
  }




}
