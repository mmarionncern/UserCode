





void
VarClassLQ::AssociateAddresses(TChain* tree) {

  vector<size_t> off;
  //off.push_back(0);
  off.push_back( varmI.size() );
  off.push_back( varmS.size() );
  off.push_back( varmB.size() );
  off.push_back( varmD.size() );
  off.push_back( varmVI.size() );
  off.push_back( varmVS.size() );
  off.push_back( varmVB.size() );
  off.push_back( varmVD.size() );

  size_t iv=0;
  
  //for(size_t iv=0;iv<branchs.size();iv++) {
    for(itI = varmI.begin(); itI != varmI.end(); itI++) { 
      tree->SetBranchStatus( ((*itI).first).c_str() , 1);
      tree->SetBranchAddress( ((*itI).first).c_str() , &((*itI).second) );
      iv++;
    }
    
   for(itS = varmS.begin(); itS != varmS.end(); itS++) { 
     tree->SetBranchStatus( ((*itS).first).c_str() , 1);
     tree->SetBranchAddress( ((*itS).first).c_str() , &((*itS).second) );
      iv++;
    }
   
   for(itB = varmB.begin(); itB != varmB.end(); itB++) { 
     tree->SetBranchStatus( ((*itB).first).c_str() , 1);
     tree->SetBranchAddress( ((*itB).first).c_str() , &((*itB).second) );
     iv++;
   }
    
   for(itD = varmD.begin(); itD != varmD.end(); itD++) { 
     tree->SetBranchStatus( ((*itD).first).c_str() , 1);
     tree->SetBranchAddress( ((*itD).first).c_str() , &((*itD).second) );
     iv++;
   }
   

  for(itVI = varmVI.begin(); itVI != varmVI.end(); itVI++) { 
    tree->SetBranchStatus( ((*itVI).first).c_str() , 1);
    ((*itVI).second) = NULL;
      tree->SetBranchAddress( ((*itVI).first).c_str() , &((*itVI).second) );
      iv++;
    }
    
   for(itVS = varmVS.begin(); itVS != varmVS.end(); itVS++) { 
     tree->SetBranchStatus( ((*itVS).first).c_str() , 1);
     ((*itVS).second) = NULL;
      tree->SetBranchAddress( ((*itVS).first).c_str() , &((*itVS).second) );
      iv++;
    }
   
   for(itVB = varmVB.begin(); itVB != varmVB.end(); itVB++) { 
     tree->SetBranchStatus( ((*itVB).first).c_str() , 1);
     ((*itVB).second) = NULL;
     tree->SetBranchAddress( ((*itVB).first).c_str() , &((*itVB).second) );
     iv++;
   }
    
   for(itVD = varmVD.begin(); itVD != varmVD.end(); itVD++) { 
     tree->SetBranchStatus( ((*itVD).first).c_str() , 1);
     ((*itVD).second) = NULL;
     tree->SetBranchAddress( ((*itVD).first).c_str() , &((*itVD).second) );
     iv++;
   }
   
}
