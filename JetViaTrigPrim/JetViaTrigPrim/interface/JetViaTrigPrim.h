
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "DataFormats/EcalDetId/interface/EcalScDetId.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EBSrFlag.h"
#include "DataFormats/EcalDigi/interface/EESrFlag.h"
#include "DataFormats/EcalDigi/interface/EcalSrFlag.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"


// Geometry files
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"


#include "Validation/EcalDigis/src/CollHandle.h"

#include <vector>
#include <string>
#include <set>
#include <utility>

//
// class declaration
//

class JetViaTrigPrim : public edm::EDAnalyzer {
   public:
      explicit JetViaTrigPrim(const edm::ParameterSet&);
      ~JetViaTrigPrim();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;


      const std::vector< DetId > RetrieveCrystals(const EcalTrigTowerDetId& ttid );
      const EcalTrigPrimDigiCollection* getTrigPrims(const edm::Event& event);
                
      void readAllCollections(const edm::Event& event); 
      
      // ----------member data ---------------------------


      bool collNotFoundWarn_;
      
      CollHandle<EBSrFlagCollection>         ebSrFlags_;
      CollHandle<EESrFlagCollection>         eeSrFlags_;   
      CollHandle<EBDigiCollection>         ebDigis_;   
 
      std::string outputFile_;
      std::string trigPrimCollection_;
      std::string trigPrimProducer_;


      /** ECAL trigger tower mapping
       */
      const EcalTrigTowerConstituentsMap * triggerTowerMap_;

   

};
