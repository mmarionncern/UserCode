import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Geometry_cff import *
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *
from Configuration.StandardSequences.MagneticField_cff import *

from PhysicsTools.PatAlgos.patLayer0_cff import *
from PhysicsTools.PatAlgos.patLayer1_cff import *

#Geometry diverses
#load("Configuration.StandardSequences.Geometry_cff")
#load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
GlobalTag.globaltag = cms.string('IDEAL_V5::All')
#load("Configuration.StandardSequences.MagneticField_cff")



# produce Z to e e candidates
zToEE = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('selectedLayer1Electrons@+ selectedLayer1Electrons@-'),
    cut = cms.string('30.0 < mass < 20000.0'),
    name = cms.string('zToEE'),
    roles = cms.vstring('positron', 'electron')
)

# require at least one Z to EE candidate
compositeFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("zToEE"),
    minNumber = cms.uint32(1)
)

electronPtFilter = cms.EDFilter("PtMinGsfElectronCountFilter",
    src = cms.InputTag("pixelMatchGsfElectrons"),
    minNumber = cms.uint32(0),
    ptMin = cms.double(10.0)                                    
)

# PAT Layer 1
#load("PhysicsTools.PatAlgos.patLayer0_cff") # need to load this
#load("PhysicsTools.PatAlgos.patLayer1_cff") # even if we run only layer 1


#Filter on SuperCluster
mergedSuperClusters = cms.EDFilter("EgammaSuperClusterMerger",
    src = cms.VInputTag(cms.InputTag("correctedHybridSuperClusters"),
                        cms.InputTag("multi5x5SuperClustersWithPreshower"))
)



mmPatProducerZee = cms.Sequence(electronPtFilter*
                                patLayer0*
                                patLayer1* 
                                zToEE*
                                mergedSuperClusters)
