import FWCore.ParameterSet.Config as cms

#from PhysicsTools.HepMCCandAlgos.genParticles_cfi import *
#from PhysicsTools.HepMCCandAlgos.genParticleCandidates_cfi import *
#from PhysicsTools.RecoAlgos.allElectrons_cfi import *


process = cms.Process("MmZeeAnalysis")

#Message Logger
#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.INFO.limit = 10
#process.MessageLogger.cerr.threshold = "DEBUG"

# DQM services
process.load("DQMServices.Core.DQM_cfg")
process.DQM.collectorHost = ''

#Geometry diverses
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V5::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

#essai
#process.load('PhysicsTools.HepMCCandAlgos.allElectronsGenParticlesMatch_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
    'file:Zeedata.root')  
   )


#output module for event data
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *_zToEE_*_*','keep *_selectedLayer1Electrons_*_*','keep *_selectedLayer1Photons_*_*','keep *_zToEEMCTruth_*_*','keep *_ElectronsMCTruth_*_*','keep *'),
    fileName = cms.untracked.string('Ini.root')
)


#Electron Analyzer
process.ZeeAnalysis = cms.EDAnalyzer('MmZeeAnalyser',

   #Collections

   electronCollection = cms.untracked.InputTag("selectedLayer1Electrons"),
   ZCandidateCollection = cms.untracked.InputTag("zToEE"),
  # electronGSFCollection = cms.untracked.InputTag("pixelMatchGsfElectrons"),
   genParticleCollection = cms.untracked.InputTag("genParticles"),
   trackCollection = cms.untracked.InputTag("pixelMatchGsfFit"),

   matchingObjectCollection = cms.untracked.InputTag('mergedSuperClusters'),

   #Output File Name
   outputFile = cms.untracked.string('ZeeAnalysisHistos_data.root'),
   #Output Histos
   histograms = cms.untracked.vstring('all'),
)


#Electron MCTruth
process.ElectronsMCTruth = cms.EDProducer("MCTruthDeltaRMatcherNew",
   src = cms.InputTag("pixelMatchGsfElectrons"),
   matched = cms.InputTag("genParticles"),
   distMin = cms.double(0.15),
   matchPDGId = cms.vint32(11)
)

process.load('PhysicsTools.HepMCCandAlgos.allElectronsGenParticlesMatch_cfi')



## #Z to ee MCTruth
process.zToEEMCTruth = cms.EDProducer("MCTruthCompositeMatcherNew",
    src = cms.InputTag("zToEE"),
    matchMaps = cms.VInputTag("electronMatch"),
    matchPDGId = cms.vint32(),                       
  #  matchPDGId = cms.vint32(23)
)

#process.load("ElectroWeakAnalysis.ZReco.zToEEGenParticlesMatch_cfi")

#produce the electron Collection
## process.electronTag = cms.EDFilter("CandViewSelector",
##     src = cms.InputTag('selectecLayer1Electrons'),
##     name = cms.string('ElectronCollection'),
##     cut = cms.string('charge = -1')
## )

## #produce the positron Collection
## process.positronTag = cms.EDFilter("CandViewSelector",
##     src = cms.InputTag('selectecLayer1Electrons'),
##     #name = cms.string('PositronCollection'),
##     cut = cms.string('charge = +1')
## )


# produce Z to e e candidates
process.zToEE = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('selectedLayer1Electrons@+ selectedLayer1Electrons@-'),
    cut = cms.string('30.0 < mass < 20000.0'),
    name = cms.string('zToEE'),
    roles = cms.vstring('positron', 'electron')
)

# require at least one Z to EE candidate
process.compositeFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("zToEE"),
    minNumber = cms.uint32(1)
)

process.electronPtFilter = cms.EDFilter("PtMinGsfElectronCountFilter",
    src = cms.InputTag("pixelMatchGsfElectrons"),
    minNumber = cms.uint32(0),
    ptMin = cms.double(10.0)                                    
)

# PAT Layer 1
process.load("PhysicsTools.PatAlgos.patLayer0_cff") # need to load this
process.load("PhysicsTools.PatAlgos.patLayer1_cff") # even if we run only layer 1


#Filter on SuperCluster
process.mergedSuperClusters = cms.EDFilter("EgammaSuperClusterMerger",
    src = cms.VInputTag(cms.InputTag("correctedHybridSuperClusters"),
                        cms.InputTag("multi5x5SuperClustersWithPreshower"))
)

#process.outpath = cms.EndPath(process.o1)

process.p = cms.Path(process.electronPtFilter*
    process.patLayer0*
    process.patLayer1* 
    process.zToEE*
    process.mergedSuperClusters*
    process.ZeeAnalysis
    )


