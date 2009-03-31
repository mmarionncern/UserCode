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

   analyseType = cms.untracked.string("EWK"),

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

#process.outpath = cms.EndPath(process.o1)

process.load("MMarionneau.MmZeeAnalyser.mmPatProducer_cfi")


process.p = cms.Path(process.mmPatProducerZee*process.ZeeAnalysis)


