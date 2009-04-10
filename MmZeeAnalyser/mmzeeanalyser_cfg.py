import FWCore.ParameterSet.Config as cms

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
 #   '/store/mc/Summer08/Zmumu/GEN-SIM-RECO/IDEAL_V9_v1/0004/FC004934-C888-DD11-A240-0015C5E9B30F.root'
    '/store/mc/2007/11/7/CSA07-DYee_M40_filter-1194479189/0013/0A3FC9C0-2A97-DC11-BEB5-001BFCDBD130.root'
 #   '/store/mc/Summer08/PhotonJetPt30/GEN-SIM-RECO/IDEAL_V9_v1/0007/0088F7BC-7EDD-DD11-BDC7-0030489453EC.root'
    )  
   )


#output module for event data
process.o1 = cms.OutputModule("PoolOutputModule",
      outputCommands = cms.untracked.vstring('keep *_zToEE_*_*',
                                             'keep *_selectedLayer1Electrons_*_*',
                                             'keep *_selectedLayer1Photons_*_*',
                                             'keep *_zToMuMu_*_*',
                                             'keep *_selectedLayer1Muons_*_*'),
    fileName = cms.untracked.string('Ini.root')
)


#Electron Analyzer
process.ZeeAnalysis = cms.EDAnalyzer('MmZeeAnalyser',

   analyseType = cms.untracked.string("TDR"),

   #Collections
   electronCollection = cms.untracked.InputTag("selectedLayer1Electrons"),
   muonCollection = cms.untracked.InputTag("selectedLayer1Muons"),                                 
   ZeeCandidateCollection = cms.untracked.InputTag("zToEE"),
   ZmumuCandidateCollection = cms.untracked.InputTag("zToMuMu"),                                

   genParticleCollection = cms.untracked.InputTag("genParticles"),
   trackCollection = cms.untracked.InputTag("pixelMatchGsfFit"),

   matchingObjectCollection = cms.untracked.InputTag('mergedSuperClusters'),

   #Output File Name
   outputFile = cms.untracked.string('ZeeAnalysisHistos_data.root'),
   #Output Histos
   histograms = cms.untracked.vstring('all'),
)

process.outpath = cms.EndPath(process.o1)

process.load("MMarionneau.MmZeeAnalyser.mmPatProducer_cfi")
#process.load("MMarionneau.MmZeeAnalyser.mmPatProducer_Zmumu_cfi")
process.load("MMarionneau.LeptonFilter.LeptonFilter_cfi")

process.p = cms.Path( process.LeptonFilter*
                      process.mmPatProducerZee*
                     (process.zToEE + process.zToMuMu)
                     *process.mergedSuperClusters
                     *process.ZeeAnalysis)


