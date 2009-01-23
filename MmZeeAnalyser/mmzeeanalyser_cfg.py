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
        fileNames = cms.untracked.vstring(#'/store/mc/2007/12/19/ReRecoIdeal-Zee-1198082306/0002/00B862E1-8AAE-DC11-B2BE-001617C3B6AA.root',

              #  '/store/relval/CMSSW_2_1_9/RelValSingleElectronPt10/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/1EB913B8-B185-DD11-83FE-000423D992DC.root'                 )
'/store/mc/Summer08/Zee/GEN-SIM-RECO/IDEAL_V9_v1/0004/F63F0574-9888-DD11-B84A-001F29084160.root',
    '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/0A28F869-F285-DD11-AF3C-001617DBD5B2.root',
    '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/162C4B5E-F585-DD11-872A-001617C3B64C.root',
    '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/205E6CE3-F485-DD11-9D53-001617C3B76A.root',
    '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/562BAFA1-F585-DD11-B931-001617DBD224.root',
    '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/565AFE10-EF85-DD11-8353-000423D6B42C.root',
    '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/5C66302A-F185-DD11-81D3-000423D98834.root',
    '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/66F60641-F685-DD11-A493-000423D987FC.root',
    '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/6E6A6E2D-F485-DD11-B707-001617DBD472.root')
#'/store/mc/CSA08/Zee/GEN-SIM-RECO/CSA08_S156_v1/0002/165E4BA2-AC2B-DD11-8163-001A644EB7CE.root')
   # 'file:Ini.root')  
                            )

#output module for event data
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *_zToEE_*_*','keep *_selectedLayer1Electrons_*_*','keep *_selectedLayer1Photons_*_*','keep *_zToEEMCTruth_*_*','keep *_ElectronsMCTruth_*_*','keep *'),
    fileName = cms.untracked.string('Ini.root')
)


#Electron Analyzer
process.ZeeAnalysis = cms.EDAnalyzer('MmZeeAnalyser',

   #Collections
  # electronCollection = cms.untracked.InputTag("selectedLayer1Electrons"),
   ZCandidateCollection = cms.untracked.InputTag("zToEE"),
  # ElectronMCTruthMap = cms.untracked.InputTag("ElectronMCTruth"),
  # ZeeMCTruthCollection = cms.untracked.InputTag(""),
                                  
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

# Zpt Filter
#process.compositeFilter = cms.EDFilter("CandViewCountFilter",
#    src = cms.InputTag("zToEE"),
#    minNumber = cms.uint32(1)
#)


#PAT Matcher
## process.matchElectrons = cms.EDProducer("PatMCMatcher",
##     src = cms.InputTag("selectedLayer1Electrons"),
##     matched = cms.InputTag("genParticles"),
##     maxDeltaR = cms.double( 0.15 ),
##     macDPTRel = cms.double( 0.5 ),                                    
##     mcPdgId = cms.vint32(11),
##     mcStatus = cms.vint32(1)
## )

# PAT Layer 1
process.load("PhysicsTools.PatAlgos.patLayer0_cff") # need to load this
process.load("PhysicsTools.PatAlgos.patLayer1_cff") # even if we run only layer 1


#Filter on SuperCluster
process.mergedSuperClusters = cms.EDFilter("EgammaSuperClusterMerger",
    src = cms.VInputTag(cms.InputTag("correctedHybridSuperClusters"),
                        cms.InputTag("multi5x5SuperClustersWithPreshower"))
)

#process.outpath = cms.EndPath(process.o1)

process.p = cms.Path(
    process.patLayer0*
    process.patLayer1* 
    process.zToEE*
    process.compositeFilter*
    process.mergedSuperClusters*
    process.ZeeAnalysis
    )
## ## #process.mergedSuperClusters*
