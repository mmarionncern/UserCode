import FWCore.ParameterSet.Config as cms

jetViaTrigPrim = cms.EDAnalyzer("JetViaTrigPrim",

                #Collections
EbSrFlagCollection = cms.InputTag("simEcalDigis","ebSrFlags"),
EeSrFlagCollection = cms.InputTag("simEcalDigis","eeSrFlags"),
trigPrimProducer = cms.string('ecalTriggerPrimitiveDigis'),
trigPrimCollection = cms.string(''),
EBDigisCollection = cms.InputTag("simEcalDigis","ebDigis"),


outputFile = cms.untracked.string('')


)
