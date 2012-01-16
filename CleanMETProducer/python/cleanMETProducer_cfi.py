import FWCore.ParameterSet.Config as cms

cleanMETProducer = cms.EDProducer('CleanMETProducer',
   pfMETInput    = cms.InputTag("pfType1CorrectedMet"),                     
   pfJETInput = cms.InputTag("selectedPatJets"),
   pfCandidateInput = cms.InputTag("particleFlow"),
   vertexInput = cms.InputTag("offlinePrimaryVertices"),
   trackInput = cms.InputTag("generalTracks"),                               
   vtxObjectInput = cms.InputTag("ZeeCandidates"),
   refObjectInput = cms.InputTag("ZeeCandidates"),

   nRefObjectInCol = cms.untracked.int32(0),
   nVtxObjectInCol = cms.untracked.int32(0),
   vtxUserDef = cms.untracked.bool(True),

   nameCleanMET = cms.string("cleanPUMET"),
   nameMinMET   = cms.string("minCPUMET"),

### Thresholds and protections
    jetPtHSThr  = cms.untracked.double(25.),
    jetBalHSThr = cms.untracked.double(0.7),
    jetDirHSThr = cms.untracked.double(2.8),                          

    jetMatchDR  = cms.untracked.double(0.1),
    jetMatchBal = cms.untracked.double(0.8), 
                                  
)
