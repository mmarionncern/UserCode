import FWCore.ParameterSet.Config as cms

cleanMETProducer = cms.EDProducer('CleanMETProducer',
   pfMETInput    = cms.InputTag("metJESCorPFAK5"),                     
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


##
## Missing ET Type I Correction
##

from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5PFJet
metJESCorPFAK5 = metJESCorAK5PFJet.clone()
metJESCorPFAK5.inputUncorJetsLabel = "pfJets"
metJESCorPFAK5.metType = "PFMET"
metJESCorPFAK5.inputUncorMetLabel = "pfMet"
metJESCorPFAK5.useTypeII = False
metJESCorPFAK5.jetPTthreshold = cms.double(10.0)

#######################################################


###
### The PF2PAT section to get the good pfJet collection to make the pfTypeI collection
### Needed here for tricky reasons

from CommonTools.ParticleFlow.PF2PAT_cff import *
