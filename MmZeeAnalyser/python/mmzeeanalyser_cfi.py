import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('MmZeeAnalyser'
    ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
