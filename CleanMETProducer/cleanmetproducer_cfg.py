import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patTemplate_cfg import *

##
## process name is already given in patTemplate_cfg
##


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/m/mmarionn/DYToEESummer11.root' 
    ),
#    skipEvents = cms.untracked.uint32(0)                        
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )


##
## To remove the "begin job processing line " from printing
##
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('cout'),
    cout         = cms.untracked.PSet(threshold = cms.untracked.string('WARNING'))
)


##
## load the producer and configure it (if needed)
##
process.load("MMarionneau.CleanMETProducer.cleanMETProducer_cfi")

##
## load the global tag
##

from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *
process.GlobalTag.globaltag = cms.string('START42_V13::All')

##
## Simple example of objects needed for vertex recognization, here Zee
##

process.goodElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("selectedPatElectrons"),
    cut = cms.string("pt > 0"),
)

process.ZeeCandidates = cms.EDProducer("CandViewShallowClonePtrCombiner",
    cut = cms.string('30.0 < mass < 1200'),
    decay = cms.string('selectedPatElectrons@+ selectedPatElectrons@-'),
    filter = cms.bool(True)                            
)

# require at least one Z to EE candidate
process.ZeeFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("ZeeCandidates"),
    minNumber = cms.uint32(1)
)


##
## Configure PAT : remove MC matching
##
from PhysicsTools.PatAlgos.tools.coreTools import *
# turn off MC matching for the process
removeMCMatching(process, ['All'])


##
## Jet configuration (get corrected jet for Type I MET computation )
##
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *
switchToPFJets(process, cms.InputTag('pfJets'), algo='AK5', postfix ="", jetCorrections=('AK5PFchs',['L2Relative', 'L3Absolute', 'Uncertainty'] ))

## If data, add the 'L2L3Residual' correction

##
## process path
##
process.p = cms.Path(
   
    # PAT, needed in the current implementation of typeI MET
    # and in the current example, it can be changed by the user
    # but needs to change the jet format in the CleanMETProducer
    process.PF2PAT*
    process.patDefaultSequence*
    process.metJESCorPFAK5*

    ### Configurable part of the analysis
      process.goodElectrons
    * process.ZeeCandidates
    * process.ZeeFilter

    ### The producer
    * process.cleanMETProducer
)
  


##
## End Path
##

process.out = cms.OutputModule("PoolOutputModule",
 #    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('PAT')),                           
     fileName = cms.untracked.string('outFile.root')
)
process.e = cms.EndPath(process.out)
