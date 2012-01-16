import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patTemplate_cfg import *

##
## process name is already given in patTemplate_cfg
##


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/tmp/mmarionn/ZmmFall11_44.root' 
    ),
#    skipEvents = cms.untracked.uint32(0)                        
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
    )


##
## To remove the "begin job processing line " from printing
##
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger = cms.Service("MessageLogger",
#    destinations = cms.untracked.vstring('cout'),
#    cout         = cms.untracked.PSet(threshold = cms.untracked.string('WARNING'))
#)


##
## load the producer and configure it (if needed)
##
process.load("MMarionneau.CleanMETProducer.cleanMETProducer_cfi")

##
## load the global tag
##

from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *
process.GlobalTag.globaltag = cms.string('START44_V12::All')

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
## Produce the pfMET Type I corrected
##
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

## Tune for MET correction : data/MC --> ak5PFL1FastL2L3Residual // ak5PFL1FastL2L3
process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")


##
## Jet configuration (get corrected jet for Type I MET computation )
##
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process, cms.InputTag('ak5PFJets'),
                    doJTA            = False,            
                    doBTagging       = False,            
                    jetCorrLabel     =  ('AK5PF',['L1FastJet','L2Relative', 'L3Absolute','Uncertainty'] ), #L2L3Residual
                    doType1MET       = False,            
                    genJetCollection = cms.InputTag("ak5GenJets"),
                    doJetID      = False,
                    jetIdLabel   = "ak5"
                    )


## If data, add the 'L2L3Residual' correction

##-------------------- Turn-on the FastJet density calculation into PAT -----------------------
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.kt6PFJets.doRhoFastjet = True
process.patJetCorrFactors.useRho = cms.bool(True)
process.patJetCorrFactors.rho = cms.InputTag("kt6PFJets:rho:PAT")
process.patJetCorrFactors.useNPV = cms.bool(False)


##
## process path
##
process.p = cms.Path(
     # Get Rho  correction for L1FastJet
    process.kt6PFJets*
    # PAT, used in the current example to get the candidates used
    # in the producer (except pFT1MET), it can be changed by the user
    # but needs to change the jet format in the CleanMETProducer
    process.patDefaultSequence*
    # create the pfMET typeI corrected
    process.producePFMETCorrections*

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
