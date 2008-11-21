import FWCore.ParameterSet.Config as cms

process = cms.Process("MmZeeAnalysis")


#Message Logger
#process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(#'/store/mc/2007/12/19/ReRecoIdeal-Zee-1198082306/0002/00B862E1-8AAE-DC11-B2BE-001617C3B6AA.root',

                 #'/store/relval/CMSSW_2_1_9/RelValSingleElectronPt10/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/1EB913B8-B185-DD11-83FE-000423D992DC.root')
'/store/relval/CMSSW_2_1_9/RelValSingleElectronPt10/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/EC56D7D8-B185-DD11-889A-000423D98EA8.root')
                            )

#output module for event data
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('Ini.root')
)


#Electron Analyzer
process.gsfElectronAnalysis = cms.EDAnalyzer('MmZeeAnalyser',

   #Collections
   electronCollection_ = cms.untracked.InputTag("recoGsfElectrons","pixelMatchGsfElectrons"),
   matchingObjectCollection = cms.untracked.InputTag('mergedSuperClusters'),

   #Output File Name
   outputFile = cms.string('gsfElectronHistos_data.root'),
                                             
   #Variables and Cuts Initialization
##     Nbinxyz = cms.int32(50),
##     Nbineop2D = cms.int32(30),
##     Nbinp = cms.int32(50),
##     Lhitsmax = cms.double(10.0),
##     Etamin = cms.double(-2.5),
##     Nbinfhits = cms.int32(20),
##     Eopmax = cms.double(5.0),
##     Pmax = cms.double(300.0),
##     Phimax = cms.double(3.2),
##     Phimin = cms.double(-3.2),
##     Dphimin = cms.double(-0.15),
##     MaxPt = cms.double(100.0),
##     Nbinlhits = cms.int32(5),
##     Nbinpteff = cms.int32(19),
##     Nbinphi2D = cms.int32(32),
##     Nbindetamatch2D = cms.int32(50),
##     Nbineta = cms.int32(50),
##     DeltaR = cms.double(0.3),
  
##     Nbinp2D = cms.int32(50),
##     Nbindeta = cms.int32(100),
##     Nbinpt2D = cms.int32(50),
##     Nbindetamatch = cms.int32(100),
##     Fhitsmax = cms.double(20.0),
##     Nbinphi = cms.int32(64),
##     Nbineta2D = cms.int32(50),
##     Eopmaxsht = cms.double(3.0),
##     MaxAbsEta = cms.double(2.5),
##     Nbindphimatch = cms.int32(100),
##     Detamax = cms.double(0.15),
##     Nbinpt = cms.int32(50),
##     Nbindphimatch2D = cms.int32(50),
##     Etamax = cms.double(2.5),
##     Dphimax = cms.double(0.15),
##     Dphimatchmax = cms.double(0.2),
##     Detamatchmax = cms.double(0.05),
##     Nbindphi = cms.int32(100),
##     Detamatchmin = cms.double(-0.05),
##     Ptmax = cms.double(100.0),
##     Nbineop = cms.int32(50),
##     Dphimatchmin = cms.double(-0.2),
##     Detamin = cms.double(-0.15)
                                             
                                             
)

#Filter on SuperCluster
process.mergedSuperClusters = cms.EDFilter("EgammaSuperClusterMerger",
    src = cms.VInputTag(cms.InputTag("correctedHybridSuperClusters"),
                        cms.InputTag("multi5x5SuperClustersWithPreshower"))
)

#process.outpath = cms.EndPath(process.o1)

process.p = cms.Path(process.mergedSuperClusters*process.gsfElectronAnalysis)
