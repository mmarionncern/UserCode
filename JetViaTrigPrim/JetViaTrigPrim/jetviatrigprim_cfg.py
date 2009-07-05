import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
         '/../user/m/mmarionn/QCD_pt100/QCD_pt_100_highlum_Prod1_1.root',
    )
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'IDEAL_V11::All' 


#output module for event data
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('sortie.root')
)

# MessageLogger:
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.INFO.limit = 10
process.MessageLogger.cerr.threshold = "DEBUG"
#process.MessageLogger.suppressWarning = ['ecalSelectiveReadoutValidation']

# ECAL Geometry:
process.load("Geometry.EcalCommonData.EcalOnly_cfi")
#  Calo geometry service model
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")

#Magnetic field:
process.load("Configuration.StandardSequences.MagneticField_cff")


process.load("CalibCalorimetry.Configuration.Ecal_FakeConditions_cff")
process.load("CalibCalorimetry.EcalTrivialCondModules.EcalTrivialCondRetriever_cfi")

#*****************************

#SR Emulation
process.simEcalDigis = cms.EDProducer("EcalSelectiveReadoutProducer",

    # ADC to GeV conversion factor used in ZS filter for EB
    ebDccAdcToGeV = cms.double(0.035),    
    # ADC to GeV conversion factor used in ZS filter for EE
    eeDccAdcToGeV = cms.double(0.06),

    # Instance name of output EB digis collection
    EBSRPdigiCollection = cms.string('ebDigis'),
    # Label of input EB and EE  collections
    digiProducer = cms.string('simEcalUnsuppressedDigis'),
    # Instance name of output EE SR flags collection
    EESrFlagCollection = cms.string('eeSrFlags'),
    # Instance name of input EE digi collections
    EEdigiCollection = cms.string('eeDigis'),
    # Instance name of input EB digi collections
    EBdigiCollection = cms.string('ebDigis'),
    # Instance name of output EB SR flags collection
    EBSrFlagCollection = cms.string('ebSrFlags'),
    # Label name of input ECAL trigger primitive collection
    trigPrimProducer = cms.string('ecalTriggerPrimitiveDigis'),
    # Instance name of output EE digis collection
    EESRPdigiCollection = cms.string('eeDigis'),
                                      
    #logical flag to write out SrFlags
    writeSrFlags = cms.untracked.bool(True),
     #logical flag to write out suppressed digis                                  
    produceDigis = cms.untracked.bool(True),
    
    # Neighbour eta range, neighborhood: (2*deltaEta+1)*(2*deltaPhi+1)
    deltaEta = cms.int32(1),
    # Neighbouring eta range, neighborhood: (2*deltaEta+1)*(2*deltaPhi+1)
    deltaPhi = cms.int32(1),
    
    # Instance name of ECAL trigger primitive collection
    trigPrimCollection = cms.string(''),

    # ZS energy threshold in GeV to apply to low interest channels of endcap
    srpEndcapLowInterestChannelZS = cms.double(0.480),
    # ZS energy threshold in GeV to apply to high interest channels of endcap
    srpEndcapHighInterestChannelZS = cms.double( -50000.480),
    # ZS energy threshold in GeV to apply to high interest channels of barrel
    srpBarrelHighInterestChannelZS = cms.double( -50000.1225 ),
    # ZS energy threshold in GeV to apply to low interest channels of barrel
    srpBarrelLowInterestChannelZS = cms.double( 0.1225),


    #for debug mode only
    trigPrimBypassWithPeakFinder = cms.bool(False),
    # Index of time sample (staring from 1) the first DCC weights is implied
    ecalDccZs1stSample = cms.int32(2),
    #DCC ZS FIR weights: weights are rounded in such way that in Hw
    #representation (weigth*1024 rounded to nearest integer) the sum is null:
    dccNormalizedWeights = cms.vdouble(-0.374, -0.374, -0.3629, 0.2721, 0.4681, 0.3707),
    #dccNormalizedWeights = cms.vdouble(-1.1865,0.0195,0.2900,0.3476,0.3008,0.2266),
    #dccNormalizedWeights = cms.vdouble(-0.3464,-0.1934,-0.6018,0.4012,0.3921,0.3483),
    # dccNormalizedWeights = cms.vdouble(-0.334,-0.333,-0.333,0.,1.,0.),

    #dccNormalizedWeights = cms.vdouble( 0.0682,0.0203,-1.026,-0.1416,1.4527,-0.374),

    #symetricZS
    symetricZS = cms.bool(False),
                                      
    #number of events whose TT and SR flags must be dumped (for debug purpose):
    dumpFlags = cms.untracked.int32(0),
    
    
    #switch to run w/o trigger primitive. For debug use only
    trigPrimBypass = cms.bool(False),
    #for debug mode only:
    trigPrimBypassHTH = cms.double(1.5),
    #for debug mode only:
    trigPrimBypassLTH = cms.double(1.5),
                                      
    #value to substitute if a TTF is missing:
    defaultTtf = cms.int32(6),
    # switch to enable/disable EB processing:
    ebOn = cms.bool(True),
    # switch to enable/disable EE processing:
    eeOn = cms.bool(True)

                                      
)


#***********************
# Ecal digi production:
process.load("SimCalorimetry.EcalSimProducers.ecaldigi_cfi")

process.simEcalUnsuppressedDigis.EBdigiCollection = 'ebDigis'
process.simEcalUnsuppressedDigis.EEdigiCollection = 'eeDigis'
process.simEcalUnsuppressedDigis.ESdigiCollection = 'esDigis'


process.tpparams = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGLinearizationConstRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.tpparams2 = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGPedestalsRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.tpparams3 = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGSlidingWindowRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.tpparams4 = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGWeightIdMapRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.tpparams5 = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGWeightGroupRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.tpparams6 = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGLutGroupRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.tpparams7 = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGLutIdMapRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.tpparams8 = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGFineGrainEBIdMapRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.tpparams9 = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGFineGrainEBGroupRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.tpparams10 = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGFineGrainStripEERcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.tpparams11 = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGFineGrainTowerEERcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.tpparams12 = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGPhysicsConstRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.EcalTrigPrimESProducer = cms.ESProducer("EcalTrigPrimESProducer",
    DatabaseFile = cms.untracked.string('TPG_startup_threshold_4.txt')
)

process.ecalTriggerPrimitiveDigis = cms.EDProducer("EcalTrigPrimProducer",
    InstanceEB = cms.string('ebDigis'),
    InstanceEE = cms.string('eeDigis'),
    Label = cms.string('simEcalUnsuppressedDigis'),

    BarrelOnly = cms.bool(False),
    Famos = cms.bool(False),
    TcpOutput = cms.bool(False),
    Debug = cms.bool(False),

    binOfMaximum = cms.int32(6), ## optional from release (100) on, from 1-10
                                                   
    TTFHighEnergyEB = cms.double(1.0),
    TTFHighEnergyEE = cms.double(1.0),
    TTFLowEnergyEB = cms.double(1.0), ## this + the following is added from 140_pre4 on
    TTFLowEnergyEE = cms.double(1.0)
)



process.load("MMarionneau.JetViaTrigPrim.JetViaTrigPrim_cfi")

process.p = cms.Path(
                     process.ecalTriggerPrimitiveDigis
                     *process.simEcalDigis
                     *process.jetViaTrigPrim
                     )

process.outpath = cms.EndPath(process.o1)

