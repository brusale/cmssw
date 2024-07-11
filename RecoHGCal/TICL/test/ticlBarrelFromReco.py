# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 -s RAW2DIGI,RECO,RECOSIM,PAT,VALIDATION:@phase2Validation+@miniAODValidation,DQM:@phase2+@miniAODDQM --conditions auto:phase2_realistic_T21 --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO -n 4 --eventcontent FEVTDEBUGHLT,MINIAODSIM,DQM --geometry Extended2026D88 --era Run3 --filein file:step2.root --fileout /ceph/abrusamolino/CLUEBarrel/SinglePion0PU/step3.root --no_exec
import os

import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('TICLBarrelFromRECO2',Phase2C17I13M9)


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#import FWCore.ParameterSet.VarParsing as VarParsing
#options = VarParsing.VarParsing('analysis')

#options.register('doSharing',
#		 False,
#		 VarParsing.VarParsing.multiplicity.singleton,
#		 VarParsing.VarParsing.varType.bool,
#		 "boolean for CLUE energy sharing")
#options.parseArguments()

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
inputFileEOS  = 'file:/eos/user/a/abrusamo/SinglePion0PU_RECO.root'
inputFileCEPH = 'file:/ceph/abrusamolino/CLUEBarrel/SinglePion0PU/step3.root'
if  (os.path.exists(inputFileEOS .split('file:')[-1])):
    inputFile = inputFileEOS
elif(os.path.exists(inputFileCEPH.split('file:')[-1])):
    inputFile = inputFileCEPH
else:
    raise RuntimeError('Could not open either "%s" or "%s"' %(inputFileEOS, inputFileCEPH))
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(inputFile),
    secondaryFileNames = cms.untracked.vstring()
)

#process.MessageLogger = cms.Service("MessageLogger",
#  destinations = cms.untracked.vstring("clue_debug_2"),
#  clue_debug_2 = cms.untracked.PSet(threshold=cms.untracked.string("DEBUG")),
#  debugModules = cms.untracked.vstring("barrelLCAssocByEnergyScoreProducerPFCluster")
#)

process.options = cms.untracked.PSet(
    #FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring('StdException'),
    Rethrow = cms.untracked.vstring(),
    #kipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:4'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(
      'TICLBarrel.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.FEVTDEBUGHLToutput.outputCommands.extend([#'keep *_ebLayerClusters_*_*',
						  'drop *',
						  'keep *_barrelLayerClusters_*_*',
						  'keep *_simBarrelLayerClusters_*_*',
						  'keep *_*_MergedCaloTruth_*',
						  'keep recoPFCluster_*_*_*',
						  'keep recoPFRecHit_*_*_*',
						  'keep *_simBarrelLayerClusterCaloParticleAssociationProducer_*_*',
						  'keep *_barrelLayerClusterCaloParticleAssociationProducer_*_*',
						  'keep *_barrelLayerClusterCaloParticleAssociationProducerPFCluster_*_*', 
						  'keep *_barrelLayerClusterSimClusterAssociationProducerPFCluster_*_*', 
						  'keep *_lcFromPFClusterProducer_*_*',
                          'keep *_barrelHcalPatternRecognition_*_*',
                          'keep *_ticlSimTracksters_*_*'])

from Validation.RecoParticleFlow.customize_pfanalysis import *
process = customize_step3(process)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
	fileName = cms.untracked.string('TICLBarrel_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '131X_mcRun4_realistic_v6', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun4_realistic_v4', '')



from RecoHGCal.TICL.customiseTICLFromReco import customiseTICLForLCDumper, customiseTICLBarrelFromReco
customiseTICLForLCDumper(process)
customiseTICLBarrelFromReco(process)

# from RecoHGCal.TICL.customiseTICLFromReco import customiseTICLBarrelFromReco
# process = customiseTICLBarrelFromReco(process)

# process.barrelValidator.doCaloParticleSelection = cms.untracked.bool(False)
# process.simBarrelValidator.doCaloParticleSelection = cms.untracked.bool(False)

# process.simBarrelLayerClusters.ebplugin.deltac = cms.vdouble(1.8*0.0175, 5*0.087, 5*0.087)
# process.simBarrelLayerClusters.ebplugin.kappa = cms.double(3.5)
# process.simBarrelLayerClusters.ebplugin.outlierDeltaFactor = cms.double(2.7*0.0175)
# process.simBarrelLayerClusters.ebplugin.fractionCutoff = cms.double(0.01)

process.barrelLayerClusters.ebplugin.deltac = cms.vdouble(1.8*0.0175, 3*0.087, 5*0.087)#1.8
process.barrelLayerClusters.ebplugin.kappa = cms.double(3.5) #3.5
process.barrelLayerClusters.ebplugin.outlierDeltaFactor = cms.double(2.7*0.0175)#2.7
process.barrelLayerClusters.ebplugin.fractionCutoff = cms.double(0.01)

# process.simBarrelLayerClusters.hbplugin.deltac = cms.vdouble(1.8*0.0175, 3*0.087, 3*0.087)
# process.simBarrelLayerClusters.hbplugin.kappa = cms.double(1.)
# process.simBarrelLayerClusters.hbplugin.outlierDeltaFactor = cms.double(7*0.087)
# process.simBarrelLayerClusters.hbplugin.fractionCutoff = cms.double(0.01)
# process.simBarrelLayerClusters.hbplugin.maxLayerIndex = cms.int32(4)

process.barrelLayerClusters.hbplugin.deltac = cms.vdouble(1.8*0.0175, 3*0.087, 3*0.087)
process.barrelLayerClusters.hbplugin.kappa = cms.double(0.)
process.barrelLayerClusters.hbplugin.outlierDeltaFactor = cms.double(5*0.087)
process.barrelLayerClusters.hbplugin.fractionCutoff = cms.double(0.01)
process.barrelLayerClusters.hbplugin.maxLayerIndex = cms.int32(4)

# process.simBarrelLayerClusters.hoplugin.deltac = cms.vdouble(1.8*0.0175, 3*0.087, 3*0.087)
# process.simBarrelLayerClusters.hoplugin.kappa = cms.double(1.)
# process.simBarrelLayerClusters.hoplugin.outlierDeltaFactor = cms.double(7*0.087)
# process.simBarrelLayerClusters.hoplugin.fractionCutoff = cms.double(0.01)
# process.simBarrelLayerClusters.hoplugin.maxLayerIndex = cms.int32(5)

process.barrelLayerClusters.hoplugin.deltac = cms.vdouble(1.8*0.0175, 3*0.087, 3*0.087)
process.barrelLayerClusters.hoplugin.kappa = cms.double(0.)
process.barrelLayerClusters.hoplugin.outlierDeltaFactor = cms.double(5*0.087)
process.barrelLayerClusters.hoplugin.fractionCutoff = cms.double(0.01)
process.barrelLayerClusters.hoplugin.maxLayerIndex = cms.int32(5)

# process.barrelValidator.histoAlgoBlock = cms.int32(20)

# from RecoHGCal.TICL.customiseTICLFromReco import customiseTICLForLCDumper
# process = customiseTICLForLCDumper(process)

# process.consumer = cms.EDAnalyzer("GenericConsumer",
#   eventProducts = cms.untracked.vstring('lcFromPFClusterProducer')
# )
# process.consumer2 = cms.EDAnalyzer("GenericConsumer",
#   eventProducts = cms.untracked.vstring('barrelLCAssocByEnergyScoreProducerPFCluster')
# )
# process.consumer3 = cms.EDAnalyzer("GenericConsumer",
#   eventProducts = cms.untracked.vstring("barrelSCAssocByEnergyScoreProducerPFCluster")
# )
# process.consumer4 = cms.EDAnalyzer("GenericConsumer",
#   eventProducts = cms.untracked.vstring("barrelLayerClusterCaloParticleAssociationProducerPFCluster")
# )
# process.consumer5 = cms.EDAnalyzer("GenericConsumer",
#   eventProducts = cms.untracked.vstring("barrelLayerClusterSimClusterAssociationProducerPFCluster")
# )

# process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput + process.consumer + process.consumer2 + process.consumer3 + process.consumer4 + process.consumer5 + process.lcDumper + process.lcDumperPF)
# process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# #process.lcDumper.layerclusters = cms.InputTag("lcFromPFClusterProducer")
# #process.lcDumper.simToRecoCollection = cms.InputTag("barrelLayerClusterCaloParticleAssociationProducerPFCluster")
# #process.lcDumper.recoToSimCollection = cms.InputTag("barrelLayerClusterCaloParticleAssociationProducerPFCluster")

# process.schedule = cms.Schedule(process.TICLBarrel,
#     			    process.TICLBarrel_Validation,
#     			    process.FEVTDEBUGHLToutput_step,
#     			    process.DQMoutput_step)


# process.TFileService = cms.Service("TFileService", 
# 				   fileName = cms.string("histo.root"))
# # Schedule definition

#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions
# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

