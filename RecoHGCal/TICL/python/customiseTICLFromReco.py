# Reconstruction
from RecoHGCal.TICL.iterativeTICL_cff import *
from RecoLocalCalo.HGCalRecProducers.hgcalMergeLayerClusters_cfi import hgcalMergeLayerClusters
from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff import hgcalLayerClustersEE, hgcalLayerClustersHSi, hgcalLayerClustersHSci
from RecoLocalCalo.HGCalRecProducers.barrelLayerClusters_cfi import barrelLayerClusters
from RecoLocalCalo.HGCalRecProducers.simBarrelLayerClusters_cfi import simBarrelLayerClusters
from RecoHGCal.TICL.ticlDumper_cfi import ticlDumper
from RecoHGCal.TICL.layerClusterDumper_cfi import layerClusterDumper
# Validation
from Validation.HGCalValidation.HGCalValidator_cfi import *
from Validation.HGCalValidation.BarrelValidator_cfi import *
from Validation.HGCalValidation.SimBarrelValidator_cfi import *
from Validation.HGCalValidation.BarrelValidatorPFCluster_cfi import *

from RecoLocalCalo.HGCalRecProducers.hgcalRecHitMapProducer_cfi import hgcalRecHitMapProducer
from RecoLocalCalo.HGCalRecProducers.barrelRecHitMapProducer_cfi import barrelRecHitMapProducer
from RecoLocalCalo.HGCalRecProducers.simBarrelRecHitMapProducer_cfi import simBarrelRecHitMapProducer
# Load DNN ESSource
from RecoTracker.IterativeTracking.iterativeTk_cff import trackdnn_source

# Automatic addition of the customisation function from RecoHGCal.Configuration.RecoHGCal_EventContent_cff
from RecoHGCal.Configuration.RecoHGCal_EventContent_cff import customiseHGCalOnlyEventContent
from SimCalorimetry.HGCalAssociatorProducers.simTracksterAssociatorByEnergyScore_cfi import simTracksterAssociatorByEnergyScore as simTsAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationLinking, tracksterSimTracksterAssociationPR,tracksterSimTracksterAssociationLinkingbyCLUE3D, tracksterSimTracksterAssociationPRbyCLUE3D

# Barrel associators
from SimCalorimetry.HGCalAssociatorProducers.barrelLayerClusterAssociatorByEnergyScore_cfi import barrelLayerClusterAssociatorByEnergyScore as barrelLCAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.barrelSimClusterAssociatorByEnergyScore_cfi import barrelSimClusterAssociatorByEnergyScore as barrelSCAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.simBarrelLayerClusterAssociatorByEnergyScore_cfi import simBarrelLayerClusterAssociatorByEnergyScore as simBarrelLCAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.simBarrelSimClusterAssociatorByEnergyScore_cfi import simBarrelSimClusterAssociatorByEnergyScore as simBarrelSCAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.simBarrelSimClusterAssociatorByEnergyScore_cfi import simBarrelSimClusterAssociatorByEnergyScore as simBarrelSCAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.BarrelLCToCPAssociation_cfi import barrelLayerClusterCaloParticleAssociation as barrelLayerClusterCaloParticleAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.SimBarrelLCToCPAssociation_cfi import simBarrelLayerClusterCaloParticleAssociation as simBarrelLayerClusterCaloParticleAssociationProdcuer
from SimCalorimetry.HGCalAssociatorProducers.BarrelLCToSCAssociation_cfi import barrelLayerClusterSimClusterAssociation as barrelLayerClusterSimClusterAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.SimBarrelLCToSCAssociation_cfi import simBarrelLayerClusterSimClusterAssociation as simBarrelSimClusterAssociationProducer
from RecoHGCal.TICL.lcFromPFClusterProducer_cfi import lcFromPFClusterProducer

def customiseTICLFromReco(process):
# TensorFlow ESSource
    process.TFESSource = cms.Task(process.trackdnn_source)
# Reconstruction
    process.hgcalLayerClustersTask = cms.Task(process.hgcalLayerClustersEE,
                                              process.hgcalLayerClustersHSi,
                                              process.hgcalLayerClustersHSci,
                                              process.hgcalMergeLayerClusters)


    process.TICL = cms.Path(process.hgcalLayerClustersTask,
                            process.TFESSource,
                            process.ticlLayerTileTask,
                            process.ticlIterationsTask,
                            process.ticlTracksterMergeTask)
# Validation
    process.TICL_ValidationProducers = cms.Task(process.hgcalRecHitMapProducer,
                                                process.lcAssocByEnergyScoreProducer,
                                                process.layerClusterCaloParticleAssociationProducer,
                                                process.scAssocByEnergyScoreProducer,
                                                process.layerClusterSimClusterAssociationProducer,
                                                process.simTsAssocByEnergyScoreProducer,  process.simTracksterHitLCAssociatorByEnergyScoreProducer, process.tracksterSimTracksterAssociationLinking, process.tracksterSimTracksterAssociationPR, process.tracksterSimTracksterAssociationLinkingbyCLUE3D, process.tracksterSimTracksterAssociationPRbyCLUE3D
                                               )

    process.TICL_Validator = cms.Task(process.hgcalValidator)
    process.TICL_Validation = cms.Path(process.TICL_ValidationProducers,
                                       process.TICL_Validator
                                      )
# Path and EndPath definitions
    process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
    process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
    process.schedule = cms.Schedule(process.TICL,
                                    process.TICL_Validation,
                                    process.FEVTDEBUGHLToutput_step,
                                    process.DQMoutput_step)
#call to customisation function customiseHGCalOnlyEventContent imported from RecoHGCal.Configuration.RecoHGCal_EventContent_cff
    process = customiseHGCalOnlyEventContent(process)

    return process

def customiseTICLBarrelFromReco(process):

    #process.lcFromPFClusterProducer = lcFromPFClusterProducer.clone()

    process.barrelLayerClustersTask = cms.Task(process.simBarrelLayerClusters,
					       process.barrelLayerClusters,
					       #process.particleFlowClusterECALUncorrected,
					       process.lcFromPFClusterProducer
					       )
    
    process.TICLBarrel = cms.Path(process.barrelLayerClustersTask)


    process.barrelLCAssocByEnergyScoreProducerPFCluster = barrelLCAssocByEnergyScoreProducer.clone()
    process.barrelSCAssocByEnergyScoreProducerPFCluster = barrelSCAssocByEnergyScoreProducer.clone()
    process.barrelLayerClusterCaloParticleAssociationProducerPFCluster = barrelLayerClusterCaloParticleAssociationProducer.clone()
    process.barrelLayerClusterSimClusterAssociationProducerPFCluster = barrelLayerClusterSimClusterAssociationProducer.clone()
    
    process.barrelLayerClusterCaloParticleAssociationProducerPFCluster.label_lc = cms.InputTag("lcFromPFClusterProducer")
    process.barrelLayerClusterSimClusterAssociationProducerPFCluster.label_lcl = cms.InputTag("lcFromPFClusterProducer")
    process.barrelLayerClusterCaloParticleAssociationProducerPFCluster.associator = cms.InputTag("barrelLCAssocByEnergyScoreProducerPFCluster")
    process.barrelLayerClusterSimClusterAssociationProducerPFCluster.label_lc = cms.InputTag("lcFromPFClusterProducer")
    process.barrelValidatorPFCluster = barrelValidatorPFCluster.clone()
    process.barrelValidatorPFCluster.associator = cms.untracked.InputTag("barrelLayerClusterCaloParticleAssociationProducerPFCluster")
    process.barrelValidatorPFCluster.associatorSim = cms.untracked.InputTag("barrelLayerClusterSimClusterAssociationProducerPFCluster")
    process.barrelValidatorPFCluster.label_layerClusterPlots = cms.InputTag("lcFromPFClusterProducer")
    process.barrelValidatorPFCluster.label_lcl = cms.InputTag("lcFromPFClusterProducer")

    process.barrelValidatorPFCluster.dirName = cms.string('PFCluster/PFClusterValidator/')

    process.TICLBarrel_ValidationProducers = cms.Task(process.barrelRecHitMapProducer,
						      process.simBarrelRecHitMapProducer,
						      process.barrelLCAssocByEnergyScoreProducerPFCluster,
						      process.barrelSCAssocByEnergyScoreProducerPFCluster,
						      process.barrelLayerClusterCaloParticleAssociationProducerPFCluster,
						      process.barrelLayerClusterSimClusterAssociationProducerPFCluster,
						      process.barrelLCAssocByEnergyScoreProducer,
						      process.barrelLayerClusterCaloParticleAssociationProducer,
						      process.barrelSCAssocByEnergyScoreProducer,
						      process.barrelLayerClusterSimClusterAssociationProducer,
						      process.simBarrelLCAssocByEnergyScoreProducer,
						      process.simBarrelLayerClusterCaloParticleAssociationProducer,
						      process.simBarrelSCAssocByEnergyScoreProducer,
						      process.simBarrelLayerClusterSimClusterAssociationProducer)
    process.TICLBarrel_Validator = cms.Task(process.barrelValidator,
					   process.simBarrelValidator, 
					   process.barrelValidatorPFCluster
					   )
    process.TICLBarrel_Validation = cms.Path(process.TICLBarrel_ValidationProducers,
					     process.TICLBarrel_Validator)

    process.consumer = cms.EDAnalyzer("GenericConsumer",
      eventProducers = cms.untracked.vstring('lcFromPFClusterProducer')
    )
    process.consumer2 = cms.EDAnalyzer("GenericConsumer",
      eventProducers = cms.untracked.vstring('barrelLCAssocByEnergyScoreProducerPFCluster')
    )
    process.consumer3 = cms.EDAnalyzer("GenericConsumer",
      eventProducers = cms.untracked.vstring("barrelSCAssocByEnergyScoreProducerPFCluster")
    )
    process.consumer4 = cms.EDAnalyzer("GenericConsumer",
      eventProducers = cms.untracked.vstring("barrelLayerClusterCaloParticleAssociationProducerPFCluster")
    )
    process.consumer5 = cms.EDAnalyzer("GenericConsumer",
      eventProducers = cms.untracked.vstring("barrelLayerClusterSimClusterAssociationProducerPFCluster")
    )
    process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput + process.consumer + process.consumer2 + process.consumer3 + process.consumer4 + process.consumer5)
    process.DQMoutput_step = cms.EndPath(process.DQMoutput)

    process.schedule = cms.Schedule(process.TICLBarrel,
				    process.TICLBarrel_Validation,
				    process.FEVTDEBUGHLToutput_step,
				    process.DQMoutput_step)
    return process

def customiseTICLForDumper(process):

				process.ticlDumper = ticlDumper.clone()
				process.TFileService = cms.Service("TFileService",
												fileName = cms.string("histo.root")
												)
				process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput + process.ticlDumper)
				return process
def customiseTICLForLCDumper(process):
				process.lcDumper = layerClusterDumper.clone()
				process.lcDumper.layerclusters = cms.InputTag("barrelLayerClusters")
				process.lcDumper.simlayerclusters = cms.InputTag("simBarrelLayerClusters")
				process.lcDumperPF = layerClusterDumper.clone()
				process.lcDumperPF.layerclusters = cms.InputTag("lcFromPFClusterProducer")
				process.lcDumperPF.simToRecoCollection = cms.InputTag("barrelLayerClusterCaloParticleAssociationProducerPFCluster")
				process.lcDumperPF.recoToSimCollection = cms.InputTag("barrelLayerClusterCaloParticleAssociationProducerPFCluster")
				process.TFileService = cms.Service("TFileService", 
								   fileName = cms.string("histo.root"))
				process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput + process.lcDumper + process.lcDumperPF)
				return process
