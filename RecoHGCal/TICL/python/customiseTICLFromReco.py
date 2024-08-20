# Reconstruction
from RecoHGCal.TICL.iterativeTICL_cff import *
from RecoLocalCalo.HGCalRecProducers.hgcalMergeLayerClusters_cfi import hgcalMergeLayerClusters
from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff import hgcalLayerClustersEE, hgcalLayerClustersHSi, hgcalLayerClustersHSci
from RecoHGCal.TICL.ticlDumper_cfi import ticlDumper
from RecoHGCal.TICL.layerClusterDumper_cfi import layerClusterDumper
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
# Validation
from Validation.HGCalValidation.HGCalValidator_cfi import *
from RecoLocalCalo.HGCalRecProducers.recHitMapProducer_cfi import recHitMapProducer

#from RecoLocalCalo.HGCalRecProducers.hgcalRecHitMapProducer_cfi import hgcalRecHitMapProducer
from RecoLocalCalo.HGCalRecProducers.recHitMapProducer_cfi import recHitMapProducer
#from RecoLocalCalo.HGCalRecProducers.simBarrelRecHitMapProducer_cfi import simBarrelRecHitMapProducer
# Load DNN ESSource
from RecoTracker.IterativeTracking.iterativeTk_cff import trackdnn_source

# Automatic addition of the customisation function from RecoHGCal.Configuration.RecoHGCal_EventContent_cff
from RecoHGCal.Configuration.RecoHGCal_EventContent_cff import customiseHGCalOnlyEventContent
from SimCalorimetry.HGCalAssociatorProducers.simTracksterAssociatorByEnergyScore_cfi import simTracksterAssociatorByEnergyScore as simTsAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationLinking, tracksterSimTracksterAssociationPR,tracksterSimTracksterAssociationLinkingbyCLUE3D, tracksterSimTracksterAssociationPRbyCLUE3D

# Barrel associators
from RecoParticleFlow.PFClusterProducer.barrelLayerClusters_cff  import barrelLayerClusters
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer
from SimCalorimetry.HGCalAssociatorProducers.barrelLCToCPAssociatorByEnergyScoreProducer_cfi import barrelLCToCPAssociatorByEnergyScoreProducer as barrelLCToCPAssociatorByEnergyScoreProducer

from SimCalorimetry.HGCalAssociatorProducers.barrelLCToSCAssociatorByEnergyScoreProducer_cfi import barrelLCToSCAssociatorByEnergyScoreProducer as barrelLCToSCAssociatorByEnergyScoreProducer

from SimCalorimetry.HGCalAssociatorProducers.LCToCPAssociation_cfi import barrelLayerClusterCaloParticleAssociation as barrelLayerClusterCaloParticleAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cfi import barrelLayerClusterSimClusterAssociation as barrelLayerClusterSimClusterAssociationProducer

from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import barrelTracksterSimTracksterAssociationPR, barrelTracksterSimTracksterAssociationLinkingPR
from SimCalorimetry.HGCalAssociatorProducers.barrelSimTracksterHitLCAssociatorByEnergyScore_cfi import barrelSimTracksterHitLCAssociatorByEnergyScore
from RecoHGCal.TICL.lcFromPFClusterProducer_cfi import lcFromPFClusterProducer  
from RecoHGCal.TICL.recHitDumper_cfi import recHitDumper

# For the moment ignore
# from SimCalorimetry.HGCalAssociatorProducers.barrelSimClusterAssociatorByEnergyScore_cfi import barrelSimClusterAssociatorByEnergyScore as barrelSCAssocByEnergyScoreProducer
# from SimCalorimetry.HGCalAssociatorProducers.simBarrelLayerClusterAssociatorByEnergyScore_cfi import simBarrelLayerClusterAssociatorByEnergyScore as simBarrelLCAssocByEnergyScoreProducer
# from SimCalorimetry.HGCalAssociatorProducers.simBarrelSimClusterAssociatorByEnergyScore_cfi import simBarrelSimClusterAssociatorByEnergyScore as simBarrelSCAssocByEnergyScoreProducer
# from SimCalorimetry.HGCalAssociatorProducers.simBarrelSimClusterAssociatorByEnergyScore_cfi import simBarrelSimClusterAssociatorByEnergyScore as simBarrelSCAssocByEnergyScoreProducer
# from SimCalorimetry.HGCalAssociatorProducers.BarrelLCToCPAssociation_cfi import barrelLayerClusterCaloParticleAssociation as barrelLayerClusterCaloParticleAssociationProducer
# from SimCalorimetry.HGCalAssociatorProducers.SimBarrelLCToCPAssociation_cfi import simBarrelLayerClusterCaloParticleAssociation as simBarrelLayerClusterCaloParticleAssociationProdcuer
# from SimCalorimetry.HGCalAssociatorProducers.BarrelLCToSCAssociation_cfi import barrelLayerClusterSimClusterAssociation as barrelLayerClusterSimClusterAssociationProducer
# from SimCalorimetry.HGCalAssociatorProducers.SimBarrelLCToSCAssociation_cfi import simBarrelLayerClusterSimClusterAssociation as simBarrelSimClusterAssociationProducer


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
    process.TICL_ValidationProducers = cms.Task(process.recHitMapProducer,
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

    process.lcFromPFClusterProducer = lcFromPFClusterProducer.clone()
    
    
    process.barrelLayerClustersTask = cms.Task(process.barrelLayerClusters
                                               ,process.lcFromPFClusterProducer
    )

    process.barrelPatternRecognition = _trackstersProducer.clone(
        detector = "Barrel",
        layer_clusters = "barrelLayerClusters",
        time_layerclusters = "barrelLayerClusters:timeLayerCluster",
        filtered_mask = "barrelLayerClusters:InitialLayerClustersMask",
        original_mask = "barrelLayerClusters:InitialLayerClustersMask",
        seeding_regions = "ticlSeedingGlobal",
        itername = "CA",
        patternRecognitionBy = "CA",
        doRegression = cms.bool(False),
        pluginPatternRecognitionByFastJet = dict (
            antikt_radius = 0.1,
            minNumLayerCluster = 0,
            algo_verbosity = 2
        ),
        pluginPatternRecognitionByCA = dict (
            min_layers_per_trackster = 1,
            root_doublet_max_distance_from_seed_squared = 9998,
            min_cos_theta = 0.99,
            siblings_maxRSquared = cms.vdouble(3*0.0175, 5*0.087, 5*0.087)
        ),
        pluginPatternRecognitionByCLUE3D = dict (
            criticalDensity = cms.vdouble(0.5, 0.5, 0.5),
            criticalSelfDensity = cms.vdouble(0., 0., 0.),
            criticalEtaPhiDistance = cms.vdouble(7*0.087, 7*0.087, 7*0.087, 7*0.087),
            densityEtaPhiDistanceSqr = cms.vdouble(0.37, 0.37, 0.37, 0.37), # 7*0.087 squared
            #nearestHigherOnSameLayer = cms.bool(True),
            #densityOnSameLayer = cms.bool(True),
            minNumLayerCluster = cms.vint32(1, 1, 1),
            useAbsoluteProjectiveScale = cms.bool(False),
            densitySiblingLayers = cms.vint32(4,4,4)
        )
    )

    process.barrelLayerClusters.hbplugin.kappa = cms.double(1.5)

    process.ticlSimTracksters.layer_clusters = cms.InputTag("barrelLayerClusters") 
    process.ticlSimTracksters.filtered_mask = cms.InputTag("barrelLayerClusters:InitialLayerClustersMask")
    process.ticlSimTracksters.time_layerclusters = cms.InputTag("barrelLayerClusters:timeLayerCluster")
    process.ticlSimTracksters.layerClusterCaloParticleAssociator = cms.InputTag("barrelLayerClusterCaloParticleAssociationProducer")
    process.ticlSimTracksters.layerClusterSimClusterAssociator = cms.InputTag("barrelLayerClusterSimClusterAssociationProducer")

    process.barrelSimTracksterHitLCAssociatorByEnergyScoreProducer = barrelSimTracksterHitLCAssociatorByEnergyScore.clone() 
    process.barrelTracksterSimTracksterAssociationPR = barrelTracksterSimTracksterAssociationPR.clone()
    process.barrelTracksterSimTracksterAssociationLinkingPR = barrelTracksterSimTracksterAssociationLinkingPR.clone() 
    process.barrellPatternRecognitionSeq = cms.Sequence(process.barrelPatternRecognition)
    process.ticlLayerTileProducer.detector = cms.string('Barrel')
    process.TICLBarrel = cms.Path(process.ticlSimTracksters 
                                  + process.ticlLayerTileProducer
                                  + process.ticlSeedingGlobal 
                                  + process.barrelPatternRecognition
                                  , process.barrelLayerClustersTask)

    # We want to run CLUE on not -cleaned Rechit collections
    process.recHitMapProducer.EBInput = cms.InputTag("particleFlowRecHitECAL")
    process.recHitMapProducer.HBInput = cms.InputTag("particleFlowRecHitHBHE")
    process.recHitMapProducer.HOInput = cms.InputTag("particleFlowRecHitHO")

    process.barrelLayerClusters.EBInput = cms.InputTag("particleFlowRecHitECAL")
    process.barrelLayerClusters.HBInput = cms.InputTag("particleFlowRecHitHBHE")
    process.barrelLayerClusters.HOInput = cms.InputTag("particleFlowRecHitHO")
       
    process.barrelLCAssocByEnergyScoreProducerPFCluster = barrelLCToCPAssociatorByEnergyScoreProducer.clone()
    process.barrelSCAssocByEnergyScoreProducerPFCluster = barrelLCToSCAssociatorByEnergyScoreProducer.clone()
    process.barrelLayerClusterCaloParticleAssociationProducerPFCluster =  barrelLayerClusterCaloParticleAssociationProducer.clone()
    process.barrelLayerClusterSimClusterAssociationProducerPFCluster = barrelLayerClusterSimClusterAssociationProducer.clone()
  
    process.barrelLayerClusterCaloParticleAssociationProducer.label_lc = cms.InputTag("barrelLayerClusters") 
    process.barrelLayerClusterCaloParticleAssociationProducerPFCluster.label_lc = cms.InputTag("lcFromPFClusterProducer")
    process.barrelLayerClusterSimClusterAssociationProducerPFCluster.label_lcl = cms.InputTag("lcFromPFClusterProducer")
    process.barrelLayerClusterCaloParticleAssociationProducerPFCluster.associator = cms.InputTag("barrelLCAssocByEnergyScoreProducerPFCluster")
    process.barrelLayerClusterSimClusterAssociationProducerPFCluster.label_lc = cms.InputTag("lcFromPFClusterProducer")

    process.barrelLCToCPAssociatorByEnergyScoreProducer.hits = cms.VInputTag("particleFlowRecHitECAL","particleFlowRecHitHBHE")
    process.barrelLCToSCAssociatorByEnergyScoreProducer.hits = cms.VInputTag("particleFlowRecHitECAL","particleFlowRecHitHBHE")

    # Uncorrected PFClusters are build without "cleaned "
    process.barrelLCAssocByEnergyScoreProducerPFCluster.hits = cms.VInputTag("particleFlowRecHitECAL","particleFlowRecHitHBHE")
    process.barrelSCAssocByEnergyScoreProducerPFCluster.hits = cms.VInputTag("particleFlowRecHitECAL","particleFlowRecHitHBHE")
    
    
    # Make the plots in CMSSW
    # process.barrelValidatorPFCluster = barrelValidatorPFCluster.clone()
    # process.barrelValidatorPFCluster.associator = cms.untracked.InputTag("barrelLayerClusterCaloParticleAssociationProducerPFCluster")
    # process.barrelValidatorPFCluster.associatorSim = cms.untracked.InputTag("barrelLayerClusterSimClusterAssociationProducerPFCluster")
    # process.barrelValidatorPFCluster.label_layerClusterPlots = cms.InputTag("lcFromPFClusterProducer")
    # process.barrelValidatorPFCluster.label_lcl = cms.InputTag("lcFromPFClusterProducer")

    # process.barrelValidatorPFCluster.dirName = cms.string('PFCluster/PFClusterValidator/')

    process.TICLBarrel_ValidationProducers = cms.Task(process.recHitMapProducer,
         					                          process.barrelLCAssocByEnergyScoreProducerPFCluster,
        					                          process.barrelSCAssocByEnergyScoreProducerPFCluster,
         					                          #process.barrelLayerClusterCaloParticleAssociationProducerPFCluster,
         					                          #process.barrelLayerClusterSimClusterAssociationProducerPFCluster,
                                                      process.barrelLCToCPAssociatorByEnergyScoreProducer,
                                                      process.barrelLCToSCAssociatorByEnergyScoreProducer,
                                                      process.barrelLayerClusterCaloParticleAssociationProducer,
                                                      process.barrelLayerClusterSimClusterAssociationProducer,
                                                      process.barrelSimTracksterHitLCAssociatorByEnergyScoreProducer,
                                                      process.barrelTracksterSimTracksterAssociationPR,
                                                      process.barrelTracksterSimTracksterAssociationLinkingPR,
    )
    #process.TICLBarrel_Validator = cms.Task(process.barrelValidator,
    #     				                    process.barrelValidatorPFCluster
    #)
    process.TICLBarrel_Validation = cms.Path(process.TICLBarrel_ValidationProducers
         				                     #+process.TICLBarrel_Validator
    )

    process.ticlDumper = ticlDumper.clone(
        saveLCs = True,
        saveCLUE3DTracksters = True,
        saveTrackstersMerged = False,
        saveSimTrackstersSC = True,
        saveSimTrackstersCP = True,
        saveTICLCandidate = False,
        saveSimTICLCandidate = False,
        saveTracks = False,
        saveAssociations = True,
        trackstersclue3d = "barrelPatternRecognition",
        layerClusters = "barrelLayerClusters",
        layer_clustersTime = "barrelLayerClusters:timeLayerCluster",
        recoToSimAssociatorSC = "barrelTracksterSimTracksterAssociationPR:recoToSim",
        simToRecoAssociatorSC = "barrelTracksterSimTracksterAssociationPR:simToReco",
        recoToSimAssociatorCP = "barrelTracksterSimTracksterAssociationLinkingPR:recoToSim",
        simToRecoAssociatorCP = "barrelTracksterSimTracksterAssociationLinkingPR:simToReco",
    )

    process.consumer = cms.EDAnalyzer("GenericConsumer",
      eventProducts = cms.untracked.vstring(['lcFromPFClusterProducer',
                                            'barrelLCAssocByEnergyScoreProducerPFCluster',
                                            'barrelSCAssocByEnergyScoreProducerPFCluster',
                                            'barrelLayerClusterCaloParticleAssociationProducerPFCluster',
                                            'barrelLayerClusterSimClusterAssociationProducerPFCluster',
                                            'barrelPatternRecognition',
                                            'ticlSimTrackster',
                                            'barrelSimTracksterHitLCAssociatorByEnergyScoreProducer',
                                            'tracksterSimTracksterAssociationPRbyCLUE3D',
      ])
    )

    process.FastTimerService.writeJSONSummary = cms.untracked.bool(True)
    process.FastTimerService.jsonFileName = cms.untracked.string('resourcesCA.json')

    process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput
                                                  +process.consumer
                                                  #+process.lcDumper
                                                  #+process.lcDumperPF 
                                                  +process.recHitDumper
                                                  +process.ticlDumper
    )
  
    process.DQMoutput_step = cms.EndPath(process.DQMoutput + process.fastTimerServiceClient)

    process.schedule = cms.Schedule(process.TICLBarrel,
				                    process.TICLBarrel_Validation,
				                    process.FEVTDEBUGHLToutput_step,
				                    process.DQMoutput_step
    )
    return process

def customiseTICLForDumper(process):

    process.ticlDumper = ticlDumper.clone(
        saveLCs=True,
        saveCLUE3DTracksters=True,
        saveTrackstersMerged=True,
        saveSimTrackstersSC=True,
        saveSimTrackstersCP=True,
        saveTICLCandidate=True,
        saveSimTICLCandidate=True,
        saveTracks=True,
        saveAssociations=True,
    )

    from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
    ticl_v5.toModify(process.ticlDumper,
                     # trackstersclue3d = cms.InputTag('mergedTrackstersProducer'), # For future separate iterations
                     trackstersclue3d = cms.InputTag('ticlTrackstersCLUE3DHigh'),
                     ticlcandidates = cms.InputTag("ticlCandidate"),
                     trackstersmerged = cms.InputTag("ticlCandidate"),
                     trackstersInCand = cms.InputTag("ticlCandidate"))

    process.TFileService = cms.Service("TFileService",
                                       fileName=cms.string("histo.root")
                                       )
    process.FEVTDEBUGHLToutput_step = cms.EndPath(
        process.FEVTDEBUGHLToutput + process.ticlDumper)
    return process

def customiseTICLForLCDumper(process):
    process.recHitDumper = recHitDumper.clone()
    process.lcDumper = layerClusterDumper.clone()
    process.lcDumper.layerclusters = cms.InputTag("barrelLayerClusters")
    #process.lcDumper.pfrechits = cms.InputTag("particleFlowRecHitECAL", "Cleaned")
    #process.lcDumper.simlayerclusters = cms.InputTag("simBarrelLayerClusters")
    process.lcDumperPF = layerClusterDumper.clone()
    process.lcDumperPF.layerclusters = cms.InputTag("lcFromPFClusterProducer")
    process.lcDumperPF.simToRecoCollection = cms.InputTag("barrelLayerClusterCaloParticleAssociationProducerPFCluster")
    process.lcDumperPF.recoToSimCollection = cms.InputTag("barrelLayerClusterCaloParticleAssociationProducerPFCluster")
    #process.lcDumperPF.pfrechits = cms.InputTag("particleFlowRecHitECAL")# without cleaned
    process.TFileService = cms.Service("TFileService",
				       fileName = cms.string("histo.root"))
    #process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput_step + process.lcDumper + process.lcDumperPF)
    return process
