import FWCore.ParameterSet.Config as cms

from Validation.HGCalValidation.CaloParticleSelectionForEfficiency_cfi import *
from Validation.HGCalValidation.HGVHistoProducerAlgoBlock_cfi import *

from SimCalorimetry.HGCalAssociatorProducers.BarrelLCToCPAssociation_cfi import barrelLayerClusterCaloParticleAssociation
from SimCalorimetry.HGCalAssociatorProducers.BarrelLCToSCAssociation_cfi import barrelLayerClusterSimClusterAssociation

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer

from RecoHGCal.TICL.iterativeTICL_cff import ticlIterLabels, ticlIterLabelsMerge

labelTst = [cms.InputTag("ticlTracksters"+iteration) for iteration in ticlIterLabelsMerge]
labelTst.extend([cms.InputTag("ticlSimTracksters", "fromCPs"), cms.InputTag("ticlSimTracksters")])
lcInputMask = [cms.InputTag("ticlTracksters"+iteration) for iteration in ticlIterLabels]
lcInputMask.extend([cms.InputTag("ticlSimTracksters", "fromCPs"), cms.InputTag("ticlSimTracksters")])
barrelValidator = DQMEDAnalyzer(
    "BarrelValidator",

    ### general settings ###
    # selection of CP for evaluation of efficiency #
    CaloParticleSelectionForEfficiency,

    ### reco input configuration ###
    #2DLayerClusters, PFClusters, Tracksters
    label_lcl = barrelLayerClusterCaloParticleAssociation.label_lc,
    label_tst = cms.VInputTag(labelTst),
    label_simTS = cms.InputTag("ticlSimTracksters"),
    label_simTSFromCP = cms.InputTag("ticlSimTracksters", "fromCPs"),

    associator = cms.untracked.InputTag("barrelLayerClusterCaloParticleAssociationProducer"),

    associatorSim = cms.untracked.InputTag("barrelLayerClusterSimClusterAssociationProducer"),

    #General info on layers etc.
    SaveGeneralInfo = cms.untracked.bool(True),
    #CaloParticle related plots
    doCaloParticlePlots = cms.untracked.bool(True),
    #Select caloParticles for efficiency or pass through
    doCaloParticleSelection = cms.untracked.bool(True),
    #SimCluster related plots
    doSimClustersPlots = cms.untracked.bool(True),
    label_SimClusters = cms.InputTag("SimClusters"),
    label_SimClustersLevel = cms.InputTag("ClusterLevel"),
    #Layer Cluster related plots
    doLayerClustersPlots = cms.untracked.bool(True),
    label_layerClusterPlots = cms.InputTag("barrelLayerClusters"),
    label_LCToCPLinking = cms.InputTag("LCToCP_association"),
    #Trackster related plots
    doTrackstersPlots = cms.untracked.bool(True),
    label_TS = cms.string("Morphology"),
    label_TSToCPLinking = cms.string("TSToCP_linking"),
    label_TSToSTSPR = cms.string("TSToSTS_patternRecognition"),

    #The cumulative material budget in front of each layer. To be more specific, it
    #is the material budget just in front of the active material (not including it).
    #This file is created using the official material budget code.
    cummatbudinxo = cms.FileInPath('Validation/HGCalValidation/data/D41.cumulative.xo'),

    ### sim input configuration ###
    label_cp_effic = barrelLayerClusterCaloParticleAssociation.label_cp,
    label_cp_fake = cms.InputTag("mix","MergedCaloTruth"),
    #simClusters
    label_scl = barrelLayerClusterSimClusterAssociation.label_scl,

    simVertices = cms.InputTag("g4SimHits"),

    LayerClustersInputMask = cms.VInputTag(lcInputMask),

    #Total number of layers of HGCal that we want to monitor
    #Could get this also from HGCalImagingAlgo::maxlayer but better to get it from here
    totallayers_to_monitor = cms.int32(6),
    #Thicknesses we want to monitor. -1 is for scintillator
    thicknesses_to_monitor = cms.vint32(120,200,300,-1),

    # HistoProducerAlgo. Defines the set of plots to be booked and filled
    histoProducerAlgoBlock = HGVHistoProducerAlgoBlock,

    ### output configuration
    dirName = cms.string('CLUE-barrel/BarrelValidator/')

)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(barrelValidator,
    label_cp_fake = "mixData:MergedCaloTruth"
)

#from Configuration.Eras.Modifier_phase2_hgcalV10_cff import phase2_hgcalV10
#phase2_hgcalV10.toModify(barrelValidator, totallayers_to_monitor = cms.int32(50))

#from Configuration.Eras.Modifier_phase2_hgcalV16_cff import phase2_hgcalV16
#phase2_hgcalV16.toModify(barrelValidator, totallayers_to_monitor = cms.int32(47))
