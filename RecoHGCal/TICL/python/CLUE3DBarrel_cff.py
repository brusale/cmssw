import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer

## links
from RecoHGCal.TICL.tracksterLinksProducer_cfi import tracksterLinksProducer as _tracksterLinksProducer

## candidates
from RecoHGCal.TICL.ticlCandidateProducer_cfi import ticlCandidateProducer as _ticlCandidateProducer

filteredLayerClustersCLUE3DBarrel = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSizeBarrel",
    min_cluster_size = 2,
    algo_number = [10, 11],
    iteration_label = "CLUE3DBarrel",
    max_layerId = 5
)

# PATTERN RECOGNITION
ticlTrackstersCLUE3DBarrel = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersCLUE3DBarrel:CLUE3DBarrel",
    seeding_regions = "ticlSeedingGlobal",
    itername = "CLUE3DBarrel",
    detector = "Barrel",
    patternRecognitionBy = "CLUE3D",
    layer_clusters_barrel_tiles = "ticlLayerTileBarrel:ticlLayerTilesBarrel",
    pluginPatternRecognitionByCLUE3D = dict(
        algo_verbosity = 9999,
        criticalDensity = [0.5, 0.5, 0.5],
        criticalSelfDensity = [0., 0., 0.],
        criticalEtaPhiDistance = [3 * 0.0175, 3 * 0.087, 3 * 0.087, 3 * 0.087],
        nearestHigherOnSameLayer = False,
        densityOnSameLayer = False,
        minNumLayerCluster = [1, 1, 1],
        useAbsoluteProjectiveScale = False,
        densitySiblingLayers = [2, 4, 4]
    )
)

ticlTrackstersLinksBarrel = _tracksterLinksProducer.clone(
    layer_clusters = ["hgcalMergeLayerClusters"],
    layer_clustersTime = ["hgcalMergeLayerClusters:timeLayerCluster"],
    original_masks = ["filteredLayerClustersCLUE3DBarrel:CLUE3DBarrel"],
    tracksters_collections = ["ticlTrackstersCLUE3DBarrel"],
    inferenceAlgo = "TracksterBarrelPIDbyDNN",
    regressionAndPid = True
)

ticlCandidateBarrel = _ticlCandidateProducer.clone(
    cutTk = '0 < abs(eta) < 1.48 && pt > 1. && quality("highPurity") && hitPattern().numberOfLostHits("MISSING_OUTER_HITS") < 5',
    detector = "Barrel",
    egamma_tracksterlinks_collections = ["ticlTrackstersLinksBarrel"],
    egamma_tracksters_collections = ["ticlTrackstersLinksBarrel"],
    general_tracksterlinks_collections = ["ticlTrackstersLinksBarrel"],
    general_tracksters_collections = ["ticlTrackstersLinksBarrel"],
    interpretationDescPSet = cms.PSet(
        algo_verbosity = cms.int32(0),
        cutTk = cms.string('0 < abs(eta) < 1.48 && pt > 1. && quality("highPurity") && hitPattern().numberOfLostHits("MISSING_OUTER_HITS") < 5'),
        delta_tk_ts_interface = cms.double(3*0.087),
        delta_tk_ts_layer1 = cms.double(3*0.0175),
        timing_quality_threshold = cms.double(0.5),
        type = cms.string('General'),
    ),
    layer_clusters = "hgcalMergeLayerClusters",
    layer_clustersTime = "hgcalMergeLayerClusters:timeLayerCluster",
    muons = "muons1stStep",
    original_masks = ["filteredLayerClustersCLUE3DBarrel:CLUE3DBarrel"],
    propagator = "PropagatorWithMaterial",
    timingQualityThreshold = 0.5,
    timingSoA = "mtdSoA",
    tracks = "generalTracks",
    useMTDTiming = False,
    useTimingAverage = False
)

ticlCLUE3DBarrelTask = cms.Task(ticlSeedingGlobal
    ,ticlTrackstersCLUE3DBarrel
)

ticlTrackstersLinksBarrelTask = cms.Task(ticlTrackstersLinksBarrel)

ticlCandidateBarrelTask = cms.Task(ticlCandidateBarrel)

 
