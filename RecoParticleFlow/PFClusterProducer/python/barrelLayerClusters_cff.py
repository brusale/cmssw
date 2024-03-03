import FWCore.ParameterSet.Config as cms

from RecoParticleFlow.PFClusterProducer.barrelLayerClusters_cfi import barrelLayerClusters

# HCAL
barrelLayerClusters.hbplugin.kappa = cms.double(1)
barrelLayerClusters.hbplugin.deltac = cms.vdouble(1.8 * 0.0175, 3 * 0.087, 3 * 0.087)
barrelLayerClusters.hbplugin.outlierDeltaFactor = cms.double(5 * 0.087)
barrelLayerClusters.hbplugin.maxLayerIndex = cms.int32(4)


barrelLayerClusters.hoplugin.kappa = cms.double(1)
barrelLayerClusters.hoplugin.deltac = cms.vdouble(1.8 * 0.0175, 3 * 0.087, 3 * 0.087)
barrelLayerClusters.hoplugin.outlierDeltaFactor = cms.double(5 * 0.087)
barrelLayerClusters.hoplugin.maxLayerIndex = cms.int32(5)
