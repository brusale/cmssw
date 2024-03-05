import FWCore.ParameterSet.Config as cms

from RecoParticleFlow.PFClusterProducer.barrelLayerClusters_cfi import barrelLayerClusters


# ECAL
barrelLayerClusters.ebplugin.kappa = cms.double(3.5)
barrelLayerClusters.ebplugin.deltac = cms.vdouble(1.8*0.0175, 5*0.087, 5*0.087)
barrelLayerClusters.ebplugin.outlierDeltaFactor = cms.double(2.7*0.0175)
barrelLayerClusters.ebplugin.maxLayerIndex = cms.int32(0)
#barrelLayerClusters.ebplugin.fractionCutoff = cms.double(0.01)

# HCAL
barrelLayerClusters.hbplugin.kappa = cms.double(1)
barrelLayerClusters.hbplugin.deltac = cms.vdouble(1.8 * 0.0175, 3 * 0.087, 3 * 0.087)
barrelLayerClusters.hbplugin.outlierDeltaFactor = cms.double(5 * 0.087)
barrelLayerClusters.hbplugin.maxLayerIndex = cms.int32(4)


barrelLayerClusters.hoplugin.kappa = cms.double(1)
barrelLayerClusters.hoplugin.deltac = cms.vdouble(1.8 * 0.0175, 3 * 0.087, 3 * 0.087)
barrelLayerClusters.hoplugin.outlierDeltaFactor = cms.double(5 * 0.087)
barrelLayerClusters.hoplugin.maxLayerIndex = cms.int32(5)
