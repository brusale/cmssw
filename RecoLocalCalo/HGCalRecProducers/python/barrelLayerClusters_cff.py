import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.barrelLayerClusters_cfi import barrelLayerClusters 

barrelLayerClusters.ebplugin.maxLayerIndex = cms.int32(0)
barrelLayerClusters.ebplugin.rhoc = cms.double(30)

barrelLayerClusters.hbplugin.maxLayerIndex = cms.int32(4)
barrelLayerClusters.hbplugin.rhoc = cms.double(30)

barrelLayerClusters.hoplugin.maxLayerIndex = cms.int32(5)
barrelLayerClusters.hoplugin.rhoc = cms.double(30)


