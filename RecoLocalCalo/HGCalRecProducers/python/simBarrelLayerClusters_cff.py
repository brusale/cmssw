import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.simBarrelLayerClusters_cfi import simBarrelLayerClusters 

simBarrelLayerClusters.ebplugin.maxLayerIndex = cms.int32(0)
simBarrelLayerClusters.ebplugin.rhoc = cms.double(30)

simBarrelLayerClusters.hbplugin.maxLayerIndex = cms.int32(3)
simBarrelLayerClusters.hbplugin.rhoc = cms.double(30)

simBarrelLayerClusters.hoplugin.maxLayerIndex = cms.int32(0)
simBarrelLayerClusters.hoplugin.rhoc = cms.double(30)


