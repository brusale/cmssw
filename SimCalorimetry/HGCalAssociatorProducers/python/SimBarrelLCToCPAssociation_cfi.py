import FWCore.ParameterSet.Config as cms

simBarrelLayerClusterCaloParticleAssociation = cms.EDProducer("SimBarrelLCToCPAssociatorEDProducer",
    associator = cms.InputTag('simBarrelLCAssocByEnergyScoreProducer'),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
    label_lc = cms.InputTag("simBarrelLayerClusters")
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(simBarrelLayerClusterCaloParticleAssociation,
    label_cp = "mixData:MergedCaloTruth"
)

