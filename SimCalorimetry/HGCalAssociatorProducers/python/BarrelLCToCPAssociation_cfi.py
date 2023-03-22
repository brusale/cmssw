import FWCore.ParameterSet.Config as cms

barrelLayerClusterCaloParticleAssociation = cms.EDProducer("LCToCPAssociatorEDProducer",
    associator = cms.InputTag('barrelLCAssocByEnergyScoreProducer'),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
    label_lc = cms.InputTag("barrekLayerClusters")
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(barrelLayerClusterCaloParticleAssociation,
    label_cp = "mixData:MergedCaloTruth"
)

