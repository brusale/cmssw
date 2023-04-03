import FWCore.ParameterSet.Config as cms

barrelLayerClusterSimClusterAssociation = cms.EDProducer("BarrelLCToSCAssociatorEDProducer",
  associator = cms.InputTag('barrelSCAssocByEnergyScoreProducer'),
  label_scl = cms.InputTag("mix", "MergedCaloTruth"),
  label_lcl = cms.InputTag("barrelLayerClusters")
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(barrelLayerClusterSimClusterAssociation,
  label_scl = "mixData:MergedCaloTruth"
)

