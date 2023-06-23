import FWCore.ParameterSet.Config as cms

simBarrelLayerClusterSimClusterAssociation = cms.EDProducer("SimBarrelLCToSCAssociatorEDProducer",
  associator = cms.InputTag('simBarrelSCAssocByEnergyScoreProducer'),
  label_scl = cms.InputTag("mix", "MergedCaloTruth"),
  label_lcl = cms.InputTag("simBarrelLayerClusters")
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(simBarrelLayerClusterSimClusterAssociation,
  label_scl = "mixData:MergedCaloTruth"
)

