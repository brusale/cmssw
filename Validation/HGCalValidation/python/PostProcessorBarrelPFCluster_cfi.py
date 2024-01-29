import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
from RecoHGCal.TICL.iterativeTICL_cff import ticlIterLabelsMerge
from Validation.HGCalValidation.BarrelValidatorPFCluster_cfi import barrelValidatorPFCluster

prefix = "PFCluster/PFClusterValidator/"
maxLayers = 6

eff_layers = ["effic_eta_layer{:02d} 'LayerCluster Efficiency vs #eta Layer{:02d}' Num_CaloParticle_Eta_perlayer{:02d} Denom_CaloParticle_Eta_perlayer{:02d}".format(i, i, i, i)  for i in range(maxLayers)]
eff_layers.extend(["effic_phi_layer{:02d} 'LayerCluster Efficiency vs #phi Layer{:02d}' Num_CaloParticle_Phi_perlayer{:02d} Denom_CaloParticle_Phi_perLayer{:02d}".format(i, i, i, i) for i in range(maxLayers)])
eff_layers.extend(["duplicate_eta_layer{:02d} 'LayerCluster Duplicate(Split) Rate vs #eta Layer{:02d}' NumDup_CaloParticle_Eta_perlayer{:02d} Denom_CaloParticle_Eta_perlayer{:02d}".format(i, i, i, i) for i in range(maxLayers)])
eff_layers.extend(["duplicate_phi_layer{:02d} 'LayerCluster Duplicate(Split) Rate vs #phi Layer{:02d}' NumDup_CaloParticle_Phi_perlayer{:02d} Denom_CaloParticle_Phi_perlayer{:02d}".format(i, i, i, i) for i in range(maxLayers)])
eff_layers.extend(["fake_eta_layer{:02d} 'LayerCluster Fake Rate vs #eta Layer{:02d}' Num_LayerCluster_Eta_perlayer{:02d} Denom_LayerCluster_Eta_perlayer{:02d} fake".format(i, i, i, i) for i in range(maxLayers)])
eff_layers.extend(["fake_phi_layer{:02d} 'LayerCluster Fake Rate vs #phi Layer{:02d}' Num_LayerCluster_Phi{:02d} Denom_LayerCluster_Phi_perlayer{:02d} fake".format(i, i, i, i) for i in range(maxLayers)])
eff_layers.extend(["merge_eta_layer{:02d} 'LayerCluster Merge Rate vs #eta Layer{:02d}' NumMerge_LayerCluster_Eta_perlayer{:02d} Denom_LayerCluster_Phi_perlayer{:02d}".format(i, i, i, i) for i in range(maxLayers)])
eff_layers.extend(["merge_phi_layer{:02d} 'LayerCluster Merge Rate vs #phi Layer{:02d}' NumMerge_LayerCluster_Phi_perlayer{:02d} Denom_LayerCluster_Phi_perlayer{:02d}".format(i, i, i, i) for i in range(maxLayers)])
lcToCP_linking = barrelValidatorPFCluster.label_LCToCPLinking._InputTag__moduleLabel

postProcessorBarrelPFClusters = DQMEDHarvester('DQMGenericClient',
    subDirs = cms.untracked.vstring(prefix + barrelValidatorPFCluster.label_layerClusterPlots._InputTag__moduleLabel + '/' + lcToCP_linking),
    efficiency = cms.vstring(eff_layers),
    resolution = cms.vstring(),
    cumulativeDists = cms.untracked.vstring(),
    noFlowDists = cms.untracked.vstring(),
    outputFileName = cms.untracked.string(""),
    verbose = cms.untracked.uint32(4))

