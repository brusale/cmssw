import FWCore.ParameterSet.Config as cms

from Validation.HGCalValidation.HGCalSimHitsClient_cff import *
from Validation.HGCalValidation.HGCalDigiClient_cff    import *
from Validation.HGCalValidation.HGCalRecHitsClient_cff import *
from Validation.HGCalValidation.PostProcessorHGCAL_cfi import postProcessorHGCALlayerclusters,postProcessorHGCALsimclusters,postProcessorHGCALTracksters
from Validation.HGCalValidation.PostProcessorBarrel_cfi import postProcessorBarrellayerclusters
from Validation.HGCalValidation.SimPostProcessorBarrel_cfi import simPostProcessorBarrellayerclusters 
from Validation.HGCalValidation.PostProcessorBarrelPFCluster_cfi import postProcessorBarrelPFClusters
hgcalPostProcessor = cms.Sequence(hgcalSimHitClientEE
    + hgcalSimHitClientHEF
    + hgcalSimHitClientHEB
    + hgcalDigiClientEE
    + hgcalDigiClientHEF
    + hgcalDigiClientHEB
    + hgcalRecHitClientEE
    + hgcalRecHitClientHEF
    + hgcalRecHitClientHEB)

hgcalValidatorPostProcessor = cms.Sequence(
    postProcessorHGCALlayerclusters+
    postProcessorBarrellayerclusters+
    simPostProcessorBarrellayerclusters+
    postProcessorBarrelPFClusters+
    postProcessorHGCALsimclusters+
    postProcessorHGCALTracksters)
