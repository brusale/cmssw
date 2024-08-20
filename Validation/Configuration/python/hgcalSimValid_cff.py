import FWCore.ParameterSet.Config as cms

from SimCalorimetry.HGCalSimProducers.hgcHitAssociation_cfi import lcAssocByEnergyScoreProducer, scAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.simTracksterAssociatorByEnergyScore_cfi import simTracksterAssociatorByEnergyScore as simTsAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.layerClusterSimTracksterAssociatorByEnergyScore_cfi import layerClusterSimTracksterAssociatorByEnergyScore as lcSimTSAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToCPAssociation_cfi import layerClusterCaloParticleAssociation as layerClusterCaloParticleAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToCPAssociation_cfi import barrelLayerClusterCaloParticleAssociation as barrelLayerClusterCaloParticleAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.hgcalSimTracksterHitLCAssociatorByEnergyScore_cfi import hgcalSimTracksterHitLCAssociatorByEnergyScore as simTracksterHitLCAssociatorByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cfi import layerClusterSimClusterAssociation as layerClusterSimClusterAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cfi import barrelLayerClusterSimClusterAssociation as barrelLayerClusterSimClusterAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToSimTSAssociation_cfi import layerClusterSimTracksterAssociation as layerClusterSimTracksterAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToCPAssociation_cfi import layerClusterCaloParticleAssociationHFNose as layerClusterCaloParticleAssociationProducerHFNose
from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cfi import layerClusterSimClusterAssociationHFNose as layerClusterSimClusterAssociationProducerHFNose
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationLinking, tracksterSimTracksterAssociationPR,tracksterSimTracksterAssociationLinkingbyCLUE3D, tracksterSimTracksterAssociationPRbyCLUE3D, tracksterSimTracksterAssociationLinkingPU, tracksterSimTracksterAssociationPRPU #, tracksterSimTracksterAssociationLinkingbyCLUE3DEM, tracksterSimTracksterAssociationLinkingbyCLUE3DHAD, tracksterSimTracksterAssociationPRbyCLUE3DEM, tracksterSimTracksterAssociationPRbyCLUE3DHAD
from RecoHGCal.TICL.mergedTrackstersProducer_cfi import mergedTrackstersProducer as _mergedTrackstersProducer
from SimCalorimetry.HGCalAssociatorProducers.SimTauProducer_cfi import *

from SimCalorimetry.HGCalAssociatorProducers.barrelLCToCPAssociatorByEnergyScoreProducer_cfi import barrelLCToCPAssociatorByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.barrelLCToSCAssociatorByEnergyScoreProducer_cfi import barrelLCToSCAssociatorByEnergyScoreProducer
from Validation.HGCalValidation.simhitValidation_cff    import *
from Validation.HGCalValidation.digiValidation_cff      import *
from Validation.HGCalValidation.rechitValidation_cff    import *
from Validation.HGCalValidation.hgcalHitValidation_cff  import *
from RecoHGCal.TICL.SimTracksters_cff import *


from Validation.HGCalValidation.HGCalValidator_cfi import hgcalValidator
from Validation.RecoParticleFlow.PFJetValidation_cff import pfJetValidation1 as _hgcalPFJetValidation

from Validation.HGCalValidation.ticlPFValidation_cfi import ticlPFValidation
hgcalTiclPFValidation = cms.Sequence(ticlPFValidation)

from Validation.HGCalValidation.ticlTrackstersEdgesValidation_cfi import ticlTrackstersEdgesValidation
hgcalTiclTrackstersEdgesValidationSequence = cms.Sequence(ticlTrackstersEdgesValidation)

hgcalValidatorSequence = cms.Sequence(hgcalValidator)
hgcalPFJetValidation = _hgcalPFJetValidation.clone(BenchmarkLabel = 'PFJetValidation/HGCAlCompWithGenJet',
    VariablePtBins=[10., 30., 80., 120., 250., 600.],
    DeltaPtOvPtHistoParameter = dict(EROn=True,EREtaMax=3.0, EREtaMin=1.6, slicingOn=True))

hgcalAssociators = cms.Task(lcAssocByEnergyScoreProducer, layerClusterCaloParticleAssociationProducer,
                            scAssocByEnergyScoreProducer, layerClusterSimClusterAssociationProducer,
                            barrelLCToCPAssociatorByEnergyScoreProducer, barrelLayerClusterCaloParticleAssociationProducer,
                            barrelLCToSCAssociatorByEnergyScoreProducer, barrelLayerClusterSimClusterAssociationProducer, 
                            lcSimTSAssocByEnergyScoreProducer, layerClusterSimTracksterAssociationProducer,
                            simTsAssocByEnergyScoreProducer,  simTracksterHitLCAssociatorByEnergyScoreProducer,
                            tracksterSimTracksterAssociationLinking, tracksterSimTracksterAssociationPR,
                            tracksterSimTracksterAssociationLinkingbyCLUE3D, tracksterSimTracksterAssociationPRbyCLUE3D,
                            tracksterSimTracksterAssociationLinkingPU, tracksterSimTracksterAssociationPRPU,
                            SimTauProducer
                            )

from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
''' For future separate iterations
mergedTrackstersProducer = _mergedTrackstersProducer.clone()
ticl_v5.toModify(hgcalAssociators, lambda x: x.add(mergedTrackstersProducer, tracksterSimTracksterAssociationLinkingbyCLUE3DEM, tracksterSimTracksterAssociationLinkingbyCLUE3DHAD, tracksterSimTracksterAssociationPRbyCLUE3DEM, tracksterSimTracksterAssociationPRbyCLUE3DHAD))
'''

hgcalValidation = cms.Sequence(hgcalSimHitValidationEE
                               + hgcalSimHitValidationHEF
                               + hgcalSimHitValidationHEB
                               + hgcalDigiValidationEE
                               + hgcalDigiValidationHEF
                               + hgcalDigiValidationHEB
                               + hgcalRecHitValidationEE
                               + hgcalRecHitValidationHEF
                               + hgcalRecHitValidationHEB
                               + hgcalHitValidationSequence
                               + hgcalValidatorSequence
                               + hgcalTiclPFValidation
                               #Currently commented out until trackster edges are saved
#                               + hgcalTiclTrackstersEdgesValidationSequence
                               + hgcalPFJetValidation)

_hfnose_hgcalAssociatorsTask = hgcalAssociators.copy()
_hfnose_hgcalAssociatorsTask.add(layerClusterCaloParticleAssociationProducerHFNose, layerClusterSimClusterAssociationProducerHFNose)
