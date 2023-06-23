import FWCore.ParameterSet.Config as cms

from SimCalorimetry.HGCalSimProducers.hgcHitAssociation_cfi import lcAssocByEnergyScoreProducer, scAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.simTracksterAssociatorByEnergyScore_cfi import simTracksterAssociatorByEnergyScore as simTsAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.layerClusterSimTracksterAssociatorByEnergyScore_cfi import layerClusterSimTracksterAssociatorByEnergyScore as lcSimTSAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToCPAssociation_cfi import layerClusterCaloParticleAssociation as layerClusterCaloParticleAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.simTracksterHitLCAssociatorByEnergyScore_cfi import simTracksterHitLCAssociatorByEnergyScore as simTracksterHitLCAssociatorByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cfi import layerClusterSimClusterAssociation as layerClusterSimClusterAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToSimTSAssociation_cfi import layerClusterSimTracksterAssociation as layerClusterSimTracksterAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToCPAssociation_cfi import layerClusterCaloParticleAssociationHFNose as layerClusterCaloParticleAssociationProducerHFNose
from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cfi import layerClusterSimClusterAssociationHFNose as layerClusterSimClusterAssociationProducerHFNose
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationLinking, tracksterSimTracksterAssociationPR,tracksterSimTracksterAssociationLinkingbyCLUE3D, tracksterSimTracksterAssociationPRbyCLUE3D
from SimCalorimetry.HGCalAssociatorProducers.BarrelLCToCPAssociation_cfi import barrelLayerClusterCaloParticleAssociation as barrelLayerClusterCaloParticleAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.barrelLayerClusterAssociatorByEnergyScore_cfi import barrelLayerClusterAssociatorByEnergyScore as barrelLCAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.BarrelLCToSCAssociation_cfi import barrelLayerClusterSimClusterAssociation as barrelLayerClusterSimClusterAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.barrelSimClusterAssociatorByEnergyScore_cfi import barrelSimClusterAssociatorByEnergyScore as barrelSCAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.SimBarrelLCToCPAssociation_cfi import simBarrelLayerClusterCaloParticleAssociation as simBarrelLayerClusterCaloParticleAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.simBarrelLayerClusterAssociatorByEnergyScore_cfi import simBarrelLayerClusterAssociatorByEnergyScore as simBarrelLCAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.SimBarrelLCToSCAssociation_cfi import simBarrelLayerClusterSimClusterAssociation as simBarrelLayerClusterSimClusterAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.simBarrelSimClusterAssociatorByEnergyScore_cfi import simBarrelSimClusterAssociatorByEnergyScore as simBarrelSCAssocByEnergyScoreProducer

from Validation.HGCalValidation.simhitValidation_cff    import *
from Validation.HGCalValidation.digiValidation_cff      import *
from Validation.HGCalValidation.rechitValidation_cff    import *
from Validation.HGCalValidation.hgcalHitValidation_cff  import *
from RecoHGCal.TICL.SimTracksters_cff import *


from Validation.HGCalValidation.HGCalValidator_cfi import hgcalValidator
from Validation.HGCalValidation.BarrelValidator_cfi import barrelValidator
from Validation.HGCalValidation.SimBarrelValidator_cfi import simBarrelValidator
from Validation.RecoParticleFlow.PFJetValidation_cff import pfJetValidation1 as _hgcalPFJetValidation

from Validation.HGCalValidation.ticlPFValidation_cfi import ticlPFValidation
hgcalTiclPFValidation = cms.Sequence(ticlPFValidation)

from Validation.HGCalValidation.ticlTrackstersEdgesValidation_cfi import ticlTrackstersEdgesValidation
hgcalTiclTrackstersEdgesValidationSequence = cms.Sequence(ticlTrackstersEdgesValidation)

hgcalValidatorSequence = cms.Sequence(hgcalValidator
				      + barrelValidator
				      + simBarrelValidator
				      )
hgcalPFJetValidation = _hgcalPFJetValidation.clone(BenchmarkLabel = 'PFJetValidation/HGCAlCompWithGenJet',
    VariablePtBins=[10., 30., 80., 120., 250., 600.],
    DeltaPtOvPtHistoParameter = dict(EROn=True,EREtaMax=3.0, EREtaMin=1.6, slicingOn=True))

hgcalAssociators = cms.Task(lcAssocByEnergyScoreProducer, layerClusterCaloParticleAssociationProducer,
                            scAssocByEnergyScoreProducer, layerClusterSimClusterAssociationProducer,
                            lcSimTSAssocByEnergyScoreProducer, layerClusterSimTracksterAssociationProducer,
                            simTsAssocByEnergyScoreProducer,  simTracksterHitLCAssociatorByEnergyScoreProducer,
                            tracksterSimTracksterAssociationLinking, tracksterSimTracksterAssociationPR,
                            tracksterSimTracksterAssociationLinkingbyCLUE3D, tracksterSimTracksterAssociationPRbyCLUE3D,
			    barrelLCAssocByEnergyScoreProducer, barrelLayerClusterCaloParticleAssociationProducer,
			    barrelSCAssocByEnergyScoreProducer, barrelLayerClusterSimClusterAssociationProducer,
			    simBarrelLCAssocByEnergyScoreProducer, simBarrelLayerClusterCaloParticleAssociationProducer,
			    simBarrelSCAssocByEnergyScoreProducer, simBarrelLayerClusterSimClusterAssociationProducer
		  )

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
