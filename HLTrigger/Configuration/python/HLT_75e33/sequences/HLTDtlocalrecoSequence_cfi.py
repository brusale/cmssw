import FWCore.ParameterSet.Config as cms

from ..modules.hltDt1DRecHits_cfi import *
from ..modules.hltDt4DSegments_cfi import *

HLTDtlocalrecoSequence = cms.Sequence(hltDt1DRecHits+hltDt4DSegments)
