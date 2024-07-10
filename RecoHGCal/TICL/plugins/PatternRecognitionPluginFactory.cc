#include "RecoHGCal/TICL/plugins/PatternRecognitionPluginFactory.h"
#include "PatternRecognitionbyCA.h"
#include "PatternRecognitionbyCLUE3D.h"
#include "PatternRecognitionbyFastJet.h"
#include "PatternRecognitionbyPassthrough.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"

EDM_REGISTER_VALIDATED_PLUGINFACTORY(PatternRecognitionFactory, "PatternRecognitionFactory");
EDM_REGISTER_VALIDATED_PLUGINFACTORY(PatternRecognitionHFNoseFactory, "PatternRecognitionHFNoseFactory");
EDM_REGISTER_VALIDATED_PLUGINFACTORY(PatternRecognitionHCALFactory, "PatternRecognitionHCALFactory");
EDM_REGISTER_VALIDATED_PLUGINFACTORY(PatternRecognitionECALFactory, "PatternRecognitionECALFactory");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactory, ticl::PatternRecognitionbyCA<TICLLayerTiles>, "CA");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactory, ticl::PatternRecognitionbyCLUE3D<TICLLayerTiles>, "CLUE3D");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactory, ticl::PatternRecognitionbyFastJet<TICLLayerTiles>, "FastJet");
// DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionHCALFactory, ticl::PatternRecognitionbyCA<TICLLayerTilesHCAL>, "CA");
// DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionHCALFactory, ticl::PatternRecognitionbyCLUE3D<TICLLayerTilesHCAL>, "CLUE3D");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionHCALFactory, ticl::PatternRecognitionbyFastJet<TICLLayerTilesHCAL>, "FastJet");
// DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionECALFactory, ticl::PatternRecognitionbyCA<TICLLayerTilesECAL>, "CA");
// DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionECALFactory, ticl::PatternRecognitionbyCLUE3D<TICLLayerTilesECAL>, "CLUE3D");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionECALFactory, ticl::PatternRecognitionbyFastJet<TICLLayerTilesECAL>, "FastJet");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactory,
                            ticl::PatternRecognitionbyPassthrough<TICLLayerTiles>,
                            "Passthrough");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionHFNoseFactory, ticl::PatternRecognitionbyCA<TICLLayerTilesHFNose>, "CA");
