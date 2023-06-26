
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalLayerClusterAlgoFactory.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalClusteringAlgoBase.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalImagingAlgo.h"
#include "RecoLocalCalo/HGCalRecProducers/plugins/HGCalCLUEAlgo.h"
#include "RecoLocalCalo/HGCalRecProducers/plugins/BarrelCLUEAlgo.h"
#include "RecoLocalCalo/HGCalRecProducers/plugins/SimBarrelCLUEAlgo.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"

EDM_REGISTER_VALIDATED_PLUGINFACTORY(HGCalLayerClusterAlgoFactory, "HGCalLayerClusterAlgoFactory");
DEFINE_EDM_VALIDATED_PLUGIN(HGCalLayerClusterAlgoFactory, HGCalImagingAlgo, "Imaging");
DEFINE_EDM_VALIDATED_PLUGIN(HGCalLayerClusterAlgoFactory, HGCalSiCLUEAlgo, "SiCLUE");
DEFINE_EDM_VALIDATED_PLUGIN(HGCalLayerClusterAlgoFactory, HGCalSciCLUEAlgo, "SciCLUE");
DEFINE_EDM_VALIDATED_PLUGIN(HGCalLayerClusterAlgoFactory, HFNoseCLUEAlgo, "HFNoseCLUE");
DEFINE_EDM_VALIDATED_PLUGIN(HGCalLayerClusterAlgoFactory, EBCLUEAlgo, "EBCLUE");
DEFINE_EDM_VALIDATED_PLUGIN(HGCalLayerClusterAlgoFactory, HBCLUEAlgo, "HBCLUE");
DEFINE_EDM_VALIDATED_PLUGIN(HGCalLayerClusterAlgoFactory, HOCLUEAlgo, "HOCLUE");
DEFINE_EDM_VALIDATED_PLUGIN(HGCalLayerClusterAlgoFactory, SimEBCLUEAlgo, "SimEBCLUE");
DEFINE_EDM_VALIDATED_PLUGIN(HGCalLayerClusterAlgoFactory, SimHBCLUEAlgo, "SimHBCLUE");
DEFINE_EDM_VALIDATED_PLUGIN(HGCalLayerClusterAlgoFactory, SimHOCLUEAlgo, "SimHOCLUE");
