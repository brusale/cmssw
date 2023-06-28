#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Validation/HGCalValidation/interface/HGCalValidator.h"
#include "Validation/HGCalValidation/interface/BarrelValidator.h"
#include "Validation/HGCalValidation/interface/SimBarrelValidator.h"

DEFINE_FWK_MODULE(HGCalValidator);
DEFINE_FWK_MODULE(BarrelValidator);
DEFINE_FWK_MODULE(SimBarrelValidator);
