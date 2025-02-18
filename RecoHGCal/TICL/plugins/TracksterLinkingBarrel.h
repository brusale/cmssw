#ifndef RecoHGCal_TICL_TracksterLinkingBarrel_H
#define RecoHGCal_TICL_TracksterLinKingBarrel_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include <array>

namespace ticl {
  class TracksterLinkingBarrel : public TracksterLinkingAlgoBase {
    public:
      TracksterLinkingBarrel(const edm::ParameterSet& conf, 
                             edm::ConsumesCollector iC,
                             cms::Ort::ONNXRuntime const* onnxRuntime = nullptr);

      ~TracksterLinkingBarrel() override {}

      void linkTracksters(const Inputs& input,
                          std::vector<Trackster>& resultTracksters,
                          std::vector<std::vector<unsigned int>>& linkedResultTracksters,
                          std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId) override;

      void initialize(const HGCalDDDConstants* hgcons,
                      const hgcal::RecHitTools rhtools,
                      const edm::ESHandle<MagneticField> bfieldH,
                      const edm::ESHandle<Propagator> propH) override;

      static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {
        iDesc.add<double>("min_cos_theta", 0.95);
        TracksterLinkingAlgoBase::fillPSetDescription(iDesc);
      }
    private:
      const HGCalDDDConstants* hgcons_;
      hgcal::RecHitTools rhtools_;
      edm::ESHandle<MagneticField> bfield_;
      edm::ESHandle<Propagator> propagator_;
      std::array<float, TileConstants::nEtaBins> eta_windows_;
      std::array<float, TileConstants::nPhiBins> phi_windows_;

      float min_cos_theta_;
  };
}

#endif

