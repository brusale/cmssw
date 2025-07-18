#ifndef RecoHGCal_TICL_TracksterBarrelPIDbyDNN_H__
#define RecoHGCal_TICL_TracksterBarrelPIDbyDNN_H__

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "RecoHGCal/TICL/interface/TracksterInferenceAlgoBase.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace ticl {
  class TracksterBarrelPIDbyDNN : public TracksterInferenceAlgoBase {
    public:
      explicit TracksterBarrelPIDbyDNN(const edm::ParameterSet& conf);
      void inputData(const std::vector<reco::CaloCluster>& layerClusters, std::vector<Trackster>& tracksters) override;
      void runInference(std::vector<Trackster>& tracksters) override;
      static void fillPSetDescription(edm::ParameterSetDescription& iDesc);

    private:
      const std::unique_ptr<cms::Ort::ONNXRuntime> onnxPIDRuntimeInstance_;
      const cms::Ort::ONNXRuntime* onnxPIDSession_;

      const std::vector<std::string> inputNames_;
      const std::vector<std::string> output_id_;
      static constexpr int eidNFeatures_ = 8; 
      hgcal::RecHitTools rhtools_;
      std::vector<std::vector<int64_t>> input_shapes_;
      std::vector<int> tracksterIndices_;
      std::vector<std::vector<float>> input_Data_;
      int batchSize_;
  };
}

#endif 
