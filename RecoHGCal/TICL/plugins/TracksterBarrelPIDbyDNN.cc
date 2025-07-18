#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "RecoHGCal/TICL/interface/TracksterBarrelPIDbyDNN.h"
#include "RecoHGCal/TICL/interface/TracksterInferenceAlgoFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"

namespace ticl {
  using namespace cms::Ort;

  TracksterBarrelPIDbyDNN::TracksterBarrelPIDbyDNN(const edm::ParameterSet& conf)
    : TracksterInferenceAlgoBase(conf),
      onnxPIDRuntimeInstance_(std::make_unique<cms::Ort::ONNXRuntime>(
        conf.getParameter<edm::FileInPath>("onnxPIDModelPath").fullPath().c_str())),
      inputNames_(conf.getParameter<std::vector<std::string>>("inputNames")),
      output_id_(conf.getParameter<std::vector<std::string>>("output_id")) {
    onnxPIDSession_ = onnxPIDRuntimeInstance_.get();
  }

  void TracksterBarrelPIDbyDNN::inputData(const std::vector<reco::CaloCluster>& layerClusters,
                                          std::vector<Trackster>& tracksters) {
    // process input data
    batchSize_ = static_cast<int>(tracksters.size());
    if (batchSize_ == 0)
      return;
    
    std::vector<int64_t> inputShape = {batchSize_, eidNFeatures_};
    input_shapes_ = {inputShape};

    input_Data_.clear();
    input_Data_.emplace_back(batchSize_ * eidNFeatures_ , 0);

    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    std::cout << "inputData.size(): " << input_Data_.size() << std::endl;
    std::cout << "inputNames.size(): " << inputNames_.size() << std::endl;
    for (const auto& data : input_Data_) {
      std::cout << "input_Data_[ ]: " << data.size() << std::endl;
    }
    for (int i = 0; i < batchSize_; i++) {
      const Trackster& trackster = tracksters[i];
      auto index = eidNFeatures_*i;
      auto tst_barycenter = trackster.barycenter();
      std::cout << "inputData_[0].size(): " << input_Data_[0].size() << std::endl;
      input_Data_[0][index] = sqrt(tst_barycenter.x()*tst_barycenter.x() + tst_barycenter.y()*tst_barycenter.y());
      input_Data_[0][index + 1] = trackster.raw_energy();
      input_Data_[0][index + 2] = trackster.raw_em_energy();
      input_Data_[0][index + 3] = tst_barycenter.z();
      input_Data_[0][index + 4] = tst_barycenter.eta();
      input_Data_[0][index + 5] = tst_barycenter.phi();
      std::vector<int> vertices_layerids;
      for (auto idx : trackster.vertices()) {
        auto firstHitDetId = layerClusters[idx].hitsAndFractions()[0].first;
        int layerId = rhtools_.getLayerWithOffset(firstHitDetId);
        vertices_layerids.push_back(layerId);
      } 
      int minLayer = *std::min_element(vertices_layerids.begin(), vertices_layerids.end());
      int maxLayer = *std::max_element(vertices_layerids.begin(), vertices_layerids.end());
      input_Data_[0][index + 6] = minLayer;
      input_Data_[0][index + 7] = maxLayer;
    }  
  }

  void TracksterBarrelPIDbyDNN::runInference(std::vector<Trackster>& tracksters) {
    if (batchSize_ == 0)
      return;

    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    auto pidOutput = onnxPIDSession_->run(inputNames_, input_Data_, input_shapes_, output_id_, batchSize_);
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    auto pidOutputTensor = pidOutput[0];
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    float *probs = pidOutputTensor.data();
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    if (!output_id_.empty()) {
      for (int i = 0; i < batchSize_; i++) {
        tracksters[tracksterIndices_[i]].setProbabilities(probs);
        probs += tracksters[tracksterIndices_[i]].id_probabilities().size();
      }
    }
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
  }

  void TracksterBarrelPIDbyDNN::fillPSetDescription(edm::ParameterSetDescription& iDesc) {
    iDesc.add<int>("algo_verbosity", 0);
    iDesc.add<edm::FileInPath>("onnxPIDModelPath", edm::FileInPath("test_2025.onnx"));
    iDesc.add<std::vector<std::string>>("inputNames", {"input"});
    iDesc.add<std::vector<std::string>>("output_id", {"output"});
  }
}
