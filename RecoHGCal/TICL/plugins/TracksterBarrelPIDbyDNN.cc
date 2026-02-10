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
      onnxEnergyRuntimeInstance_(std::make_unique<cms::Ort::ONNXRuntime>(
        conf.getParameter<edm::FileInPath>("energyModelPath").fullPath().c_str())),
      inputNames_(conf.getParameter<std::vector<std::string>>("inputNames")),
      output_id_(conf.getParameter<std::vector<std::string>>("output_id")),
      output_energy_(conf.getParameter<std::vector<std::string>>("output_energy")) {
    onnxPIDSession_ = onnxPIDRuntimeInstance_.get();
    onnxEnergySession_ = onnxEnergyRuntimeInstance_.get();
  }

  void TracksterBarrelPIDbyDNN::inputData(const std::vector<reco::CaloCluster>& layerClusters,
                                          std::vector<Trackster>& tracksters, 
                                          const hgcal::RecHitTools& rhtools) {
    // process input data
    batchSize_ = static_cast<int>(tracksters.size());
    if (batchSize_ == 0)
      return;
    
    std::vector<int64_t> inputShape = {batchSize_, eidNFeatures_};
    input_shapes_ = {inputShape};

    input_Data_.clear();
    input_Data_.emplace_back(batchSize_ * eidNFeatures_ , 0);

    for (int i = 0; i < batchSize_; i++) {
      const Trackster& trackster = tracksters[i];
      auto index = eidNFeatures_*i;
      auto tst_barycenter = trackster.barycenter();
      input_Data_[0][index] = static_cast<float>(sqrt(tst_barycenter.x()*tst_barycenter.x() + tst_barycenter.y()*tst_barycenter.y()));
      input_Data_[0][index + 1] = static_cast<float>(trackster.raw_energy());
      input_Data_[0][index + 2] = static_cast<float>(trackster.raw_em_energy());
      input_Data_[0][index + 3] = static_cast<float>(tst_barycenter.z());
      input_Data_[0][index + 4] = static_cast<float>(tst_barycenter.eta());
      input_Data_[0][index + 5] = static_cast<float>(tst_barycenter.phi());
      std::vector<int> vertices_layerids;
      for (auto idx : trackster.vertices()) {
        auto firstHitDetId = layerClusters[idx].hitsAndFractions()[0].first;
        int layerId = rhtools_.getLayerWithOffset(firstHitDetId);
        vertices_layerids.push_back(layerId);
      } 
      int minLayer = *std::min_element(vertices_layerids.begin(), vertices_layerids.end());
      int maxLayer = *std::max_element(vertices_layerids.begin(), vertices_layerids.end());
      input_Data_[0][index + 6] = static_cast<float>(minLayer);
      input_Data_[0][index + 7] = static_cast<float>(maxLayer);
    }  
  }

  void TracksterBarrelPIDbyDNN::runInference(std::vector<Trackster>& tracksters) {
    if (batchSize_ == 0)
      return;

    auto pidOutput = onnxPIDSession_->run(inputNames_, input_Data_, input_shapes_, output_id_, batchSize_);
    auto pidOutputTensor = pidOutput[0];
    float *probs = pidOutputTensor.data();
    if (!output_id_.empty()) {
      for (int i = 0; i < batchSize_; i++) {
        tracksters[i].setIdProbability(Trackster::ParticleType::photon, probs[i]);
      }
    }

    auto energyOutput = onnxEnergySession_->run(inputNames_, input_Data_, input_shapes_, output_energy_, batchSize_);
    auto energyOutputTensor = energyOutput[0];
    float *energies = energyOutputTensor.data();
    if (!output_energy_.empty()) {
      for (int i = 0; i < batchSize_; i++) {
        tracksters[i].setRegressedEnergy(energies[i]);
      }
    }
  }

  void TracksterBarrelPIDbyDNN::fillPSetDescription(edm::ParameterSetDescription& iDesc) {
    iDesc.add<int>("algo_verbosity", 0);
    iDesc.add<edm::FileInPath>("onnxPIDModelPath", edm::FileInPath("model_pid.onnx"));
    iDesc.add<edm::FileInPath>("energyModelPath", edm::FileInPath("model_energy.onnx"));
    iDesc.add<std::vector<std::string>>("inputNames", {"input"});
    iDesc.add<std::vector<std::string>>("output_id", {"output_id"});
    iDesc.add<std::vector<std::string>>("output_energy", {"output_energy"});
  }
}
