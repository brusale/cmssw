#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoHGCal/TICL/plugins/TracksterLinkingBarrel.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "RecoHGCal/TICL/plugins/TICLGraph.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include <numbers>

using namespace ticl;

TracksterLinkingBarrel::TracksterLinkingBarrel(const edm::ParameterSet &conf, 
                                               edm::ConsumesCollector iC,
                                               cms::Ort::ONNXRuntime const *onnxRuntime)
  : TracksterLinkingAlgoBase(conf, iC),
    min_cos_theta_(conf.getParameter<double>("min_cos_theta")),
    eta_window_(conf.getParameter<unsigned int>("eta_window")),
    phi_window_(conf.getParameter<unsigned int>("phi_window")),
    min_span_(conf.getParameter<unsigned int>("min_span")) {}

void TracksterLinkingBarrel::initialize(const HGCalDDDConstants *hgcons,
                                        const hgcal::RecHitTools rhtools,
                                        const edm::ESHandle<MagneticField> bfield,
                                        const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;
  bfield_ = bfield;
  propagator_ = propH;

  eta_width_ = (TileConstantsBarrel::maxEta - TileConstantsBarrel::minEta) / TileConstantsBarrel::nEtaBins;
  phi_width_ = 2*std::numbers::pi_v<float> / TileConstantsBarrel::nPhiBins;
}

void TracksterLinkingBarrel::linkTracksters(
    const Inputs &input,
    std::vector<Trackster> &resultTracksters,
    std::vector<std::vector<unsigned int>> &linkedResultTracksters,
    std::vector<std::vector<unsigned int>> &linkedTracksterIdToInputTracksterId) {
  const auto &tracksters = input.tracksters;
  const auto &layerClusters = input.layerClusters;

  // vector of trackster indices sorted by energy
  std::vector<unsigned int> sortedTracksters(tracksters.size());
  std::iota(sortedTracksters.begin(), sortedTracksters.end(), 0);
  std::sort(sortedTracksters.begin(), sortedTracksters.end(), [&tracksters](unsigned int i, unsigned int j) {
    return tracksters[i].raw_energy() > tracksters[j].raw_energy();
  });

  // mask for tracksters that have already been linked
  std::vector<int> maskedTracksters(tracksters.size(), 0);
  std::vector<int> isRootTrackster(tracksters.size(), 1);
  // project barycenter on first layer
  TICLLayerTileBarrel tracksterTile;
  std::vector<ticl::Vector> tracksterBarycenter(tracksters.size());
  for (auto const t_idx : sortedTracksters) {
    tracksterBarycenter[t_idx] = tracksters[t_idx].barycenter();
    tracksterTile.fill(tracksters[t_idx].barycenter().eta(), tracksters[t_idx].barycenter().phi(), t_idx);
  }

  std::vector<ticl::Node> allNodes;
  for (size_t it = 0; it < tracksters.size(); ++it) {
    allNodes.emplace_back(it);
  }


  auto getTracksterSpan = [&](const auto& trackster) {
    std::vector<int> verticesLayerId;
    for (size_t lcId : trackster.vertices()) {
      const auto& layerCluster = layerClusters[lcId];
      auto firstHitDetId = layerCluster.hitsAndFractions()[0].first;
      auto layerId = rhtools_.getLayerWithOffset(firstHitDetId);
      verticesLayerId.push_back(layerId);
    }
    int maxLayer = *std::max_element(verticesLayerId.begin(), verticesLayerId.end());
    int minLayer = *std::min_element(verticesLayerId.begin(), verticesLayerId.end());
    return (maxLayer - minLayer) + 1;
  };

  for (auto const &t_idx : sortedTracksters) {
    //if (maskedTracksters[t_idx]) continue;
    linkedResultTracksters.reserve(t_idx);
    auto const& trackster = tracksters[t_idx];
    unsigned int trackster_span = getTracksterSpan(trackster);
    if (trackster_span < min_span_) continue;
    // if trackster spans only one layer use barycenter
    // otherwise use PCA axis
    auto const& barycenter = (trackster_span == 1) ? trackster.barycenter() : trackster.eigenvectors(0);
    //auto const& barycenter = trackster.barycenter();
    float t_norm = std::sqrt(barycenter.mag2());

    float eta_min = std::max(abs(barycenter.eta()) - eta_window_ * eta_width_, TileConstantsBarrel::minEta);
    float eta_max = std::min(abs(barycenter.eta()) + eta_window_ * eta_width_, TileConstantsBarrel::maxEta);
    float phi_min = std::max(abs(barycenter.phi()) - phi_window_ * phi_width_, -std::numbers::pi_v<float>);
    float phi_max = std::min(abs(barycenter.phi()) + phi_window_ * phi_width_, std::numbers::pi_v<float>);
 
    std::array<int, 4> searchBox = tracksterTile.searchBoxEtaPhi(eta_min, eta_max, phi_min, phi_max);
    if (searchBox[2] > searchBox[3]) searchBox[3] += TileConstantsBarrel::nPhiBins; 

    for (int ieta = searchBox[0]; ieta <= searchBox[1]; ++ieta) {
      for (int iphi = searchBox[2]; iphi <= searchBox[3]; ++iphi) {
        auto& neighbours = tracksterTile[tracksterTile.globalBin(ieta, (iphi % TileConstantsBarrel::nPhiBins))];
        for (auto n : neighbours) {
          if (t_idx == n || maskedTracksters[n]) continue;
          auto& other_trackster = tracksters[n];
          auto other_trackster_span = getTracksterSpan(other_trackster);
          auto& other_barycenter = (other_trackster_span == 0) ? other_trackster.barycenter() : other_trackster.eigenvectors(0);
          float dot_product = other_barycenter.Dot(barycenter);
          float cos_alpha = dot_product / (t_norm * std::sqrt(other_barycenter.mag2()));
          if (abs(cos_alpha) > min_cos_theta_) {
            if (trackster.raw_energy() >= other_trackster.raw_energy()) {
              allNodes[t_idx].addOuterNeighbour(n);
              allNodes[n].addInnerNeighbour(t_idx);
              isRootTrackster[n] = 0;
            } else {
              allNodes[n].addOuterNeighbour(t_idx);
              allNodes[t_idx].addInnerNeighbour(n);
              isRootTrackster[t_idx] = 0;
            }
            //maskedTracksters[n] = 1; 
          }
        }
      }
    }
  }

  TICLGraph graph(allNodes);  
  auto sortedRootNodes = graph.getRootNodes();
  std::sort(sortedRootNodes.begin(), sortedRootNodes.end(), [&tracksters](const ticl::Node &n1, const ticl::Node &n2) {
    unsigned int n1Id = n1.getId();
    unsigned int n2Id = n2.getId();
    return tracksters[n1Id].raw_energy() > tracksters[n2Id].raw_energy();
  });

  int ic = 0;
  auto const &components = graph.findSubComponents(sortedRootNodes);
  linkedTracksterIdToInputTracksterId.resize(components.size());
  for (auto const &comp : components) {
    std::vector<unsigned int> linkedTracksters;
    Trackster outTrackster;
    //if (comp.size() == 1) continue;
    for (auto const &node : comp) {
      linkedTracksterIdToInputTracksterId[ic].push_back(node);
      outTrackster.mergeTracksters(input.tracksters[node]);
    }
    linkedTracksters.push_back(resultTracksters.size());
    resultTracksters.push_back(outTrackster);
    ++ic;
  }
}

   
     
