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
    min_cos_theta_(conf.getParameter<double>("min_cos_theta")) {}

void TracksterLinkingBarrel::initialize(const HGCalDDDConstants *hgcons,
                                        const hgcal::RecHitTools rhtools,
                                        const edm::ESHandle<MagneticField> bfield,
                                        const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;
  bfield_ = bfield;
  propagator_ = propH;

  float etaWidth = (TileConstantsBarrel::maxEta - TileConstantsBarrel::minEta) / TileConstantsBarrel::nEtaBins;
  float phiWidth = 2*std::numbers::pi_v<float> / TileConstants::nPhiBins;
 
  for (int i = 0; i < TileConstantsBarrel::nEtaBins; ++i) 
    eta_windows_[i] = TileConstantsBarrel::minEta + i*etaWidth;
  for (int i = 0; i < TileConstantsBarrel::nPhiBins; ++i)
    phi_windows_[i] = -std::numbers::pi_v<float> + i*phiWidth;
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

  for (auto const &t_idx : sortedTracksters) {
    //if (maskedTracksters[t_idx]) continue;
    linkedResultTracksters.reserve(t_idx);
    auto const& trackster = tracksters[t_idx];
    auto const& barycenter = (trackster.vertices().size() == 1) ? trackster.barycenter() : trackster.eigenvectors(0);
    //auto const& barycenter = trackster.barycenter();
    float t_norm = std::sqrt(barycenter.mag2());

    std::cout << "=== TracksterLinkingBarrel:" << __LINE__ << " ===" << std::endl;
    std::cout << "Opening window for trackster " << t_idx << std::endl;
    std::cout << "barycenter_eta: " << barycenter.eta() << std::endl;
    std::cout << "barycenter_phi: " << barycenter.phi() << std::endl;
    std::cout << "n_vertices: " << trackster.vertices().size() << std::endl;
    float const eta_window = eta_windows_[tracksterTile.etaBin(barycenter.eta())];
    float const phi_window = phi_windows_[tracksterTile.etaBin(barycenter.phi())];
    std::cout << "eta_window: " << eta_windows_[tracksterTile.etaBin(barycenter.eta())] << std::endl;
    std::cout << "phi_window: " << eta_windows_[tracksterTile.etaBin(barycenter.phi())] << std::endl;
    
    float eta_min = std::max(abs(barycenter.eta()) - eta_window, TileConstantsBarrel::minEta);
    float eta_max = std::min(abs(barycenter.eta()) + eta_window, TileConstantsBarrel::maxEta);
    float phi_min = std::max(abs(barycenter.phi()) - phi_window, -std::numbers::pi_v<float>);
    float phi_max = std::min(abs(barycenter.phi()) + phi_window, std::numbers::pi_v<float>);
 
    std::array<int, 4> searchBox = tracksterTile.searchBoxEtaPhi(eta_min, eta_max, phi_min, phi_max);
    if (searchBox[2] > searchBox[3]) searchBox[3] += TileConstantsBarrel::nPhiBins; 

    for (int ieta = searchBox[0]; ieta <= searchBox[1]; ++ieta) {
      std::cout << "ieta bin: " << ieta << std::endl;
      for (int iphi = searchBox[2]; iphi <= searchBox[3]; ++iphi) {
        std::cout << "iphi bin: " << iphi << std::endl;
        auto& neighbours = tracksterTile[tracksterTile.globalBin(ieta, (iphi % TileConstantsBarrel::nPhiBins))];
        for (auto n : neighbours) {
          if (t_idx == n || maskedTracksters[n]) continue;
          auto& other_trackster = tracksters[n];
          float dot_product = other_trackster.eigenvectors(0).Dot(barycenter);
          //float dot_product = other_trackster.barycenter().Dot(barycenter);
          float cos_alpha = dot_product / (t_norm * std::sqrt(other_trackster.eigenvectors(0).mag2()));
          //float cos_alpha = dot_product / (t_norm * std::sqrt(other_trackster.barycenter().mag2()));
          if (abs(cos_alpha) > min_cos_theta_) {
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            std::cout << "t barycenter: (" << trackster.barycenter().eta() << "," << trackster.barycenter().phi() << ")" << std::endl;
            std::cout << "other_t barycenter: (" << other_trackster.barycenter().eta() << "," << other_trackster.barycenter().phi() << ")" << std::endl;
            if (trackster.raw_energy() >= other_trackster.raw_energy()) {
              allNodes[t_idx].addOuterNeighbour(n);
              allNodes[n].addInnerNeighbour(t_idx);
              isRootTrackster[n] = 0;
              std::cout << __FILE__ << " " << __LINE__ << std::endl;
              std::cout << "Linking trackster " << n << " to trackster " << t_idx << std::endl;
            } else {
              allNodes[n].addOuterNeighbour(t_idx);
              allNodes[t_idx].addInnerNeighbour(n);
              isRootTrackster[t_idx] = 0;
              std::cout << __FILE__ << " " << __LINE__ << std::endl;
              std::cout << "Linking trackster " << t_idx << " to trackster " << n << std::endl;
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

   
     
