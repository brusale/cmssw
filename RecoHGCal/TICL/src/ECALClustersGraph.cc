#include "RecoHGCal/TICL/interface/ECALClustersGraph.h"

namespace ticl {

  ECALClustersGraph::ECALClustersGraph(const std::vector<reco::CaloCluster>& clusters,
				       float minEtSeed): clusters_(clusters), minEtSeed_(minEtSeed){
    //Loop over clusters and fill the tiles
    for (size_t i = 0; i < clusters_.size(); ++i) {
      auto& cluster = clusters_[i];
      tiles_.fill(0, static_cast<float>(cluster.eta()),
		  static_cast<float>(cluster.phi()),
		  i);
    } 
  }

  TICLGraph ECALClustersGraph::buildGraph() {
    // Now looking for neighboard clusters with lower energy
    // using the tiles
    std::vector<ticl::Node> nodes;
    std::vector<int> isRootNode;
    for (size_t i = 0; i < clusters_.size(); ++i){
      auto neighbours = getClustersInSeedWindow(i);
      ticl::Node node(i, false/*isTrackster*/);

      if (neighbours.size()>0){
	for (const auto& k_cl : neighbours)
	  node.addInnerNeighbour(k_cl);
	isRootNode.push_back(1);
      }
      else{
	isRootNode.push_back(0);
      }
      nodes.push_back(node);
    }
    return TICLGraph(nodes, isRootNode);
  }
  
  std::vector<std::vector<unsigned int>> ECALClustersGraph::getConnectedGraphs(const TICLGraph& graph) {
    return graph.getConnectedComponents(true /*outgoing*/, true /*ingoing*/);
  }
  
  std::vector<unsigned int> ECALClustersGraph::getClustersInSeedWindow(unsigned int clIndex) {
    std::vector<unsigned int> out;
    auto& cluster = clusters_[clIndex];
    auto seed_eta = cluster.eta();
    auto seed_et = cluster.energy()/std::cosh(seed_eta);
    if (seed_et < minEtSeed_) return out;
    std::array<float, 3> window = getEtaPhiWindow(seed_eta);
    // Now we use the tiles to get the clusters in the window
    float etaMin, etaMax;
    float phiMin = cluster.phi() - window[2];
    float phiMax = cluster.phi() + window[2];
    if (seed_eta>0.) {
      etaMax = seed_eta + window[0];
      etaMin = seed_eta + window[1];
    }
    else {
      etaMax = seed_eta - window[0];
      etaMin = seed_eta - window[1];
    }
    auto tileIndices = tiles_[0].searchBoxEtaPhi(etaMin, etaMax, phiMin, phiMax);
    for (int i = tileIndices[0]; i <= tileIndices[1]; ++i) {
      for (int j = tileIndices[2]; j <= tileIndices[3]; ++j) {
      	auto& tile = tiles_[0][tiles_[0].globalBin(i, j)];
	for (const auto& k_cl : tile){
	  // Now I check the distance
	  auto& neighbour = clusters_[k_cl];

	  auto neig_eta = neighbour.eta();
	  auto neig_phi = neighbour.phi();

	  // Do not link to higher energy clusters
	  if (neighbour.energy()/std::cosh(neighbour.eta()) > seed_et) continue;
          // Geometric check
	  if (neig_eta > etaMax || neig_eta < etaMin ||
	      neig_phi > phiMax || neig_phi < phiMin) continue;
	  
	  out.push_back(k_cl);
	}	
      }
    }
    return out;
  }  
};
