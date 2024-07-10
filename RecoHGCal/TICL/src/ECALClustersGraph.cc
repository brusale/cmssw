 #include "RecoHGCal/TICL/interface/ECALClustersGraph.h"

#include <iostream>

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
      std::cout << "Cluster " << i << " has " << neighbours.size() << " neighbours" << std::endl;
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
    std::cout << "Looking for neighbours of cluster " << clIndex << std::endl;
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
      etaMin = seed_eta - window[0];
      etaMax = seed_eta - window[1];
    }
    std::cout << "Eta-phi window: " << etaMin << " " << etaMax << " " << phiMin << " " << phiMax << std::endl;
    std::cout << "Seed et: " << seed_et << std::endl;
    
    auto tileIndices = tiles_[0].searchBoxEtaPhi(etaMin, etaMax, phiMin, phiMax);
    std::cout << "Looking in tiles: ";
    for (const auto& i : tileIndices) std::cout << i << ", ";
    std::cout << std::endl;
    
    for (int i = tileIndices[0]; i <= tileIndices[1]; ++i) {
      for (int j = tileIndices[2]; j <= tileIndices[3]; ++j) {
      	auto& tile = tiles_[0][tiles_[0].globalBin(i, j)];
	for (const auto& k_cl : tile){
	  if (k_cl == clIndex) continue; // Exclude the seed
	  
	  // Now I check the distance
	  auto& neighbour = clusters_[k_cl];

	  auto neig_eta = neighbour.eta();
	  auto neig_phi = neighbour.phi();

	  std::cout << "Found neighbour " << k_cl 
		    << " - Neighbour eta-phi-et: " << neig_eta
		    << " " << neig_phi
		    << " " << neighbour.energy()/std::cosh(neighbour.eta())
		    << std::endl;

	  
	  // Do not link to higher energy clusters
	  if (neighbour.energy()/std::cosh(neighbour.eta()) > seed_et) continue;
          // Geometric check
	
	  if (neig_eta > etaMax || neig_eta < etaMin ||
	      neig_phi > phiMax || neig_phi < phiMin) continue;

	  std::cout << "Connected!" << std::endl;
	  out.push_back(k_cl);
	}	
      }
    }
    return out;
  }  
};
