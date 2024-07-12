#ifndef RecoHGCal_TICL_ECALClustersGraph_h
#define RecoHGCal_TICL_ECALClustersGraph_h

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "RecoHGCal/TICL/interface/TICLGraph.h"

#include <vector>

namespace ticl {

  class ECALClustersGraph {


  public:
    ECALClustersGraph(const std::vector<reco::CaloCluster>& clusters,
		      float minEtSeed);
    ~ECALClustersGraph(){};

    TICLGraph buildGraph();
    std::vector<std::vector<unsigned int>> getConnectedGraphs(const  TICLGraph& graph);
    
    std::vector<unsigned int> getClustersInSeedWindow(unsigned int clIndex);
    
  private:


    std::array<float, 3> getEtaPhiWindow(float eta) {
      /// Returns the eta and phi window for a SuperCluster region
      auto aeta = std::abs(eta);
      std::array<float, 3> out; //deta_up, deta_down, dphi
      if (aeta <= 1.5){
	out[0] = (0.1/1.5)*aeta + 0.1;
	out[1] = -0.1;
      }
      else{
	out[0] = (0.1/1.5)*(aeta-1.5) + 0.2;
	out[1] = -0.1 + (-0.2/1.5)*(aeta-1.5);
      }
      out[2] = 0.7 + (-0.1/3)*aeta;
      return out;
    }

    const std::vector<reco::CaloCluster>& clusters_;
    float minEtSeed_;

    TICLLayerTilesBarrel tiles_;

  };
};



#endif
