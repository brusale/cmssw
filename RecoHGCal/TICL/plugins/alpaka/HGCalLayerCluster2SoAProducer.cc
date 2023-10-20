#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

//new LayerClusters
#include "DataFormats/CaloRecHit/interface/LayerClusters.h"
#include "DataFormats/CaloRecHit/interface/CaloClustersHostCollection.h"
#include "DataFormats/CaloRecHit/interface/LayerClustersHostCollection.h"
#include "DataFormats/CaloRecHit/interface/alpaka/LayerClustersDeviceCollection.h"

#include "HeterogeneousCore/AlpakaInterface/interface/OneToManyAssoc.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/Math/interface/Point3D.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalLayerCluster2SoAProducer : public stream::EDProducer<> {
    public: 
      HGCalLayerCluster2SoAProducer(const edm::ParameterSet&);
      ~HGCalLayerCluster2SoAProducer() override {}
      
      struct hitAndFraction {
	float fraction;
	uint32_t hit;
      }

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      void produce(device::Event& evt, device::EventSetup const&) override;
      static void fillErrorEntry(reco::LayerClustersHostCollection::View::element& soa_element, reco::CaloCluster& cluster, const CaloGeometry* geom, hgcal::RecHitTools rhtools);
      static void fillSoA(reco::LayerClustersHostCollection::View::element soa_element, reco::CaloCluster cluster, const CaloGeometry* geom, hgcal::RecHitTools rhtools);
      static void fillOneToManyAssociation(reco::LayerClustersHostCollection::View::element& soa_element, reco::CaloCluster cluster);

    private:
      hgcal::RecHitTools rhtools_; //let's put this here for now
      const CaloGeometry* geom_;
     
      edm::EDGetTokenT<reco::CaloClusterCollection> input_token_;
      //edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> time_token_;
      const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
      const device::EDPutToken<LayerClustersDeviceCollection> device_token_;

      void beginRun(const edm::Run&, const edm::EventSetup&) override;
  };

  void HGCalLayerCluster2SoAProducer::beginRun(const edm::Run&, edm::EventSetup const& es) {
    geom_ = &(es.getData(geometry_token_));
    rhtools_.setGeometry(*geom_);
  }

  HGCalLayerCluster2SoAProducer::HGCalLayerCluster2SoAProducer(const edm::ParameterSet& config)
    : input_token_(consumes(config.getParameter<edm::InputTag>("input_clusters"))),
      //time_token_(consumes(config.getParameter<edm::InputTag>("time_map"))),
      geometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()),
      device_token_(produces()) {}

  void HGCalLayerCluster2SoAProducer::fillErrorEntry(reco::LayerClustersHostCollection::View::element& soa_element, reco::CaloCluster& cluster, const CaloGeometry* geom, hgcal::RecHitTools rhtools) {
    float sum_x = 0.;
    float sum_y = 0.;
    float sum_sqr_x = 0.;
    float sum_sqr_y = 0.;
    float ref_x = cluster.x();
    float ref_y = cluster.y();
    auto& hits_and_fractions = cluster.hitsAndFractions();
    float invClsize = 1. / hits_and_fractions.size();
    for (auto const& haf : hits_and_fractions) {
      auto position = geom->getPosition(haf.first);
      auto diff_x = position.x() - ref_x;
      auto diff_y = position.y() - ref_y;
      sum_x += diff_x;
      sum_sqr_x += diff_x * diff_x;
      sum_y += diff_y;
      sum_sqr_y = diff_y * diff_y;
    }	
    float radius_x = sqrt((sum_sqr_x - (sum_x * sum_x) * invClsize) * invClsize);
    float radius_y = sqrt((sum_sqr_y - (sum_y * sum_y) * invClsize) * invClsize);

    auto detId = hits_and_fractions[0].first;
    if (invClsize == 1.) {
      if (rhtools.isSilicon(detId)) {
        //Silicon case
        radius_x = radius_y = rhtools.getRadiusToSide(detId);
      } else {
	  auto const& position = rhtools.getPosition(detId);
	  auto const& eta_phi_window = rhtools.getScintDEtaDPhi(detId);
	  radius_x = radius_y = position.perp() * eta_phi_window.second;
      }
    }
    soa_element.radius() = radius_x + radius_y;
  }

void HGCalLayerCluster2SoAProducer::fillOneToManyAssociation(reco::LayerClustersHostCollection::View::element& soa_element, reco::CaloCluster cluster) {
  auto hits_and_fractions = cluster.hitsAndFractions();
  std::vector<hitAndFraction> haf;
  int nh = hits_and_fractions.size();
  haf.reserve(nh);
  
  for (int i = 0; i < nh; ++i) {
    haf[i].hit = hits_and_fractions[i].first;
    haf[i].fraction = hits_and_fractions[i].second;
  }

  soa_element.hitsAndFractions().bulkFill(apc, haf, nh + 1); 
  
}

  void HGCalLayerCluster2SoAProducer::fillSoA(reco::LayerClustersHostCollection::View::element soa_element, reco::CaloCluster cluster, const CaloGeometry* geom, hgcal::RecHitTools rhtools) {
    soa_element.x() = cluster.x();
    soa_element.y() = cluster.y();
    soa_element.z() = cluster.z();
    soa_element.eta() = cluster.eta();
    soa_element.phi() = cluster.phi();
    soa_element.energy() = cluster.energy();
    fillErrorEntry(soa_element, cluster, geom, rhtools);
    fillOneToManyAssociation(soa_element, cluster);
  }

  void HGCalLayerCluster2SoAProducer::produce(device::Event& evt, device::EventSetup const& es) {
    //edm::Handle<reco::CaloClusterCollection> input_h;
    //evt.getByToken(input_token_, input_h);
    auto& input_h = evt.get(input_token_);
    //geom_ = &(es.getData(geometry_token_));
    //auto& time = evt.get(time_token_);
    //reco::CaloClusterHostCollection calo_clusters{(*input_h).size(), event.queue()};
    //auto& view = calo_clusters.view();
    int32_t num_clusters = input_h.size();
										       
    reco::LayerClustersHostCollection h_layer_clusters{num_clusters, evt.queue()};
    auto& h_layer_clusters_view = h_layer_clusters.view();
										       
    for (int i = 0; i < num_clusters; ++i) {
      fillSoA(h_layer_clusters_view[i], (input_h)[i], geom_, rhtools_);
    }
										       
    LayerClustersDeviceCollection d_layer_clusters{num_clusters, evt.queue()};
    alpaka::memcpy(evt.queue(), d_layer_clusters.buffer(), h_layer_clusters.buffer());
    //synchronize
    //alpaka::wait(event.queue());
										       
    evt.emplace(device_token_, std::move(d_layer_clusters));
    //event.put(std::move(layer_clusters));
  }

  void HGCalLayerCluster2SoAProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("input_clusters", edm::InputTag("hgcalMergeLayerClusters"));
    desc.add<edm::InputTag>("time_map", edm::InputTag("hgcalMergeLayerClusters", "timeLayerClusters"));

    descriptions.addWithDefaultLabel(desc);
  }
} // ALPAKA_ACC_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HGCalLayerCluster2SoAProducer);
