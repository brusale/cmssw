#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/TICLLayerClusters.h"
#include "DataFormats/CaloRecHit/interface/TICLLayerClustersHostCollection.h"

#include "TTree.h"
#include "TFile.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

class LayerClustersSoAnalyzer : public edm::global::EDAnalyzer<> {
  public:
    LayerClustersSoAnalyzer(const edm::ParameterSet&);
    ~LayerClustersSoAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    void analyze(edm::StreamID, const edm::Event&, const edm::EventSetup&) const override;

  private:
    const edm::EDGetTokenT<std::vector<reco::CaloCluster>> old_clusters_token_;
    const edm::EDGetTokenT<reco::TICLLayerClustersHostCollection> new_clusters_token_;

    void beginJob() override;
    void endJob() override {};

    //void clearVariables();
    TTree* layerclusters_soa_;
   
    /*std::vector<float> old_x;
    std::vector<float> new_x; 
    std::vector<float> old_y; 
    std::vector<float> new_y; 
    std::vector<float> old_z; 
    std::vector<float> new_z; 
    std::vector<float> old_eta; 
    std::vector<float> new_eta; 
    std::vector<float> old_phi; 
    std::vector<float> new_phi;
    std::vector<float> old_energy; 
    std::vector<float> new_energy;*/
};

/*void LayerClustersSoAnalyzer::clearVariables() {
  old_x.clear();
  old_y.clear();
  old_z.clear();
  old_eta.clear();
  old_phi.clear();
  old_energy.clear();
  new_x.clear();
  new_y.clear();
  new_z.clear();
  new_eta.clear();
  new_phi.clear();
  new_energy.clear();
}*/
  
LayerClustersSoAnalyzer::LayerClustersSoAnalyzer(const edm::ParameterSet& ps) 
  : old_clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("old_clusters"))), 
    new_clusters_token_(consumes<reco::TICLLayerClustersHostCollection>(ps.getParameter<edm::InputTag>("new_clusters"))) {}

LayerClustersSoAnalyzer::~LayerClustersSoAnalyzer() { /*clearVariables();*/ }

void LayerClustersSoAnalyzer::beginJob() {
  //edm::Service<TFileService> fs;
  //layerclusters_soa_ = fs->make<TTree>("layerclusters_soa", "layerclusters_soa");

  /*layerclusters_soa_->Branch("old_x", &old_x);
  layerclusters_soa_->Branch("old_y", &old_y);
  layerclusters_soa_->Branch("old_z", &old_z);
  layerclusters_soa_->Branch("old_eta", &old_eta);
  layerclusters_soa_->Branch("old_phi", &old_phi);
  layerclusters_soa_->Branch("old_energy", &old_energy);
  layerclusters_soa_->Branch("new_x", &new_x);
  layerclusters_soa_->Branch("new_y", &new_y);
  layerclusters_soa_->Branch("new_z", &new_z);
  layerclusters_soa_->Branch("new_eta", &new_eta);
  layerclusters_soa_->Branch("new_phi", &new_phi);
  layerclusters_soa_->Branch("new_energy", &new_energy);*/
}

void LayerClustersSoAnalyzer::analyze(edm::StreamID, const edm::Event& event, const edm::EventSetup&) const {
  //clearVariables();

  edm::Handle<std::vector<reco::CaloCluster>> old_clusters_h;
  event.getByToken(old_clusters_token_, old_clusters_h);
  auto& old_clusters = *old_clusters_h;

  //edm::Handle<reco::LayerClustersHostCollection> new_clusters_h;
  //event.getByToken(new_clusters_token_, new_clusters_h);
  auto& new_clusters = event.get(new_clusters_token_);
  //const auto& new_clusters = *new_clusters_h;
  int nclusters = old_clusters.size();
  /*for (const auto& cluster : old_clusters) {
    old_x.push_back(cluster.x());
    old_y.push_back(cluster.y());
    old_z.push_back(cluster.z());
    old_eta.push_back(cluster.eta());
    old_phi.push_back(cluster.phi());
    old_energy.push_back(cluster.energy());
  }*/


  const auto& lc_view = new_clusters.view();
  /*for (int i = 0;i < nclusters; ++i) {
    auto cluster = lc_view[i];
    new_x.push_back(cluster.x());
    new_y.push_back(cluster.y());
    new_z.push_back(cluster.z());
    new_eta.push_back(cluster.eta());
    new_phi.push_back(cluster.phi());
    new_energy.push_back(cluster.energy());
  }*/
  for (int i = 0; i < nclusters; ++i) {
    std::cout << "============\n";
    std::cout << "CaloCluster\n  energy: " << old_clusters[i].energy() << " x: " << old_clusters[i].x() << " y: " << old_clusters[i].y() << " y: " << old_clusters[i].z() << "\n";
    std::cout << "LayerCluster\n  energy: " << lc_view[i].energy() << " x: " << lc_view[i].x() << " y: " << lc_view[i].y() << " y: " << lc_view[i].z() << "\n";
  }
}

void LayerClustersSoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("old_clusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("new_clusters", edm::InputTag("caloclusters2layerclusters"));
  descriptions.add("layerClustersSoAnalyzer", desc);
}

DEFINE_FWK_MODULE(LayerClustersSoAnalyzer);
