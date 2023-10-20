#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/LayerClusters.h"
#include "DataFormats/CaloRecHit/interface/LayerClustersHostCollection.h"

#include "TTree.h"
#include "TFile.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

class LayerClustersSoAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
  public:
    explicit LayerClustersSoAnalyzer(const edm::ParameterSet&);
    ~LayerClustersSoAnalyzer();

  private:
    const edm::EDGetTokenT<std::vector<reco::CaloCluster>> old_clusters_token_;
    const edm::EDGetTokenT<reco::LayerClustersHostCollection> new_clusters_token_;

    void beginJob() override;
    void beginRun(const edm::Run&, const edm::EventSetup&) override {};
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endRun(edm::Run const& iEvent, edm::EventSetup const&) override {};
    void endJob() override {};

    void clearVariables();
    TTree* layerclusters_soa_;
   
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    std::vector<float> old_x, new_x, old_y, new_y, old_z, new_z, old_eta, new_eta, old_phi, new_phi;
    std::vector<float> old_energy, new_energy;
};

void LayerClustersSoAnalyzer::clearVariables() {
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
}

LayerClustersSoAnalyzer::LayerClustersSoAnalyzer(const edm::ParameterSet& ps) 
  : old_clusters_token_(consumes<reco::CaloCluster>(ps.getParameter<edm::InputTag>("old_clusters"), 
    new_clusters_token_(consumes<reco::LayerClustersHostCollection>(ps.getParameter<edm::InputTag>("new_clusters"))) {}

LayerClustersSoAnalyzer::~LayerClustersSoAnalyzer() { clearVariables(); }

void LayerClustersSoAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  layerclusters_soa_ = fs->make<TTree>("layerclusters_soa", "layerclusters_soa");

  layerclusters_soa_->Branch("old_x", &old_x);
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
  layerclusters_soa_->Branch("new_energy", &new_energy);
}

void LayerClustersSoAnalyzer::analyze(const edm::Event& event, const edm::EventSetup&) {
  clearVariables();

  edm::Handle<std::vector<reco::CaloCluster>> old_clusters_h;
  event.getByToken(old_clusters_token_, old_clusters_h);
  const auto& old_clusters = *old_clusters_h;

  edm::Handle<reco::LayerClustersHostCollection> new_clusters_h;
  event.getByToken(new_clusters_token_, new_clusters_h);
  const auto& new_clusters = *new_clusters_h;

  for (const auto& cluster : old_clusters) {
    old_x.push_back(cluster.x());
    old_y.push_back(cluster.y());
    old_z.push_back(cluster.z());
    old_eta.push_back(cluster.eta());
    old_phi.push_back(cluster.phi());
    old_energy.push_back(cluster.energy());
  }


  const auto& lc_view = new_clusters::view;
  for (int i = 0;i < lc_view.size(); ++i) {
    auto cluster = lc_view[i];
    new_x.push_back(cluster.x());
    new_y.push_back(cluster.y());
    new_z.push_back(cluster.z());
    new_eta.push_back(cluster.eta());
    new_phi.push_back(cluster.phi());
    new_energy.push_back(cluster.energy());
  }
}

void LayerClustersSoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("old_clusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("new_clusters", edm::InputTag("caloclusters2layerclusters"));
  descriptions.add("layerClustersSoAnalyzer", desc);
}

DEFINE_FWK_MODULE(LayerClustersSoAnalyzer);
