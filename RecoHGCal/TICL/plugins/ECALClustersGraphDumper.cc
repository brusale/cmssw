// -*- C++ -*-
//
// Package:    RecoHGCal/TICL
// Class:      ECALClustersGraphDumper
//
/**\class ECALClustersGraphDumper ECALClustersGraphDumper.cc RecoHGCal/TICL/plugins/ECALClustersGraphDumper.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Davide Valsecchi
//         Created:  Wed, 10 Jul 2024 09:58:08 GMT
//
//

// system include files
#include <memory>
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include "RecoHGCal/TICL/interface/ECALClustersGraph.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class ECALClustersGraphDumper : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ECALClustersGraphDumper(const edm::ParameterSet&);
  ~ECALClustersGraphDumper() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::CaloClusterCollection> caloClustersToken_;
   // Output tree
  TTree* tree_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ECALClustersGraphDumper::ECALClustersGraphDumper(const edm::ParameterSet& iConfig)
  : caloClustersToken_(consumes<reco::CaloClusterCollection>(iConfig.getUntrackedParameter<edm::InputTag>("layer_clusters"))) {
}

ECALClustersGraphDumper::~ECALClustersGraphDumper() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void ECALClustersGraphDumper::analyze(const edm::Event& event, const edm::EventSetup& iSetup) {
  using namespace edm;
  edm::Service<TFileService> fs;

  tree_ = fs->make<TTree>("graphs", "TICL tracksters");

  //get all the tracksters
  edm::Handle<reco::CaloClusterCollection> caloClusterHandle;
  event.getByToken(caloClustersToken_, caloClusterHandle);
  const auto& clusters = *caloClusterHandle;

  std::cout << "Running ECAL GRAPh tools" << std::endl;
  ticl::ECALClustersGraph ecalGraphTool(clusters, 5.0);
  auto graph = ecalGraphTool.buildGraph();

  auto connected = graph.getConnectedComponents(true, true);
  //Print the groups
  for (size_t i = 0; i < connected.size(); ++i) {
    std::cout << "Group " << i << " has "
	      << connected[i].size() << " clusters: " ;
    for (const auto& cl : connected[i]) {
      std::cout << cl << " ";
    }
    std::cout << std::endl;
  }
  
  
  

}

// ------------ method called once each job just before starting event loop  ------------
void ECALClustersGraphDumper::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void ECALClustersGraphDumper::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ECALClustersGraphDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  desc.addUntracked<edm::InputTag>("layer_clusters", edm::InputTag("ecalClusters"));
  descriptions.addDefault(desc);
  
  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ECALClustersGraphDumper);
