#include "RecoHGCal/TICL/interface/TICLInterpretationAlgoBase.h"
#include "RecoHGCal/TICL/plugins/GeneralInterpretationAlgo.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

using namespace ticl;

using Vector = ticl::Trackster::Vector;

GeneralInterpretationAlgo::~GeneralInterpretationAlgo() {}

GeneralInterpretationAlgo::GeneralInterpretationAlgo(const edm::ParameterSet &conf, edm::ConsumesCollector cc)
    : TICLInterpretationAlgoBase(conf, cc),
      del_tk_ts_layer1_(conf.getParameter<double>("delta_tk_ts_layer1")),
      del_tk_ts_int_(conf.getParameter<double>("delta_tk_ts_interface")),
      timing_quality_threshold_(conf.getParameter<double>("timing_quality_threshold")) {}

void GeneralInterpretationAlgo::initialize(const HGCalDDDConstants *hgcons,
                                           const hgcal::RecHitTools rhtools,
                                           const edm::ESHandle<MagneticField> bfieldH,
                                           const edm::ESHandle<Propagator> propH,
                                           const std::string detector) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;
  detector_ = detector;
  buildLayers();

  bfield_ = bfieldH;
  propagator_ = propH;
}

void GeneralInterpretationAlgo::buildLayers() {
  // build disks at HGCal front & EM-Had interface for track propagation

  if (detector_ == "HGCAL") {
    float zVal = hgcons_->waferZ(1, true);
    std::pair<float, float> rMinMax = hgcons_->rangeR(zVal, true);
  
    float zVal_interface = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z();
    std::pair<float, float> rMinMax_interface = hgcons_->rangeR(zVal_interface, true);
  
    for (int iSide = 0; iSide < 2; ++iSide) {
      float zSide = (iSide == 0) ? (-1. * zVal) : zVal;
      firstDisk_[iSide] =
          std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
                                                Disk::RotationType(),
                                                SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5, zSide + 0.5))
                                        .get());
  
      zSide = (iSide == 0) ? (-1. * zVal_interface) : zVal_interface;
      interfaceDisk_[iSide] = std::make_unique<GeomDet>(
          Disk::build(Disk::PositionType(0, 0, zSide),
                      Disk::RotationType(),
                      SimpleDiskBounds(rMinMax_interface.first, rMinMax_interface.second, zSide - 0.5, zSide + 0.5))
              .get());
    }
  } else {
    GlobalPoint ecalSurface = rhtools_.getPositionLayer(0, false, true);
    ecalCylinder_ = Cylinder::build(ecalSurface.mag(),
                                           Surface::PositionType(0., 0., 0.),
                                           Surface::RotationType());
    
    GlobalPoint hcalSurface = rhtools_.getPositionLayer(1, false, true);
    hcalCylinder_ = Cylinder::build(hcalSurface.mag(),
                                    Surface::PositionType(0., 0., 0.),
                                    Surface::RotationType());
  }
}

Vector GeneralInterpretationAlgo::propagateTrackster(const Trackster &t,
                                                     const unsigned idx,
                                                     float zVal,
                                                     std::array<TICLLayerTile, 2> &tracksterTiles,
                                                     TICLLayerTileBarrel &tracksterTilesBarrel) {
  // needs only the positive Z co-ordinate of the surface to propagate to
  // the correct sign is calculated inside according to the barycenter of trackster
  Vector const &baryc = t.barycenter();
  Vector directnv = t.eigenvectors(0);

  // barycenter as direction for tracksters w/ poor PCA
  // propagation still done to get the cartesian coords
  // which are anyway converted to eta, phi in linking
  // -> can be simplified later

  //FP: disable PCA propagation for the moment and fallback to barycenter position
  // if (t.eigenvalues()[0] / t.eigenvalues()[1] < 20)
  directnv = baryc.unit();
  if (detector_ == "HGCAL") {
    zVal *= (baryc.Z() > 0) ? 1 : -1;
    float par = (zVal - baryc.Z()) / directnv.Z();
    float xOnSurface = par * directnv.X() + baryc.X();
    float yOnSurface = par * directnv.Y() + baryc.Y();
    Vector tPoint(xOnSurface, yOnSurface, zVal);
    if (tPoint.Eta() > 0) {
      tracksterTiles[1].fill(tPoint.Eta(), tPoint.Phi(), idx);
    } else if (tPoint.Eta() < 0) {
      tracksterTiles[0].fill(tPoint.Eta(), tPoint.Phi(), idx);
    }

    return tPoint;
  } else {
    zVal *= (baryc.Eta() > 0) ? 1 : -1;
    //float 
    float par = (zVal - baryc.R()) / directnv.R();
    float xOnSurface = par * directnv.X() + baryc.X();
    float yOnSurface = par * directnv.Y() + baryc.Y();
    Vector tPoint(xOnSurface, yOnSurface, zVal);
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    std::cout << "xOnSurface: " << xOnSurface << std::endl;
    std::cout << "yOnSurface: " << yOnSurface << std::endl;
    std::cout << "zVal: " << zVal << std::endl;
    std::cout << "tPoint.Eta(): " << tPoint.Eta() << std::endl;
    std::cout << "tPoint.Phi(): " << tPoint.Phi() << std::endl;
    std::cout << "idx: " << idx << std::endl;
    tracksterTilesBarrel.fill(tPoint.Eta(), tPoint.Phi(), idx);
    return tPoint;
  }
}

void GeneralInterpretationAlgo::findTrackstersInWindow(const edm::MultiSpan<Trackster> &tracksters,
                                                       const std::vector<std::pair<Vector, unsigned>> &seedingCollection,
                                                       const std::array<TICLLayerTile, 2> &tracksterTiles,
                                                       const TICLLayerTileBarrel &tracksterTilesBarrel,
                                                       const std::vector<Vector> &tracksterPropPoints,
                                                       const float delta,
                                                       unsigned trackstersSize,
                                                       std::vector<std::vector<unsigned>> &resultCollection,
                                                       bool useMask = false) {
  // Finds tracksters in tracksterTiles within an eta-phi window
  // (given by delta) of the objects (track/trackster) in the seedingCollection.
  // Element i in resultCollection is the vector of trackster
  // indices found close to the i-th object in the seedingCollection.
  // If specified, Tracksters are masked once found as close to an object.
  std::vector<int> mask(trackstersSize, 0);
  const float delta2 = delta * delta;

  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  std::cout << "seedingCollection.size(): " << seedingCollection.size() << std::endl;
  std::cout << "detector: " << detector_ << std::endl;
  for (auto &i : seedingCollection) {
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    float seed_eta = i.first.Eta();
    float seed_phi = i.first.Phi();
    unsigned seedId = i.second;
    std::array<int, 4> search_box;
    std::vector<unsigned> in_delta;
    std::vector<float> energies;
    // std::vector<float> distances2;
    if (detector_ == "HGCAL") {
      auto sideZ = seed_eta > 0;  //forward or backward region
      const TICLLayerTile &tile = tracksterTiles[sideZ];
      float eta_min = std::max(std::fabs(seed_eta) - delta, (float)TileConstants::minEta);
      float eta_max = std::min(std::fabs(seed_eta) + delta, (float)TileConstants::maxEta);
      search_box = tile.searchBoxEtaPhi(eta_min, eta_max, seed_phi - delta, seed_phi + delta);
    
      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          const auto &in_tile = tile[tile.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
          for (const unsigned &t_i : in_tile) {
            // calculate actual distances of tracksters to the seed for a more accurate cut
            auto sep2 = (tracksterPropPoints[t_i].Eta() - seed_eta) * (tracksterPropPoints[t_i].Eta() - seed_eta) +
                        (tracksterPropPoints[t_i].Phi() - seed_phi) * (tracksterPropPoints[t_i].Phi() - seed_phi);
            if (sep2 < delta2) {
              in_delta.push_back(t_i);
              // distances2.push_back(sep2);
              energies.push_back(tracksters[t_i].raw_energy());
            }
          }
        }
      }
    } else {
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      float eta_min = std::max(std::fabs(seed_eta) - delta, (float)TileConstantsBarrel::minEta);
      float eta_max = std::min(std::fabs(seed_eta) + delta, (float)TileConstantsBarrel::maxEta);
      std::cout << "eta_min: " << eta_min << std::endl;
      std::cout << "eta_max: " << eta_max << std::endl;
      search_box = tracksterTilesBarrel.searchBoxEtaPhi(eta_min, eta_max, seed_phi - delta, seed_phi + delta);
      // get range of bins touched by delta
      
      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          const auto &in_tile = tracksterTilesBarrel[tracksterTilesBarrel.globalBin(eta_i, (phi_i % TileConstantsBarrel::nPhiBins))];
          std::cout << "seed_phi: " << seed_phi << std::endl;
          std::cout << "delta: " << delta <<  std::endl;
          std::cout << "eta_i: " << eta_i << std::endl;
          std::cout << "phi_i: " << phi_i << std::endl;
          std::cout << "globalBin: " << tracksterTilesBarrel.globalBin(eta_i, (phi_i % TileConstantsBarrel::nPhiBins)) << std::endl;
          std::cout << "t_i: " << in_tile.size() << std::endl;
          for (const unsigned &t_i : in_tile) {
            // calculate actual distances of tracksters to the seed for a more accurate cut
            auto sep2 = (tracksterPropPoints[t_i].Eta() - seed_eta) * (tracksterPropPoints[t_i].Eta() - seed_eta) +
                        (tracksterPropPoints[t_i].Phi() - seed_phi) * (tracksterPropPoints[t_i].Phi() - seed_phi);
            if (sep2 < delta2) {
              in_delta.push_back(t_i);
              // distances2.push_back(sep2);
              energies.push_back(tracksters[t_i].raw_energy());
            }
          }
        }
      }
    }

    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    // sort tracksters found in ascending order of their distances from the seed
    std::vector<unsigned> indices(in_delta.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](unsigned i, unsigned j) { return energies[i] > energies[j]; });

    // push back sorted tracksters in the result collection
    for (const unsigned &index : indices) {
      const auto &t_i = in_delta[index];
      if (!mask[t_i]) {
        resultCollection[seedId].push_back(t_i);
        if (useMask)
          mask[t_i] = 1;
      }
    }
  }  // seeding collection loop
}

bool GeneralInterpretationAlgo::timeAndEnergyCompatible(float &total_raw_energy,
                                                        const reco::Track &track,
                                                        const Trackster &trackster,
                                                        const float &tkT,
                                                        const float &tkTErr,
                                                        const float &tkQual,
                                                        const float &tkBeta,
                                                        const GlobalPoint &tkMtdPos,
                                                        bool useMTDTiming) {
  float threshold = std::min(0.2 * trackster.raw_energy(), 10.0);
  bool energyCompatible = (total_raw_energy + trackster.raw_energy() < track.p() + threshold);

  if (!useMTDTiming)
    return energyCompatible;

  // compatible if trackster time is within 3sigma of
  // track time; compatible if either: no time assigned
  // to trackster or track

  float tsT = trackster.time();
  float tsTErr = trackster.timeError();

  bool timeCompatible = false;
  if (tsT == -99. or tkTErr == -1 or tkQual < timing_quality_threshold_)
    timeCompatible = true;
  else {
    const auto &barycenter = trackster.barycenter();

    const auto deltaSoverV = std::sqrt((barycenter.x() - tkMtdPos.x()) * (barycenter.x() - tkMtdPos.x()) +
                                       (barycenter.y() - tkMtdPos.y()) * (barycenter.y() - tkMtdPos.y()) +
                                       (barycenter.z() - tkMtdPos.z()) * (barycenter.z() - tkMtdPos.z())) /
                             (tkBeta * 29.9792458);

    const auto deltaT = tsT - tkT;

    //  timeCompatible = (std::abs(deltaSoverV - deltaT) < maxDeltaT_ * sqrt(tsTErr * tsTErr + tkTErr * tkTErr));
    // use sqrt(2) * error on the track for the total error, because the time of the trackster is too small
    timeCompatible = std::abs(deltaSoverV - deltaT) < maxDeltaT_ * std::sqrt(tsTErr * tsTErr + tkTErr * tkTErr);
  }

  if (TICLInterpretationAlgoBase::algo_verbosity_ > VerbosityLevel::Advanced) {
    if (!(energyCompatible))
      LogDebug("GeneralInterpretationAlgo")
          << "energy incompatible : track p " << track.p() << " trackster energy " << trackster.raw_energy() << "\n"
          << "                      total_raw_energy " << total_raw_energy << " greater than track p + threshold "
          << track.p() + threshold << "\n";
    if (!(timeCompatible))
      LogDebug("GeneralInterpretationAlgo") << "time incompatible : track time " << tkT << " +/- " << tkTErr
                                            << " trackster time " << tsT << " +/- " << tsTErr << "\n";
  }

  return energyCompatible && timeCompatible;
}

void GeneralInterpretationAlgo::makeCandidates(const Inputs &input,
                                               edm::Handle<MtdHostCollection> inputTiming_h,
                                               std::vector<Trackster> &resultTracksters,
                                               std::vector<int> &resultCandidate) {
  bool useMTDTiming = inputTiming_h.isValid();
  const auto tkH = input.tracksHandle;
  const auto maskTracks = input.maskedTracks;
  const auto &tracks = *tkH;
  const auto &tracksters = input.tracksters;

  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);

  // propagated point collections
  // elements in the propagated points collecions are used
  // to look for potential linkages in the appropriate tiles
  std::vector<std::pair<Vector, unsigned>> trackPColl;     // propagated track points and index of track in collection
  std::vector<std::pair<Vector, unsigned>> tkPropIntColl;  // tracks propagated to lastLayerEE

  trackPColl.reserve(tracks.size());
  tkPropIntColl.reserve(tracks.size());

  std::array<TICLLayerTile, 2> tracksterPropTiles = {};  // all Tracksters, propagated to layer 1
  std::array<TICLLayerTile, 2> tsPropIntTiles = {};      // all Tracksters, propagated to lastLayerEE
  TICLLayerTileBarrel tracksterPropTilesBarrel;
  TICLLayerTileBarrel tsPropIntTilesBarrel;

  if (TICLInterpretationAlgoBase::algo_verbosity_ > VerbosityLevel::Advanced)
    LogDebug("GeneralInterpretationAlgo") << "------- Geometric Linking ------- \n";

  // Propagate tracks
  std::vector<unsigned> candidateTrackIds;
  candidateTrackIds.reserve(tracks.size());
  for (unsigned i = 0; i < tracks.size(); ++i) {
    if (!maskTracks[i])
      continue;
    candidateTrackIds.push_back(i);
  }

  std::sort(candidateTrackIds.begin(), candidateTrackIds.end(), [&](unsigned i, unsigned j) {
    return tracks[i].p() > tracks[j].p();
  });

  for (auto const i : candidateTrackIds) {
    const auto &tk = tracks[i];
    const auto &fts = trajectoryStateTransform::outerFreeState((tk), bFieldProd);
    if (detector_ == "HGCAL") {
      int iSide = int(tk.eta() > 0);
      // to the HGCal front
      const auto &tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
      if (tsos.isValid()) {
        Vector trackP(tsos.globalPosition().x(), tsos.globalPosition().y(), tsos.globalPosition().z());
        trackPColl.emplace_back(trackP, i);
      }
      // to lastLayerEE
      const auto &tsos_int = prop.propagate(fts, interfaceDisk_[iSide]->surface());
      if (tsos_int.isValid()) {
        Vector trackP(tsos_int.globalPosition().x(), tsos_int.globalPosition().y(), tsos_int.globalPosition().z());
        tkPropIntColl.emplace_back(trackP, i);
      }
    } else {
      GlobalPoint barrelSurfacePoint = rhtools_.getPositionLayer(0, false, true);
      const auto &tsos = prop.propagate(fts, ecalCylinder_->fastTangent(barrelSurfacePoint));
      if (tsos.isValid()) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        Vector trackP(tsos.globalPosition().x(), tsos.globalPosition().y(), tsos.globalPosition().z());
        trackPColl.emplace_back(trackP, i);
      }
      GlobalPoint barrelInterfacePoint = rhtools_.getPositionLayer(1, false, true);
      const auto &tsos_int = prop.propagate(fts, hcalCylinder_->fastTangent(barrelInterfacePoint));
      if (tsos_int.isValid()) {
        Vector trackP(tsos_int.globalPosition().x(), tsos_int.globalPosition().y(), tsos_int.globalPosition().z());
        tkPropIntColl.emplace_back(trackP, i);
      }
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      std::cout << "trackPColl: " << trackPColl.size() << std::endl;
      std::cout << "tkPropIntColl: " << trackPColl.size() << std::endl;
    }
  }  // Tracks
  tkPropIntColl.shrink_to_fit();
  trackPColl.shrink_to_fit();
  candidateTrackIds.shrink_to_fit();

  // Propagate tracksters

  // Record postions of all tracksters propagated to layer 1 and lastLayerEE,
  // to be used later for distance calculation in the link finding stage
  // indexed by trackster index in event collection
  std::vector<Vector> tsAllProp;
  std::vector<Vector> tsAllPropInt;
  tsAllProp.reserve(tracksters.size());
  tsAllPropInt.reserve(tracksters.size());
  // Propagate tracksters

  if (detector_ == "Barrel") {
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    std::cout << "nTracksters: " << tracksters.size() << std::endl;
    std::cout << "nCandidateTracks: " << candidateTrackIds.size() << std::endl;
  }
  for (unsigned i = 0; i < tracksters.size(); ++i) {
    const auto &t = tracksters[i];
    if (TICLInterpretationAlgoBase::algo_verbosity_ > VerbosityLevel::Advanced)
      LogDebug("GeneralInterpretationAlgo")
          << "trackster " << i << " - eta " << t.barycenter().eta() << " phi " << t.barycenter().phi() << " time "
          << t.time() << " energy " << t.raw_energy() << "\n";

    if (detector_ == "HGCAL") {
      // to HGCal front
      float zVal = hgcons_->waferZ(1, true);
      auto tsP = propagateTrackster(t, i, zVal, tracksterPropTiles, tracksterPropTilesBarrel);
      tsAllProp.emplace_back(tsP);

      // to lastLayerEE
      zVal = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z();
      tsP = propagateTrackster(t, i, zVal, tsPropIntTiles, tsPropIntTilesBarrel);
      tsAllPropInt.emplace_back(tsP);
    } else {
      float ecalR = rhtools_.getPositionLayer(0, false, true).mag();
      auto tsP = propagateTrackster(t, i, ecalR, tracksterPropTiles, tracksterPropTilesBarrel);
      tsAllProp.emplace_back(tsP);

      float hcalR = rhtools_.getPositionLayer(1, false, true).mag();
      tsP = propagateTrackster(t, i, hcalR, tsPropIntTiles, tsPropIntTilesBarrel);
      tsAllPropInt.emplace_back(tsP);
    }   
  }  // TS

  // step 1: tracks -> all tracksters, at firstLayerEE
  std::vector<std::vector<unsigned>> tsNearTk(tracks.size());
  findTrackstersInWindow(
      tracksters, trackPColl, tracksterPropTiles, tracksterPropTilesBarrel, tsAllProp, del_tk_ts_layer1_, tracksters.size(), tsNearTk);

  // step 2: tracks -> all tracksters, at lastLayerEE
  std::vector<std::vector<unsigned>> tsNearTkAtInt(tracks.size());
  findTrackstersInWindow(
      tracksters, tkPropIntColl, tsPropIntTiles, tsPropIntTilesBarrel, tsAllPropInt, del_tk_ts_int_, tracksters.size(), tsNearTkAtInt);

  std::vector<unsigned int> chargedHadronsFromTk;
  std::vector<std::vector<unsigned int>> trackstersInTrackIndices;
  trackstersInTrackIndices.resize(tracks.size());

  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  std::vector<bool> chargedMask(tracksters.size(), true);
  for (unsigned &i : candidateTrackIds) {
    std::cout << "tsNearTk[i].size(): " << tsNearTk[i].size() << std::endl;
    std::cout << "tsNearTkAtInt[i].size(): " << tsNearTkAtInt[i].size() << std::endl;
    if (tsNearTk[i].empty() && tsNearTkAtInt[i].empty()) {  // nothing linked to track, make charged hadrons
      continue;
    }

    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    std::vector<unsigned int> chargedCandidate;
    float total_raw_energy = 0.f;

    float track_time = 0.f;
    float track_timeErr = 0.f;
    float track_quality = 0.f;
    float track_beta = 0.f;
    GlobalPoint track_MtdPos{0.f, 0.f, 0.f};
    if (useMTDTiming) {
      auto const &inputTimingView = (*inputTiming_h).const_view();
      track_time = inputTimingView.time()[i];
      track_timeErr = inputTimingView.timeErr()[i];
      track_quality = inputTimingView.MVAquality()[i];
      track_beta = inputTimingView.beta()[i];
      track_MtdPos = {
          inputTimingView.posInMTD_x()[i], inputTimingView.posInMTD_y()[i], inputTimingView.posInMTD_z()[i]};
    }

    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    for (auto const tsIdx : tsNearTk[i]) {
      std::cout << "tsIdx: " << tsIdx << std::endl;
      std::cout << "chargedMask[tsIdx]: " << chargedMask[tsIdx] << std::endl;
      bool isCompatible = (chargedMask[tsIdx] && timeAndEnergyCompatible(total_raw_energy,
                                                        tracks[i],
                                                        tracksters[tsIdx],
                                                        track_time,
                                                        track_timeErr,
                                                        track_quality,
                                                        track_beta,
                                                        track_MtdPos,
                                                        useMTDTiming));
        if (chargedMask[tsIdx]) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        chargedCandidate.push_back(tsIdx);
        chargedMask[tsIdx] = false;
        total_raw_energy += tracksters[tsIdx].raw_energy();
      }
    }
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    for (const unsigned tsIdx : tsNearTkAtInt[i]) {  // do the same for tk -> ts links at the interface
      if (chargedMask[tsIdx] && timeAndEnergyCompatible(total_raw_energy,
                                                        tracks[i],
                                                        tracksters[tsIdx],
                                                        track_time,
                                                        track_timeErr,
                                                        track_quality,
                                                        track_beta,
                                                        track_MtdPos,
                                                        useMTDTiming)) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        chargedCandidate.push_back(tsIdx);
        chargedMask[tsIdx] = false;
        total_raw_energy += tracksters[tsIdx].raw_energy();
      }
    }
    trackstersInTrackIndices[i] = chargedCandidate;
  }

  for (size_t iTrack = 0; iTrack < trackstersInTrackIndices.size(); iTrack++) {
    if (!trackstersInTrackIndices[iTrack].empty()) {
      if (trackstersInTrackIndices[iTrack].size() == 1) {
        auto tracksterId = trackstersInTrackIndices[iTrack][0];
        resultCandidate[iTrack] = resultTracksters.size();
        resultTracksters.push_back(input.tracksters[tracksterId]);
      } else {
        // in this case mergeTracksters() clears the pid probabilities and the regressed energy is not set
        // TODO: fix probabilities when CNN will be splitted
        Trackster outTrackster;
        bool isHadron = false;
        for (auto const tracksterId : trackstersInTrackIndices[iTrack]) {
          //maskTracksters[tracksterId] = 0;
          outTrackster.mergeTracksters(input.tracksters[tracksterId]);
          if (input.tracksters[tracksterId].isHadronic())
            isHadron = true;
        }
        resultCandidate[iTrack] = resultTracksters.size();
        resultTracksters.push_back(outTrackster);
        // since a track has been linked it can only be electron or charged hadron
        if (isHadron)
          resultTracksters.back().setIdProbability(ticl::Trackster::ParticleType::charged_hadron, 1.f);
        else
          resultTracksters.back().setIdProbability(ticl::Trackster::ParticleType::electron, 1.f);
      }
    }
  }

  for (size_t iTrackster = 0; iTrackster < input.tracksters.size(); iTrackster++) {
    if (chargedMask[iTrackster]) {
      resultTracksters.push_back(input.tracksters[iTrackster]);
    }
  }
};

void GeneralInterpretationAlgo::fillPSetDescription(edm::ParameterSetDescription &desc) {
  desc.add<double>("delta_tk_ts_layer1", 0.02);
  desc.add<double>("delta_tk_ts_interface", 0.03);
  desc.add<double>("timing_quality_threshold", 0.5);
  TICLInterpretationAlgoBase::fillPSetDescription(desc);
}
