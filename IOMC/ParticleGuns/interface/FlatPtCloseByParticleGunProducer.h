#ifndef IOMC_ParticleGun_FlatPtCloseByParticleGunProducer_H
#define IOMC_ParticleGun_FlatPtCloseByParticleGunProducer_H

#include "IOMC/ParticleGuns/interface/BaseFlatGunProducer.h"

namespace edm {

  class FlatPtCloseByParticleGunProducer : public BaseFlatGunProducer {
  public:
    FlatPtCloseByParticleGunProducer(const ParameterSet&);
    ~FlatPtCloseByParticleGunProducer() override;

  private:
    void produce(Event& e, const EventSetup& es) override;

  protected:
    // data members
    bool fControlledByEta;
    double fPtMin, fPtMax, fEtaMin, fEtaMax, fRMin, fRMax, fZMin, fZMax, fDelta, fPhiMin, fPhiMax;
    int fNParticles;
    bool fMaxPtSpread = false;
    bool fPointing = false;
    bool fOverlapping = false;
    bool fRandomShoot = false;
    std::vector<int> fPartIDs;
  };
}  // namespace edm

#endif
