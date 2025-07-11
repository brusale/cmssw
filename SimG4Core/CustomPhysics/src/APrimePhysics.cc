#include "SimG4Core/CustomPhysics/interface/APrimePhysics.h"
#include "SimG4Core/CustomPhysics/interface/CMSAPrime.h"
#include "SimG4Core/CustomPhysics/interface/CMSmuDarkBremsstrahlung.h"
// Geant 4
#include "G4Electron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4ProcessManager.hh"
#include <CLHEP/Units/SystemOfUnits.h>

APrimePhysics::APrimePhysics(double APMass, const G4String& scalefile, const G4double cxBias, const G4String& name)
    : G4VPhysicsConstructor(name), aprimeDef_(nullptr) {
  apmass = APMass;
  mgfile = scalefile;
  biasFactor = cxBias;
}

void APrimePhysics::ConstructParticle() {
  /**
        * Insert A-prime into the Geant4 particle table.
        * For now we flag it as stable.
        */
  aprimeDef_ = CMSAPrime::APrime(apmass);
}

void APrimePhysics::ConstructProcess() {
  G4ParticleDefinition* muonminus = G4MuonMinus::MuonMinusDefinition();
  G4ParticleDefinition* muonplus = G4MuonPlus::MuonPlusDefinition();
  G4ProcessManager* pmplus = muonplus->GetProcessManager();
  G4ProcessManager* pmminus = muonminus->GetProcessManager();
  pmplus->AddDiscreteProcess(new CMSmuDarkBremsstrahlung(mgfile, biasFactor), 6);
  pmminus->AddDiscreteProcess(new CMSmuDarkBremsstrahlung(mgfile, biasFactor), 6);
}
