// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4SingleParticleSource> : std::false_type { };
  template<> struct DefaultConstructible<G4SingleParticleSource> : std::false_type { };
template<> struct SuperType<G4SingleParticleSource> { typedef G4VPrimaryGenerator type; };
}

// Class generating the wrapper for type G4SingleParticleSource
// signature to use in the veto file: G4SingleParticleSource
struct JlG4SingleParticleSource: public Wrapper {

  JlG4SingleParticleSource(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4SingleParticleSource (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:70:7
    jlcxx::TypeWrapper<G4SingleParticleSource>  t = jlModule.add_type<G4SingleParticleSource>("G4SingleParticleSource",
      jlcxx::julia_base_type<G4VPrimaryGenerator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4SingleParticleSource>>(new jlcxx::TypeWrapper<G4SingleParticleSource>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for void G4SingleParticleSource::GeneratePrimaryVertex(G4Event *) (" __HERE__ ")");
    // signature to use in the veto list: void G4SingleParticleSource::GeneratePrimaryVertex(G4Event *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:81:10
    t.method("GeneratePrimaryVertex", static_cast<void (G4SingleParticleSource::*)(G4Event *) >(&G4SingleParticleSource::GeneratePrimaryVertex));

    DEBUG_MSG("Adding wrapper for G4SPSPosDistribution * G4SingleParticleSource::GetPosDist() (" __HERE__ ")");
    // signature to use in the veto list: G4SPSPosDistribution * G4SingleParticleSource::GetPosDist()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:84:34
    t.method("GetPosDist", static_cast<G4SPSPosDistribution * (G4SingleParticleSource::*)()  const>(&G4SingleParticleSource::GetPosDist));

    DEBUG_MSG("Adding wrapper for G4SPSAngDistribution * G4SingleParticleSource::GetAngDist() (" __HERE__ ")");
    // signature to use in the veto list: G4SPSAngDistribution * G4SingleParticleSource::GetAngDist()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:87:34
    t.method("GetAngDist", static_cast<G4SPSAngDistribution * (G4SingleParticleSource::*)()  const>(&G4SingleParticleSource::GetAngDist));

    DEBUG_MSG("Adding wrapper for G4SPSEneDistribution * G4SingleParticleSource::GetEneDist() (" __HERE__ ")");
    // signature to use in the veto list: G4SPSEneDistribution * G4SingleParticleSource::GetEneDist()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:90:34
    t.method("GetEneDist", static_cast<G4SPSEneDistribution * (G4SingleParticleSource::*)()  const>(&G4SingleParticleSource::GetEneDist));

    DEBUG_MSG("Adding wrapper for G4SPSRandomGenerator * G4SingleParticleSource::GetBiasRndm() (" __HERE__ ")");
    // signature to use in the veto list: G4SPSRandomGenerator * G4SingleParticleSource::GetBiasRndm()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:93:34
    t.method("GetBiasRndm", static_cast<G4SPSRandomGenerator * (G4SingleParticleSource::*)()  const>(&G4SingleParticleSource::GetBiasRndm));

    DEBUG_MSG("Adding wrapper for void G4SingleParticleSource::SetVerbosity(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4SingleParticleSource::SetVerbosity(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:96:10
    t.method("SetVerbosity", static_cast<void (G4SingleParticleSource::*)(G4int) >(&G4SingleParticleSource::SetVerbosity));

    DEBUG_MSG("Adding wrapper for void G4SingleParticleSource::SetParticleDefinition(G4ParticleDefinition *) (" __HERE__ ")");
    // signature to use in the veto list: void G4SingleParticleSource::SetParticleDefinition(G4ParticleDefinition *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:99:10
    t.method("SetParticleDefinition", static_cast<void (G4SingleParticleSource::*)(G4ParticleDefinition *) >(&G4SingleParticleSource::SetParticleDefinition));

    DEBUG_MSG("Adding wrapper for G4ParticleDefinition * G4SingleParticleSource::GetParticleDefinition() (" __HERE__ ")");
    // signature to use in the veto list: G4ParticleDefinition * G4SingleParticleSource::GetParticleDefinition()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:100:34
    t.method("GetParticleDefinition", static_cast<G4ParticleDefinition * (G4SingleParticleSource::*)()  const>(&G4SingleParticleSource::GetParticleDefinition));

    DEBUG_MSG("Adding wrapper for void G4SingleParticleSource::SetParticleCharge(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SingleParticleSource::SetParticleCharge(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:104:17
    t.method("SetParticleCharge", static_cast<void (G4SingleParticleSource::*)(G4double) >(&G4SingleParticleSource::SetParticleCharge));

    DEBUG_MSG("Adding wrapper for void G4SingleParticleSource::SetParticlePolarization(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4SingleParticleSource::SetParticlePolarization(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:107:17
    t.method("SetParticlePolarization", static_cast<void (G4SingleParticleSource::*)(const G4ThreeVector &) >(&G4SingleParticleSource::SetParticlePolarization));

    DEBUG_MSG("Adding wrapper for const G4ThreeVector & G4SingleParticleSource::GetParticlePolarization() (" __HERE__ ")");
    // signature to use in the veto list: const G4ThreeVector & G4SingleParticleSource::GetParticlePolarization()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:109:33
    t.method("GetParticlePolarization", static_cast<const G4ThreeVector & (G4SingleParticleSource::*)()  const>(&G4SingleParticleSource::GetParticlePolarization));

    DEBUG_MSG("Adding wrapper for void G4SingleParticleSource::SetParticleTime(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SingleParticleSource::SetParticleTime(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:113:17
    t.method("SetParticleTime", static_cast<void (G4SingleParticleSource::*)(G4double) >(&G4SingleParticleSource::SetParticleTime));

    DEBUG_MSG("Adding wrapper for G4double G4SingleParticleSource::GetParticleTime() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4SingleParticleSource::GetParticleTime()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:114:21
    t.method("GetParticleTime", static_cast<G4double (G4SingleParticleSource::*)()  const>(&G4SingleParticleSource::GetParticleTime));

    DEBUG_MSG("Adding wrapper for void G4SingleParticleSource::SetNumberOfParticles(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4SingleParticleSource::SetNumberOfParticles(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:117:17
    t.method("SetNumberOfParticles", static_cast<void (G4SingleParticleSource::*)(G4int) >(&G4SingleParticleSource::SetNumberOfParticles));

    DEBUG_MSG("Adding wrapper for G4int G4SingleParticleSource::GetNumberOfParticles() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4SingleParticleSource::GetNumberOfParticles()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:119:18
    t.method("GetNumberOfParticles", static_cast<G4int (G4SingleParticleSource::*)()  const>(&G4SingleParticleSource::GetNumberOfParticles));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4SingleParticleSource::GetParticlePosition() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4SingleParticleSource::GetParticlePosition()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:123:26
    t.method("GetParticlePosition", static_cast<G4ThreeVector (G4SingleParticleSource::*)()  const>(&G4SingleParticleSource::GetParticlePosition));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4SingleParticleSource::GetParticleMomentumDirection() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4SingleParticleSource::GetParticleMomentumDirection()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:125:26
    t.method("GetParticleMomentumDirection", static_cast<G4ThreeVector (G4SingleParticleSource::*)()  const>(&G4SingleParticleSource::GetParticleMomentumDirection));

    DEBUG_MSG("Adding wrapper for G4double G4SingleParticleSource::GetParticleEnergy() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4SingleParticleSource::GetParticleEnergy()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SingleParticleSource.hh:127:21
    t.method("GetParticleEnergy", static_cast<G4double (G4SingleParticleSource::*)()  const>(&G4SingleParticleSource::GetParticleEnergy));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4SingleParticleSource>> type_;
};
std::shared_ptr<Wrapper> newJlG4SingleParticleSource(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4SingleParticleSource(module));
}
