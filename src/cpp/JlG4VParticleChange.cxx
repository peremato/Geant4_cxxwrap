// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4VParticleChange> : std::false_type { };
  template<> struct DefaultConstructible<G4VParticleChange> : std::false_type { };
}

// Class generating the wrapper for type G4VParticleChange
// signature to use in the veto file: G4VParticleChange
struct JlG4VParticleChange: public Wrapper {

  JlG4VParticleChange(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4VParticleChange (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VParticleChange.hh:68:7
    jlcxx::TypeWrapper<G4VParticleChange>  t = jlModule.add_type<G4VParticleChange>("G4VParticleChange");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4VParticleChange>>(new jlcxx::TypeWrapper<G4VParticleChange>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4VParticleChange>> type_;
};
std::shared_ptr<Wrapper> newJlG4VParticleChange(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4VParticleChange(module));
}
