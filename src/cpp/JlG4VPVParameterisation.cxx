// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4VPVParameterisation> : std::false_type { };
  template<> struct DefaultConstructible<G4VPVParameterisation> : std::false_type { };
}

// Class generating the wrapper for type G4VPVParameterisation
// signature to use in the veto file: G4VPVParameterisation
struct JlG4VPVParameterisation: public Wrapper {

  JlG4VPVParameterisation(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4VPVParameterisation (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VPVParameterisation.hh:68:7
    jlcxx::TypeWrapper<G4VPVParameterisation>  t = jlModule.add_type<G4VPVParameterisation>("G4VPVParameterisation");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4VPVParameterisation>>(new jlcxx::TypeWrapper<G4VPVParameterisation>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4VPVParameterisation>> type_;
};
std::shared_ptr<Wrapper> newJlG4VPVParameterisation(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4VPVParameterisation(module));
}