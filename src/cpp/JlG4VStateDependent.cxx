// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4VStateDependent> : std::false_type { };
  template<> struct DefaultConstructible<G4VStateDependent> : std::false_type { };
}

// Class generating the wrapper for type G4VStateDependent
// signature to use in the veto file: G4VStateDependent
struct JlG4VStateDependent: public Wrapper {

  JlG4VStateDependent(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4VStateDependent (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VStateDependent.hh:43:7
    jlcxx::TypeWrapper<G4VStateDependent>  t = jlModule.add_type<G4VStateDependent>("G4VStateDependent");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4VStateDependent>>(new jlcxx::TypeWrapper<G4VStateDependent>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;

    module_.set_override_module(jl_base_module);

    DEBUG_MSG("Adding wrapper for G4bool G4VStateDependent::operator==(const G4VStateDependent &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VStateDependent::operator==(const G4VStateDependent &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VStateDependent.hh:48:10
    t.method("==", static_cast<G4bool (G4VStateDependent::*)(const G4VStateDependent &)  const>(&G4VStateDependent::operator==));

    DEBUG_MSG("Adding wrapper for G4bool G4VStateDependent::operator!=(const G4VStateDependent &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VStateDependent::operator!=(const G4VStateDependent &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VStateDependent.hh:49:10
    t.method("!=", static_cast<G4bool (G4VStateDependent::*)(const G4VStateDependent &)  const>(&G4VStateDependent::operator!=));

    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for G4bool G4VStateDependent::Notify(G4ApplicationState) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VStateDependent::Notify(G4ApplicationState)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VStateDependent.hh:51:18
    t.method("Notify", static_cast<G4bool (G4VStateDependent::*)(G4ApplicationState) >(&G4VStateDependent::Notify));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4VStateDependent>> type_;
};
std::shared_ptr<Wrapper> newJlG4VStateDependent(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4VStateDependent(module));
}
