// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4PolyhedraHistorical> : std::false_type { };
  template<> struct DefaultConstructible<G4PolyhedraHistorical> : std::false_type { };
}

// Class generating the wrapper for type G4PolyhedraHistorical
// signature to use in the veto file: G4PolyhedraHistorical
struct JlG4PolyhedraHistorical: public Wrapper {

  JlG4PolyhedraHistorical(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4PolyhedraHistorical (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4PolyhedraHistorical.hh:39:7
    jlcxx::TypeWrapper<G4PolyhedraHistorical>  t = jlModule.add_type<G4PolyhedraHistorical>("G4PolyhedraHistorical");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4PolyhedraHistorical>>(new jlcxx::TypeWrapper<G4PolyhedraHistorical>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4PolyhedraHistorical>> type_;
};
std::shared_ptr<Wrapper> newJlG4PolyhedraHistorical(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4PolyhedraHistorical(module));
}
