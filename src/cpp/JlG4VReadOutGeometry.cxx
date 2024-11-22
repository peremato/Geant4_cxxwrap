// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4VReadOutGeometry> : std::false_type { };
  template<> struct DefaultConstructible<G4VReadOutGeometry> : std::false_type { };
}

// Class generating the wrapper for type G4VReadOutGeometry
// signature to use in the veto file: G4VReadOutGeometry
struct JlG4VReadOutGeometry: public Wrapper {

  JlG4VReadOutGeometry(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4VReadOutGeometry (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VReadOutGeometry.hh:39:7
    jlcxx::TypeWrapper<G4VReadOutGeometry>  t = jlModule.add_type<G4VReadOutGeometry>("G4VReadOutGeometry");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4VReadOutGeometry>>(new jlcxx::TypeWrapper<G4VReadOutGeometry>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4VReadOutGeometry>> type_;
};
std::shared_ptr<Wrapper> newJlG4VReadOutGeometry(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4VReadOutGeometry(module));
}
