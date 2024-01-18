// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4PhysicsListHelper> : std::false_type { };
  template<> struct DefaultConstructible<G4PhysicsListHelper> : std::false_type { };
}

// Class generating the wrapper for type G4PhysicsListHelper
// signature to use in the veto file: G4PhysicsListHelper
struct JlG4PhysicsListHelper: public Wrapper {

  JlG4PhysicsListHelper(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4PhysicsListHelper (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PhysicsListHelper.hh:49:7
    jlcxx::TypeWrapper<G4PhysicsListHelper>  t = jlModule.add_type<G4PhysicsListHelper>("G4PhysicsListHelper");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4PhysicsListHelper>>(new jlcxx::TypeWrapper<G4PhysicsListHelper>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4PhysicsListHelper>> type_;
};
std::shared_ptr<Wrapper> newJlG4PhysicsListHelper(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4PhysicsListHelper(module));
}