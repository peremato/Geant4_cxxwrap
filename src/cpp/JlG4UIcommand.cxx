// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4UIcommand> : std::false_type { };
  template<> struct DefaultConstructible<G4UIcommand> : std::false_type { };
}

// Class generating the wrapper for type G4UIcommand
// signature to use in the veto file: G4UIcommand
struct JlG4UIcommand: public Wrapper {

  JlG4UIcommand(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4UIcommand (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4UIcommand.hh:51:7
    jlcxx::TypeWrapper<G4UIcommand>  t = jlModule.add_type<G4UIcommand>("G4UIcommand");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4UIcommand>>(new jlcxx::TypeWrapper<G4UIcommand>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4UIcommand>> type_;
};
std::shared_ptr<Wrapper> newJlG4UIcommand(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4UIcommand(module));
}
