// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4DCofThisEvent> : std::false_type { };
  template<> struct DefaultConstructible<G4DCofThisEvent> : std::false_type { };
}

// Class generating the wrapper for type G4DCofThisEvent
// signature to use in the veto file: G4DCofThisEvent
struct JlG4DCofThisEvent: public Wrapper {

  JlG4DCofThisEvent(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4DCofThisEvent (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4DCofThisEvent.hh:50:7
    jlcxx::TypeWrapper<G4DCofThisEvent>  t = jlModule.add_type<G4DCofThisEvent>("G4DCofThisEvent");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4DCofThisEvent>>(new jlcxx::TypeWrapper<G4DCofThisEvent>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes    );
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4DCofThisEvent>> type_;
};
std::shared_ptr<Wrapper> newJlG4DCofThisEvent(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4DCofThisEvent(module));
}
