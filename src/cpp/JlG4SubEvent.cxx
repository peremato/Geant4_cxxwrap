// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4SubEvent> : std::false_type { };
  template<> struct DefaultConstructible<G4SubEvent> : std::false_type { };
}

// Class generating the wrapper for type G4SubEvent
// signature to use in the veto file: G4SubEvent
struct JlG4SubEvent: public Wrapper {

  JlG4SubEvent(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4SubEvent (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SubEvent.hh:42:7
    jlcxx::TypeWrapper<G4SubEvent>  t = jlModule.add_type<G4SubEvent>("G4SubEvent");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4SubEvent>>(new jlcxx::TypeWrapper<G4SubEvent>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4SubEvent>> type_;
};
std::shared_ptr<Wrapper> newJlG4SubEvent(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4SubEvent(module));
}
