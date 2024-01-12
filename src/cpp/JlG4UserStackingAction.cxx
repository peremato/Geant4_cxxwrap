// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4UserStackingAction> : std::false_type { };
  template<> struct DefaultConstructible<G4UserStackingAction> : std::false_type { };
}

// Class generating the wrapper for type G4UserStackingAction
// signature to use in the veto file: G4UserStackingAction
struct JlG4UserStackingAction: public Wrapper {

  JlG4UserStackingAction(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4UserStackingAction (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4UserStackingAction.hh:44:7
    jlcxx::TypeWrapper<G4UserStackingAction>  t = jlModule.add_type<G4UserStackingAction>("G4UserStackingAction");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4UserStackingAction>>(new jlcxx::TypeWrapper<G4UserStackingAction>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4UserStackingAction>> type_;
};
std::shared_ptr<Wrapper> newJlG4UserStackingAction(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4UserStackingAction(module));
}
