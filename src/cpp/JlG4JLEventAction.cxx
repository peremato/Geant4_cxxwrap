// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4JLEventAction> : std::false_type { };
  template<> struct DefaultConstructible<G4JLEventAction> : std::false_type { };
template<> struct SuperType<G4JLEventAction> { typedef G4UserEventAction type; };
}

// Class generating the wrapper for type G4JLEventAction
// signature to use in the veto file: G4JLEventAction
struct JlG4JLEventAction: public Wrapper {

  JlG4JLEventAction(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4JLEventAction (" __HERE__ ")");
    // defined in ./cpp/Geant4Wrap.h:219:7
    jlcxx::TypeWrapper<G4JLEventAction>  t = jlModule.add_type<G4JLEventAction>("G4JLEventAction",
      jlcxx::julia_base_type<G4UserEventAction>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4JLEventAction>>(new jlcxx::TypeWrapper<G4JLEventAction>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void G4JLEventAction::G4JLEventAction(eventaction_f, void *, eventaction_f, void *) (" __HERE__ ")");
    // defined in ./cpp/Geant4Wrap.h:221:5
    t.constructor<eventaction_f>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<eventaction_f, void *>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<eventaction_f, void *, eventaction_f>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<eventaction_f, void *, eventaction_f, void *>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void G4JLEventAction::BeginOfEventAction(const G4Event *) (" __HERE__ ")");
    // signature to use in the veto list: void G4JLEventAction::BeginOfEventAction(const G4Event *)
    // defined in ./cpp/Geant4Wrap.h:225:18
    t.method("BeginOfEventAction", static_cast<void (G4JLEventAction::*)(const G4Event *) >(&G4JLEventAction::BeginOfEventAction));

    DEBUG_MSG("Adding wrapper for void G4JLEventAction::EndOfEventAction(const G4Event *) (" __HERE__ ")");
    // signature to use in the veto list: void G4JLEventAction::EndOfEventAction(const G4Event *)
    // defined in ./cpp/Geant4Wrap.h:226:20
    t.method("EndOfEventAction", static_cast<void (G4JLEventAction::*)(const G4Event *) >(&G4JLEventAction::EndOfEventAction));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4JLEventAction>> type_;
};
std::shared_ptr<Wrapper> newJlG4JLEventAction(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4JLEventAction(module));
}
