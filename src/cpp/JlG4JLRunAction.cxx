// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4JLRunAction> : std::false_type { };
  template<> struct DefaultConstructible<G4JLRunAction> : std::false_type { };
template<> struct SuperType<G4JLRunAction> { typedef G4UserRunAction type; };
}

// Class generating the wrapper for type G4JLRunAction
// signature to use in the veto file: G4JLRunAction
struct JlG4JLRunAction: public Wrapper {

  JlG4JLRunAction(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4JLRunAction (" __HERE__ ")");
    // defined in ./cpp/Geant4Wrap.h:194:7
    jlcxx::TypeWrapper<G4JLRunAction>  t = jlModule.add_type<G4JLRunAction>("G4JLRunAction",
      jlcxx::julia_base_type<G4UserRunAction>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4JLRunAction>>(new jlcxx::TypeWrapper<G4JLRunAction>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/false);


    DEBUG_MSG("Adding wrapper for void G4JLRunAction::G4JLRunAction(runaction_f, void *, runaction_f, void *) (" __HERE__ ")");
    // defined in ./cpp/Geant4Wrap.h:196:5
    t.constructor<runaction_f>(/*finalize=*/false);
    t.constructor<runaction_f, void *>(/*finalize=*/false);
    t.constructor<runaction_f, void *, runaction_f>(/*finalize=*/false);
    t.constructor<runaction_f, void *, runaction_f, void *>(/*finalize=*/false);

    DEBUG_MSG("Adding wrapper for void G4JLRunAction::BeginOfRunAction(const G4Run *) (" __HERE__ ")");
    // signature to use in the veto list: void G4JLRunAction::BeginOfRunAction(const G4Run *)
    // defined in ./cpp/Geant4Wrap.h:200:18
    t.method("BeginOfRunAction", static_cast<void (G4JLRunAction::*)(const G4Run *) >(&G4JLRunAction::BeginOfRunAction));

    DEBUG_MSG("Adding wrapper for void G4JLRunAction::EndOfRunAction(const G4Run *) (" __HERE__ ")");
    // signature to use in the veto list: void G4JLRunAction::EndOfRunAction(const G4Run *)
    // defined in ./cpp/Geant4Wrap.h:201:20
    t.method("EndOfRunAction", static_cast<void (G4JLRunAction::*)(const G4Run *) >(&G4JLRunAction::EndOfRunAction));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4JLRunAction>> type_;
};
std::shared_ptr<Wrapper> newJlG4JLRunAction(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4JLRunAction(module));
}
