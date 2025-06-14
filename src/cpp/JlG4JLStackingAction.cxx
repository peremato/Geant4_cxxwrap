// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4JLStackingAction> : std::false_type { };
  template<> struct DefaultConstructible<G4JLStackingAction> : std::false_type { };
template<> struct SuperType<G4JLStackingAction> { typedef G4UserStackingAction type; };
}

// Class generating the wrapper for type G4JLStackingAction
// signature to use in the veto file: G4JLStackingAction
struct JlG4JLStackingAction: public Wrapper {

  JlG4JLStackingAction(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4JLStackingAction (" __HERE__ ")");
    // defined in cpp/Geant4Wrap.h:237:7
    jlcxx::TypeWrapper<G4JLStackingAction>  t = jlModule.add_type<G4JLStackingAction>("G4JLStackingAction",
      jlcxx::julia_base_type<G4UserStackingAction>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4JLStackingAction>>(new jlcxx::TypeWrapper<G4JLStackingAction>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes    );


    DEBUG_MSG("Adding wrapper for void G4JLStackingAction::G4JLStackingAction(classify_f, void *, stackaction_f, void *, stackaction_f, void *) (" __HERE__ ")");
    // defined in cpp/Geant4Wrap.h:239:5
    t.constructor<classify_f>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("classify")    );
    t.constructor<classify_f, void *>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("classify"), jlcxx::arg("classifyd")    );
    t.constructor<classify_f, void *, stackaction_f>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("classify"), jlcxx::arg("classifyd"), jlcxx::arg("newstage")    );
    t.constructor<classify_f, void *, stackaction_f, void *>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("classify"), jlcxx::arg("classifyd"), jlcxx::arg("newstage"), jlcxx::arg("newstaged")    );
    t.constructor<classify_f, void *, stackaction_f, void *, stackaction_f>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("classify"), jlcxx::arg("classifyd"), jlcxx::arg("newstage"), jlcxx::arg("newstaged"), jlcxx::arg("newevent")    );
    t.constructor<classify_f, void *, stackaction_f, void *, stackaction_f, void *>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("classify"), jlcxx::arg("classifyd"), jlcxx::arg("newstage"), jlcxx::arg("newstaged"), jlcxx::arg("newevent"), jlcxx::arg("neweventd")    );

    DEBUG_MSG("Adding wrapper for G4ClassificationOfNewTrack G4JLStackingAction::ClassifyNewTrack(const G4Track *) (" __HERE__ ")");
    // signature to use in the veto list: G4ClassificationOfNewTrack G4JLStackingAction::ClassifyNewTrack(const G4Track *)
    // defined in cpp/Geant4Wrap.h:246:40
    t.method("ClassifyNewTrack", [](G4JLStackingAction& a, const G4Track * arg0)->G4ClassificationOfNewTrack { return a.ClassifyNewTrack(arg0); }, jlcxx::arg("this"), jlcxx::arg("track"));
    t.method("ClassifyNewTrack", [](G4JLStackingAction* a, const G4Track * arg0)->G4ClassificationOfNewTrack { return a->ClassifyNewTrack(arg0); }, jlcxx::arg("this"), jlcxx::arg("track"));

    DEBUG_MSG("Adding wrapper for void G4JLStackingAction::NewStage() (" __HERE__ ")");
    // signature to use in the veto list: void G4JLStackingAction::NewStage()
    // defined in cpp/Geant4Wrap.h:253:18
    t.method("NewStage", [](G4JLStackingAction& a)->void { a.NewStage(); }, jlcxx::arg("this"));
    t.method("NewStage", [](G4JLStackingAction* a)->void { a->NewStage(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4JLStackingAction::PrepareNewEvent() (" __HERE__ ")");
    // signature to use in the veto list: void G4JLStackingAction::PrepareNewEvent()
    // defined in cpp/Geant4Wrap.h:254:18
    t.method("PrepareNewEvent", [](G4JLStackingAction& a)->void { a.PrepareNewEvent(); }, jlcxx::arg("this"));
    t.method("PrepareNewEvent", [](G4JLStackingAction* a)->void { a->PrepareNewEvent(); }, jlcxx::arg("this"));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4JLStackingAction>> type_;
};
std::shared_ptr<Wrapper> newJlG4JLStackingAction(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4JLStackingAction(module));
}
