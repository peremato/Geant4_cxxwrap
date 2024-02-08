// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4JLSensDet> : std::false_type { };
  template<> struct DefaultConstructible<G4JLSensDet> : std::false_type { };
template<> struct SuperType<G4JLSensDet> { typedef G4VSensitiveDetector type; };
}

// Class generating the wrapper for type G4JLSensDet
// signature to use in the veto file: G4JLSensDet
struct JlG4JLSensDet: public Wrapper {

  JlG4JLSensDet(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4JLSensDet (" __HERE__ ")");
    // defined in ./cpp/Geant4Wrap.h:46:7
    jlcxx::TypeWrapper<G4JLSensDet>  t = jlModule.add_type<G4JLSensDet>("G4JLSensDet",
      jlcxx::julia_base_type<G4VSensitiveDetector>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4JLSensDet>>(new jlcxx::TypeWrapper<G4JLSensDet>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4JLSensDet::G4JLSensDet(const G4String &, processhits_f, void *) (" __HERE__ ")");
    // defined in ./cpp/Geant4Wrap.h:48:5
    t.constructor<const G4String &, processhits_f, void *>(/*finalize=*/false);

    DEBUG_MSG("Adding wrapper for void G4JLSensDet::Initialize(G4HCofThisEvent *) (" __HERE__ ")");
    // signature to use in the veto list: void G4JLSensDet::Initialize(G4HCofThisEvent *)
    // defined in ./cpp/Geant4Wrap.h:55:18
    t.method("Initialize", static_cast<void (G4JLSensDet::*)(G4HCofThisEvent *) >(&G4JLSensDet::Initialize));

    DEBUG_MSG("Adding wrapper for void G4JLSensDet::EndOfEvent(G4HCofThisEvent *) (" __HERE__ ")");
    // signature to use in the veto list: void G4JLSensDet::EndOfEvent(G4HCofThisEvent *)
    // defined in ./cpp/Geant4Wrap.h:56:18
    t.method("EndOfEvent", static_cast<void (G4JLSensDet::*)(G4HCofThisEvent *) >(&G4JLSensDet::EndOfEvent));

    DEBUG_MSG("Adding wrapper for G4bool G4JLSensDet::ProcessHits(G4Step *, G4TouchableHistory *) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4JLSensDet::ProcessHits(G4Step *, G4TouchableHistory *)
    // defined in ./cpp/Geant4Wrap.h:57:20
    t.method("ProcessHits", static_cast<G4bool (G4JLSensDet::*)(G4Step *, G4TouchableHistory *) >(&G4JLSensDet::ProcessHits));

    DEBUG_MSG("Adding wrapper for void G4JLSensDet::SetInitialize(initend_f, void *) (" __HERE__ ")");
    // signature to use in the veto list: void G4JLSensDet::SetInitialize(initend_f, void *)
    // defined in ./cpp/Geant4Wrap.h:58:10
    t.method("SetInitialize", static_cast<void (G4JLSensDet::*)(initend_f, void *) >(&G4JLSensDet::SetInitialize));

    DEBUG_MSG("Adding wrapper for void G4JLSensDet::SetEndOfEvent(initend_f, void *) (" __HERE__ ")");
    // signature to use in the veto list: void G4JLSensDet::SetEndOfEvent(initend_f, void *)
    // defined in ./cpp/Geant4Wrap.h:59:10
    t.method("SetEndOfEvent", static_cast<void (G4JLSensDet::*)(initend_f, void *) >(&G4JLSensDet::SetEndOfEvent));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4JLSensDet>> type_;
};
std::shared_ptr<Wrapper> newJlG4JLSensDet(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4JLSensDet(module));
}
