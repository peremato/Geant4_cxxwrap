// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4StateManager> : std::false_type { };
  template<> struct DefaultConstructible<G4StateManager> : std::false_type { };
}

// Class generating the wrapper for type G4StateManager
// signature to use in the veto file: G4StateManager
struct JlG4StateManager: public Wrapper {

  JlG4StateManager(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4StateManager (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:49:7
    jlcxx::TypeWrapper<G4StateManager>  t = jlModule.add_type<G4StateManager>("G4StateManager");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4StateManager>>(new jlcxx::TypeWrapper<G4StateManager>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;

    DEBUG_MSG("Adding wrapper for G4StateManager * G4StateManager::GetStateManager() (" __HERE__ ")");
    // signature to use in the veto list: G4StateManager * G4StateManager::GetStateManager()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:52:26
    t.method("G4StateManager!GetStateManager", static_cast<G4StateManager * (*)() >(&G4StateManager::GetStateManager));

    DEBUG_MSG("Adding wrapper for const G4ApplicationState & G4StateManager::GetCurrentState() (" __HERE__ ")");
    // signature to use in the veto list: const G4ApplicationState & G4StateManager::GetCurrentState()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:64:29
    t.method("GetCurrentState", static_cast<const G4ApplicationState & (G4StateManager::*)()  const>(&G4StateManager::GetCurrentState));

    DEBUG_MSG("Adding wrapper for const G4ApplicationState & G4StateManager::GetPreviousState() (" __HERE__ ")");
    // signature to use in the veto list: const G4ApplicationState & G4StateManager::GetPreviousState()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:66:29
    t.method("GetPreviousState", static_cast<const G4ApplicationState & (G4StateManager::*)()  const>(&G4StateManager::GetPreviousState));

    DEBUG_MSG("Adding wrapper for G4bool G4StateManager::SetNewState(const G4ApplicationState &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4StateManager::SetNewState(const G4ApplicationState &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:68:10
    t.method("SetNewState", static_cast<G4bool (G4StateManager::*)(const G4ApplicationState &) >(&G4StateManager::SetNewState));

    DEBUG_MSG("Adding wrapper for G4bool G4StateManager::SetNewState(const G4ApplicationState &, const char *) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4StateManager::SetNewState(const G4ApplicationState &, const char *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:72:10
    t.method("SetNewState", static_cast<G4bool (G4StateManager::*)(const G4ApplicationState &, const char *) >(&G4StateManager::SetNewState));

    DEBUG_MSG("Adding wrapper for G4bool G4StateManager::RegisterDependent(G4VStateDependent *, G4bool) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4StateManager::RegisterDependent(G4VStateDependent *, G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:77:10
    t.method("RegisterDependent", static_cast<G4bool (G4StateManager::*)(G4VStateDependent *, G4bool) >(&G4StateManager::RegisterDependent));
    t.method("RegisterDependent", [](G4StateManager& a, G4VStateDependent * arg0)->G4bool { return a.RegisterDependent(arg0); });
    t.method("RegisterDependent", [](G4StateManager* a, G4VStateDependent * arg0)->G4bool { return a->RegisterDependent(arg0); });

    DEBUG_MSG("Adding wrapper for G4bool G4StateManager::DeregisterDependent(G4VStateDependent *) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4StateManager::DeregisterDependent(G4VStateDependent *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:83:10
    t.method("DeregisterDependent", static_cast<G4bool (G4StateManager::*)(G4VStateDependent *) >(&G4StateManager::DeregisterDependent));

    DEBUG_MSG("Adding wrapper for G4VStateDependent * G4StateManager::RemoveDependent(const G4VStateDependent *) (" __HERE__ ")");
    // signature to use in the veto list: G4VStateDependent * G4StateManager::RemoveDependent(const G4VStateDependent *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:86:22
    t.method("RemoveDependent", static_cast<G4VStateDependent * (G4StateManager::*)(const G4VStateDependent *) >(&G4StateManager::RemoveDependent));

    DEBUG_MSG("Adding wrapper for G4String G4StateManager::GetStateString(const G4ApplicationState &) (" __HERE__ ")");
    // signature to use in the veto list: G4String G4StateManager::GetStateString(const G4ApplicationState &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:89:12
    t.method("GetStateString", static_cast<G4String (G4StateManager::*)(const G4ApplicationState &)  const>(&G4StateManager::GetStateString));

    DEBUG_MSG("Adding wrapper for void G4StateManager::SetSuppressAbortion(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4StateManager::SetSuppressAbortion(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:92:15
    t.method("SetSuppressAbortion", static_cast<void (G4StateManager::*)(G4int) >(&G4StateManager::SetSuppressAbortion));

    DEBUG_MSG("Adding wrapper for G4int G4StateManager::GetSuppressAbortion() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4StateManager::GetSuppressAbortion()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:93:16
    t.method("GetSuppressAbortion", static_cast<G4int (G4StateManager::*)()  const>(&G4StateManager::GetSuppressAbortion));

    DEBUG_MSG("Adding wrapper for const char * G4StateManager::GetMessage() (" __HERE__ ")");
    // signature to use in the veto list: const char * G4StateManager::GetMessage()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:94:22
    t.method("GetMessage", static_cast<const char * (G4StateManager::*)()  const>(&G4StateManager::GetMessage));

    DEBUG_MSG("Adding wrapper for void G4StateManager::SetExceptionHandler(G4VExceptionHandler *) (" __HERE__ ")");
    // signature to use in the veto list: void G4StateManager::SetExceptionHandler(G4VExceptionHandler *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:95:15
    t.method("SetExceptionHandler", static_cast<void (G4StateManager::*)(G4VExceptionHandler *) >(&G4StateManager::SetExceptionHandler));

    DEBUG_MSG("Adding wrapper for G4VExceptionHandler * G4StateManager::GetExceptionHandler() (" __HERE__ ")");
    // signature to use in the veto list: G4VExceptionHandler * G4StateManager::GetExceptionHandler()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:96:31
    t.method("GetExceptionHandler", static_cast<G4VExceptionHandler * (G4StateManager::*)()  const>(&G4StateManager::GetExceptionHandler));

    DEBUG_MSG("Adding wrapper for void G4StateManager::SetVerboseLevel(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4StateManager::SetVerboseLevel(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4StateManager.hh:97:15
    t.method("G4StateManager!SetVerboseLevel", static_cast<void (*)(G4int) >(&G4StateManager::SetVerboseLevel));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4StateManager>> type_;
};
std::shared_ptr<Wrapper> newJlG4StateManager(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4StateManager(module));
}
