// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4TransportationManager> : std::false_type { };
  template<> struct DefaultConstructible<G4TransportationManager> : std::false_type { };
}

// Class generating the wrapper for type G4TransportationManager
// signature to use in the veto file: G4TransportationManager
struct JlG4TransportationManager: public Wrapper {

  JlG4TransportationManager(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4TransportationManager (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:53:7
    jlcxx::TypeWrapper<G4TransportationManager>  t = jlModule.add_type<G4TransportationManager>("G4TransportationManager");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4TransportationManager>>(new jlcxx::TypeWrapper<G4TransportationManager>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;

    DEBUG_MSG("Adding wrapper for G4TransportationManager * G4TransportationManager::GetTransportationManager() (" __HERE__ ")");
    // signature to use in the veto list: G4TransportationManager * G4TransportationManager::GetTransportationManager()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:57:38
    module_.method("G4TransportationManager!GetTransportationManager", static_cast<G4TransportationManager * (*)() >(&G4TransportationManager::GetTransportationManager));

    DEBUG_MSG("Adding wrapper for G4TransportationManager * G4TransportationManager::GetInstanceIfExist() (" __HERE__ ")");
    // signature to use in the veto list: G4TransportationManager * G4TransportationManager::GetInstanceIfExist()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:59:38
    module_.method("G4TransportationManager!GetInstanceIfExist", static_cast<G4TransportationManager * (*)() >(&G4TransportationManager::GetInstanceIfExist));

    DEBUG_MSG("Adding wrapper for G4FieldManager * G4TransportationManager::GetFieldManager() (" __HERE__ ")");
    // signature to use in the veto list: G4FieldManager * G4TransportationManager::GetFieldManager()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:64:29
    t.method("GetFieldManager", static_cast<G4FieldManager * (G4TransportationManager::*)()  const>(&G4TransportationManager::GetFieldManager));

    DEBUG_MSG("Adding wrapper for void G4TransportationManager::SetFieldManager(G4FieldManager *) (" __HERE__ ")");
    // signature to use in the veto list: void G4TransportationManager::SetFieldManager(G4FieldManager *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:65:11
    t.method("SetFieldManager", static_cast<void (G4TransportationManager::*)(G4FieldManager *) >(&G4TransportationManager::SetFieldManager));

    DEBUG_MSG("Adding wrapper for G4Navigator * G4TransportationManager::GetNavigatorForTracking() (" __HERE__ ")");
    // signature to use in the veto list: G4Navigator * G4TransportationManager::GetNavigatorForTracking()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:68:26
    t.method("GetNavigatorForTracking", static_cast<G4Navigator * (G4TransportationManager::*)()  const>(&G4TransportationManager::GetNavigatorForTracking));

    DEBUG_MSG("Adding wrapper for void G4TransportationManager::SetNavigatorForTracking(G4Navigator *) (" __HERE__ ")");
    // signature to use in the veto list: void G4TransportationManager::SetNavigatorForTracking(G4Navigator *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:69:11
    t.method("SetNavigatorForTracking", static_cast<void (G4TransportationManager::*)(G4Navigator *) >(&G4TransportationManager::SetNavigatorForTracking));

    DEBUG_MSG("Adding wrapper for void G4TransportationManager::SetWorldForTracking(G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4TransportationManager::SetWorldForTracking(G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:72:18
    t.method("SetWorldForTracking", static_cast<void (G4TransportationManager::*)(G4VPhysicalVolume *) >(&G4TransportationManager::SetWorldForTracking));

    DEBUG_MSG("Adding wrapper for G4SafetyHelper * G4TransportationManager::GetSafetyHelper() (" __HERE__ ")");
    // signature to use in the veto list: G4SafetyHelper * G4TransportationManager::GetSafetyHelper()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:84:29
    t.method("GetSafetyHelper", static_cast<G4SafetyHelper * (G4TransportationManager::*)()  const>(&G4TransportationManager::GetSafetyHelper));

    DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4TransportationManager::GetParallelWorld(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4VPhysicalVolume * G4TransportationManager::GetParallelWorld(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:87:25
    t.method("GetParallelWorld", static_cast<G4VPhysicalVolume * (G4TransportationManager::*)(const G4String &) >(&G4TransportationManager::GetParallelWorld));

    DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4TransportationManager::IsWorldExisting(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4VPhysicalVolume * G4TransportationManager::IsWorldExisting(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:91:25
    t.method("IsWorldExisting", static_cast<G4VPhysicalVolume * (G4TransportationManager::*)(const G4String &) >(&G4TransportationManager::IsWorldExisting));

    DEBUG_MSG("Adding wrapper for G4Navigator * G4TransportationManager::GetNavigator(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4Navigator * G4TransportationManager::GetNavigator(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:95:19
    t.method("GetNavigator", static_cast<G4Navigator * (G4TransportationManager::*)(const G4String &) >(&G4TransportationManager::GetNavigator));

    DEBUG_MSG("Adding wrapper for G4Navigator * G4TransportationManager::GetNavigator(G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: G4Navigator * G4TransportationManager::GetNavigator(G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:96:19
    t.method("GetNavigator", static_cast<G4Navigator * (G4TransportationManager::*)(G4VPhysicalVolume *) >(&G4TransportationManager::GetNavigator));

    DEBUG_MSG("Adding wrapper for G4bool G4TransportationManager::RegisterWorld(G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4TransportationManager::RegisterWorld(G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:101:13
    t.method("RegisterWorld", static_cast<G4bool (G4TransportationManager::*)(G4VPhysicalVolume *) >(&G4TransportationManager::RegisterWorld));

    DEBUG_MSG("Adding wrapper for void G4TransportationManager::DeRegisterNavigator(G4Navigator *) (" __HERE__ ")");
    // signature to use in the veto list: void G4TransportationManager::DeRegisterNavigator(G4Navigator *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:102:11
    t.method("DeRegisterNavigator", static_cast<void (G4TransportationManager::*)(G4Navigator *) >(&G4TransportationManager::DeRegisterNavigator));

    DEBUG_MSG("Adding wrapper for G4int G4TransportationManager::ActivateNavigator(G4Navigator *) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4TransportationManager::ActivateNavigator(G4Navigator *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:103:13
    t.method("ActivateNavigator", static_cast<G4int (G4TransportationManager::*)(G4Navigator *) >(&G4TransportationManager::ActivateNavigator));

    DEBUG_MSG("Adding wrapper for void G4TransportationManager::DeActivateNavigator(G4Navigator *) (" __HERE__ ")");
    // signature to use in the veto list: void G4TransportationManager::DeActivateNavigator(G4Navigator *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:104:11
    t.method("DeActivateNavigator", static_cast<void (G4TransportationManager::*)(G4Navigator *) >(&G4TransportationManager::DeActivateNavigator));

    DEBUG_MSG("Adding wrapper for void G4TransportationManager::InactivateAll() (" __HERE__ ")");
    // signature to use in the veto list: void G4TransportationManager::InactivateAll()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:105:11
    t.method("InactivateAll", static_cast<void (G4TransportationManager::*)() >(&G4TransportationManager::InactivateAll));

    DEBUG_MSG("Adding wrapper for G4Navigator * G4TransportationManager::GetFirstTrackingNavigator() (" __HERE__ ")");
    // signature to use in the veto list: G4Navigator * G4TransportationManager::GetFirstTrackingNavigator()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:109:26
    module_.method("G4TransportationManager!GetFirstTrackingNavigator", static_cast<G4Navigator * (*)() >(&G4TransportationManager::GetFirstTrackingNavigator));

    DEBUG_MSG("Adding wrapper for void G4TransportationManager::SetFirstTrackingNavigator(G4Navigator *) (" __HERE__ ")");
    // signature to use in the veto list: void G4TransportationManager::SetFirstTrackingNavigator(G4Navigator *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:110:18
    module_.method("G4TransportationManager!SetFirstTrackingNavigator", static_cast<void (*)(G4Navigator *) >(&G4TransportationManager::SetFirstTrackingNavigator));

    DEBUG_MSG("Adding wrapper for void G4TransportationManager::ClearParallelWorlds() (" __HERE__ ")");
    // signature to use in the veto list: void G4TransportationManager::ClearParallelWorlds()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4TransportationManager.hh:118:11
    t.method("ClearParallelWorlds", static_cast<void (G4TransportationManager::*)() >(&G4TransportationManager::ClearParallelWorlds));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4TransportationManager>> type_;
};
std::shared_ptr<Wrapper> newJlG4TransportationManager(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4TransportationManager(module));
}
