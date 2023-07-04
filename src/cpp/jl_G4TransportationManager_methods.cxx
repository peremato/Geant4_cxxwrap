// this file was auto-generated by wrapit v0.1.0-54-g4322429
#include <type_traits>
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

#include "cpp/jlGeant4.h"


#ifdef VERBOSE_IMPORT
#  define DEBUG_MSG(a) std::cerr << a << "\n"
#else
#  define DEBUG_MSG(a)
#endif
#define __HERE__  __FILE__ ":" QUOTE2(__LINE__)
#define QUOTE(arg) #arg
#define QUOTE2(arg) QUOTE(arg)
void add_methods_for_G4TransportationManager(jlcxx::Module& types, jlcxx::TypeWrapper<G4TransportationManager>& t154) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4TransportationManager
   */

  DEBUG_MSG("Adding wrapper for G4TransportationManager * G4TransportationManager::GetTransportationManager() (" __HERE__ ")");
  // signature to use in the veto list: G4TransportationManager * G4TransportationManager::GetTransportationManager()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:57:38
  t154.method("G4TransportationManager!GetTransportationManager", static_cast<G4TransportationManager * (*)() >(&G4TransportationManager::GetTransportationManager));

  DEBUG_MSG("Adding wrapper for G4TransportationManager * G4TransportationManager::GetInstanceIfExist() (" __HERE__ ")");
  // signature to use in the veto list: G4TransportationManager * G4TransportationManager::GetInstanceIfExist()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:59:38
  t154.method("G4TransportationManager!GetInstanceIfExist", static_cast<G4TransportationManager * (*)() >(&G4TransportationManager::GetInstanceIfExist));

  DEBUG_MSG("Adding wrapper for G4FieldManager * G4TransportationManager::GetFieldManager() (" __HERE__ ")");
  // signature to use in the veto list: G4FieldManager * G4TransportationManager::GetFieldManager()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:64:29
  t154.method("GetFieldManager", static_cast<G4FieldManager * (G4TransportationManager::*)()  const>(&G4TransportationManager::GetFieldManager));

  DEBUG_MSG("Adding wrapper for void G4TransportationManager::SetFieldManager(G4FieldManager *) (" __HERE__ ")");
  // signature to use in the veto list: void G4TransportationManager::SetFieldManager(G4FieldManager *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:65:11
  t154.method("SetFieldManager", static_cast<void (G4TransportationManager::*)(G4FieldManager *) >(&G4TransportationManager::SetFieldManager));

  DEBUG_MSG("Adding wrapper for G4Navigator * G4TransportationManager::GetNavigatorForTracking() (" __HERE__ ")");
  // signature to use in the veto list: G4Navigator * G4TransportationManager::GetNavigatorForTracking()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:68:26
  t154.method("GetNavigatorForTracking", static_cast<G4Navigator * (G4TransportationManager::*)()  const>(&G4TransportationManager::GetNavigatorForTracking));

  DEBUG_MSG("Adding wrapper for void G4TransportationManager::SetNavigatorForTracking(G4Navigator *) (" __HERE__ ")");
  // signature to use in the veto list: void G4TransportationManager::SetNavigatorForTracking(G4Navigator *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:69:11
  t154.method("SetNavigatorForTracking", static_cast<void (G4TransportationManager::*)(G4Navigator *) >(&G4TransportationManager::SetNavigatorForTracking));

  DEBUG_MSG("Adding wrapper for void G4TransportationManager::SetWorldForTracking(G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4TransportationManager::SetWorldForTracking(G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:72:18
  t154.method("SetWorldForTracking", static_cast<void (G4TransportationManager::*)(G4VPhysicalVolume *) >(&G4TransportationManager::SetWorldForTracking));

  DEBUG_MSG("Adding wrapper for size_t G4TransportationManager::GetNoActiveNavigators() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4TransportationManager::GetNoActiveNavigators()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:76:25
  t154.method("GetNoActiveNavigators", static_cast<size_t (G4TransportationManager::*)()  const>(&G4TransportationManager::GetNoActiveNavigators));

  DEBUG_MSG("Adding wrapper for size_t G4TransportationManager::GetNoWorlds() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4TransportationManager::GetNoWorlds()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:80:25
  t154.method("GetNoWorlds", static_cast<size_t (G4TransportationManager::*)()  const>(&G4TransportationManager::GetNoWorlds));

  DEBUG_MSG("Adding wrapper for G4SafetyHelper * G4TransportationManager::GetSafetyHelper() (" __HERE__ ")");
  // signature to use in the veto list: G4SafetyHelper * G4TransportationManager::GetSafetyHelper()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:84:29
  t154.method("GetSafetyHelper", static_cast<G4SafetyHelper * (G4TransportationManager::*)()  const>(&G4TransportationManager::GetSafetyHelper));

  DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4TransportationManager::GetParallelWorld(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4VPhysicalVolume * G4TransportationManager::GetParallelWorld(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:87:25
  t154.method("GetParallelWorld", static_cast<G4VPhysicalVolume * (G4TransportationManager::*)(const G4String &) >(&G4TransportationManager::GetParallelWorld));

  DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4TransportationManager::IsWorldExisting(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4VPhysicalVolume * G4TransportationManager::IsWorldExisting(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:91:25
  t154.method("IsWorldExisting", static_cast<G4VPhysicalVolume * (G4TransportationManager::*)(const G4String &) >(&G4TransportationManager::IsWorldExisting));

  DEBUG_MSG("Adding wrapper for G4Navigator * G4TransportationManager::GetNavigator(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4Navigator * G4TransportationManager::GetNavigator(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:95:19
  t154.method("GetNavigator", static_cast<G4Navigator * (G4TransportationManager::*)(const G4String &) >(&G4TransportationManager::GetNavigator));

  DEBUG_MSG("Adding wrapper for G4Navigator * G4TransportationManager::GetNavigator(G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: G4Navigator * G4TransportationManager::GetNavigator(G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:96:19
  t154.method("GetNavigator", static_cast<G4Navigator * (G4TransportationManager::*)(G4VPhysicalVolume *) >(&G4TransportationManager::GetNavigator));

  DEBUG_MSG("Adding wrapper for G4bool G4TransportationManager::RegisterWorld(G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4TransportationManager::RegisterWorld(G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:101:13
  t154.method("RegisterWorld", static_cast<G4bool (G4TransportationManager::*)(G4VPhysicalVolume *) >(&G4TransportationManager::RegisterWorld));

  DEBUG_MSG("Adding wrapper for void G4TransportationManager::DeRegisterNavigator(G4Navigator *) (" __HERE__ ")");
  // signature to use in the veto list: void G4TransportationManager::DeRegisterNavigator(G4Navigator *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:102:11
  t154.method("DeRegisterNavigator", static_cast<void (G4TransportationManager::*)(G4Navigator *) >(&G4TransportationManager::DeRegisterNavigator));

  DEBUG_MSG("Adding wrapper for G4int G4TransportationManager::ActivateNavigator(G4Navigator *) (" __HERE__ ")");
  // signature to use in the veto list: G4int G4TransportationManager::ActivateNavigator(G4Navigator *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:103:13
  t154.method("ActivateNavigator", static_cast<G4int (G4TransportationManager::*)(G4Navigator *) >(&G4TransportationManager::ActivateNavigator));

  DEBUG_MSG("Adding wrapper for void G4TransportationManager::DeActivateNavigator(G4Navigator *) (" __HERE__ ")");
  // signature to use in the veto list: void G4TransportationManager::DeActivateNavigator(G4Navigator *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:104:11
  t154.method("DeActivateNavigator", static_cast<void (G4TransportationManager::*)(G4Navigator *) >(&G4TransportationManager::DeActivateNavigator));

  DEBUG_MSG("Adding wrapper for void G4TransportationManager::InactivateAll() (" __HERE__ ")");
  // signature to use in the veto list: void G4TransportationManager::InactivateAll()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:105:11
  t154.method("InactivateAll", static_cast<void (G4TransportationManager::*)() >(&G4TransportationManager::InactivateAll));

  DEBUG_MSG("Adding wrapper for G4Navigator * G4TransportationManager::GetFirstTrackingNavigator() (" __HERE__ ")");
  // signature to use in the veto list: G4Navigator * G4TransportationManager::GetFirstTrackingNavigator()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:109:26
  t154.method("G4TransportationManager!GetFirstTrackingNavigator", static_cast<G4Navigator * (*)() >(&G4TransportationManager::GetFirstTrackingNavigator));

  DEBUG_MSG("Adding wrapper for void G4TransportationManager::SetFirstTrackingNavigator(G4Navigator *) (" __HERE__ ")");
  // signature to use in the veto list: void G4TransportationManager::SetFirstTrackingNavigator(G4Navigator *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:110:18
  t154.method("G4TransportationManager!SetFirstTrackingNavigator", static_cast<void (*)(G4Navigator *) >(&G4TransportationManager::SetFirstTrackingNavigator));

  DEBUG_MSG("Adding wrapper for void G4TransportationManager::ClearParallelWorlds() (" __HERE__ ")");
  // signature to use in the veto list: void G4TransportationManager::ClearParallelWorlds()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TransportationManager.hh:118:11
  t154.method("ClearParallelWorlds", static_cast<void (G4TransportationManager::*)() >(&G4TransportationManager::ClearParallelWorlds));

  /* End of G4TransportationManager class method wrappers
   **********************************************************************/

}