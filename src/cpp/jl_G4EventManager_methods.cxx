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
void add_methods_for_G4EventManager(jlcxx::Module& types, jlcxx::TypeWrapper<G4EventManager>& t159) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4EventManager
   */

  DEBUG_MSG("Adding wrapper for G4EventManager * G4EventManager::GetEventManager() (" __HERE__ ")");
  // signature to use in the veto list: G4EventManager * G4EventManager::GetEventManager()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:62:28
  t159.method("G4EventManager!GetEventManager", static_cast<G4EventManager * (*)() >(&G4EventManager::GetEventManager));

  DEBUG_MSG("Adding wrapper for void G4EventManager::ProcessOneEvent(G4Event *) (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::ProcessOneEvent(G4Event *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:71:10
  t159.method("ProcessOneEvent", static_cast<void (G4EventManager::*)(G4Event *) >(&G4EventManager::ProcessOneEvent));

  DEBUG_MSG("Adding wrapper for void G4EventManager::ProcessOneEvent(G4TrackVector *, G4Event *) (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::ProcessOneEvent(G4TrackVector *, G4Event *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:74:10
  t159.method("ProcessOneEvent", static_cast<void (G4EventManager::*)(G4TrackVector *, G4Event *) >(&G4EventManager::ProcessOneEvent));
  t159.method("ProcessOneEvent", [](G4EventManager& a, G4TrackVector * arg0)->void{ a.ProcessOneEvent(arg0); });
  t159.method("ProcessOneEvent", [](G4EventManager* a, G4TrackVector * arg0)->void{ a->ProcessOneEvent(arg0); });

  DEBUG_MSG("Adding wrapper for void G4EventManager::StackTracks(G4TrackVector *, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::StackTracks(G4TrackVector *, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:91:10
  t159.method("StackTracks", static_cast<void (G4EventManager::*)(G4TrackVector *, G4bool) >(&G4EventManager::StackTracks));
  t159.method("StackTracks", [](G4EventManager& a, G4TrackVector * arg0)->void{ a.StackTracks(arg0); });
  t159.method("StackTracks", [](G4EventManager* a, G4TrackVector * arg0)->void{ a->StackTracks(arg0); });

  DEBUG_MSG("Adding wrapper for const G4Event * G4EventManager::GetConstCurrentEvent() (" __HERE__ ")");
  // signature to use in the veto list: const G4Event * G4EventManager::GetConstCurrentEvent()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:95:27
  t159.method("GetConstCurrentEvent", static_cast<const G4Event * (G4EventManager::*)() >(&G4EventManager::GetConstCurrentEvent));

  DEBUG_MSG("Adding wrapper for G4Event * G4EventManager::GetNonconstCurrentEvent() (" __HERE__ ")");
  // signature to use in the veto list: G4Event * G4EventManager::GetNonconstCurrentEvent()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:97:21
  t159.method("GetNonconstCurrentEvent", static_cast<G4Event * (G4EventManager::*)() >(&G4EventManager::GetNonconstCurrentEvent));

  DEBUG_MSG("Adding wrapper for void G4EventManager::AbortCurrentEvent() (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::AbortCurrentEvent()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:103:10
  t159.method("AbortCurrentEvent", static_cast<void (G4EventManager::*)() >(&G4EventManager::AbortCurrentEvent));

  DEBUG_MSG("Adding wrapper for void G4EventManager::SetUserAction(G4UserEventAction *) (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::SetUserAction(G4UserEventAction *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:109:10
  t159.method("SetUserAction", static_cast<void (G4EventManager::*)(G4UserEventAction *) >(&G4EventManager::SetUserAction));

  DEBUG_MSG("Adding wrapper for void G4EventManager::SetUserAction(G4UserStackingAction *) (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::SetUserAction(G4UserStackingAction *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:110:10
  t159.method("SetUserAction", static_cast<void (G4EventManager::*)(G4UserStackingAction *) >(&G4EventManager::SetUserAction));

  DEBUG_MSG("Adding wrapper for void G4EventManager::SetUserAction(G4UserTrackingAction *) (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::SetUserAction(G4UserTrackingAction *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:111:10
  t159.method("SetUserAction", static_cast<void (G4EventManager::*)(G4UserTrackingAction *) >(&G4EventManager::SetUserAction));

  DEBUG_MSG("Adding wrapper for void G4EventManager::SetUserAction(G4UserSteppingAction *) (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::SetUserAction(G4UserSteppingAction *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:112:10
  t159.method("SetUserAction", static_cast<void (G4EventManager::*)(G4UserSteppingAction *) >(&G4EventManager::SetUserAction));

  DEBUG_MSG("Adding wrapper for G4UserEventAction * G4EventManager::GetUserEventAction() (" __HERE__ ")");
  // signature to use in the veto list: G4UserEventAction * G4EventManager::GetUserEventAction()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:113:31
  t159.method("GetUserEventAction", static_cast<G4UserEventAction * (G4EventManager::*)() >(&G4EventManager::GetUserEventAction));

  DEBUG_MSG("Adding wrapper for G4UserStackingAction * G4EventManager::GetUserStackingAction() (" __HERE__ ")");
  // signature to use in the veto list: G4UserStackingAction * G4EventManager::GetUserStackingAction()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:115:34
  t159.method("GetUserStackingAction", static_cast<G4UserStackingAction * (G4EventManager::*)() >(&G4EventManager::GetUserStackingAction));

  DEBUG_MSG("Adding wrapper for G4UserTrackingAction * G4EventManager::GetUserTrackingAction() (" __HERE__ ")");
  // signature to use in the veto list: G4UserTrackingAction * G4EventManager::GetUserTrackingAction()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:117:34
  t159.method("GetUserTrackingAction", static_cast<G4UserTrackingAction * (G4EventManager::*)() >(&G4EventManager::GetUserTrackingAction));

  DEBUG_MSG("Adding wrapper for G4UserSteppingAction * G4EventManager::GetUserSteppingAction() (" __HERE__ ")");
  // signature to use in the veto list: G4UserSteppingAction * G4EventManager::GetUserSteppingAction()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:119:34
  t159.method("GetUserSteppingAction", static_cast<G4UserSteppingAction * (G4EventManager::*)() >(&G4EventManager::GetUserSteppingAction));

  DEBUG_MSG("Adding wrapper for void G4EventManager::SetNumberOfAdditionalWaitingStacks(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::SetNumberOfAdditionalWaitingStacks(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:125:10
  t159.method("SetNumberOfAdditionalWaitingStacks", static_cast<void (G4EventManager::*)(G4int) >(&G4EventManager::SetNumberOfAdditionalWaitingStacks));

  DEBUG_MSG("Adding wrapper for void G4EventManager::KeepTheCurrentEvent() (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::KeepTheCurrentEvent()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:128:10
  t159.method("KeepTheCurrentEvent", static_cast<void (G4EventManager::*)() >(&G4EventManager::KeepTheCurrentEvent));

  DEBUG_MSG("Adding wrapper for G4StackManager * G4EventManager::GetStackManager() (" __HERE__ ")");
  // signature to use in the veto list: G4StackManager * G4EventManager::GetStackManager()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:132:28
  t159.method("GetStackManager", static_cast<G4StackManager * (G4EventManager::*)()  const>(&G4EventManager::GetStackManager));

  DEBUG_MSG("Adding wrapper for G4TrackingManager * G4EventManager::GetTrackingManager() (" __HERE__ ")");
  // signature to use in the veto list: G4TrackingManager * G4EventManager::GetTrackingManager()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:134:31
  t159.method("GetTrackingManager", static_cast<G4TrackingManager * (G4EventManager::*)()  const>(&G4EventManager::GetTrackingManager));

  DEBUG_MSG("Adding wrapper for G4int G4EventManager::GetVerboseLevel() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4EventManager::GetVerboseLevel()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:137:18
  t159.method("GetVerboseLevel", static_cast<G4int (G4EventManager::*)() >(&G4EventManager::GetVerboseLevel));

  DEBUG_MSG("Adding wrapper for void G4EventManager::SetVerboseLevel(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::SetVerboseLevel(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:139:17
  t159.method("SetVerboseLevel", static_cast<void (G4EventManager::*)(G4int) >(&G4EventManager::SetVerboseLevel));

  DEBUG_MSG("Adding wrapper for void G4EventManager::SetUserInformation(G4VUserEventInformation *) (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::SetUserInformation(G4VUserEventInformation *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:147:10
  t159.method("SetUserInformation", static_cast<void (G4EventManager::*)(G4VUserEventInformation *) >(&G4EventManager::SetUserInformation));

  DEBUG_MSG("Adding wrapper for G4VUserEventInformation * G4EventManager::GetUserInformation() (" __HERE__ ")");
  // signature to use in the veto list: G4VUserEventInformation * G4EventManager::GetUserInformation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:148:30
  t159.method("GetUserInformation", static_cast<G4VUserEventInformation * (G4EventManager::*)() >(&G4EventManager::GetUserInformation));

  DEBUG_MSG("Adding wrapper for G4PrimaryTransformer * G4EventManager::GetPrimaryTransformer() (" __HERE__ ")");
  // signature to use in the veto list: G4PrimaryTransformer * G4EventManager::GetPrimaryTransformer()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:153:34
  t159.method("GetPrimaryTransformer", static_cast<G4PrimaryTransformer * (G4EventManager::*)()  const>(&G4EventManager::GetPrimaryTransformer));

  DEBUG_MSG("Adding wrapper for void G4EventManager::SetPrimaryTransformer(G4PrimaryTransformer *) (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::SetPrimaryTransformer(G4PrimaryTransformer *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:155:17
  t159.method("SetPrimaryTransformer", static_cast<void (G4EventManager::*)(G4PrimaryTransformer *) >(&G4EventManager::SetPrimaryTransformer));

  DEBUG_MSG("Adding wrapper for void G4EventManager::StoreRandomNumberStatusToG4Event(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4EventManager::StoreRandomNumberStatusToG4Event(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EventManager.hh:157:17
  t159.method("StoreRandomNumberStatusToG4Event", static_cast<void (G4EventManager::*)(G4int) >(&G4EventManager::StoreRandomNumberStatusToG4Event));

  /* End of G4EventManager class method wrappers
   **********************************************************************/

}