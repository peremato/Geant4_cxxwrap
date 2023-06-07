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
void add_methods_for_G4Run(jlcxx::Module& types, jlcxx::TypeWrapper<G4Run>& t114) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4Run
   */

  DEBUG_MSG("Adding wrapper for void G4Run::RecordEvent(const G4Event *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Run::RecordEvent(const G4Event *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Run.hh:61:18
  t114.method("RecordEvent", static_cast<void (G4Run::*)(const G4Event *) >(&G4Run::RecordEvent));

  DEBUG_MSG("Adding wrapper for void G4Run::Merge(const G4Run *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Run::Merge(const G4Run *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Run.hh:66:18
  t114.method("Merge", static_cast<void (G4Run::*)(const G4Run *) >(&G4Run::Merge));

  DEBUG_MSG("Adding wrapper for void G4Run::StoreEvent(G4Event *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Run::StoreEvent(G4Event *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Run.hh:69:10
  t114.method("StoreEvent", static_cast<void (G4Run::*)(G4Event *) >(&G4Run::StoreEvent));

  DEBUG_MSG("Adding wrapper for G4int G4Run::GetRunID() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4Run::GetRunID()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Run.hh:78:18
  t114.method("GetRunID", static_cast<G4int (G4Run::*)()  const>(&G4Run::GetRunID));

  DEBUG_MSG("Adding wrapper for G4int G4Run::GetNumberOfEvent() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4Run::GetNumberOfEvent()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Run.hh:80:18
  t114.method("GetNumberOfEvent", static_cast<G4int (G4Run::*)()  const>(&G4Run::GetNumberOfEvent));

  DEBUG_MSG("Adding wrapper for G4int G4Run::GetNumberOfEventToBeProcessed() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4Run::GetNumberOfEventToBeProcessed()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Run.hh:83:18
  t114.method("GetNumberOfEventToBeProcessed", static_cast<G4int (G4Run::*)()  const>(&G4Run::GetNumberOfEventToBeProcessed));

  DEBUG_MSG("Adding wrapper for const G4HCtable * G4Run::GetHCtable() (" __HERE__ ")");
  // signature to use in the veto list: const G4HCtable * G4Run::GetHCtable()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Run.hh:87:29
  t114.method("GetHCtable", static_cast<const G4HCtable * (G4Run::*)()  const>(&G4Run::GetHCtable));

  DEBUG_MSG("Adding wrapper for const G4String & G4Run::GetRandomNumberStatus() (" __HERE__ ")");
  // signature to use in the veto list: const G4String & G4Run::GetRandomNumberStatus()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Run.hh:91:28
  t114.method("GetRandomNumberStatus", static_cast<const G4String & (G4Run::*)()  const>(&G4Run::GetRandomNumberStatus));

  DEBUG_MSG("Adding wrapper for void G4Run::SetRunID(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4Run::SetRunID(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Run.hh:102:17
  t114.method("SetRunID", static_cast<void (G4Run::*)(G4int) >(&G4Run::SetRunID));

  DEBUG_MSG("Adding wrapper for void G4Run::SetNumberOfEventToBeProcessed(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4Run::SetNumberOfEventToBeProcessed(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Run.hh:103:17
  t114.method("SetNumberOfEventToBeProcessed", static_cast<void (G4Run::*)(G4int) >(&G4Run::SetNumberOfEventToBeProcessed));

  DEBUG_MSG("Adding wrapper for void G4Run::SetHCtable(G4HCtable *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Run::SetHCtable(G4HCtable *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Run.hh:107:17
  t114.method("SetHCtable", static_cast<void (G4Run::*)(G4HCtable *) >(&G4Run::SetHCtable));

  DEBUG_MSG("Adding wrapper for void G4Run::SetRandomNumberStatus(G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4Run::SetRandomNumberStatus(G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Run.hh:109:17
  t114.method("SetRandomNumberStatus", static_cast<void (G4Run::*)(G4String &) >(&G4Run::SetRandomNumberStatus));

  /* End of G4Run class method wrappers
   **********************************************************************/

}
