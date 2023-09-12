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
void add_methods_for_G4VModularPhysicsList(jlcxx::Module& types, jlcxx::TypeWrapper<G4VModularPhysicsList>& t210) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4VModularPhysicsList
   */

  DEBUG_MSG("Adding wrapper for void G4VModularPhysicsList::ConstructParticle() (" __HERE__ ")");
  // signature to use in the veto list: void G4VModularPhysicsList::ConstructParticle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:93:18
  t210.method("ConstructParticle", static_cast<void (G4VModularPhysicsList::*)() >(&G4VModularPhysicsList::ConstructParticle));

  DEBUG_MSG("Adding wrapper for void G4VModularPhysicsList::ConstructProcess() (" __HERE__ ")");
  // signature to use in the veto list: void G4VModularPhysicsList::ConstructProcess()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:97:18
  t210.method("ConstructProcess", static_cast<void (G4VModularPhysicsList::*)() >(&G4VModularPhysicsList::ConstructProcess));

  DEBUG_MSG("Adding wrapper for void G4VModularPhysicsList::RegisterPhysics(G4VPhysicsConstructor *) (" __HERE__ ")");
  // signature to use in the veto list: void G4VModularPhysicsList::RegisterPhysics(G4VPhysicsConstructor *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:102:10
  t210.method("RegisterPhysics", static_cast<void (G4VModularPhysicsList::*)(G4VPhysicsConstructor *) >(&G4VModularPhysicsList::RegisterPhysics));

  DEBUG_MSG("Adding wrapper for const G4VPhysicsConstructor * G4VModularPhysicsList::GetPhysics(G4int) (" __HERE__ ")");
  // signature to use in the veto list: const G4VPhysicsConstructor * G4VModularPhysicsList::GetPhysics(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:105:34
  t210.method("GetPhysics", static_cast<const G4VPhysicsConstructor * (G4VModularPhysicsList::*)(G4int)  const>(&G4VModularPhysicsList::GetPhysics));

  DEBUG_MSG("Adding wrapper for const G4VPhysicsConstructor * G4VModularPhysicsList::GetPhysics(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: const G4VPhysicsConstructor * G4VModularPhysicsList::GetPhysics(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:106:34
  t210.method("GetPhysics", static_cast<const G4VPhysicsConstructor * (G4VModularPhysicsList::*)(const G4String &)  const>(&G4VModularPhysicsList::GetPhysics));

  DEBUG_MSG("Adding wrapper for const G4VPhysicsConstructor * G4VModularPhysicsList::GetPhysicsWithType(G4int) (" __HERE__ ")");
  // signature to use in the veto list: const G4VPhysicsConstructor * G4VModularPhysicsList::GetPhysicsWithType(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:107:34
  t210.method("GetPhysicsWithType", static_cast<const G4VPhysicsConstructor * (G4VModularPhysicsList::*)(G4int)  const>(&G4VModularPhysicsList::GetPhysicsWithType));

  DEBUG_MSG("Adding wrapper for void G4VModularPhysicsList::ReplacePhysics(G4VPhysicsConstructor *) (" __HERE__ ")");
  // signature to use in the veto list: void G4VModularPhysicsList::ReplacePhysics(G4VPhysicsConstructor *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:109:10
  t210.method("ReplacePhysics", static_cast<void (G4VModularPhysicsList::*)(G4VPhysicsConstructor *) >(&G4VModularPhysicsList::ReplacePhysics));

  DEBUG_MSG("Adding wrapper for void G4VModularPhysicsList::RemovePhysics(G4VPhysicsConstructor *) (" __HERE__ ")");
  // signature to use in the veto list: void G4VModularPhysicsList::RemovePhysics(G4VPhysicsConstructor *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:116:10
  t210.method("RemovePhysics", static_cast<void (G4VModularPhysicsList::*)(G4VPhysicsConstructor *) >(&G4VModularPhysicsList::RemovePhysics));

  DEBUG_MSG("Adding wrapper for void G4VModularPhysicsList::RemovePhysics(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4VModularPhysicsList::RemovePhysics(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:117:10
  t210.method("RemovePhysics", static_cast<void (G4VModularPhysicsList::*)(G4int) >(&G4VModularPhysicsList::RemovePhysics));

  DEBUG_MSG("Adding wrapper for void G4VModularPhysicsList::RemovePhysics(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4VModularPhysicsList::RemovePhysics(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:118:10
  t210.method("RemovePhysics", static_cast<void (G4VModularPhysicsList::*)(const G4String &) >(&G4VModularPhysicsList::RemovePhysics));

  DEBUG_MSG("Adding wrapper for G4int G4VModularPhysicsList::GetInstanceID() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VModularPhysicsList::GetInstanceID()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:121:18
  t210.method("GetInstanceID", static_cast<G4int (G4VModularPhysicsList::*)()  const>(&G4VModularPhysicsList::GetInstanceID));

  DEBUG_MSG("Adding wrapper for void G4VModularPhysicsList::TerminateWorker() (" __HERE__ ")");
  // signature to use in the veto list: void G4VModularPhysicsList::TerminateWorker()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:123:18
  t210.method("TerminateWorker", static_cast<void (G4VModularPhysicsList::*)() >(&G4VModularPhysicsList::TerminateWorker));

  DEBUG_MSG("Adding wrapper for void G4VModularPhysicsList::SetVerboseLevel(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4VModularPhysicsList::SetVerboseLevel(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:125:10
  t210.method("SetVerboseLevel", static_cast<void (G4VModularPhysicsList::*)(G4int) >(&G4VModularPhysicsList::SetVerboseLevel));

  DEBUG_MSG("Adding wrapper for G4int G4VModularPhysicsList::GetVerboseLevel() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VModularPhysicsList::GetVerboseLevel()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VModularPhysicsList.hh:126:11
  t210.method("GetVerboseLevel", static_cast<G4int (G4VModularPhysicsList::*)()  const>(&G4VModularPhysicsList::GetVerboseLevel));

  /* End of G4VModularPhysicsList class method wrappers
   **********************************************************************/

}
