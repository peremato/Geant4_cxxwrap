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
void add_methods_for_G4ScoringManager(jlcxx::Module& types, jlcxx::TypeWrapper<G4ScoringManager>& t159) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4ScoringManager
   */

  DEBUG_MSG("Adding wrapper for G4ScoringManager * G4ScoringManager::GetScoringManager() (" __HERE__ ")");
  // signature to use in the veto list: G4ScoringManager * G4ScoringManager::GetScoringManager()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:64:28
  t159.method("G4ScoringManager!GetScoringManager", static_cast<G4ScoringManager * (*)() >(&G4ScoringManager::GetScoringManager));

  DEBUG_MSG("Adding wrapper for G4ScoringManager * G4ScoringManager::GetScoringManagerIfExist() (" __HERE__ ")");
  // signature to use in the veto list: G4ScoringManager * G4ScoringManager::GetScoringManagerIfExist()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:67:28
  t159.method("G4ScoringManager!GetScoringManagerIfExist", static_cast<G4ScoringManager * (*)() >(&G4ScoringManager::GetScoringManagerIfExist));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::SetReplicaLevel(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::SetReplicaLevel(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:72:15
  t159.method("G4ScoringManager!SetReplicaLevel", static_cast<void (*)(G4int) >(&G4ScoringManager::SetReplicaLevel));

  DEBUG_MSG("Adding wrapper for G4int G4ScoringManager::GetReplicaLevel() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4ScoringManager::GetReplicaLevel()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:73:16
  t159.method("G4ScoringManager!GetReplicaLevel", static_cast<G4int (*)() >(&G4ScoringManager::GetReplicaLevel));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::Accumulate(G4VHitsCollection *) (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::Accumulate(G4VHitsCollection *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:79:8
  t159.method("Accumulate", static_cast<void (G4ScoringManager::*)(G4VHitsCollection *) >(&G4ScoringManager::Accumulate));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::Merge(const G4ScoringManager *) (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::Merge(const G4ScoringManager *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:80:8
  t159.method("Merge", static_cast<void (G4ScoringManager::*)(const G4ScoringManager *) >(&G4ScoringManager::Merge));

  DEBUG_MSG("Adding wrapper for G4VScoringMesh * G4ScoringManager::FindMesh(G4VHitsCollection *) (" __HERE__ ")");
  // signature to use in the veto list: G4VScoringMesh * G4ScoringManager::FindMesh(G4VHitsCollection *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:81:19
  t159.method("FindMesh", static_cast<G4VScoringMesh * (G4ScoringManager::*)(G4VHitsCollection *) >(&G4ScoringManager::FindMesh));

  DEBUG_MSG("Adding wrapper for G4VScoringMesh * G4ScoringManager::FindMesh(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4VScoringMesh * G4ScoringManager::FindMesh(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:82:19
  t159.method("FindMesh", static_cast<G4VScoringMesh * (G4ScoringManager::*)(const G4String &) >(&G4ScoringManager::FindMesh));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::List() (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::List()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:83:8
  t159.method("List", static_cast<void (G4ScoringManager::*)()  const>(&G4ScoringManager::List));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::Dump() (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::Dump()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:84:8
  t159.method("Dump", static_cast<void (G4ScoringManager::*)()  const>(&G4ScoringManager::Dump));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::DrawMesh(const G4String &, const G4String &, const G4String &, G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::DrawMesh(const G4String &, const G4String &, const G4String &, G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:85:8
  t159.method("DrawMesh", static_cast<void (G4ScoringManager::*)(const G4String &, const G4String &, const G4String &, G4int) >(&G4ScoringManager::DrawMesh));
  t159.method("DrawMesh", [](G4ScoringManager& a, const G4String & arg0, const G4String & arg1, const G4String & arg2)->void{ a.DrawMesh(arg0, arg1, arg2); });
  t159.method("DrawMesh", [](G4ScoringManager* a, const G4String & arg0, const G4String & arg1, const G4String & arg2)->void{ a->DrawMesh(arg0, arg1, arg2); });

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::DrawMesh(const G4String &, const G4String &, G4int, G4int, const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::DrawMesh(const G4String &, const G4String &, G4int, G4int, const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:87:8
  t159.method("DrawMesh", static_cast<void (G4ScoringManager::*)(const G4String &, const G4String &, G4int, G4int, const G4String &) >(&G4ScoringManager::DrawMesh));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::DumpQuantityToFile(const G4String &, const G4String &, const G4String &, const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::DumpQuantityToFile(const G4String &, const G4String &, const G4String &, const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:89:8
  t159.method("DumpQuantityToFile", static_cast<void (G4ScoringManager::*)(const G4String &, const G4String &, const G4String &, const G4String &) >(&G4ScoringManager::DumpQuantityToFile));
  t159.method("DumpQuantityToFile", [](G4ScoringManager& a, const G4String & arg0, const G4String & arg1, const G4String & arg2)->void{ a.DumpQuantityToFile(arg0, arg1, arg2); });
  t159.method("DumpQuantityToFile", [](G4ScoringManager* a, const G4String & arg0, const G4String & arg1, const G4String & arg2)->void{ a->DumpQuantityToFile(arg0, arg1, arg2); });

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::DumpAllQuantitiesToFile(const G4String &, const G4String &, const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::DumpAllQuantitiesToFile(const G4String &, const G4String &, const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:92:8
  t159.method("DumpAllQuantitiesToFile", static_cast<void (G4ScoringManager::*)(const G4String &, const G4String &, const G4String &) >(&G4ScoringManager::DumpAllQuantitiesToFile));
  t159.method("DumpAllQuantitiesToFile", [](G4ScoringManager& a, const G4String & arg0, const G4String & arg1)->void{ a.DumpAllQuantitiesToFile(arg0, arg1); });
  t159.method("DumpAllQuantitiesToFile", [](G4ScoringManager* a, const G4String & arg0, const G4String & arg1)->void{ a->DumpAllQuantitiesToFile(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::ListScoreColorMaps() (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::ListScoreColorMaps()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:96:8
  t159.method("ListScoreColorMaps", static_cast<void (G4ScoringManager::*)() >(&G4ScoringManager::ListScoreColorMaps));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::SetCurrentMesh(G4VScoringMesh *) (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::SetCurrentMesh(G4VScoringMesh *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:98:15
  t159.method("SetCurrentMesh", static_cast<void (G4ScoringManager::*)(G4VScoringMesh *) >(&G4ScoringManager::SetCurrentMesh));

  DEBUG_MSG("Adding wrapper for G4VScoringMesh * G4ScoringManager::GetCurrentMesh() (" __HERE__ ")");
  // signature to use in the veto list: G4VScoringMesh * G4ScoringManager::GetCurrentMesh()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:99:26
  t159.method("GetCurrentMesh", static_cast<G4VScoringMesh * (G4ScoringManager::*)()  const>(&G4ScoringManager::GetCurrentMesh));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::CloseCurrentMesh() (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::CloseCurrentMesh()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:100:15
  t159.method("CloseCurrentMesh", static_cast<void (G4ScoringManager::*)() >(&G4ScoringManager::CloseCurrentMesh));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::SetVerboseLevel(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::SetVerboseLevel(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:101:15
  t159.method("SetVerboseLevel", static_cast<void (G4ScoringManager::*)(G4int) >(&G4ScoringManager::SetVerboseLevel));

  DEBUG_MSG("Adding wrapper for G4int G4ScoringManager::GetVerboseLevel() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4ScoringManager::GetVerboseLevel()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:111:16
  t159.method("GetVerboseLevel", static_cast<G4int (G4ScoringManager::*)()  const>(&G4ScoringManager::GetVerboseLevel));

  DEBUG_MSG("Adding wrapper for size_t G4ScoringManager::GetNumberOfMesh() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4ScoringManager::GetNumberOfMesh()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:112:17
  t159.method("GetNumberOfMesh", static_cast<size_t (G4ScoringManager::*)()  const>(&G4ScoringManager::GetNumberOfMesh));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::RegisterScoringMesh(G4VScoringMesh *) (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::RegisterScoringMesh(G4VScoringMesh *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:113:15
  t159.method("RegisterScoringMesh", static_cast<void (G4ScoringManager::*)(G4VScoringMesh *) >(&G4ScoringManager::RegisterScoringMesh));

  DEBUG_MSG("Adding wrapper for G4VScoringMesh * G4ScoringManager::GetMesh(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4VScoringMesh * G4ScoringManager::GetMesh(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:119:26
  t159.method("GetMesh", static_cast<G4VScoringMesh * (G4ScoringManager::*)(G4int)  const>(&G4ScoringManager::GetMesh));

  DEBUG_MSG("Adding wrapper for G4String G4ScoringManager::GetWorldName(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4String G4ScoringManager::GetWorldName(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:120:19
  t159.method("GetWorldName", static_cast<G4String (G4ScoringManager::*)(G4int)  const>(&G4ScoringManager::GetWorldName));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::SetScoreWriter(G4VScoreWriter *) (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::SetScoreWriter(G4VScoreWriter *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:126:15
  t159.method("SetScoreWriter", static_cast<void (G4ScoringManager::*)(G4VScoreWriter *) >(&G4ScoringManager::SetScoreWriter));

  DEBUG_MSG("Adding wrapper for void G4ScoringManager::SetFactor(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4ScoringManager::SetFactor(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:136:15
  t159.method("SetFactor", static_cast<void (G4ScoringManager::*)(G4double) >(&G4ScoringManager::SetFactor));
  t159.method("SetFactor", [](G4ScoringManager& a)->void{ a.SetFactor(); });
  t159.method("SetFactor", [](G4ScoringManager* a)->void{ a->SetFactor(); });

  DEBUG_MSG("Adding wrapper for G4double G4ScoringManager::GetFactor() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4ScoringManager::GetFactor()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ScoringManager.hh:141:19
  t159.method("GetFactor", static_cast<G4double (G4ScoringManager::*)()  const>(&G4ScoringManager::GetFactor));

  /* End of G4ScoringManager class method wrappers
   **********************************************************************/

}
