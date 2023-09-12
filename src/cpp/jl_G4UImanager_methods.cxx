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
void add_methods_for_G4UImanager(jlcxx::Module& types, jlcxx::TypeWrapper<G4UImanager>& t156) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4UImanager
   */

  DEBUG_MSG("Adding wrapper for G4UImanager * G4UImanager::GetUIpointer() (" __HERE__ ")");
  // signature to use in the veto list: G4UImanager * G4UImanager::GetUIpointer()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:61:25
  t156.method("G4UImanager!GetUIpointer", static_cast<G4UImanager * (*)() >(&G4UImanager::GetUIpointer));

  DEBUG_MSG("Adding wrapper for G4UImanager * G4UImanager::GetMasterUIpointer() (" __HERE__ ")");
  // signature to use in the veto list: G4UImanager * G4UImanager::GetMasterUIpointer()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:62:25
  t156.method("G4UImanager!GetMasterUIpointer", static_cast<G4UImanager * (*)() >(&G4UImanager::GetMasterUIpointer));

  DEBUG_MSG("Adding wrapper for G4String G4UImanager::GetCurrentValues(const char *) (" __HERE__ ")");
  // signature to use in the veto list: G4String G4UImanager::GetCurrentValues(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:73:14
  t156.method("GetCurrentValues", static_cast<G4String (G4UImanager::*)(const char *) >(&G4UImanager::GetCurrentValues));

  DEBUG_MSG("Adding wrapper for void G4UImanager::AddNewCommand(G4UIcommand *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::AddNewCommand(G4UIcommand *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:79:10
  t156.method("AddNewCommand", static_cast<void (G4UImanager::*)(G4UIcommand *) >(&G4UImanager::AddNewCommand));

  DEBUG_MSG("Adding wrapper for void G4UImanager::RemoveCommand(G4UIcommand *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::RemoveCommand(G4UIcommand *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:82:10
  t156.method("RemoveCommand", static_cast<void (G4UImanager::*)(G4UIcommand *) >(&G4UImanager::RemoveCommand));

  DEBUG_MSG("Adding wrapper for void G4UImanager::ExecuteMacroFile(const char *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::ExecuteMacroFile(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:86:10
  t156.method("ExecuteMacroFile", static_cast<void (G4UImanager::*)(const char *) >(&G4UImanager::ExecuteMacroFile));

  DEBUG_MSG("Adding wrapper for void G4UImanager::Loop(const char *, const char *, G4double, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::Loop(const char *, const char *, G4double, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:89:10
  t156.method("Loop", static_cast<void (G4UImanager::*)(const char *, const char *, G4double, G4double, G4double) >(&G4UImanager::Loop));
  t156.method("Loop", [](G4UImanager& a, const char * arg0, const char * arg1, G4double arg2, G4double arg3)->void{ a.Loop(arg0, arg1, arg2, arg3); });
  t156.method("Loop", [](G4UImanager* a, const char * arg0, const char * arg1, G4double arg2, G4double arg3)->void{ a->Loop(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for void G4UImanager::Foreach(const char *, const char *, const char *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::Foreach(const char *, const char *, const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:94:10
  t156.method("Foreach", static_cast<void (G4UImanager::*)(const char *, const char *, const char *) >(&G4UImanager::Foreach));

  DEBUG_MSG("Adding wrapper for G4int G4UImanager::ApplyCommand(const char *) (" __HERE__ ")");
  // signature to use in the veto list: G4int G4UImanager::ApplyCommand(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:99:11
  t156.method("ApplyCommand", static_cast<G4int (G4UImanager::*)(const char *) >(&G4UImanager::ApplyCommand));

  DEBUG_MSG("Adding wrapper for G4UIcommand * G4UImanager::FindCommand(const char *) (" __HERE__ ")");
  // signature to use in the veto list: G4UIcommand * G4UImanager::FindCommand(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:110:18
  t156.method("FindCommand", static_cast<G4UIcommand * (G4UImanager::*)(const char *) >(&G4UImanager::FindCommand));

  DEBUG_MSG("Adding wrapper for G4UIcommand * G4UImanager::FindCommand(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4UIcommand * G4UImanager::FindCommand(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:111:18
  t156.method("FindCommand", static_cast<G4UIcommand * (G4UImanager::*)(const G4String &) >(&G4UImanager::FindCommand));

  DEBUG_MSG("Adding wrapper for void G4UImanager::StoreHistory(const char *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::StoreHistory(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:117:10
  t156.method("StoreHistory", static_cast<void (G4UImanager::*)(const char *) >(&G4UImanager::StoreHistory));
  t156.method("StoreHistory", [](G4UImanager& a)->void{ a.StoreHistory(); });
  t156.method("StoreHistory", [](G4UImanager* a)->void{ a->StoreHistory(); });

  DEBUG_MSG("Adding wrapper for void G4UImanager::StoreHistory(G4bool, const char *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::StoreHistory(G4bool, const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:118:10
  t156.method("StoreHistory", static_cast<void (G4UImanager::*)(G4bool, const char *) >(&G4UImanager::StoreHistory));
  t156.method("StoreHistory", [](G4UImanager& a, G4bool arg0)->void{ a.StoreHistory(arg0); });
  t156.method("StoreHistory", [](G4UImanager* a, G4bool arg0)->void{ a->StoreHistory(arg0); });

  DEBUG_MSG("Adding wrapper for void G4UImanager::ListCommands(const char *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::ListCommands(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:123:10
  t156.method("ListCommands", static_cast<void (G4UImanager::*)(const char *) >(&G4UImanager::ListCommands));

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetAlias(const char *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetAlias(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:127:10
  t156.method("SetAlias", static_cast<void (G4UImanager::*)(const char *) >(&G4UImanager::SetAlias));

  DEBUG_MSG("Adding wrapper for void G4UImanager::RemoveAlias(const char *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::RemoveAlias(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:132:10
  t156.method("RemoveAlias", static_cast<void (G4UImanager::*)(const char *) >(&G4UImanager::RemoveAlias));

  DEBUG_MSG("Adding wrapper for void G4UImanager::ListAlias() (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::ListAlias()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:135:10
  t156.method("ListAlias", static_cast<void (G4UImanager::*)() >(&G4UImanager::ListAlias));

  DEBUG_MSG("Adding wrapper for G4String G4UImanager::SolveAlias(const char *) (" __HERE__ ")");
  // signature to use in the veto list: G4String G4UImanager::SolveAlias(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:138:14
  t156.method("SolveAlias", static_cast<G4String (G4UImanager::*)(const char *) >(&G4UImanager::SolveAlias));

  DEBUG_MSG("Adding wrapper for void G4UImanager::CreateHTML(const char *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::CreateHTML(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:141:10
  t156.method("CreateHTML", static_cast<void (G4UImanager::*)(const char *) >(&G4UImanager::CreateHTML));
  t156.method("CreateHTML", [](G4UImanager& a)->void{ a.CreateHTML(); });
  t156.method("CreateHTML", [](G4UImanager* a)->void{ a->CreateHTML(); });

  DEBUG_MSG("Adding wrapper for void G4UImanager::LoopS(const char *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::LoopS(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:144:10
  t156.method("LoopS", static_cast<void (G4UImanager::*)(const char *) >(&G4UImanager::LoopS));

  DEBUG_MSG("Adding wrapper for void G4UImanager::ForeachS(const char *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::ForeachS(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:145:10
  t156.method("ForeachS", static_cast<void (G4UImanager::*)(const char *) >(&G4UImanager::ForeachS));

  DEBUG_MSG("Adding wrapper for G4bool G4UImanager::Notify(G4ApplicationState) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4UImanager::Notify(G4ApplicationState)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:149:12
  t156.method("Notify", static_cast<G4bool (G4UImanager::*)(G4ApplicationState) >(&G4UImanager::Notify));

  DEBUG_MSG("Adding wrapper for G4String G4UImanager::GetCurrentStringValue(const char *, G4int, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4String G4UImanager::GetCurrentStringValue(const char *, G4int, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:152:14
  t156.method("GetCurrentStringValue", static_cast<G4String (G4UImanager::*)(const char *, G4int, G4bool) >(&G4UImanager::GetCurrentStringValue));
  t156.method("GetCurrentStringValue", [](G4UImanager& a, const char * arg0)->G4String{ return a.GetCurrentStringValue(arg0); });
  t156.method("GetCurrentStringValue", [](G4UImanager& a, const char * arg0, G4int arg1)->G4String{ return a.GetCurrentStringValue(arg0, arg1); });
  t156.method("GetCurrentStringValue", [](G4UImanager* a, const char * arg0)->G4String{ return a->GetCurrentStringValue(arg0); });
  t156.method("GetCurrentStringValue", [](G4UImanager* a, const char * arg0, G4int arg1)->G4String{ return a->GetCurrentStringValue(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for G4int G4UImanager::GetCurrentIntValue(const char *, G4int, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4int G4UImanager::GetCurrentIntValue(const char *, G4int, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:155:11
  t156.method("GetCurrentIntValue", static_cast<G4int (G4UImanager::*)(const char *, G4int, G4bool) >(&G4UImanager::GetCurrentIntValue));
  t156.method("GetCurrentIntValue", [](G4UImanager& a, const char * arg0)->G4int{ return a.GetCurrentIntValue(arg0); });
  t156.method("GetCurrentIntValue", [](G4UImanager& a, const char * arg0, G4int arg1)->G4int{ return a.GetCurrentIntValue(arg0, arg1); });
  t156.method("GetCurrentIntValue", [](G4UImanager* a, const char * arg0)->G4int{ return a->GetCurrentIntValue(arg0); });
  t156.method("GetCurrentIntValue", [](G4UImanager* a, const char * arg0, G4int arg1)->G4int{ return a->GetCurrentIntValue(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for G4double G4UImanager::GetCurrentDoubleValue(const char *, G4int, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4UImanager::GetCurrentDoubleValue(const char *, G4int, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:157:14
  t156.method("GetCurrentDoubleValue", static_cast<G4double (G4UImanager::*)(const char *, G4int, G4bool) >(&G4UImanager::GetCurrentDoubleValue));
  t156.method("GetCurrentDoubleValue", [](G4UImanager& a, const char * arg0)->G4double{ return a.GetCurrentDoubleValue(arg0); });
  t156.method("GetCurrentDoubleValue", [](G4UImanager& a, const char * arg0, G4int arg1)->G4double{ return a.GetCurrentDoubleValue(arg0, arg1); });
  t156.method("GetCurrentDoubleValue", [](G4UImanager* a, const char * arg0)->G4double{ return a->GetCurrentDoubleValue(arg0); });
  t156.method("GetCurrentDoubleValue", [](G4UImanager* a, const char * arg0, G4int arg1)->G4double{ return a->GetCurrentDoubleValue(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for G4String G4UImanager::GetCurrentStringValue(const char *, const char *, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4String G4UImanager::GetCurrentStringValue(const char *, const char *, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:160:14
  t156.method("GetCurrentStringValue", static_cast<G4String (G4UImanager::*)(const char *, const char *, G4bool) >(&G4UImanager::GetCurrentStringValue));
  t156.method("GetCurrentStringValue", [](G4UImanager& a, const char * arg0, const char * arg1)->G4String{ return a.GetCurrentStringValue(arg0, arg1); });
  t156.method("GetCurrentStringValue", [](G4UImanager* a, const char * arg0, const char * arg1)->G4String{ return a->GetCurrentStringValue(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for G4int G4UImanager::GetCurrentIntValue(const char *, const char *, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4int G4UImanager::GetCurrentIntValue(const char *, const char *, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:163:11
  t156.method("GetCurrentIntValue", static_cast<G4int (G4UImanager::*)(const char *, const char *, G4bool) >(&G4UImanager::GetCurrentIntValue));
  t156.method("GetCurrentIntValue", [](G4UImanager& a, const char * arg0, const char * arg1)->G4int{ return a.GetCurrentIntValue(arg0, arg1); });
  t156.method("GetCurrentIntValue", [](G4UImanager* a, const char * arg0, const char * arg1)->G4int{ return a->GetCurrentIntValue(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for G4double G4UImanager::GetCurrentDoubleValue(const char *, const char *, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4UImanager::GetCurrentDoubleValue(const char *, const char *, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:165:14
  t156.method("GetCurrentDoubleValue", static_cast<G4double (G4UImanager::*)(const char *, const char *, G4bool) >(&G4UImanager::GetCurrentDoubleValue));
  t156.method("GetCurrentDoubleValue", [](G4UImanager& a, const char * arg0, const char * arg1)->G4double{ return a.GetCurrentDoubleValue(arg0, arg1); });
  t156.method("GetCurrentDoubleValue", [](G4UImanager* a, const char * arg0, const char * arg1)->G4double{ return a->GetCurrentDoubleValue(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetPauseAtBeginOfEvent(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetPauseAtBeginOfEvent(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:177:17
  t156.method("SetPauseAtBeginOfEvent", static_cast<void (G4UImanager::*)(G4bool) >(&G4UImanager::SetPauseAtBeginOfEvent));

  DEBUG_MSG("Adding wrapper for G4bool G4UImanager::GetPauseAtBeginOfEvent() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4UImanager::GetPauseAtBeginOfEvent()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:178:19
  t156.method("GetPauseAtBeginOfEvent", static_cast<G4bool (G4UImanager::*)()  const>(&G4UImanager::GetPauseAtBeginOfEvent));

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetPauseAtEndOfEvent(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetPauseAtEndOfEvent(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:179:17
  t156.method("SetPauseAtEndOfEvent", static_cast<void (G4UImanager::*)(G4bool) >(&G4UImanager::SetPauseAtEndOfEvent));

  DEBUG_MSG("Adding wrapper for G4bool G4UImanager::GetPauseAtEndOfEvent() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4UImanager::GetPauseAtEndOfEvent()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:180:19
  t156.method("GetPauseAtEndOfEvent", static_cast<G4bool (G4UImanager::*)()  const>(&G4UImanager::GetPauseAtEndOfEvent));

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetVerboseLevel(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetVerboseLevel(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:199:17
  t156.method("SetVerboseLevel", static_cast<void (G4UImanager::*)(G4int) >(&G4UImanager::SetVerboseLevel));

  DEBUG_MSG("Adding wrapper for G4int G4UImanager::GetVerboseLevel() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4UImanager::GetVerboseLevel()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:200:18
  t156.method("GetVerboseLevel", static_cast<G4int (G4UImanager::*)()  const>(&G4UImanager::GetVerboseLevel));

  DEBUG_MSG("Adding wrapper for G4int G4UImanager::GetNumberOfHistory() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4UImanager::GetNumberOfHistory()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:201:18
  t156.method("GetNumberOfHistory", static_cast<G4int (G4UImanager::*)()  const>(&G4UImanager::GetNumberOfHistory));

  DEBUG_MSG("Adding wrapper for G4String G4UImanager::GetPreviousCommand(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4String G4UImanager::GetPreviousCommand(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:202:21
  t156.method("GetPreviousCommand", static_cast<G4String (G4UImanager::*)(G4int)  const>(&G4UImanager::GetPreviousCommand));

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetMaxHistSize(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetMaxHistSize(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:211:17
  t156.method("SetMaxHistSize", static_cast<void (G4UImanager::*)(G4int) >(&G4UImanager::SetMaxHistSize));

  DEBUG_MSG("Adding wrapper for G4int G4UImanager::GetMaxHistSize() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4UImanager::GetMaxHistSize()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:212:18
  t156.method("GetMaxHistSize", static_cast<G4int (G4UImanager::*)()  const>(&G4UImanager::GetMaxHistSize));

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetMacroSearchPath(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetMacroSearchPath(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:214:17
  t156.method("SetMacroSearchPath", static_cast<void (G4UImanager::*)(const G4String &) >(&G4UImanager::SetMacroSearchPath));

  DEBUG_MSG("Adding wrapper for const G4String & G4UImanager::GetMacroSearchPath() (" __HERE__ ")");
  // signature to use in the veto list: const G4String & G4UImanager::GetMacroSearchPath()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:215:28
  t156.method("GetMacroSearchPath", static_cast<const G4String & (G4UImanager::*)()  const>(&G4UImanager::GetMacroSearchPath));

  DEBUG_MSG("Adding wrapper for void G4UImanager::ParseMacroSearchPath() (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::ParseMacroSearchPath()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:216:10
  t156.method("ParseMacroSearchPath", static_cast<void (G4UImanager::*)() >(&G4UImanager::ParseMacroSearchPath));

  DEBUG_MSG("Adding wrapper for G4String G4UImanager::FindMacroPath(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4String G4UImanager::FindMacroPath(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:217:14
  t156.method("FindMacroPath", static_cast<G4String (G4UImanager::*)(const G4String &)  const>(&G4UImanager::FindMacroPath));

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetMasterUIManager(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetMasterUIManager(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:219:17
  t156.method("SetMasterUIManager", static_cast<void (G4UImanager::*)(G4bool) >(&G4UImanager::SetMasterUIManager));

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetIgnoreCmdNotFound(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetIgnoreCmdNotFound(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:229:17
  t156.method("SetIgnoreCmdNotFound", static_cast<void (G4UImanager::*)(G4bool) >(&G4UImanager::SetIgnoreCmdNotFound));

  DEBUG_MSG("Adding wrapper for std::vector<G4String> * G4UImanager::GetCommandStack() (" __HERE__ ")");
  // signature to use in the veto list: std::vector<G4String> * G4UImanager::GetCommandStack()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:231:28
  t156.method("GetCommandStack", static_cast<std::vector<G4String> * (G4UImanager::*)() >(&G4UImanager::GetCommandStack));

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetUpForAThread(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetUpForAThread(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:234:10
  t156.method("SetUpForAThread", static_cast<void (G4UImanager::*)(G4int) >(&G4UImanager::SetUpForAThread));

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetUpForSpecialThread(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetUpForSpecialThread(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:237:10
  t156.method("SetUpForSpecialThread", static_cast<void (G4UImanager::*)(const G4String &) >(&G4UImanager::SetUpForSpecialThread));

  DEBUG_MSG("Adding wrapper for G4int G4UImanager::GetThreadID() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4UImanager::GetThreadID()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:239:18
  t156.method("GetThreadID", static_cast<G4int (G4UImanager::*)()  const>(&G4UImanager::GetThreadID));

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetCoutFileName(const G4String &, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetCoutFileName(const G4String &, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:241:10
  t156.method("SetCoutFileName", static_cast<void (G4UImanager::*)(const G4String &, G4bool) >(&G4UImanager::SetCoutFileName));
  t156.method("SetCoutFileName", [](G4UImanager& a)->void{ a.SetCoutFileName(); });
  t156.method("SetCoutFileName", [](G4UImanager& a, const G4String & arg0)->void{ a.SetCoutFileName(arg0); });
  t156.method("SetCoutFileName", [](G4UImanager* a)->void{ a->SetCoutFileName(); });
  t156.method("SetCoutFileName", [](G4UImanager* a, const G4String & arg0)->void{ a->SetCoutFileName(arg0); });

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetCerrFileName(const G4String &, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetCerrFileName(const G4String &, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:243:10
  t156.method("SetCerrFileName", static_cast<void (G4UImanager::*)(const G4String &, G4bool) >(&G4UImanager::SetCerrFileName));
  t156.method("SetCerrFileName", [](G4UImanager& a)->void{ a.SetCerrFileName(); });
  t156.method("SetCerrFileName", [](G4UImanager& a, const G4String & arg0)->void{ a.SetCerrFileName(arg0); });
  t156.method("SetCerrFileName", [](G4UImanager* a)->void{ a->SetCerrFileName(); });
  t156.method("SetCerrFileName", [](G4UImanager* a, const G4String & arg0)->void{ a->SetCerrFileName(arg0); });

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetThreadPrefixString(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetThreadPrefixString(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:245:10
  t156.method("SetThreadPrefixString", static_cast<void (G4UImanager::*)(const G4String &) >(&G4UImanager::SetThreadPrefixString));
  t156.method("SetThreadPrefixString", [](G4UImanager& a)->void{ a.SetThreadPrefixString(); });
  t156.method("SetThreadPrefixString", [](G4UImanager* a)->void{ a->SetThreadPrefixString(); });

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetThreadUseBuffer(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetThreadUseBuffer(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:246:10
  t156.method("SetThreadUseBuffer", static_cast<void (G4UImanager::*)(G4bool) >(&G4UImanager::SetThreadUseBuffer));
  t156.method("SetThreadUseBuffer", [](G4UImanager& a)->void{ a.SetThreadUseBuffer(); });
  t156.method("SetThreadUseBuffer", [](G4UImanager* a)->void{ a->SetThreadUseBuffer(); });

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetThreadIgnore(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetThreadIgnore(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:247:10
  t156.method("SetThreadIgnore", static_cast<void (G4UImanager::*)(G4int) >(&G4UImanager::SetThreadIgnore));
  t156.method("SetThreadIgnore", [](G4UImanager& a)->void{ a.SetThreadIgnore(); });
  t156.method("SetThreadIgnore", [](G4UImanager* a)->void{ a->SetThreadIgnore(); });

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetThreadIgnoreInit(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetThreadIgnoreInit(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:248:10
  t156.method("SetThreadIgnoreInit", static_cast<void (G4UImanager::*)(G4bool) >(&G4UImanager::SetThreadIgnoreInit));
  t156.method("SetThreadIgnoreInit", [](G4UImanager& a)->void{ a.SetThreadIgnoreInit(); });
  t156.method("SetThreadIgnoreInit", [](G4UImanager* a)->void{ a->SetThreadIgnoreInit(); });

  DEBUG_MSG("Adding wrapper for void G4UImanager::UseDoublePrecisionStr(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::UseDoublePrecisionStr(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:251:17
  t156.method("G4UImanager!UseDoublePrecisionStr", static_cast<void (*)(G4bool) >(&G4UImanager::UseDoublePrecisionStr));

  DEBUG_MSG("Adding wrapper for G4bool G4UImanager::DoublePrecisionStr() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4UImanager::DoublePrecisionStr()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:252:19
  t156.method("G4UImanager!DoublePrecisionStr", static_cast<G4bool (*)() >(&G4UImanager::DoublePrecisionStr));

  DEBUG_MSG("Adding wrapper for G4int G4UImanager::GetLastReturnCode() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4UImanager::GetLastReturnCode()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:254:18
  t156.method("GetLastReturnCode", static_cast<G4int (G4UImanager::*)()  const>(&G4UImanager::GetLastReturnCode));

  DEBUG_MSG("Adding wrapper for bool G4UImanager::IsLastCommandOutputTreated() (" __HERE__ ")");
  // signature to use in the veto list: bool G4UImanager::IsLastCommandOutputTreated()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:256:17
  t156.method("IsLastCommandOutputTreated", static_cast<bool (G4UImanager::*)() >(&G4UImanager::IsLastCommandOutputTreated));

  DEBUG_MSG("Adding wrapper for void G4UImanager::SetLastCommandOutputTreated() (" __HERE__ ")");
  // signature to use in the veto list: void G4UImanager::SetLastCommandOutputTreated()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UImanager.hh:257:17
  t156.method("SetLastCommandOutputTreated", static_cast<void (G4UImanager::*)() >(&G4UImanager::SetLastCommandOutputTreated));

  /* End of G4UImanager class method wrappers
   **********************************************************************/

}
