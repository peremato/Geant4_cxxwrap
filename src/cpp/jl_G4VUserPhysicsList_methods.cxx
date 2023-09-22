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
void add_methods_for_G4VUserPhysicsList(jlcxx::Module& types, jlcxx::TypeWrapper<G4VUserPhysicsList>& t155) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4VUserPhysicsList
   */


  DEBUG_MSG("Adding wrapper for G4VUserPhysicsList & G4VUserPhysicsList::operator=(const G4VUserPhysicsList &) (" __HERE__ ")");
  // signature to use in the veto list: G4VUserPhysicsList & G4VUserPhysicsList::operator=(const G4VUserPhysicsList &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:114:25
  t155.method("assign", static_cast<G4VUserPhysicsList & (G4VUserPhysicsList::*)(const G4VUserPhysicsList &) >(&G4VUserPhysicsList::operator=));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::ConstructParticle() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::ConstructParticle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:117:18
  t155.method("ConstructParticle", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::ConstructParticle));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::Construct() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::Construct()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:121:10
  t155.method("Construct", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::Construct));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::ConstructProcess() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::ConstructProcess()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:125:18
  t155.method("ConstructProcess", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::ConstructProcess));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetCuts() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::SetCuts()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:130:18
  t155.method("SetCuts", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::SetCuts));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetDefaultCutValue(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::SetDefaultCutValue(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:133:10
  t155.method("SetDefaultCutValue", static_cast<void (G4VUserPhysicsList::*)(G4double) >(&G4VUserPhysicsList::SetDefaultCutValue));

  DEBUG_MSG("Adding wrapper for G4double G4VUserPhysicsList::GetDefaultCutValue() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VUserPhysicsList::GetDefaultCutValue()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:134:14
  t155.method("GetDefaultCutValue", static_cast<G4double (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::GetDefaultCutValue));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::BuildPhysicsTable() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::BuildPhysicsTable()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:139:10
  t155.method("BuildPhysicsTable", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::BuildPhysicsTable));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::PreparePhysicsTable(G4ParticleDefinition *) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::PreparePhysicsTable(G4ParticleDefinition *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:144:10
  t155.method("PreparePhysicsTable", static_cast<void (G4VUserPhysicsList::*)(G4ParticleDefinition *) >(&G4VUserPhysicsList::PreparePhysicsTable));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::BuildPhysicsTable(G4ParticleDefinition *) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::BuildPhysicsTable(G4ParticleDefinition *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:147:10
  t155.method("BuildPhysicsTable", static_cast<void (G4VUserPhysicsList::*)(G4ParticleDefinition *) >(&G4VUserPhysicsList::BuildPhysicsTable));

  DEBUG_MSG("Adding wrapper for G4bool G4VUserPhysicsList::StorePhysicsTable(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VUserPhysicsList::StorePhysicsTable(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:150:12
  t155.method("StorePhysicsTable", static_cast<G4bool (G4VUserPhysicsList::*)(const G4String &) >(&G4VUserPhysicsList::StorePhysicsTable));
  t155.method("StorePhysicsTable", [](G4VUserPhysicsList& a)->G4bool{ return a.StorePhysicsTable(); });
  t155.method("StorePhysicsTable", [](G4VUserPhysicsList* a)->G4bool{ return a->StorePhysicsTable(); });

  DEBUG_MSG("Adding wrapper for G4bool G4VUserPhysicsList::IsPhysicsTableRetrieved() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VUserPhysicsList::IsPhysicsTableRetrieved()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:155:12
  t155.method("IsPhysicsTableRetrieved", static_cast<G4bool (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::IsPhysicsTableRetrieved));

  DEBUG_MSG("Adding wrapper for G4bool G4VUserPhysicsList::IsStoredInAscii() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VUserPhysicsList::IsStoredInAscii()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:156:12
  t155.method("IsStoredInAscii", static_cast<G4bool (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::IsStoredInAscii));

  DEBUG_MSG("Adding wrapper for const G4String & G4VUserPhysicsList::GetPhysicsTableDirectory() (" __HERE__ ")");
  // signature to use in the veto list: const G4String & G4VUserPhysicsList::GetPhysicsTableDirectory()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:160:21
  t155.method("GetPhysicsTableDirectory", static_cast<const G4String & (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::GetPhysicsTableDirectory));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetPhysicsTableRetrieved(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::SetPhysicsTableRetrieved(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:163:10
  t155.method("SetPhysicsTableRetrieved", static_cast<void (G4VUserPhysicsList::*)(const G4String &) >(&G4VUserPhysicsList::SetPhysicsTableRetrieved));
  t155.method("SetPhysicsTableRetrieved", [](G4VUserPhysicsList& a)->void{ a.SetPhysicsTableRetrieved(); });
  t155.method("SetPhysicsTableRetrieved", [](G4VUserPhysicsList* a)->void{ a->SetPhysicsTableRetrieved(); });

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetStoredInAscii() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::SetStoredInAscii()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:164:10
  t155.method("SetStoredInAscii", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::SetStoredInAscii));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::ResetPhysicsTableRetrieved() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::ResetPhysicsTableRetrieved()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:169:10
  t155.method("ResetPhysicsTableRetrieved", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::ResetPhysicsTableRetrieved));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::ResetStoredInAscii() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::ResetStoredInAscii()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:170:10
  t155.method("ResetStoredInAscii", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::ResetStoredInAscii));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::DumpList() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::DumpList()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:173:10
  t155.method("DumpList", static_cast<void (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::DumpList));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::DumpCutValuesTable(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::DumpCutValuesTable(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:176:10
  t155.method("DumpCutValuesTable", static_cast<void (G4VUserPhysicsList::*)(G4int) >(&G4VUserPhysicsList::DumpCutValuesTable));
  t155.method("DumpCutValuesTable", [](G4VUserPhysicsList& a)->void{ a.DumpCutValuesTable(); });
  t155.method("DumpCutValuesTable", [](G4VUserPhysicsList* a)->void{ a->DumpCutValuesTable(); });

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::DumpCutValuesTableIfRequested() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::DumpCutValuesTableIfRequested()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:180:10
  t155.method("DumpCutValuesTableIfRequested", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::DumpCutValuesTableIfRequested));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetVerboseLevel(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::SetVerboseLevel(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:184:10
  t155.method("SetVerboseLevel", static_cast<void (G4VUserPhysicsList::*)(G4int) >(&G4VUserPhysicsList::SetVerboseLevel));

  DEBUG_MSG("Adding wrapper for G4int G4VUserPhysicsList::GetVerboseLevel() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VUserPhysicsList::GetVerboseLevel()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:185:11
  t155.method("GetVerboseLevel", static_cast<G4int (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::GetVerboseLevel));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::UseCoupledTransportation(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::UseCoupledTransportation(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:191:10
  t155.method("UseCoupledTransportation", static_cast<void (G4VUserPhysicsList::*)(G4bool) >(&G4VUserPhysicsList::UseCoupledTransportation));
  t155.method("UseCoupledTransportation", [](G4VUserPhysicsList& a)->void{ a.UseCoupledTransportation(); });
  t155.method("UseCoupledTransportation", [](G4VUserPhysicsList* a)->void{ a->UseCoupledTransportation(); });

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetCutsWithDefault() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::SetCutsWithDefault()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:193:10
  t155.method("SetCutsWithDefault", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::SetCutsWithDefault));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetCutValue(G4double, const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::SetCutValue(G4double, const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:198:10
  t155.method("SetCutValue", static_cast<void (G4VUserPhysicsList::*)(G4double, const G4String &) >(&G4VUserPhysicsList::SetCutValue));

  DEBUG_MSG("Adding wrapper for G4double G4VUserPhysicsList::GetCutValue(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VUserPhysicsList::GetCutValue(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:201:14
  t155.method("GetCutValue", static_cast<G4double (G4VUserPhysicsList::*)(const G4String &)  const>(&G4VUserPhysicsList::GetCutValue));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetCutValue(G4double, const G4String &, const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::SetCutValue(G4double, const G4String &, const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:204:10
  t155.method("SetCutValue", static_cast<void (G4VUserPhysicsList::*)(G4double, const G4String &, const G4String &) >(&G4VUserPhysicsList::SetCutValue));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetParticleCuts(G4double, G4ParticleDefinition *, G4Region *) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::SetParticleCuts(G4double, G4ParticleDefinition *, G4Region *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:208:10
  t155.method("SetParticleCuts", static_cast<void (G4VUserPhysicsList::*)(G4double, G4ParticleDefinition *, G4Region *) >(&G4VUserPhysicsList::SetParticleCuts));
  t155.method("SetParticleCuts", [](G4VUserPhysicsList& a, G4double arg0, G4ParticleDefinition * arg1)->void{ a.SetParticleCuts(arg0, arg1); });
  t155.method("SetParticleCuts", [](G4VUserPhysicsList* a, G4double arg0, G4ParticleDefinition * arg1)->void{ a->SetParticleCuts(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetParticleCuts(G4double, const G4String &, G4Region *) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::SetParticleCuts(G4double, const G4String &, G4Region *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:210:10
  t155.method("SetParticleCuts", static_cast<void (G4VUserPhysicsList::*)(G4double, const G4String &, G4Region *) >(&G4VUserPhysicsList::SetParticleCuts));
  t155.method("SetParticleCuts", [](G4VUserPhysicsList& a, G4double arg0, const G4String & arg1)->void{ a.SetParticleCuts(arg0, arg1); });
  t155.method("SetParticleCuts", [](G4VUserPhysicsList* a, G4double arg0, const G4String & arg1)->void{ a->SetParticleCuts(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetCutsForRegion(G4double, const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::SetCutsForRegion(G4double, const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:217:10
  t155.method("SetCutsForRegion", static_cast<void (G4VUserPhysicsList::*)(G4double, const G4String &) >(&G4VUserPhysicsList::SetCutsForRegion));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetApplyCuts(G4bool, const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::SetApplyCuts(G4bool, const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:220:10
  t155.method("SetApplyCuts", static_cast<void (G4VUserPhysicsList::*)(G4bool, const G4String &) >(&G4VUserPhysicsList::SetApplyCuts));

  DEBUG_MSG("Adding wrapper for G4bool G4VUserPhysicsList::GetApplyCuts(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VUserPhysicsList::GetApplyCuts(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:221:12
  t155.method("GetApplyCuts", static_cast<G4bool (G4VUserPhysicsList::*)(const G4String &)  const>(&G4VUserPhysicsList::GetApplyCuts));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::RemoveProcessManager() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::RemoveProcessManager()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:224:10
  t155.method("RemoveProcessManager", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::RemoveProcessManager));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::RemoveTrackingManager() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::RemoveTrackingManager()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:228:10
  t155.method("RemoveTrackingManager", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::RemoveTrackingManager));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::AddProcessManager(G4ParticleDefinition *, G4ProcessManager *) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::AddProcessManager(G4ParticleDefinition *, G4ProcessManager *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:232:10
  t155.method("AddProcessManager", static_cast<void (G4VUserPhysicsList::*)(G4ParticleDefinition *, G4ProcessManager *) >(&G4VUserPhysicsList::AddProcessManager));
  t155.method("AddProcessManager", [](G4VUserPhysicsList& a, G4ParticleDefinition * arg0)->void{ a.AddProcessManager(arg0); });
  t155.method("AddProcessManager", [](G4VUserPhysicsList* a, G4ParticleDefinition * arg0)->void{ a->AddProcessManager(arg0); });

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::CheckParticleList() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::CheckParticleList()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:236:10
  t155.method("CheckParticleList", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::CheckParticleList));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::DisableCheckParticleList() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::DisableCheckParticleList()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:239:10
  t155.method("DisableCheckParticleList", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::DisableCheckParticleList));

  DEBUG_MSG("Adding wrapper for G4int G4VUserPhysicsList::GetInstanceID() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VUserPhysicsList::GetInstanceID()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:241:18
  t155.method("GetInstanceID", static_cast<G4int (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::GetInstanceID));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::InitializeWorker() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::InitializeWorker()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:246:18
  t155.method("InitializeWorker", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::InitializeWorker));

  DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::TerminateWorker() (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPhysicsList::TerminateWorker()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPhysicsList.hh:249:18
  t155.method("TerminateWorker", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::TerminateWorker));

  /* End of G4VUserPhysicsList class method wrappers
   **********************************************************************/

}
