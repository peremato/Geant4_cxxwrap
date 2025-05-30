// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4VUserPhysicsList> : std::false_type { };
  template<> struct DefaultConstructible<G4VUserPhysicsList> : std::false_type { };
}

// Class generating the wrapper for type G4VUserPhysicsList
// signature to use in the veto file: G4VUserPhysicsList
struct JlG4VUserPhysicsList: public Wrapper {

  JlG4VUserPhysicsList(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4VUserPhysicsList (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:103:7
    jlcxx::TypeWrapper<G4VUserPhysicsList>  t = jlModule.add_type<G4VUserPhysicsList>("G4VUserPhysicsList");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4VUserPhysicsList>>(new jlcxx::TypeWrapper<G4VUserPhysicsList>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for G4VUserPhysicsList & G4VUserPhysicsList::operator=(const G4VUserPhysicsList &) (" __HERE__ ")");
    // signature to use in the veto list: G4VUserPhysicsList & G4VUserPhysicsList::operator=(const G4VUserPhysicsList &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:109:25
    t.method("assign", static_cast<G4VUserPhysicsList & (G4VUserPhysicsList::*)(const G4VUserPhysicsList &) >(&G4VUserPhysicsList::operator=));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::ConstructParticle() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::ConstructParticle()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:113:18
    t.method("ConstructParticle", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::ConstructParticle));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::Construct() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::Construct()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:117:10
    t.method("Construct", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::Construct));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::ConstructProcess() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::ConstructProcess()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:122:18
    t.method("ConstructProcess", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::ConstructProcess));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetCuts() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::SetCuts()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:125:18
    t.method("SetCuts", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::SetCuts));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetDefaultCutValue(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::SetDefaultCutValue(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:130:10
    t.method("SetDefaultCutValue", static_cast<void (G4VUserPhysicsList::*)(G4double) >(&G4VUserPhysicsList::SetDefaultCutValue));

    DEBUG_MSG("Adding wrapper for G4double G4VUserPhysicsList::GetDefaultCutValue() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VUserPhysicsList::GetDefaultCutValue()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:131:14
    t.method("GetDefaultCutValue", static_cast<G4double (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::GetDefaultCutValue));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::BuildPhysicsTable() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::BuildPhysicsTable()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:136:10
    t.method("BuildPhysicsTable", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::BuildPhysicsTable));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::PreparePhysicsTable(G4ParticleDefinition *) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::PreparePhysicsTable(G4ParticleDefinition *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:139:10
    t.method("PreparePhysicsTable", static_cast<void (G4VUserPhysicsList::*)(G4ParticleDefinition *) >(&G4VUserPhysicsList::PreparePhysicsTable));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::BuildPhysicsTable(G4ParticleDefinition *) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::BuildPhysicsTable(G4ParticleDefinition *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:142:10
    t.method("BuildPhysicsTable", static_cast<void (G4VUserPhysicsList::*)(G4ParticleDefinition *) >(&G4VUserPhysicsList::BuildPhysicsTable));

    DEBUG_MSG("Adding wrapper for G4bool G4VUserPhysicsList::StorePhysicsTable(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VUserPhysicsList::StorePhysicsTable(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:147:12
    t.method("StorePhysicsTable", static_cast<G4bool (G4VUserPhysicsList::*)(const G4String &) >(&G4VUserPhysicsList::StorePhysicsTable));
    t.method("StorePhysicsTable", [](G4VUserPhysicsList& a)->G4bool { return a.StorePhysicsTable(); });
    t.method("StorePhysicsTable", [](G4VUserPhysicsList* a)->G4bool { return a->StorePhysicsTable(); });

    DEBUG_MSG("Adding wrapper for G4bool G4VUserPhysicsList::IsPhysicsTableRetrieved() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VUserPhysicsList::IsPhysicsTableRetrieved()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:151:12
    t.method("IsPhysicsTableRetrieved", static_cast<G4bool (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::IsPhysicsTableRetrieved));

    DEBUG_MSG("Adding wrapper for G4bool G4VUserPhysicsList::IsStoredInAscii() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VUserPhysicsList::IsStoredInAscii()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:152:12
    t.method("IsStoredInAscii", static_cast<G4bool (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::IsStoredInAscii));

    DEBUG_MSG("Adding wrapper for const G4String & G4VUserPhysicsList::GetPhysicsTableDirectory() (" __HERE__ ")");
    // signature to use in the veto list: const G4String & G4VUserPhysicsList::GetPhysicsTableDirectory()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:155:21
    t.method("GetPhysicsTableDirectory", static_cast<const G4String & (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::GetPhysicsTableDirectory));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetPhysicsTableRetrieved(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::SetPhysicsTableRetrieved(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:160:10
    t.method("SetPhysicsTableRetrieved", static_cast<void (G4VUserPhysicsList::*)(const G4String &) >(&G4VUserPhysicsList::SetPhysicsTableRetrieved));
    t.method("SetPhysicsTableRetrieved", [](G4VUserPhysicsList& a)->void { a.SetPhysicsTableRetrieved(); });
    t.method("SetPhysicsTableRetrieved", [](G4VUserPhysicsList* a)->void { a->SetPhysicsTableRetrieved(); });

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetStoredInAscii() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::SetStoredInAscii()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:161:10
    t.method("SetStoredInAscii", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::SetStoredInAscii));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::ResetPhysicsTableRetrieved() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::ResetPhysicsTableRetrieved()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:164:10
    t.method("ResetPhysicsTableRetrieved", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::ResetPhysicsTableRetrieved));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::ResetStoredInAscii() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::ResetStoredInAscii()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:165:10
    t.method("ResetStoredInAscii", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::ResetStoredInAscii));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::DumpList() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::DumpList()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:168:10
    t.method("DumpList", static_cast<void (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::DumpList));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::DumpCutValuesTable(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::DumpCutValuesTable(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:172:10
    t.method("DumpCutValuesTable", static_cast<void (G4VUserPhysicsList::*)(G4int) >(&G4VUserPhysicsList::DumpCutValuesTable));
    t.method("DumpCutValuesTable", [](G4VUserPhysicsList& a)->void { a.DumpCutValuesTable(); });
    t.method("DumpCutValuesTable", [](G4VUserPhysicsList* a)->void { a->DumpCutValuesTable(); });

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::DumpCutValuesTableIfRequested() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::DumpCutValuesTableIfRequested()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:176:10
    t.method("DumpCutValuesTableIfRequested", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::DumpCutValuesTableIfRequested));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetVerboseLevel(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::SetVerboseLevel(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:182:10
    t.method("SetVerboseLevel", static_cast<void (G4VUserPhysicsList::*)(G4int) >(&G4VUserPhysicsList::SetVerboseLevel));

    DEBUG_MSG("Adding wrapper for G4int G4VUserPhysicsList::GetVerboseLevel() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4VUserPhysicsList::GetVerboseLevel()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:183:11
    t.method("GetVerboseLevel", static_cast<G4int (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::GetVerboseLevel));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::UseCoupledTransportation(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::UseCoupledTransportation(G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:185:10
    t.method("UseCoupledTransportation", static_cast<void (G4VUserPhysicsList::*)(G4bool) >(&G4VUserPhysicsList::UseCoupledTransportation));
    t.method("UseCoupledTransportation", [](G4VUserPhysicsList& a)->void { a.UseCoupledTransportation(); });
    t.method("UseCoupledTransportation", [](G4VUserPhysicsList* a)->void { a->UseCoupledTransportation(); });

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetCutsWithDefault() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::SetCutsWithDefault()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:190:10
    t.method("SetCutsWithDefault", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::SetCutsWithDefault));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetCutValue(G4double, const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::SetCutValue(G4double, const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:193:10
    t.method("SetCutValue", static_cast<void (G4VUserPhysicsList::*)(G4double, const G4String &) >(&G4VUserPhysicsList::SetCutValue));

    DEBUG_MSG("Adding wrapper for G4double G4VUserPhysicsList::GetCutValue(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VUserPhysicsList::GetCutValue(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:196:14
    t.method("GetCutValue", static_cast<G4double (G4VUserPhysicsList::*)(const G4String &)  const>(&G4VUserPhysicsList::GetCutValue));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetCutValue(G4double, const G4String &, const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::SetCutValue(G4double, const G4String &, const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:199:10
    t.method("SetCutValue", static_cast<void (G4VUserPhysicsList::*)(G4double, const G4String &, const G4String &) >(&G4VUserPhysicsList::SetCutValue));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetParticleCuts(G4double, G4ParticleDefinition *, G4Region *) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::SetParticleCuts(G4double, G4ParticleDefinition *, G4Region *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:205:10
    t.method("SetParticleCuts", static_cast<void (G4VUserPhysicsList::*)(G4double, G4ParticleDefinition *, G4Region *) >(&G4VUserPhysicsList::SetParticleCuts));
    t.method("SetParticleCuts", [](G4VUserPhysicsList& a, G4double arg0, G4ParticleDefinition * arg1)->void { a.SetParticleCuts(arg0, arg1); });
    t.method("SetParticleCuts", [](G4VUserPhysicsList* a, G4double arg0, G4ParticleDefinition * arg1)->void { a->SetParticleCuts(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetParticleCuts(G4double, const G4String &, G4Region *) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::SetParticleCuts(G4double, const G4String &, G4Region *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:206:10
    t.method("SetParticleCuts", static_cast<void (G4VUserPhysicsList::*)(G4double, const G4String &, G4Region *) >(&G4VUserPhysicsList::SetParticleCuts));
    t.method("SetParticleCuts", [](G4VUserPhysicsList& a, G4double arg0, const G4String & arg1)->void { a.SetParticleCuts(arg0, arg1); });
    t.method("SetParticleCuts", [](G4VUserPhysicsList* a, G4double arg0, const G4String & arg1)->void { a->SetParticleCuts(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetCutsForRegion(G4double, const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::SetCutsForRegion(G4double, const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:209:10
    t.method("SetCutsForRegion", static_cast<void (G4VUserPhysicsList::*)(G4double, const G4String &) >(&G4VUserPhysicsList::SetCutsForRegion));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::SetApplyCuts(G4bool, const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::SetApplyCuts(G4bool, const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:212:10
    t.method("SetApplyCuts", static_cast<void (G4VUserPhysicsList::*)(G4bool, const G4String &) >(&G4VUserPhysicsList::SetApplyCuts));

    DEBUG_MSG("Adding wrapper for G4bool G4VUserPhysicsList::GetApplyCuts(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VUserPhysicsList::GetApplyCuts(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:213:12
    t.method("GetApplyCuts", static_cast<G4bool (G4VUserPhysicsList::*)(const G4String &)  const>(&G4VUserPhysicsList::GetApplyCuts));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::RemoveProcessManager() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::RemoveProcessManager()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:217:10
    t.method("RemoveProcessManager", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::RemoveProcessManager));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::RemoveTrackingManager() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::RemoveTrackingManager()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:221:10
    t.method("RemoveTrackingManager", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::RemoveTrackingManager));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::AddProcessManager(G4ParticleDefinition *, G4ProcessManager *) (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::AddProcessManager(G4ParticleDefinition *, G4ProcessManager *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:224:10
    t.method("AddProcessManager", static_cast<void (G4VUserPhysicsList::*)(G4ParticleDefinition *, G4ProcessManager *) >(&G4VUserPhysicsList::AddProcessManager));
    t.method("AddProcessManager", [](G4VUserPhysicsList& a, G4ParticleDefinition * arg0)->void { a.AddProcessManager(arg0); });
    t.method("AddProcessManager", [](G4VUserPhysicsList* a, G4ParticleDefinition * arg0)->void { a->AddProcessManager(arg0); });

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::CheckParticleList() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::CheckParticleList()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:228:10
    t.method("CheckParticleList", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::CheckParticleList));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::DisableCheckParticleList() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::DisableCheckParticleList()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:230:10
    t.method("DisableCheckParticleList", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::DisableCheckParticleList));

    DEBUG_MSG("Adding wrapper for G4int G4VUserPhysicsList::GetInstanceID() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4VUserPhysicsList::GetInstanceID()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:232:18
    t.method("GetInstanceID", static_cast<G4int (G4VUserPhysicsList::*)()  const>(&G4VUserPhysicsList::GetInstanceID));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::InitializeWorker() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::InitializeWorker()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:238:18
    t.method("InitializeWorker", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::InitializeWorker));

    DEBUG_MSG("Adding wrapper for void G4VUserPhysicsList::TerminateWorker() (" __HERE__ ")");
    // signature to use in the veto list: void G4VUserPhysicsList::TerminateWorker()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VUserPhysicsList.hh:242:18
    t.method("TerminateWorker", static_cast<void (G4VUserPhysicsList::*)() >(&G4VUserPhysicsList::TerminateWorker));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4VUserPhysicsList>> type_;
};
std::shared_ptr<Wrapper> newJlG4VUserPhysicsList(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4VUserPhysicsList(module));
}
