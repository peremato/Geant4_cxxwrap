// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4VProcess> : std::false_type { };
  template<> struct DefaultConstructible<G4VProcess> : std::false_type { };
}

// Class generating the wrapper for type G4VProcess
// signature to use in the veto file: G4VProcess
struct JlG4VProcess: public Wrapper {

  JlG4VProcess(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4VProcess (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:60:7
    jlcxx::TypeWrapper<G4VProcess>  t = jlModule.add_type<G4VProcess>("G4VProcess");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4VProcess>>(new jlcxx::TypeWrapper<G4VProcess>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    module_.set_override_module(jl_base_module);

    DEBUG_MSG("Adding wrapper for G4bool G4VProcess::operator==(const G4VProcess &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VProcess::operator==(const G4VProcess &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:78:12
    t.method("==", static_cast<G4bool (G4VProcess::*)(const G4VProcess &)  const>(&G4VProcess::operator==));

    DEBUG_MSG("Adding wrapper for G4bool G4VProcess::operator!=(const G4VProcess &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VProcess::operator!=(const G4VProcess &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:79:12
    t.method("!=", static_cast<G4bool (G4VProcess::*)(const G4VProcess &)  const>(&G4VProcess::operator!=));

    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for G4VParticleChange * G4VProcess::PostStepDoIt(const G4Track &, const G4Step &) (" __HERE__ ")");
    // signature to use in the veto list: G4VParticleChange * G4VProcess::PostStepDoIt(const G4Track &, const G4Step &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:86:32
    t.method("PostStepDoIt", static_cast<G4VParticleChange * (G4VProcess::*)(const G4Track &, const G4Step &) >(&G4VProcess::PostStepDoIt));

    DEBUG_MSG("Adding wrapper for G4VParticleChange * G4VProcess::AlongStepDoIt(const G4Track &, const G4Step &) (" __HERE__ ")");
    // signature to use in the veto list: G4VParticleChange * G4VProcess::AlongStepDoIt(const G4Track &, const G4Step &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:91:32
    t.method("AlongStepDoIt", static_cast<G4VParticleChange * (G4VProcess::*)(const G4Track &, const G4Step &) >(&G4VProcess::AlongStepDoIt));

    DEBUG_MSG("Adding wrapper for G4VParticleChange * G4VProcess::AtRestDoIt(const G4Track &, const G4Step &) (" __HERE__ ")");
    // signature to use in the veto list: G4VParticleChange * G4VProcess::AtRestDoIt(const G4Track &, const G4Step &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:95:32
    t.method("AtRestDoIt", static_cast<G4VParticleChange * (G4VProcess::*)(const G4Track &, const G4Step &) >(&G4VProcess::AtRestDoIt));

    DEBUG_MSG("Adding wrapper for G4double G4VProcess::AlongStepGetPhysicalInteractionLength(const G4Track &, G4double, G4double, G4double &, G4GPILSelection *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VProcess::AlongStepGetPhysicalInteractionLength(const G4Track &, G4double, G4double, G4double &, G4GPILSelection *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:113:22
    t.method("AlongStepGetPhysicalInteractionLength", static_cast<G4double (G4VProcess::*)(const G4Track &, G4double, G4double, G4double &, G4GPILSelection *) >(&G4VProcess::AlongStepGetPhysicalInteractionLength));

    DEBUG_MSG("Adding wrapper for G4double G4VProcess::AtRestGetPhysicalInteractionLength(const G4Track &, G4ForceCondition *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VProcess::AtRestGetPhysicalInteractionLength(const G4Track &, G4ForceCondition *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:120:22
    t.method("AtRestGetPhysicalInteractionLength", static_cast<G4double (G4VProcess::*)(const G4Track &, G4ForceCondition *) >(&G4VProcess::AtRestGetPhysicalInteractionLength));

    DEBUG_MSG("Adding wrapper for G4double G4VProcess::PostStepGetPhysicalInteractionLength(const G4Track &, G4double, G4ForceCondition *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VProcess::PostStepGetPhysicalInteractionLength(const G4Track &, G4double, G4ForceCondition *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:124:22
    t.method("PostStepGetPhysicalInteractionLength", static_cast<G4double (G4VProcess::*)(const G4Track &, G4double, G4ForceCondition *) >(&G4VProcess::PostStepGetPhysicalInteractionLength));

    DEBUG_MSG("Adding wrapper for G4double G4VProcess::GetCurrentInteractionLength() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VProcess::GetCurrentInteractionLength()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:154:21
    t.method("GetCurrentInteractionLength", static_cast<G4double (G4VProcess::*)()  const>(&G4VProcess::GetCurrentInteractionLength));

    DEBUG_MSG("Adding wrapper for void G4VProcess::SetPILfactor(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::SetPILfactor(G4double)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:159:17
    t.method("SetPILfactor", static_cast<void (G4VProcess::*)(G4double) >(&G4VProcess::SetPILfactor));

    DEBUG_MSG("Adding wrapper for G4double G4VProcess::GetPILfactor() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VProcess::GetPILfactor()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:160:21
    t.method("GetPILfactor", static_cast<G4double (G4VProcess::*)()  const>(&G4VProcess::GetPILfactor));

    DEBUG_MSG("Adding wrapper for G4double G4VProcess::AlongStepGPIL(const G4Track &, G4double, G4double, G4double &, G4GPILSelection *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VProcess::AlongStepGPIL(const G4Track &, G4double, G4double, G4double &, G4GPILSelection *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:169:21
    t.method("AlongStepGPIL", static_cast<G4double (G4VProcess::*)(const G4Track &, G4double, G4double, G4double &, G4GPILSelection *) >(&G4VProcess::AlongStepGPIL));

    DEBUG_MSG("Adding wrapper for G4double G4VProcess::AtRestGPIL(const G4Track &, G4ForceCondition *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VProcess::AtRestGPIL(const G4Track &, G4ForceCondition *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:175:21
    t.method("AtRestGPIL", static_cast<G4double (G4VProcess::*)(const G4Track &, G4ForceCondition *) >(&G4VProcess::AtRestGPIL));

    DEBUG_MSG("Adding wrapper for G4double G4VProcess::PostStepGPIL(const G4Track &, G4double, G4ForceCondition *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VProcess::PostStepGPIL(const G4Track &, G4double, G4ForceCondition *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:178:21
    t.method("PostStepGPIL", static_cast<G4double (G4VProcess::*)(const G4Track &, G4double, G4ForceCondition *) >(&G4VProcess::PostStepGPIL));

    DEBUG_MSG("Adding wrapper for G4bool G4VProcess::IsApplicable(const G4ParticleDefinition &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VProcess::IsApplicable(const G4ParticleDefinition &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:182:20
    t.method("IsApplicable", static_cast<G4bool (G4VProcess::*)(const G4ParticleDefinition &) >(&G4VProcess::IsApplicable));

    DEBUG_MSG("Adding wrapper for void G4VProcess::BuildPhysicsTable(const G4ParticleDefinition &) (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::BuildPhysicsTable(const G4ParticleDefinition &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:187:18
    t.method("BuildPhysicsTable", static_cast<void (G4VProcess::*)(const G4ParticleDefinition &) >(&G4VProcess::BuildPhysicsTable));

    DEBUG_MSG("Adding wrapper for void G4VProcess::PreparePhysicsTable(const G4ParticleDefinition &) (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::PreparePhysicsTable(const G4ParticleDefinition &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:194:18
    t.method("PreparePhysicsTable", static_cast<void (G4VProcess::*)(const G4ParticleDefinition &) >(&G4VProcess::PreparePhysicsTable));

    DEBUG_MSG("Adding wrapper for G4bool G4VProcess::StorePhysicsTable(const G4ParticleDefinition *, const G4String &, G4bool) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VProcess::StorePhysicsTable(const G4ParticleDefinition *, const G4String &, G4bool)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:206:20
    t.method("StorePhysicsTable", static_cast<G4bool (G4VProcess::*)(const G4ParticleDefinition *, const G4String &, G4bool) >(&G4VProcess::StorePhysicsTable));

    DEBUG_MSG("Adding wrapper for G4bool G4VProcess::RetrievePhysicsTable(const G4ParticleDefinition *, const G4String &, G4bool) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VProcess::RetrievePhysicsTable(const G4ParticleDefinition *, const G4String &, G4bool)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:211:20
    t.method("RetrievePhysicsTable", static_cast<G4bool (G4VProcess::*)(const G4ParticleDefinition *, const G4String &, G4bool) >(&G4VProcess::RetrievePhysicsTable));

    DEBUG_MSG("Adding wrapper for const G4String & G4VProcess::GetPhysicsTableFileName(const G4ParticleDefinition *, const G4String &, const G4String &, G4bool) (" __HERE__ ")");
    // signature to use in the veto list: const G4String & G4VProcess::GetPhysicsTableFileName(const G4ParticleDefinition *, const G4String &, const G4String &, G4bool)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:219:21
    t.method("GetPhysicsTableFileName", static_cast<const G4String & (G4VProcess::*)(const G4ParticleDefinition *, const G4String &, const G4String &, G4bool) >(&G4VProcess::GetPhysicsTableFileName));
    t.method("GetPhysicsTableFileName", [](G4VProcess& a, const G4ParticleDefinition * arg0, const G4String & arg1, const G4String & arg2)->const G4String & { return a.GetPhysicsTableFileName(arg0, arg1, arg2); });
    t.method("GetPhysicsTableFileName", [](G4VProcess* a, const G4ParticleDefinition * arg0, const G4String & arg1, const G4String & arg2)->const G4String & { return a->GetPhysicsTableFileName(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for const G4String & G4VProcess::GetProcessName() (" __HERE__ ")");
    // signature to use in the veto list: const G4String & G4VProcess::GetProcessName()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:225:28
    t.method("GetProcessName", static_cast<const G4String & (G4VProcess::*)()  const>(&G4VProcess::GetProcessName));

    DEBUG_MSG("Adding wrapper for G4ProcessType G4VProcess::GetProcessType() (" __HERE__ ")");
    // signature to use in the veto list: G4ProcessType G4VProcess::GetProcessType()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:228:26
    t.method("GetProcessType", static_cast<G4ProcessType (G4VProcess::*)()  const>(&G4VProcess::GetProcessType));

    DEBUG_MSG("Adding wrapper for void G4VProcess::SetProcessType(G4ProcessType) (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::SetProcessType(G4ProcessType)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:231:17
    t.method("SetProcessType", static_cast<void (G4VProcess::*)(G4ProcessType) >(&G4VProcess::SetProcessType));

    DEBUG_MSG("Adding wrapper for G4int G4VProcess::GetProcessSubType() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4VProcess::GetProcessSubType()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:234:18
    t.method("GetProcessSubType", static_cast<G4int (G4VProcess::*)()  const>(&G4VProcess::GetProcessSubType));

    DEBUG_MSG("Adding wrapper for void G4VProcess::SetProcessSubType(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::SetProcessSubType(G4int)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:237:17
    t.method("SetProcessSubType", static_cast<void (G4VProcess::*)(G4int) >(&G4VProcess::SetProcessSubType));

    DEBUG_MSG("Adding wrapper for const G4String & G4VProcess::GetProcessTypeName(G4ProcessType) (" __HERE__ ")");
    // signature to use in the veto list: const G4String & G4VProcess::GetProcessTypeName(G4ProcessType)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:240:28
    module_.method("G4VProcess!GetProcessTypeName", static_cast<const G4String & (*)(G4ProcessType) >(&G4VProcess::GetProcessTypeName));

    DEBUG_MSG("Adding wrapper for const G4VProcess * G4VProcess::GetCreatorProcess() (" __HERE__ ")");
    // signature to use in the veto list: const G4VProcess * G4VProcess::GetCreatorProcess()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:243:31
    t.method("GetCreatorProcess", static_cast<const G4VProcess * (G4VProcess::*)()  const>(&G4VProcess::GetCreatorProcess));

    DEBUG_MSG("Adding wrapper for void G4VProcess::StartTracking(G4Track *) (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::StartTracking(G4Track *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:247:18
    t.method("StartTracking", static_cast<void (G4VProcess::*)(G4Track *) >(&G4VProcess::StartTracking));

    DEBUG_MSG("Adding wrapper for void G4VProcess::EndTracking() (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::EndTracking()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:248:18
    t.method("EndTracking", static_cast<void (G4VProcess::*)() >(&G4VProcess::EndTracking));

    DEBUG_MSG("Adding wrapper for void G4VProcess::SetProcessManager(const G4ProcessManager *) (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::SetProcessManager(const G4ProcessManager *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:251:18
    t.method("SetProcessManager", static_cast<void (G4VProcess::*)(const G4ProcessManager *) >(&G4VProcess::SetProcessManager));

    DEBUG_MSG("Adding wrapper for const G4ProcessManager * G4VProcess::GetProcessManager() (" __HERE__ ")");
    // signature to use in the veto list: const G4ProcessManager * G4VProcess::GetProcessManager()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:254:37
    t.method("GetProcessManager", static_cast<const G4ProcessManager * (G4VProcess::*)() >(&G4VProcess::GetProcessManager));

    DEBUG_MSG("Adding wrapper for void G4VProcess::ResetNumberOfInteractionLengthLeft() (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::ResetNumberOfInteractionLengthLeft()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:257:18
    t.method("ResetNumberOfInteractionLengthLeft", static_cast<void (G4VProcess::*)() >(&G4VProcess::ResetNumberOfInteractionLengthLeft));

    DEBUG_MSG("Adding wrapper for G4double G4VProcess::GetNumberOfInteractionLengthLeft() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VProcess::GetNumberOfInteractionLengthLeft()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:260:21
    t.method("GetNumberOfInteractionLengthLeft", static_cast<G4double (G4VProcess::*)()  const>(&G4VProcess::GetNumberOfInteractionLengthLeft));

    DEBUG_MSG("Adding wrapper for G4double G4VProcess::GetTotalNumberOfInteractionLengthTraversed() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VProcess::GetTotalNumberOfInteractionLengthTraversed()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:263:21
    t.method("GetTotalNumberOfInteractionLengthTraversed", static_cast<G4double (G4VProcess::*)()  const>(&G4VProcess::GetTotalNumberOfInteractionLengthTraversed));

    DEBUG_MSG("Adding wrapper for G4bool G4VProcess::isAtRestDoItIsEnabled() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VProcess::isAtRestDoItIsEnabled()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:267:19
    t.method("isAtRestDoItIsEnabled", static_cast<G4bool (G4VProcess::*)()  const>(&G4VProcess::isAtRestDoItIsEnabled));

    DEBUG_MSG("Adding wrapper for G4bool G4VProcess::isAlongStepDoItIsEnabled() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VProcess::isAlongStepDoItIsEnabled()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:268:19
    t.method("isAlongStepDoItIsEnabled", static_cast<G4bool (G4VProcess::*)()  const>(&G4VProcess::isAlongStepDoItIsEnabled));

    DEBUG_MSG("Adding wrapper for G4bool G4VProcess::isPostStepDoItIsEnabled() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4VProcess::isPostStepDoItIsEnabled()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:269:19
    t.method("isPostStepDoItIsEnabled", static_cast<G4bool (G4VProcess::*)()  const>(&G4VProcess::isPostStepDoItIsEnabled));

    DEBUG_MSG("Adding wrapper for void G4VProcess::DumpInfo() (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::DumpInfo()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:274:19
    t.method("DumpInfo", static_cast<void (G4VProcess::*)()  const>(&G4VProcess::DumpInfo));

    DEBUG_MSG("Adding wrapper for void G4VProcess::SetVerboseLevel(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::SetVerboseLevel(G4int)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:280:18
    t.method("SetVerboseLevel", static_cast<void (G4VProcess::*)(G4int) >(&G4VProcess::SetVerboseLevel));

    DEBUG_MSG("Adding wrapper for G4int G4VProcess::GetVerboseLevel() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4VProcess::GetVerboseLevel()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:281:18
    t.method("GetVerboseLevel", static_cast<G4int (G4VProcess::*)()  const>(&G4VProcess::GetVerboseLevel));

    DEBUG_MSG("Adding wrapper for void G4VProcess::SetMasterProcess(G4VProcess *) (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::SetMasterProcess(G4VProcess *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:287:18
    t.method("SetMasterProcess", static_cast<void (G4VProcess::*)(G4VProcess *) >(&G4VProcess::SetMasterProcess));

    DEBUG_MSG("Adding wrapper for const G4VProcess * G4VProcess::GetMasterProcess() (" __HERE__ ")");
    // signature to use in the veto list: const G4VProcess * G4VProcess::GetMasterProcess()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:289:30
    t.method("GetMasterProcess", static_cast<const G4VProcess * (G4VProcess::*)()  const>(&G4VProcess::GetMasterProcess));

    DEBUG_MSG("Adding wrapper for void G4VProcess::BuildWorkerPhysicsTable(const G4ParticleDefinition &) (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::BuildWorkerPhysicsTable(const G4ParticleDefinition &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:296:18
    t.method("BuildWorkerPhysicsTable", static_cast<void (G4VProcess::*)(const G4ParticleDefinition &) >(&G4VProcess::BuildWorkerPhysicsTable));

    DEBUG_MSG("Adding wrapper for void G4VProcess::PrepareWorkerPhysicsTable(const G4ParticleDefinition &) (" __HERE__ ")");
    // signature to use in the veto list: void G4VProcess::PrepareWorkerPhysicsTable(const G4ParticleDefinition &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4VProcess.hh:304:18
    t.method("PrepareWorkerPhysicsTable", static_cast<void (G4VProcess::*)(const G4ParticleDefinition &) >(&G4VProcess::PrepareWorkerPhysicsTable));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4VProcess>> type_;
};
std::shared_ptr<Wrapper> newJlG4VProcess(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4VProcess(module));
}
