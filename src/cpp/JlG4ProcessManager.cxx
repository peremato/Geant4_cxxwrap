// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4ProcessManager> : std::false_type { };
  template<> struct DefaultConstructible<G4ProcessManager> : std::false_type { };
}

// Class generating the wrapper for type G4ProcessManager
// signature to use in the veto file: G4ProcessManager
struct JlG4ProcessManager: public Wrapper {

  JlG4ProcessManager(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4ProcessManager (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:94:7
    jlcxx::TypeWrapper<G4ProcessManager>  t = jlModule.add_type<G4ProcessManager>("G4ProcessManager");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4ProcessManager>>(new jlcxx::TypeWrapper<G4ProcessManager>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4ProcessManager::G4ProcessManager(const G4ParticleDefinition *) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:100:5
    t.constructor<const G4ParticleDefinition *>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4ProcessManager::G4ProcessManager(G4ProcessManager &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:103:5
    t.constructor<G4ProcessManager &>(/*finalize=*/true);
    module_.set_override_module(jl_base_module);

    DEBUG_MSG("Adding wrapper for G4bool G4ProcessManager::operator==(const G4ProcessManager &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4ProcessManager::operator==(const G4ProcessManager &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:113:12
    t.method("==", static_cast<G4bool (G4ProcessManager::*)(const G4ProcessManager &)  const>(&G4ProcessManager::operator==));

    DEBUG_MSG("Adding wrapper for G4bool G4ProcessManager::operator!=(const G4ProcessManager &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4ProcessManager::operator!=(const G4ProcessManager &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:114:12
    t.method("!=", static_cast<G4bool (G4ProcessManager::*)(const G4ProcessManager &)  const>(&G4ProcessManager::operator!=));

    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for G4ProcessVector * G4ProcessManager::GetProcessList() (" __HERE__ ")");
    // signature to use in the veto list: G4ProcessVector * G4ProcessManager::GetProcessList()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:116:29
    t.method("GetProcessList", static_cast<G4ProcessVector * (G4ProcessManager::*)()  const>(&G4ProcessManager::GetProcessList));

    DEBUG_MSG("Adding wrapper for G4int G4ProcessManager::GetProcessListLength() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ProcessManager::GetProcessListLength()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:119:18
    t.method("GetProcessListLength", static_cast<G4int (G4ProcessManager::*)()  const>(&G4ProcessManager::GetProcessListLength));

    DEBUG_MSG("Adding wrapper for G4int G4ProcessManager::GetProcessIndex(G4VProcess *) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ProcessManager::GetProcessIndex(G4VProcess *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:122:18
    t.method("GetProcessIndex", static_cast<G4int (G4ProcessManager::*)(G4VProcess *)  const>(&G4ProcessManager::GetProcessIndex));

    DEBUG_MSG("Adding wrapper for G4ProcessVector * G4ProcessManager::GetProcessVector(G4ProcessVectorDoItIndex, G4ProcessVectorTypeIndex) (" __HERE__ ")");
    // signature to use in the veto list: G4ProcessVector * G4ProcessManager::GetProcessVector(G4ProcessVectorDoItIndex, G4ProcessVectorTypeIndex)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:125:29
    t.method("GetProcessVector", static_cast<G4ProcessVector * (G4ProcessManager::*)(G4ProcessVectorDoItIndex, G4ProcessVectorTypeIndex)  const>(&G4ProcessManager::GetProcessVector));
    t.method("GetProcessVector", [](G4ProcessManager const& a, G4ProcessVectorDoItIndex arg0)->G4ProcessVector * { return a.GetProcessVector(arg0); });
    t.method("GetProcessVector", [](G4ProcessManager const* a, G4ProcessVectorDoItIndex arg0)->G4ProcessVector * { return a->GetProcessVector(arg0); });

    DEBUG_MSG("Adding wrapper for G4ProcessVector * G4ProcessManager::GetAtRestProcessVector(G4ProcessVectorTypeIndex) (" __HERE__ ")");
    // signature to use in the veto list: G4ProcessVector * G4ProcessManager::GetAtRestProcessVector(G4ProcessVectorTypeIndex)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:131:29
    t.method("GetAtRestProcessVector", static_cast<G4ProcessVector * (G4ProcessManager::*)(G4ProcessVectorTypeIndex)  const>(&G4ProcessManager::GetAtRestProcessVector));
    t.method("GetAtRestProcessVector", [](G4ProcessManager const& a)->G4ProcessVector * { return a.GetAtRestProcessVector(); });
    t.method("GetAtRestProcessVector", [](G4ProcessManager const* a)->G4ProcessVector * { return a->GetAtRestProcessVector(); });

    DEBUG_MSG("Adding wrapper for G4ProcessVector * G4ProcessManager::GetAlongStepProcessVector(G4ProcessVectorTypeIndex) (" __HERE__ ")");
    // signature to use in the veto list: G4ProcessVector * G4ProcessManager::GetAlongStepProcessVector(G4ProcessVectorTypeIndex)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:138:29
    t.method("GetAlongStepProcessVector", static_cast<G4ProcessVector * (G4ProcessManager::*)(G4ProcessVectorTypeIndex)  const>(&G4ProcessManager::GetAlongStepProcessVector));
    t.method("GetAlongStepProcessVector", [](G4ProcessManager const& a)->G4ProcessVector * { return a.GetAlongStepProcessVector(); });
    t.method("GetAlongStepProcessVector", [](G4ProcessManager const* a)->G4ProcessVector * { return a->GetAlongStepProcessVector(); });

    DEBUG_MSG("Adding wrapper for G4ProcessVector * G4ProcessManager::GetPostStepProcessVector(G4ProcessVectorTypeIndex) (" __HERE__ ")");
    // signature to use in the veto list: G4ProcessVector * G4ProcessManager::GetPostStepProcessVector(G4ProcessVectorTypeIndex)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:145:29
    t.method("GetPostStepProcessVector", static_cast<G4ProcessVector * (G4ProcessManager::*)(G4ProcessVectorTypeIndex)  const>(&G4ProcessManager::GetPostStepProcessVector));
    t.method("GetPostStepProcessVector", [](G4ProcessManager const& a)->G4ProcessVector * { return a.GetPostStepProcessVector(); });
    t.method("GetPostStepProcessVector", [](G4ProcessManager const* a)->G4ProcessVector * { return a->GetPostStepProcessVector(); });

    DEBUG_MSG("Adding wrapper for G4int G4ProcessManager::GetProcessVectorIndex(G4VProcess *, G4ProcessVectorDoItIndex, G4ProcessVectorTypeIndex) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ProcessManager::GetProcessVectorIndex(G4VProcess *, G4ProcessVectorDoItIndex, G4ProcessVectorTypeIndex)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:152:11
    t.method("GetProcessVectorIndex", static_cast<G4int (G4ProcessManager::*)(G4VProcess *, G4ProcessVectorDoItIndex, G4ProcessVectorTypeIndex)  const>(&G4ProcessManager::GetProcessVectorIndex));
    t.method("GetProcessVectorIndex", [](G4ProcessManager const& a, G4VProcess * arg0, G4ProcessVectorDoItIndex arg1)->G4int { return a.GetProcessVectorIndex(arg0, arg1); });
    t.method("GetProcessVectorIndex", [](G4ProcessManager const* a, G4VProcess * arg0, G4ProcessVectorDoItIndex arg1)->G4int { return a->GetProcessVectorIndex(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for G4int G4ProcessManager::GetAtRestIndex(G4VProcess *, G4ProcessVectorTypeIndex) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ProcessManager::GetAtRestIndex(G4VProcess *, G4ProcessVectorTypeIndex)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:157:18
    t.method("GetAtRestIndex", static_cast<G4int (G4ProcessManager::*)(G4VProcess *, G4ProcessVectorTypeIndex)  const>(&G4ProcessManager::GetAtRestIndex));
    t.method("GetAtRestIndex", [](G4ProcessManager const& a, G4VProcess * arg0)->G4int { return a.GetAtRestIndex(arg0); });
    t.method("GetAtRestIndex", [](G4ProcessManager const* a, G4VProcess * arg0)->G4int { return a->GetAtRestIndex(arg0); });

    DEBUG_MSG("Adding wrapper for G4int G4ProcessManager::GetAlongStepIndex(G4VProcess *, G4ProcessVectorTypeIndex) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ProcessManager::GetAlongStepIndex(G4VProcess *, G4ProcessVectorTypeIndex)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:161:18
    t.method("GetAlongStepIndex", static_cast<G4int (G4ProcessManager::*)(G4VProcess *, G4ProcessVectorTypeIndex)  const>(&G4ProcessManager::GetAlongStepIndex));
    t.method("GetAlongStepIndex", [](G4ProcessManager const& a, G4VProcess * arg0)->G4int { return a.GetAlongStepIndex(arg0); });
    t.method("GetAlongStepIndex", [](G4ProcessManager const* a, G4VProcess * arg0)->G4int { return a->GetAlongStepIndex(arg0); });

    DEBUG_MSG("Adding wrapper for G4int G4ProcessManager::GetPostStepIndex(G4VProcess *, G4ProcessVectorTypeIndex) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ProcessManager::GetPostStepIndex(G4VProcess *, G4ProcessVectorTypeIndex)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:165:18
    t.method("GetPostStepIndex", static_cast<G4int (G4ProcessManager::*)(G4VProcess *, G4ProcessVectorTypeIndex)  const>(&G4ProcessManager::GetPostStepIndex));
    t.method("GetPostStepIndex", [](G4ProcessManager const& a, G4VProcess * arg0)->G4int { return a.GetPostStepIndex(arg0); });
    t.method("GetPostStepIndex", [](G4ProcessManager const* a, G4VProcess * arg0)->G4int { return a->GetPostStepIndex(arg0); });

    DEBUG_MSG("Adding wrapper for G4int G4ProcessManager::AddProcess(G4VProcess *, G4int, G4int, G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ProcessManager::AddProcess(G4VProcess *, G4int, G4int, G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:171:11
    t.method("AddProcess", static_cast<G4int (G4ProcessManager::*)(G4VProcess *, G4int, G4int, G4int) >(&G4ProcessManager::AddProcess));
    t.method("AddProcess", [](G4ProcessManager& a, G4VProcess * arg0)->G4int { return a.AddProcess(arg0); });
    t.method("AddProcess", [](G4ProcessManager& a, G4VProcess * arg0, G4int arg1)->G4int { return a.AddProcess(arg0, arg1); });
    t.method("AddProcess", [](G4ProcessManager& a, G4VProcess * arg0, G4int arg1, G4int arg2)->G4int { return a.AddProcess(arg0, arg1, arg2); });
    t.method("AddProcess", [](G4ProcessManager* a, G4VProcess * arg0)->G4int { return a->AddProcess(arg0); });
    t.method("AddProcess", [](G4ProcessManager* a, G4VProcess * arg0, G4int arg1)->G4int { return a->AddProcess(arg0, arg1); });
    t.method("AddProcess", [](G4ProcessManager* a, G4VProcess * arg0, G4int arg1, G4int arg2)->G4int { return a->AddProcess(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for G4int G4ProcessManager::AddRestProcess(G4VProcess *, G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ProcessManager::AddRestProcess(G4VProcess *, G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:196:18
    t.method("AddRestProcess", static_cast<G4int (G4ProcessManager::*)(G4VProcess *, G4int) >(&G4ProcessManager::AddRestProcess));
    t.method("AddRestProcess", [](G4ProcessManager& a, G4VProcess * arg0)->G4int { return a.AddRestProcess(arg0); });
    t.method("AddRestProcess", [](G4ProcessManager* a, G4VProcess * arg0)->G4int { return a->AddRestProcess(arg0); });

    DEBUG_MSG("Adding wrapper for G4int G4ProcessManager::AddDiscreteProcess(G4VProcess *, G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ProcessManager::AddDiscreteProcess(G4VProcess *, G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:197:18
    t.method("AddDiscreteProcess", static_cast<G4int (G4ProcessManager::*)(G4VProcess *, G4int) >(&G4ProcessManager::AddDiscreteProcess));
    t.method("AddDiscreteProcess", [](G4ProcessManager& a, G4VProcess * arg0)->G4int { return a.AddDiscreteProcess(arg0); });
    t.method("AddDiscreteProcess", [](G4ProcessManager* a, G4VProcess * arg0)->G4int { return a->AddDiscreteProcess(arg0); });

    DEBUG_MSG("Adding wrapper for G4int G4ProcessManager::AddContinuousProcess(G4VProcess *, G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ProcessManager::AddContinuousProcess(G4VProcess *, G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:198:18
    t.method("AddContinuousProcess", static_cast<G4int (G4ProcessManager::*)(G4VProcess *, G4int) >(&G4ProcessManager::AddContinuousProcess));
    t.method("AddContinuousProcess", [](G4ProcessManager& a, G4VProcess * arg0)->G4int { return a.AddContinuousProcess(arg0); });
    t.method("AddContinuousProcess", [](G4ProcessManager* a, G4VProcess * arg0)->G4int { return a->AddContinuousProcess(arg0); });

    DEBUG_MSG("Adding wrapper for G4int G4ProcessManager::GetProcessOrdering(G4VProcess *, G4ProcessVectorDoItIndex) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ProcessManager::GetProcessOrdering(G4VProcess *, G4ProcessVectorDoItIndex)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:205:11
    t.method("GetProcessOrdering", static_cast<G4int (G4ProcessManager::*)(G4VProcess *, G4ProcessVectorDoItIndex) >(&G4ProcessManager::GetProcessOrdering));

    DEBUG_MSG("Adding wrapper for void G4ProcessManager::SetProcessOrdering(G4VProcess *, G4ProcessVectorDoItIndex, G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4ProcessManager::SetProcessOrdering(G4VProcess *, G4ProcessVectorDoItIndex, G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:210:10
    t.method("SetProcessOrdering", static_cast<void (G4ProcessManager::*)(G4VProcess *, G4ProcessVectorDoItIndex, G4int) >(&G4ProcessManager::SetProcessOrdering));
    t.method("SetProcessOrdering", [](G4ProcessManager& a, G4VProcess * arg0, G4ProcessVectorDoItIndex arg1)->void { a.SetProcessOrdering(arg0, arg1); });
    t.method("SetProcessOrdering", [](G4ProcessManager* a, G4VProcess * arg0, G4ProcessVectorDoItIndex arg1)->void { a->SetProcessOrdering(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void G4ProcessManager::SetProcessOrderingToFirst(G4VProcess *, G4ProcessVectorDoItIndex) (" __HERE__ ")");
    // signature to use in the veto list: void G4ProcessManager::SetProcessOrderingToFirst(G4VProcess *, G4ProcessVectorDoItIndex)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:222:10
    t.method("SetProcessOrderingToFirst", static_cast<void (G4ProcessManager::*)(G4VProcess *, G4ProcessVectorDoItIndex) >(&G4ProcessManager::SetProcessOrderingToFirst));

    DEBUG_MSG("Adding wrapper for void G4ProcessManager::SetProcessOrderingToSecond(G4VProcess *, G4ProcessVectorDoItIndex) (" __HERE__ ")");
    // signature to use in the veto list: void G4ProcessManager::SetProcessOrderingToSecond(G4VProcess *, G4ProcessVectorDoItIndex)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:231:10
    t.method("SetProcessOrderingToSecond", static_cast<void (G4ProcessManager::*)(G4VProcess *, G4ProcessVectorDoItIndex) >(&G4ProcessManager::SetProcessOrderingToSecond));

    DEBUG_MSG("Adding wrapper for void G4ProcessManager::SetProcessOrderingToLast(G4VProcess *, G4ProcessVectorDoItIndex) (" __HERE__ ")");
    // signature to use in the veto list: void G4ProcessManager::SetProcessOrderingToLast(G4VProcess *, G4ProcessVectorDoItIndex)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:241:10
    t.method("SetProcessOrderingToLast", static_cast<void (G4ProcessManager::*)(G4VProcess *, G4ProcessVectorDoItIndex) >(&G4ProcessManager::SetProcessOrderingToLast));

    DEBUG_MSG("Adding wrapper for G4VProcess * G4ProcessManager::RemoveProcess(G4VProcess *) (" __HERE__ ")");
    // signature to use in the veto list: G4VProcess * G4ProcessManager::RemoveProcess(G4VProcess *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:252:17
    t.method("RemoveProcess", static_cast<G4VProcess * (G4ProcessManager::*)(G4VProcess *) >(&G4ProcessManager::RemoveProcess));

    DEBUG_MSG("Adding wrapper for G4VProcess * G4ProcessManager::RemoveProcess(G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4VProcess * G4ProcessManager::RemoveProcess(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:253:17
    t.method("RemoveProcess", static_cast<G4VProcess * (G4ProcessManager::*)(G4int) >(&G4ProcessManager::RemoveProcess));

    DEBUG_MSG("Adding wrapper for G4VProcess * G4ProcessManager::SetProcessActivation(G4VProcess *, G4bool) (" __HERE__ ")");
    // signature to use in the veto list: G4VProcess * G4ProcessManager::SetProcessActivation(G4VProcess *, G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:258:17
    t.method("SetProcessActivation", static_cast<G4VProcess * (G4ProcessManager::*)(G4VProcess *, G4bool) >(&G4ProcessManager::SetProcessActivation));

    DEBUG_MSG("Adding wrapper for G4VProcess * G4ProcessManager::SetProcessActivation(G4int, G4bool) (" __HERE__ ")");
    // signature to use in the veto list: G4VProcess * G4ProcessManager::SetProcessActivation(G4int, G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:259:17
    t.method("SetProcessActivation", static_cast<G4VProcess * (G4ProcessManager::*)(G4int, G4bool) >(&G4ProcessManager::SetProcessActivation));

    DEBUG_MSG("Adding wrapper for G4bool G4ProcessManager::GetProcessActivation(G4VProcess *) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4ProcessManager::GetProcessActivation(G4VProcess *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:264:12
    t.method("GetProcessActivation", static_cast<G4bool (G4ProcessManager::*)(G4VProcess *)  const>(&G4ProcessManager::GetProcessActivation));

    DEBUG_MSG("Adding wrapper for G4bool G4ProcessManager::GetProcessActivation(G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4ProcessManager::GetProcessActivation(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:265:12
    t.method("GetProcessActivation", static_cast<G4bool (G4ProcessManager::*)(G4int)  const>(&G4ProcessManager::GetProcessActivation));

    DEBUG_MSG("Adding wrapper for G4ParticleDefinition * G4ProcessManager::GetParticleType() (" __HERE__ ")");
    // signature to use in the veto list: G4ParticleDefinition * G4ProcessManager::GetParticleType()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:268:34
    t.method("GetParticleType", static_cast<G4ParticleDefinition * (G4ProcessManager::*)()  const>(&G4ProcessManager::GetParticleType));

    DEBUG_MSG("Adding wrapper for void G4ProcessManager::SetParticleType(const G4ParticleDefinition *) (" __HERE__ ")");
    // signature to use in the veto list: void G4ProcessManager::SetParticleType(const G4ParticleDefinition *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:270:17
    t.method("SetParticleType", static_cast<void (G4ProcessManager::*)(const G4ParticleDefinition *) >(&G4ProcessManager::SetParticleType));

    DEBUG_MSG("Adding wrapper for G4VProcess * G4ProcessManager::GetProcess(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4VProcess * G4ProcessManager::GetProcess(const G4String &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:273:17
    t.method("GetProcess", static_cast<G4VProcess * (G4ProcessManager::*)(const G4String &)  const>(&G4ProcessManager::GetProcess));

    DEBUG_MSG("Adding wrapper for void G4ProcessManager::StartTracking(G4Track *) (" __HERE__ ")");
    // signature to use in the veto list: void G4ProcessManager::StartTracking(G4Track *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:276:10
    t.method("StartTracking", static_cast<void (G4ProcessManager::*)(G4Track *) >(&G4ProcessManager::StartTracking));
    t.method("StartTracking", [](G4ProcessManager& a)->void { a.StartTracking(); });
    t.method("StartTracking", [](G4ProcessManager* a)->void { a->StartTracking(); });

    DEBUG_MSG("Adding wrapper for void G4ProcessManager::EndTracking() (" __HERE__ ")");
    // signature to use in the veto list: void G4ProcessManager::EndTracking()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:277:10
    t.method("EndTracking", static_cast<void (G4ProcessManager::*)() >(&G4ProcessManager::EndTracking));

    DEBUG_MSG("Adding wrapper for void G4ProcessManager::DumpInfo() (" __HERE__ ")");
    // signature to use in the veto list: void G4ProcessManager::DumpInfo()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:282:10
    t.method("DumpInfo", static_cast<void (G4ProcessManager::*)() >(&G4ProcessManager::DumpInfo));

    DEBUG_MSG("Adding wrapper for void G4ProcessManager::SetVerboseLevel(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4ProcessManager::SetVerboseLevel(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:284:17
    t.method("SetVerboseLevel", static_cast<void (G4ProcessManager::*)(G4int) >(&G4ProcessManager::SetVerboseLevel));

    DEBUG_MSG("Adding wrapper for G4int G4ProcessManager::GetVerboseLevel() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ProcessManager::GetVerboseLevel()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ProcessManager.hh:285:18
    t.method("GetVerboseLevel", static_cast<G4int (G4ProcessManager::*)()  const>(&G4ProcessManager::GetVerboseLevel));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4ProcessManager>> type_;
};
std::shared_ptr<Wrapper> newJlG4ProcessManager(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4ProcessManager(module));
}