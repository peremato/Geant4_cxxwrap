// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4FastSimulationManager> : std::false_type { };
  template<> struct DefaultConstructible<G4FastSimulationManager> : std::false_type { };
}

// Class generating the wrapper for type G4FastSimulationManager
// signature to use in the veto file: G4FastSimulationManager
struct JlG4FastSimulationManager: public Wrapper {

  JlG4FastSimulationManager(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4FastSimulationManager (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:73:7
    jlcxx::TypeWrapper<G4FastSimulationManager>  t = jlModule.add_type<G4FastSimulationManager>("G4FastSimulationManager");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4FastSimulationManager>>(new jlcxx::TypeWrapper<G4FastSimulationManager>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::G4FastSimulationManager(G4Envelope *, G4bool) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:86:5
    t.constructor<G4Envelope *>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<G4Envelope *, G4bool>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::AddFastSimulationModel(G4VFastSimulationModel *) (" __HERE__ ")");
    // signature to use in the veto list: void G4FastSimulationManager::AddFastSimulationModel(G4VFastSimulationModel *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:104:10
    t.method("AddFastSimulationModel", static_cast<void (G4FastSimulationManager::*)(G4VFastSimulationModel *) >(&G4FastSimulationManager::AddFastSimulationModel));

    DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::RemoveFastSimulationModel(G4VFastSimulationModel *) (" __HERE__ ")");
    // signature to use in the veto list: void G4FastSimulationManager::RemoveFastSimulationModel(G4VFastSimulationModel *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:107:10
    t.method("RemoveFastSimulationModel", static_cast<void (G4FastSimulationManager::*)(G4VFastSimulationModel *) >(&G4FastSimulationManager::RemoveFastSimulationModel));

    DEBUG_MSG("Adding wrapper for G4bool G4FastSimulationManager::ActivateFastSimulationModel(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4FastSimulationManager::ActivateFastSimulationModel(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:110:12
    t.method("ActivateFastSimulationModel", static_cast<G4bool (G4FastSimulationManager::*)(const G4String &) >(&G4FastSimulationManager::ActivateFastSimulationModel));

    DEBUG_MSG("Adding wrapper for G4bool G4FastSimulationManager::InActivateFastSimulationModel(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4FastSimulationManager::InActivateFastSimulationModel(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:113:12
    t.method("InActivateFastSimulationModel", static_cast<G4bool (G4FastSimulationManager::*)(const G4String &) >(&G4FastSimulationManager::InActivateFastSimulationModel));

    DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::ListTitle() (" __HERE__ ")");
    // signature to use in the veto list: void G4FastSimulationManager::ListTitle()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:116:10
    t.method("ListTitle", static_cast<void (G4FastSimulationManager::*)()  const>(&G4FastSimulationManager::ListTitle));

    DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::ListModels() (" __HERE__ ")");
    // signature to use in the veto list: void G4FastSimulationManager::ListModels()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:117:10
    t.method("ListModels", static_cast<void (G4FastSimulationManager::*)()  const>(&G4FastSimulationManager::ListModels));

    DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::ListModels(const G4ParticleDefinition *) (" __HERE__ ")");
    // signature to use in the veto list: void G4FastSimulationManager::ListModels(const G4ParticleDefinition *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:118:10
    t.method("ListModels", static_cast<void (G4FastSimulationManager::*)(const G4ParticleDefinition *)  const>(&G4FastSimulationManager::ListModels));

    DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::ListModels(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4FastSimulationManager::ListModels(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:119:10
    t.method("ListModels", static_cast<void (G4FastSimulationManager::*)(const G4String &)  const>(&G4FastSimulationManager::ListModels));

    DEBUG_MSG("Adding wrapper for const G4Envelope * G4FastSimulationManager::GetEnvelope() (" __HERE__ ")");
    // signature to use in the veto list: const G4Envelope * G4FastSimulationManager::GetEnvelope()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:120:23
    t.method("GetEnvelope", static_cast<const G4Envelope * (G4FastSimulationManager::*)()  const>(&G4FastSimulationManager::GetEnvelope));

    DEBUG_MSG("Adding wrapper for G4VFastSimulationModel * G4FastSimulationManager::GetFastSimulationModel(const G4String &, const G4VFastSimulationModel *, G4bool &) (" __HERE__ ")");
    // signature to use in the veto list: G4VFastSimulationModel * G4FastSimulationManager::GetFastSimulationModel(const G4String &, const G4VFastSimulationModel *, G4bool &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:122:29
    t.method("GetFastSimulationModel", static_cast<G4VFastSimulationModel * (G4FastSimulationManager::*)(const G4String &, const G4VFastSimulationModel *, G4bool &)  const>(&G4FastSimulationManager::GetFastSimulationModel));

    DEBUG_MSG("Adding wrapper for const std::vector<G4VFastSimulationModel *> & G4FastSimulationManager::GetFastSimulationModelList() (" __HERE__ ")");
    // signature to use in the veto list: const std::vector<G4VFastSimulationModel *> & G4FastSimulationManager::GetFastSimulationModelList()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:126:49
    t.method("GetFastSimulationModelList", static_cast<const std::vector<G4VFastSimulationModel *> & (G4FastSimulationManager::*)()  const>(&G4FastSimulationManager::GetFastSimulationModelList));

    DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::FlushModels() (" __HERE__ ")");
    // signature to use in the veto list: void G4FastSimulationManager::FlushModels()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:131:10
    t.method("FlushModels", static_cast<void (G4FastSimulationManager::*)() >(&G4FastSimulationManager::FlushModels));

    DEBUG_MSG("Adding wrapper for G4bool G4FastSimulationManager::PostStepGetFastSimulationManagerTrigger(const G4Track &, const G4Navigator *) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4FastSimulationManager::PostStepGetFastSimulationManagerTrigger(const G4Track &, const G4Navigator *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:138:12
    t.method("PostStepGetFastSimulationManagerTrigger", static_cast<G4bool (G4FastSimulationManager::*)(const G4Track &, const G4Navigator *) >(&G4FastSimulationManager::PostStepGetFastSimulationManagerTrigger));
    t.method("PostStepGetFastSimulationManagerTrigger", [](G4FastSimulationManager& a, const G4Track & arg0)->G4bool { return a.PostStepGetFastSimulationManagerTrigger(arg0); });
    t.method("PostStepGetFastSimulationManagerTrigger", [](G4FastSimulationManager* a, const G4Track & arg0)->G4bool { return a->PostStepGetFastSimulationManagerTrigger(arg0); });

    DEBUG_MSG("Adding wrapper for G4VParticleChange * G4FastSimulationManager::InvokePostStepDoIt() (" __HERE__ ")");
    // signature to use in the veto list: G4VParticleChange * G4FastSimulationManager::InvokePostStepDoIt()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:140:24
    t.method("InvokePostStepDoIt", static_cast<G4VParticleChange * (G4FastSimulationManager::*)() >(&G4FastSimulationManager::InvokePostStepDoIt));

    DEBUG_MSG("Adding wrapper for G4bool G4FastSimulationManager::AtRestGetFastSimulationManagerTrigger(const G4Track &, const G4Navigator *) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4FastSimulationManager::AtRestGetFastSimulationManagerTrigger(const G4Track &, const G4Navigator *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:143:12
    t.method("AtRestGetFastSimulationManagerTrigger", static_cast<G4bool (G4FastSimulationManager::*)(const G4Track &, const G4Navigator *) >(&G4FastSimulationManager::AtRestGetFastSimulationManagerTrigger));
    t.method("AtRestGetFastSimulationManagerTrigger", [](G4FastSimulationManager& a, const G4Track & arg0)->G4bool { return a.AtRestGetFastSimulationManagerTrigger(arg0); });
    t.method("AtRestGetFastSimulationManagerTrigger", [](G4FastSimulationManager* a, const G4Track & arg0)->G4bool { return a->AtRestGetFastSimulationManagerTrigger(arg0); });

    DEBUG_MSG("Adding wrapper for G4VParticleChange * G4FastSimulationManager::InvokeAtRestDoIt() (" __HERE__ ")");
    // signature to use in the veto list: G4VParticleChange * G4FastSimulationManager::InvokeAtRestDoIt()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:144:24
    t.method("InvokeAtRestDoIt", static_cast<G4VParticleChange * (G4FastSimulationManager::*)() >(&G4FastSimulationManager::InvokeAtRestDoIt));
    module_.set_override_module(jl_base_module);

    DEBUG_MSG("Adding wrapper for G4bool G4FastSimulationManager::operator==(const G4FastSimulationManager &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4FastSimulationManager::operator==(const G4FastSimulationManager &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4FastSimulationManager.hh:147:12
    t.method("==", static_cast<G4bool (G4FastSimulationManager::*)(const G4FastSimulationManager &)  const>(&G4FastSimulationManager::operator==));

    module_.unset_override_module();
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4FastSimulationManager>> type_;
};
std::shared_ptr<Wrapper> newJlG4FastSimulationManager(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4FastSimulationManager(module));
}
