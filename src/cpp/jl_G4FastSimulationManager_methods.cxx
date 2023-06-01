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
void add_methods_for_G4FastSimulationManager(jlcxx::Module& types, jlcxx::TypeWrapper<G4FastSimulationManager>& t49) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4FastSimulationManager
   */


  DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::G4FastSimulationManager(G4Envelope *, G4bool) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:90:3
  t49.constructor<G4Envelope *>(/*finalize=*/true);
  t49.constructor<G4Envelope *, G4bool>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::AddFastSimulationModel(G4VFastSimulationModel *) (" __HERE__ ")");
  // signature to use in the veto list: void G4FastSimulationManager::AddFastSimulationModel(G4VFastSimulationModel *)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:113:8
  t49.method("AddFastSimulationModel", static_cast<void (G4FastSimulationManager::*)(G4VFastSimulationModel *) >(&G4FastSimulationManager::AddFastSimulationModel));

  DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::RemoveFastSimulationModel(G4VFastSimulationModel *) (" __HERE__ ")");
  // signature to use in the veto list: void G4FastSimulationManager::RemoveFastSimulationModel(G4VFastSimulationModel *)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:116:8
  t49.method("RemoveFastSimulationModel", static_cast<void (G4FastSimulationManager::*)(G4VFastSimulationModel *) >(&G4FastSimulationManager::RemoveFastSimulationModel));

  DEBUG_MSG("Adding wrapper for G4bool G4FastSimulationManager::ActivateFastSimulationModel(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4FastSimulationManager::ActivateFastSimulationModel(const G4String &)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:122:10
  t49.method("ActivateFastSimulationModel", static_cast<G4bool (G4FastSimulationManager::*)(const G4String &) >(&G4FastSimulationManager::ActivateFastSimulationModel));

  DEBUG_MSG("Adding wrapper for G4bool G4FastSimulationManager::InActivateFastSimulationModel(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4FastSimulationManager::InActivateFastSimulationModel(const G4String &)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:125:10
  t49.method("InActivateFastSimulationModel", static_cast<G4bool (G4FastSimulationManager::*)(const G4String &) >(&G4FastSimulationManager::InActivateFastSimulationModel));

  DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::ListTitle() (" __HERE__ ")");
  // signature to use in the veto list: void G4FastSimulationManager::ListTitle()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:130:8
  t49.method("ListTitle", static_cast<void (G4FastSimulationManager::*)()  const>(&G4FastSimulationManager::ListTitle));

  DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::ListModels() (" __HERE__ ")");
  // signature to use in the veto list: void G4FastSimulationManager::ListModels()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:131:8
  t49.method("ListModels", static_cast<void (G4FastSimulationManager::*)()  const>(&G4FastSimulationManager::ListModels));

  DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::ListModels(const G4ParticleDefinition *) (" __HERE__ ")");
  // signature to use in the veto list: void G4FastSimulationManager::ListModels(const G4ParticleDefinition *)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:132:8
  t49.method("ListModels", static_cast<void (G4FastSimulationManager::*)(const G4ParticleDefinition *)  const>(&G4FastSimulationManager::ListModels));

  DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::ListModels(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4FastSimulationManager::ListModels(const G4String &)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:133:8
  t49.method("ListModels", static_cast<void (G4FastSimulationManager::*)(const G4String &)  const>(&G4FastSimulationManager::ListModels));

  DEBUG_MSG("Adding wrapper for const G4Envelope * G4FastSimulationManager::GetEnvelope() (" __HERE__ ")");
  // signature to use in the veto list: const G4Envelope * G4FastSimulationManager::GetEnvelope()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:134:21
  t49.method("GetEnvelope", static_cast<const G4Envelope * (G4FastSimulationManager::*)()  const>(&G4FastSimulationManager::GetEnvelope));

  DEBUG_MSG("Adding wrapper for G4VFastSimulationModel * G4FastSimulationManager::GetFastSimulationModel(const G4String &, const G4VFastSimulationModel *, bool &) (" __HERE__ ")");
  // signature to use in the veto list: G4VFastSimulationModel * G4FastSimulationManager::GetFastSimulationModel(const G4String &, const G4VFastSimulationModel *, bool &)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:136:27
  t49.method("GetFastSimulationModel", static_cast<G4VFastSimulationModel * (G4FastSimulationManager::*)(const G4String &, const G4VFastSimulationModel *, bool &)  const>(&G4FastSimulationManager::GetFastSimulationModel));

  DEBUG_MSG("Adding wrapper for const std::vector<G4VFastSimulationModel *> & G4FastSimulationManager::GetFastSimulationModelList() (" __HERE__ ")");
  // signature to use in the veto list: const std::vector<G4VFastSimulationModel *> & G4FastSimulationManager::GetFastSimulationModelList()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:140:47
  t49.method("GetFastSimulationModelList", static_cast<const std::vector<G4VFastSimulationModel *> & (G4FastSimulationManager::*)()  const>(&G4FastSimulationManager::GetFastSimulationModelList));

  DEBUG_MSG("Adding wrapper for void G4FastSimulationManager::FlushModels() (" __HERE__ ")");
  // signature to use in the veto list: void G4FastSimulationManager::FlushModels()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:143:8
  t49.method("FlushModels", static_cast<void (G4FastSimulationManager::*)() >(&G4FastSimulationManager::FlushModels));

  DEBUG_MSG("Adding wrapper for G4bool G4FastSimulationManager::PostStepGetFastSimulationManagerTrigger(const G4Track &, const G4Navigator *) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4FastSimulationManager::PostStepGetFastSimulationManagerTrigger(const G4Track &, const G4Navigator *)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:150:10
  t49.method("PostStepGetFastSimulationManagerTrigger", static_cast<G4bool (G4FastSimulationManager::*)(const G4Track &, const G4Navigator *) >(&G4FastSimulationManager::PostStepGetFastSimulationManagerTrigger));
  t49.method("PostStepGetFastSimulationManagerTrigger", [](G4FastSimulationManager& a, const G4Track & arg0)->G4bool{ return a.PostStepGetFastSimulationManagerTrigger(arg0); });
  t49.method("PostStepGetFastSimulationManagerTrigger", [](G4FastSimulationManager* a, const G4Track & arg0)->G4bool{ return a->PostStepGetFastSimulationManagerTrigger(arg0); });

  DEBUG_MSG("Adding wrapper for G4VParticleChange * G4FastSimulationManager::InvokePostStepDoIt() (" __HERE__ ")");
  // signature to use in the veto list: G4VParticleChange * G4FastSimulationManager::InvokePostStepDoIt()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:153:22
  t49.method("InvokePostStepDoIt", static_cast<G4VParticleChange * (G4FastSimulationManager::*)() >(&G4FastSimulationManager::InvokePostStepDoIt));

  DEBUG_MSG("Adding wrapper for G4bool G4FastSimulationManager::AtRestGetFastSimulationManagerTrigger(const G4Track &, const G4Navigator *) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4FastSimulationManager::AtRestGetFastSimulationManagerTrigger(const G4Track &, const G4Navigator *)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:156:10
  t49.method("AtRestGetFastSimulationManagerTrigger", static_cast<G4bool (G4FastSimulationManager::*)(const G4Track &, const G4Navigator *) >(&G4FastSimulationManager::AtRestGetFastSimulationManagerTrigger));
  t49.method("AtRestGetFastSimulationManagerTrigger", [](G4FastSimulationManager& a, const G4Track & arg0)->G4bool{ return a.AtRestGetFastSimulationManagerTrigger(arg0); });
  t49.method("AtRestGetFastSimulationManagerTrigger", [](G4FastSimulationManager* a, const G4Track & arg0)->G4bool{ return a->AtRestGetFastSimulationManagerTrigger(arg0); });

  DEBUG_MSG("Adding wrapper for G4VParticleChange * G4FastSimulationManager::InvokeAtRestDoIt() (" __HERE__ ")");
  // signature to use in the veto list: G4VParticleChange * G4FastSimulationManager::InvokeAtRestDoIt()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:158:23
  t49.method("InvokeAtRestDoIt", static_cast<G4VParticleChange * (G4FastSimulationManager::*)() >(&G4FastSimulationManager::InvokeAtRestDoIt));
  types.set_override_module(jl_base_module);

  DEBUG_MSG("Adding wrapper for G4bool G4FastSimulationManager::operator==(const G4FastSimulationManager &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4FastSimulationManager::operator==(const G4FastSimulationManager &)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4FastSimulationManager.hh:161:10
  t49.method("==", static_cast<G4bool (G4FastSimulationManager::*)(const G4FastSimulationManager &)  const>(&G4FastSimulationManager::operator==));

  /* End of G4FastSimulationManager class method wrappers
   **********************************************************************/

}
