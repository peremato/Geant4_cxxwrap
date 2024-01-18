// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4MTRunManager> : std::false_type { };
  template<> struct DefaultConstructible<G4MTRunManager> : std::false_type { };
template<> struct SuperType<G4MTRunManager> { typedef G4RunManager type; };
}

// Class generating the wrapper for type G4MTRunManager
// signature to use in the veto file: G4MTRunManager
struct JlG4MTRunManager: public Wrapper {

  JlG4MTRunManager(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4MTRunManager (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:58:7
    jlcxx::TypeWrapper<G4MTRunManager>  t = jlModule.add_type<G4MTRunManager>("G4MTRunManager",
      jlcxx::julia_base_type<G4RunManager>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4MTRunManager>>(new jlcxx::TypeWrapper<G4MTRunManager>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetNumberOfThreads(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetNumberOfThreads(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:74:10
    t.method("SetNumberOfThreads", static_cast<void (G4MTRunManager::*)(G4int) >(&G4MTRunManager::SetNumberOfThreads));

    DEBUG_MSG("Adding wrapper for G4int G4MTRunManager::GetNumberOfThreads() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4MTRunManager::GetNumberOfThreads()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:75:11
    t.method("GetNumberOfThreads", static_cast<G4int (G4MTRunManager::*)()  const>(&G4MTRunManager::GetNumberOfThreads));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetPinAffinity(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetPinAffinity(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:76:10
    t.method("SetPinAffinity", static_cast<void (G4MTRunManager::*)(G4int) >(&G4MTRunManager::SetPinAffinity));
    t.method("SetPinAffinity", [](G4MTRunManager& a)->void { a.SetPinAffinity(); });
    t.method("SetPinAffinity", [](G4MTRunManager* a)->void { a->SetPinAffinity(); });

    DEBUG_MSG("Adding wrapper for G4int G4MTRunManager::GetPinAffinity() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4MTRunManager::GetPinAffinity()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:77:18
    t.method("GetPinAffinity", static_cast<G4int (G4MTRunManager::*)()  const>(&G4MTRunManager::GetPinAffinity));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::Initialize() (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::Initialize()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:80:10
    t.method("Initialize", static_cast<void (G4MTRunManager::*)() >(&G4MTRunManager::Initialize));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::InitializeEventLoop(G4int, const char *, G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::InitializeEventLoop(G4int, const char *, G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:81:10
    t.method("InitializeEventLoop", static_cast<void (G4MTRunManager::*)(G4int, const char *, G4int) >(&G4MTRunManager::InitializeEventLoop));
    t.method("InitializeEventLoop", [](G4MTRunManager& a, G4int arg0)->void { a.InitializeEventLoop(arg0); });
    t.method("InitializeEventLoop", [](G4MTRunManager& a, G4int arg0, const char * arg1)->void { a.InitializeEventLoop(arg0, arg1); });
    t.method("InitializeEventLoop", [](G4MTRunManager* a, G4int arg0)->void { a->InitializeEventLoop(arg0); });
    t.method("InitializeEventLoop", [](G4MTRunManager* a, G4int arg0, const char * arg1)->void { a->InitializeEventLoop(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::InitializeThreadPool() (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::InitializeThreadPool()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:83:18
    t.method("InitializeThreadPool", static_cast<void (G4MTRunManager::*)() >(&G4MTRunManager::InitializeThreadPool));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::TerminateOneEvent() (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::TerminateOneEvent()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:86:10
    t.method("TerminateOneEvent", static_cast<void (G4MTRunManager::*)() >(&G4MTRunManager::TerminateOneEvent));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::ProcessOneEvent(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::ProcessOneEvent(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:87:10
    t.method("ProcessOneEvent", static_cast<void (G4MTRunManager::*)(G4int) >(&G4MTRunManager::ProcessOneEvent));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::ConstructScoringWorlds() (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::ConstructScoringWorlds()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:88:10
    t.method("ConstructScoringWorlds", static_cast<void (G4MTRunManager::*)() >(&G4MTRunManager::ConstructScoringWorlds));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::RunTermination() (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::RunTermination()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:89:10
    t.method("RunTermination", static_cast<void (G4MTRunManager::*)() >(&G4MTRunManager::RunTermination));

    DEBUG_MSG("Adding wrapper for G4bool G4MTRunManager::SetUpAnEvent(G4Event *, G4long &, G4long &, G4long &, G4bool) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4MTRunManager::SetUpAnEvent(G4Event *, G4long &, G4long &, G4long &, G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:98:20
    t.method("SetUpAnEvent", static_cast<G4bool (G4MTRunManager::*)(G4Event *, G4long &, G4long &, G4long &, G4bool) >(&G4MTRunManager::SetUpAnEvent));
    t.method("SetUpAnEvent", [](G4MTRunManager& a, G4Event * arg0, G4long & arg1, G4long & arg2, G4long & arg3)->G4bool { return a.SetUpAnEvent(arg0, arg1, arg2, arg3); });
    t.method("SetUpAnEvent", [](G4MTRunManager* a, G4Event * arg0, G4long & arg1, G4long & arg2, G4long & arg3)->G4bool { return a->SetUpAnEvent(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for std::vector<G4String> G4MTRunManager::GetCommandStack() (" __HERE__ ")");
    // signature to use in the veto list: std::vector<G4String> G4MTRunManager::GetCommandStack()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:113:27
    t.method("GetCommandStack", static_cast<std::vector<G4String> (G4MTRunManager::*)() >(&G4MTRunManager::GetCommandStack));

    DEBUG_MSG("Adding wrapper for size_t G4MTRunManager::GetNumberActiveThreads() (" __HERE__ ")");
    // signature to use in the veto list: size_t G4MTRunManager::GetNumberActiveThreads()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:119:20
    t.method("GetNumberActiveThreads", static_cast<size_t (G4MTRunManager::*)()  const>(&G4MTRunManager::GetNumberActiveThreads));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::ThisWorkerReady() (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::ThisWorkerReady()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:126:18
    t.method("ThisWorkerReady", static_cast<void (G4MTRunManager::*)() >(&G4MTRunManager::ThisWorkerReady));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::ThisWorkerEndEventLoop() (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::ThisWorkerEndEventLoop()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:130:18
    t.method("ThisWorkerEndEventLoop", static_cast<void (G4MTRunManager::*)() >(&G4MTRunManager::ThisWorkerEndEventLoop));

    DEBUG_MSG("Adding wrapper for G4ScoringManager * G4MTRunManager::GetMasterScoringManager() (" __HERE__ ")");
    // signature to use in the veto list: G4ScoringManager * G4MTRunManager::GetMasterScoringManager()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:132:30
    t.method("G4MTRunManager!GetMasterScoringManager", static_cast<G4ScoringManager * (*)() >(&G4MTRunManager::GetMasterScoringManager));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::addWorld(G4int, G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::addWorld(G4int, G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:134:17
    t.method("G4MTRunManager!addWorld", static_cast<void (*)(G4int, G4VPhysicalVolume *) >(&G4MTRunManager::addWorld));

    DEBUG_MSG("Adding wrapper for const CLHEP::HepRandomEngine * G4MTRunManager::getMasterRandomEngine() (" __HERE__ ")");
    // signature to use in the veto list: const CLHEP::HepRandomEngine * G4MTRunManager::getMasterRandomEngine()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:136:42
    t.method("getMasterRandomEngine", static_cast<const CLHEP::HepRandomEngine * (G4MTRunManager::*)()  const>(&G4MTRunManager::getMasterRandomEngine));

    DEBUG_MSG("Adding wrapper for G4MTRunManager * G4MTRunManager::GetMasterRunManager() (" __HERE__ ")");
    // signature to use in the veto list: G4MTRunManager * G4MTRunManager::GetMasterRunManager()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:140:28
    t.method("G4MTRunManager!GetMasterRunManager", static_cast<G4MTRunManager * (*)() >(&G4MTRunManager::GetMasterRunManager));

    DEBUG_MSG("Adding wrapper for G4RunManagerKernel * G4MTRunManager::GetMasterRunManagerKernel() (" __HERE__ ")");
    // signature to use in the veto list: G4RunManagerKernel * G4MTRunManager::GetMasterRunManagerKernel()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:144:32
    t.method("G4MTRunManager!GetMasterRunManagerKernel", static_cast<G4RunManagerKernel * (*)() >(&G4MTRunManager::GetMasterRunManagerKernel));

    DEBUG_MSG("Adding wrapper for G4MTRunManagerKernel * G4MTRunManager::GetMTMasterRunManagerKernel() (" __HERE__ ")");
    // signature to use in the veto list: G4MTRunManagerKernel * G4MTRunManager::GetMTMasterRunManagerKernel()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:145:34
    t.method("G4MTRunManager!GetMTMasterRunManagerKernel", static_cast<G4MTRunManagerKernel * (*)() >(&G4MTRunManager::GetMTMasterRunManagerKernel));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetUserInitialization(G4VUserPhysicsList *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetUserInitialization(G4VUserPhysicsList *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:147:10
    t.method("SetUserInitialization", static_cast<void (G4MTRunManager::*)(G4VUserPhysicsList *) >(&G4MTRunManager::SetUserInitialization));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetUserInitialization(G4VUserDetectorConstruction *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetUserInitialization(G4VUserDetectorConstruction *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:148:10
    t.method("SetUserInitialization", static_cast<void (G4MTRunManager::*)(G4VUserDetectorConstruction *) >(&G4MTRunManager::SetUserInitialization));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetUserInitialization(G4UserWorkerInitialization *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetUserInitialization(G4UserWorkerInitialization *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:149:10
    t.method("SetUserInitialization", static_cast<void (G4MTRunManager::*)(G4UserWorkerInitialization *) >(&G4MTRunManager::SetUserInitialization));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetUserInitialization(G4VUserActionInitialization *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetUserInitialization(G4VUserActionInitialization *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:151:10
    t.method("SetUserInitialization", static_cast<void (G4MTRunManager::*)(G4VUserActionInitialization *) >(&G4MTRunManager::SetUserInitialization));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetUserAction(G4UserRunAction *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetUserAction(G4UserRunAction *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:152:10
    t.method("SetUserAction", static_cast<void (G4MTRunManager::*)(G4UserRunAction *) >(&G4MTRunManager::SetUserAction));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetUserAction(G4VUserPrimaryGeneratorAction *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetUserAction(G4VUserPrimaryGeneratorAction *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:153:10
    t.method("SetUserAction", static_cast<void (G4MTRunManager::*)(G4VUserPrimaryGeneratorAction *) >(&G4MTRunManager::SetUserAction));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetUserAction(G4UserEventAction *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetUserAction(G4UserEventAction *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:154:10
    t.method("SetUserAction", static_cast<void (G4MTRunManager::*)(G4UserEventAction *) >(&G4MTRunManager::SetUserAction));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetUserAction(G4UserStackingAction *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetUserAction(G4UserStackingAction *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:155:10
    t.method("SetUserAction", static_cast<void (G4MTRunManager::*)(G4UserStackingAction *) >(&G4MTRunManager::SetUserAction));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetUserAction(G4UserTrackingAction *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetUserAction(G4UserTrackingAction *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:156:10
    t.method("SetUserAction", static_cast<void (G4MTRunManager::*)(G4UserTrackingAction *) >(&G4MTRunManager::SetUserAction));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetUserAction(G4UserSteppingAction *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetUserAction(G4UserSteppingAction *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:157:10
    t.method("SetUserAction", static_cast<void (G4MTRunManager::*)(G4UserSteppingAction *) >(&G4MTRunManager::SetUserAction));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::MergeScores(const G4ScoringManager *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::MergeScores(const G4ScoringManager *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:160:10
    t.method("MergeScores", static_cast<void (G4MTRunManager::*)(const G4ScoringManager *) >(&G4MTRunManager::MergeScores));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::MergeRun(const G4Run *) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::MergeRun(const G4Run *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:161:10
    t.method("MergeRun", static_cast<void (G4MTRunManager::*)(const G4Run *) >(&G4MTRunManager::MergeRun));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::RequestWorkersProcessCommandsStack() (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::RequestWorkersProcessCommandsStack()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:174:18
    t.method("RequestWorkersProcessCommandsStack", static_cast<void (G4MTRunManager::*)() >(&G4MTRunManager::RequestWorkersProcessCommandsStack));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::ThisWorkerProcessCommandsStackDone() (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::ThisWorkerProcessCommandsStackDone()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:178:18
    t.method("ThisWorkerProcessCommandsStackDone", static_cast<void (G4MTRunManager::*)() >(&G4MTRunManager::ThisWorkerProcessCommandsStackDone));

    DEBUG_MSG("Adding wrapper for G4MTRunManager::WorkerActionRequest G4MTRunManager::ThisWorkerWaitForNextAction() (" __HERE__ ")");
    // signature to use in the veto list: G4MTRunManager::WorkerActionRequest G4MTRunManager::ThisWorkerWaitForNextAction()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:184:33
    t.method("ThisWorkerWaitForNextAction", static_cast<G4MTRunManager::WorkerActionRequest (G4MTRunManager::*)() >(&G4MTRunManager::ThisWorkerWaitForNextAction));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetEventModulo(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetEventModulo(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:186:17
    t.method("SetEventModulo", static_cast<void (G4MTRunManager::*)(G4int) >(&G4MTRunManager::SetEventModulo));
    t.method("SetEventModulo", [](G4MTRunManager& a)->void { a.SetEventModulo(); });
    t.method("SetEventModulo", [](G4MTRunManager* a)->void { a->SetEventModulo(); });

    DEBUG_MSG("Adding wrapper for G4int G4MTRunManager::GetEventModulo() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4MTRunManager::GetEventModulo()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:187:18
    t.method("GetEventModulo", static_cast<G4int (G4MTRunManager::*)()  const>(&G4MTRunManager::GetEventModulo));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::AbortRun(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::AbortRun(G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:189:10
    t.method("AbortRun", static_cast<void (G4MTRunManager::*)(G4bool) >(&G4MTRunManager::AbortRun));
    t.method("AbortRun", [](G4MTRunManager& a)->void { a.AbortRun(); });
    t.method("AbortRun", [](G4MTRunManager* a)->void { a->AbortRun(); });

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::AbortEvent() (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::AbortEvent()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:190:10
    t.method("AbortEvent", static_cast<void (G4MTRunManager::*)() >(&G4MTRunManager::AbortEvent));

    DEBUG_MSG("Adding wrapper for G4int G4MTRunManager::SeedOncePerCommunication() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4MTRunManager::SeedOncePerCommunication()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:192:18
    t.method("G4MTRunManager!SeedOncePerCommunication", static_cast<G4int (*)() >(&G4MTRunManager::SeedOncePerCommunication));

    DEBUG_MSG("Adding wrapper for void G4MTRunManager::SetSeedOncePerCommunication(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4MTRunManager::SetSeedOncePerCommunication(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4MTRunManager.hh:193:17
    t.method("G4MTRunManager!SetSeedOncePerCommunication", static_cast<void (*)(G4int) >(&G4MTRunManager::SetSeedOncePerCommunication));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4MTRunManager>> type_;
};
std::shared_ptr<Wrapper> newJlG4MTRunManager(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4MTRunManager(module));
}