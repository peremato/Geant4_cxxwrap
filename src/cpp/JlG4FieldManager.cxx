// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4FieldManager> : std::false_type { };
  template<> struct DefaultConstructible<G4FieldManager> : std::false_type { };
}

// Class generating the wrapper for type G4FieldManager
// signature to use in the veto file: G4FieldManager
struct JlG4FieldManager: public Wrapper {

  JlG4FieldManager(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4FieldManager (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:84:7
    jlcxx::TypeWrapper<G4FieldManager>  t = jlModule.add_type<G4FieldManager>("G4FieldManager");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4FieldManager>>(new jlcxx::TypeWrapper<G4FieldManager>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4FieldManager::G4FieldManager(G4MagneticField *) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:93:5
    t.constructor<G4MagneticField *>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for G4bool G4FieldManager::SetDetectorField(G4Field *, G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4FieldManager::SetDetectorField(G4Field *, G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:102:12
    t.method("SetDetectorField", static_cast<G4bool (G4FieldManager::*)(G4Field *, G4int) >(&G4FieldManager::SetDetectorField));
    t.method("SetDetectorField", [](G4FieldManager& a, G4Field * arg0)->G4bool { return a.SetDetectorField(arg0); });
    t.method("SetDetectorField", [](G4FieldManager* a, G4Field * arg0)->G4bool { return a->SetDetectorField(arg0); });

    DEBUG_MSG("Adding wrapper for void G4FieldManager::ProposeDetectorField(G4Field *) (" __HERE__ ")");
    // signature to use in the veto list: void G4FieldManager::ProposeDetectorField(G4Field *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:112:17
    t.method("ProposeDetectorField", static_cast<void (G4FieldManager::*)(G4Field *) >(&G4FieldManager::ProposeDetectorField));

    DEBUG_MSG("Adding wrapper for void G4FieldManager::ChangeDetectorField(G4Field *) (" __HERE__ ")");
    // signature to use in the veto list: void G4FieldManager::ChangeDetectorField(G4Field *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:120:18
    t.method("ChangeDetectorField", static_cast<void (G4FieldManager::*)(G4Field *) >(&G4FieldManager::ChangeDetectorField));

    DEBUG_MSG("Adding wrapper for const G4Field * G4FieldManager::GetDetectorField() (" __HERE__ ")");
    // signature to use in the veto list: const G4Field * G4FieldManager::GetDetectorField()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:125:28
    t.method("GetDetectorField", static_cast<const G4Field * (G4FieldManager::*)()  const>(&G4FieldManager::GetDetectorField));

    DEBUG_MSG("Adding wrapper for G4bool G4FieldManager::DoesFieldExist() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4FieldManager::DoesFieldExist()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:126:28
    t.method("DoesFieldExist", static_cast<G4bool (G4FieldManager::*)()  const>(&G4FieldManager::DoesFieldExist));

    DEBUG_MSG("Adding wrapper for void G4FieldManager::CreateChordFinder(G4MagneticField *) (" __HERE__ ")");
    // signature to use in the veto list: void G4FieldManager::CreateChordFinder(G4MagneticField *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:129:10
    t.method("CreateChordFinder", static_cast<void (G4FieldManager::*)(G4MagneticField *) >(&G4FieldManager::CreateChordFinder));

    DEBUG_MSG("Adding wrapper for void G4FieldManager::ConfigureForTrack(const G4Track *) (" __HERE__ ")");
    // signature to use in the veto list: void G4FieldManager::ConfigureForTrack(const G4Track *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:135:20
    t.method("ConfigureForTrack", static_cast<void (G4FieldManager::*)(const G4Track *) >(&G4FieldManager::ConfigureForTrack));

    DEBUG_MSG("Adding wrapper for G4double G4FieldManager::GetDeltaIntersection() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4FieldManager::GetDeltaIntersection()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:143:21
    t.method("GetDeltaIntersection", static_cast<G4double (G4FieldManager::*)()  const>(&G4FieldManager::GetDeltaIntersection));

    DEBUG_MSG("Adding wrapper for G4double G4FieldManager::GetDeltaOneStep() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4FieldManager::GetDeltaOneStep()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:146:21
    t.method("GetDeltaOneStep", static_cast<G4double (G4FieldManager::*)()  const>(&G4FieldManager::GetDeltaOneStep));

    DEBUG_MSG("Adding wrapper for void G4FieldManager::SetAccuraciesWithDeltaOneStep(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4FieldManager::SetAccuraciesWithDeltaOneStep(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:149:17
    t.method("SetAccuraciesWithDeltaOneStep", static_cast<void (G4FieldManager::*)(G4double) >(&G4FieldManager::SetAccuraciesWithDeltaOneStep));

    DEBUG_MSG("Adding wrapper for void G4FieldManager::SetDeltaOneStep(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4FieldManager::SetDeltaOneStep(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:153:21
    t.method("SetDeltaOneStep", static_cast<void (G4FieldManager::*)(G4double) >(&G4FieldManager::SetDeltaOneStep));

    DEBUG_MSG("Adding wrapper for void G4FieldManager::SetDeltaIntersection(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4FieldManager::SetDeltaIntersection(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:155:21
    t.method("SetDeltaIntersection", static_cast<void (G4FieldManager::*)(G4double) >(&G4FieldManager::SetDeltaIntersection));

    DEBUG_MSG("Adding wrapper for G4double G4FieldManager::GetMinimumEpsilonStep() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4FieldManager::GetMinimumEpsilonStep()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:158:22
    t.method("GetMinimumEpsilonStep", static_cast<G4double (G4FieldManager::*)()  const>(&G4FieldManager::GetMinimumEpsilonStep));

    DEBUG_MSG("Adding wrapper for G4bool G4FieldManager::SetMinimumEpsilonStep(G4double) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4FieldManager::SetMinimumEpsilonStep(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:159:22
    t.method("SetMinimumEpsilonStep", static_cast<G4bool (G4FieldManager::*)(G4double) >(&G4FieldManager::SetMinimumEpsilonStep));

    DEBUG_MSG("Adding wrapper for G4double G4FieldManager::GetMaximumEpsilonStep() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4FieldManager::GetMaximumEpsilonStep()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:162:22
    t.method("GetMaximumEpsilonStep", static_cast<G4double (G4FieldManager::*)()  const>(&G4FieldManager::GetMaximumEpsilonStep));

    DEBUG_MSG("Adding wrapper for G4bool G4FieldManager::SetMaximumEpsilonStep(G4double) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4FieldManager::SetMaximumEpsilonStep(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:163:22
    t.method("SetMaximumEpsilonStep", static_cast<G4bool (G4FieldManager::*)(G4double) >(&G4FieldManager::SetMaximumEpsilonStep));

    DEBUG_MSG("Adding wrapper for G4bool G4FieldManager::DoesFieldChangeEnergy() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4FieldManager::DoesFieldChangeEnergy()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:166:21
    t.method("DoesFieldChangeEnergy", static_cast<G4bool (G4FieldManager::*)()  const>(&G4FieldManager::DoesFieldChangeEnergy));

    DEBUG_MSG("Adding wrapper for void G4FieldManager::SetFieldChangesEnergy(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4FieldManager::SetFieldChangesEnergy(G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:167:21
    t.method("SetFieldChangesEnergy", static_cast<void (G4FieldManager::*)(G4bool) >(&G4FieldManager::SetFieldChangesEnergy));

    DEBUG_MSG("Adding wrapper for G4FieldManager * G4FieldManager::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4FieldManager * G4FieldManager::Clone()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:171:29
    t.method("Clone", static_cast<G4FieldManager * (G4FieldManager::*)()  const>(&G4FieldManager::Clone));

    DEBUG_MSG("Adding wrapper for G4double G4FieldManager::GetMaxAcceptedEpsilon() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4FieldManager::GetMaxAcceptedEpsilon()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:175:21
    t.method("G4FieldManager!GetMaxAcceptedEpsilon", static_cast<G4double (*)() >(&G4FieldManager::GetMaxAcceptedEpsilon));

    DEBUG_MSG("Adding wrapper for G4bool G4FieldManager::SetMaxAcceptedEpsilon(G4double, G4bool) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4FieldManager::SetMaxAcceptedEpsilon(G4double, G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4FieldManager.hh:176:21
    t.method("G4FieldManager!SetMaxAcceptedEpsilon", static_cast<G4bool (*)(G4double, G4bool) >(&G4FieldManager::SetMaxAcceptedEpsilon));
    t.method("G4FieldManager!SetMaxAcceptedEpsilon", [](G4double arg0)->G4bool { return G4FieldManager::SetMaxAcceptedEpsilon(arg0); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4FieldManager>> type_;
};
std::shared_ptr<Wrapper> newJlG4FieldManager(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4FieldManager(module));
}