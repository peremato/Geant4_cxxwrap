// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4Navigator> : std::false_type { };
  template<> struct DefaultConstructible<G4Navigator> : std::false_type { };
}

// Class generating the wrapper for type G4Navigator
// signature to use in the veto file: G4Navigator
struct JlG4Navigator: public Wrapper {

  JlG4Navigator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4Navigator (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:66:7
    jlcxx::TypeWrapper<G4Navigator>  t = jlModule.add_type<G4Navigator>("G4Navigator");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4Navigator>>(new jlcxx::TypeWrapper<G4Navigator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for G4double G4Navigator::ComputeStep(const G4ThreeVector &, const G4ThreeVector &, const G4double, G4double &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Navigator::ComputeStep(const G4ThreeVector &, const G4ThreeVector &, const G4double, G4double &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:82:22
    t.method("ComputeStep", static_cast<G4double (G4Navigator::*)(const G4ThreeVector &, const G4ThreeVector &, const G4double, G4double &) >(&G4Navigator::ComputeStep));

    DEBUG_MSG("Adding wrapper for G4double G4Navigator::CheckNextStep(const G4ThreeVector &, const G4ThreeVector &, const G4double, G4double &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Navigator::CheckNextStep(const G4ThreeVector &, const G4ThreeVector &, const G4double, G4double &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:99:14
    t.method("CheckNextStep", static_cast<G4double (G4Navigator::*)(const G4ThreeVector &, const G4ThreeVector &, const G4double, G4double &) >(&G4Navigator::CheckNextStep));

    DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4Navigator::ResetHierarchyAndLocate(const G4ThreeVector &, const G4ThreeVector &, const G4TouchableHistory &) (" __HERE__ ")");
    // signature to use in the veto list: G4VPhysicalVolume * G4Navigator::ResetHierarchyAndLocate(const G4ThreeVector &, const G4ThreeVector &, const G4TouchableHistory &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:106:24
    t.method("ResetHierarchyAndLocate", static_cast<G4VPhysicalVolume * (G4Navigator::*)(const G4ThreeVector &, const G4ThreeVector &, const G4TouchableHistory &) >(&G4Navigator::ResetHierarchyAndLocate));

    DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4Navigator::LocateGlobalPointAndSetup(const G4ThreeVector &, const G4ThreeVector *, const G4bool, const G4bool) (" __HERE__ ")");
    // signature to use in the veto list: G4VPhysicalVolume * G4Navigator::LocateGlobalPointAndSetup(const G4ThreeVector &, const G4ThreeVector *, const G4bool, const G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:120:24
    t.method("LocateGlobalPointAndSetup", static_cast<G4VPhysicalVolume * (G4Navigator::*)(const G4ThreeVector &, const G4ThreeVector *, const G4bool, const G4bool) >(&G4Navigator::LocateGlobalPointAndSetup));
    t.method("LocateGlobalPointAndSetup", [](G4Navigator& a, const G4ThreeVector & arg0)->G4VPhysicalVolume * { return a.LocateGlobalPointAndSetup(arg0); });
    t.method("LocateGlobalPointAndSetup", [](G4Navigator& a, const G4ThreeVector & arg0, const G4ThreeVector * arg1)->G4VPhysicalVolume * { return a.LocateGlobalPointAndSetup(arg0, arg1); });
    t.method("LocateGlobalPointAndSetup", [](G4Navigator& a, const G4ThreeVector & arg0, const G4ThreeVector * arg1, const G4bool arg2)->G4VPhysicalVolume * { return a.LocateGlobalPointAndSetup(arg0, arg1, arg2); });
    t.method("LocateGlobalPointAndSetup", [](G4Navigator* a, const G4ThreeVector & arg0)->G4VPhysicalVolume * { return a->LocateGlobalPointAndSetup(arg0); });
    t.method("LocateGlobalPointAndSetup", [](G4Navigator* a, const G4ThreeVector & arg0, const G4ThreeVector * arg1)->G4VPhysicalVolume * { return a->LocateGlobalPointAndSetup(arg0, arg1); });
    t.method("LocateGlobalPointAndSetup", [](G4Navigator* a, const G4ThreeVector & arg0, const G4ThreeVector * arg1, const G4bool arg2)->G4VPhysicalVolume * { return a->LocateGlobalPointAndSetup(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void G4Navigator::LocateGlobalPointWithinVolume(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::LocateGlobalPointWithinVolume(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:142:10
    t.method("LocateGlobalPointWithinVolume", static_cast<void (G4Navigator::*)(const G4ThreeVector &) >(&G4Navigator::LocateGlobalPointWithinVolume));

    DEBUG_MSG("Adding wrapper for void G4Navigator::LocateGlobalPointAndUpdateTouchableHandle(const G4ThreeVector &, const G4ThreeVector &, G4TouchableHandle &, const G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::LocateGlobalPointAndUpdateTouchableHandle(const G4ThreeVector &, const G4ThreeVector &, G4TouchableHandle &, const G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:151:17
    t.method("LocateGlobalPointAndUpdateTouchableHandle", static_cast<void (G4Navigator::*)(const G4ThreeVector &, const G4ThreeVector &, G4TouchableHandle &, const G4bool) >(&G4Navigator::LocateGlobalPointAndUpdateTouchableHandle));
    t.method("LocateGlobalPointAndUpdateTouchableHandle", [](G4Navigator& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, G4TouchableHandle & arg2)->void { a.LocateGlobalPointAndUpdateTouchableHandle(arg0, arg1, arg2); });
    t.method("LocateGlobalPointAndUpdateTouchableHandle", [](G4Navigator* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, G4TouchableHandle & arg2)->void { a->LocateGlobalPointAndUpdateTouchableHandle(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void G4Navigator::LocateGlobalPointAndUpdateTouchable(const G4ThreeVector &, const G4ThreeVector &, G4VTouchable *, const G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::LocateGlobalPointAndUpdateTouchable(const G4ThreeVector &, const G4ThreeVector &, G4VTouchable *, const G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:160:17
    t.method("LocateGlobalPointAndUpdateTouchable", static_cast<void (G4Navigator::*)(const G4ThreeVector &, const G4ThreeVector &, G4VTouchable *, const G4bool) >(&G4Navigator::LocateGlobalPointAndUpdateTouchable));
    t.method("LocateGlobalPointAndUpdateTouchable", [](G4Navigator& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, G4VTouchable * arg2)->void { a.LocateGlobalPointAndUpdateTouchable(arg0, arg1, arg2); });
    t.method("LocateGlobalPointAndUpdateTouchable", [](G4Navigator* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, G4VTouchable * arg2)->void { a->LocateGlobalPointAndUpdateTouchable(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void G4Navigator::LocateGlobalPointAndUpdateTouchable(const G4ThreeVector &, G4VTouchable *, const G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::LocateGlobalPointAndUpdateTouchable(const G4ThreeVector &, G4VTouchable *, const G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:169:17
    t.method("LocateGlobalPointAndUpdateTouchable", static_cast<void (G4Navigator::*)(const G4ThreeVector &, G4VTouchable *, const G4bool) >(&G4Navigator::LocateGlobalPointAndUpdateTouchable));
    t.method("LocateGlobalPointAndUpdateTouchable", [](G4Navigator& a, const G4ThreeVector & arg0, G4VTouchable * arg1)->void { a.LocateGlobalPointAndUpdateTouchable(arg0, arg1); });
    t.method("LocateGlobalPointAndUpdateTouchable", [](G4Navigator* a, const G4ThreeVector & arg0, G4VTouchable * arg1)->void { a->LocateGlobalPointAndUpdateTouchable(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void G4Navigator::SetGeometricallyLimitedStep() (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::SetGeometricallyLimitedStep()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:175:17
    t.method("SetGeometricallyLimitedStep", static_cast<void (G4Navigator::*)() >(&G4Navigator::SetGeometricallyLimitedStep));

    DEBUG_MSG("Adding wrapper for G4double G4Navigator::ComputeSafety(const G4ThreeVector &, const G4double, const G4bool) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Navigator::ComputeSafety(const G4ThreeVector &, const G4double, const G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:179:22
    t.method("ComputeSafety", static_cast<G4double (G4Navigator::*)(const G4ThreeVector &, const G4double, const G4bool) >(&G4Navigator::ComputeSafety));
    t.method("ComputeSafety", [](G4Navigator& a, const G4ThreeVector & arg0)->G4double { return a.ComputeSafety(arg0); });
    t.method("ComputeSafety", [](G4Navigator& a, const G4ThreeVector & arg0, const G4double arg1)->G4double { return a.ComputeSafety(arg0, arg1); });
    t.method("ComputeSafety", [](G4Navigator* a, const G4ThreeVector & arg0)->G4double { return a->ComputeSafety(arg0); });
    t.method("ComputeSafety", [](G4Navigator* a, const G4ThreeVector & arg0, const G4double arg1)->G4double { return a->ComputeSafety(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4Navigator::GetWorldVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4VPhysicalVolume * G4Navigator::GetWorldVolume()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:190:31
    t.method("GetWorldVolume", static_cast<G4VPhysicalVolume * (G4Navigator::*)()  const>(&G4Navigator::GetWorldVolume));

    DEBUG_MSG("Adding wrapper for void G4Navigator::SetWorldVolume(G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::SetWorldVolume(G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:193:17
    t.method("SetWorldVolume", static_cast<void (G4Navigator::*)(G4VPhysicalVolume *) >(&G4Navigator::SetWorldVolume));

    DEBUG_MSG("Adding wrapper for G4TouchableHistory * G4Navigator::CreateTouchableHistory() (" __HERE__ ")");
    // signature to use in the veto list: G4TouchableHistory * G4Navigator::CreateTouchableHistory()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:197:32
    t.method("CreateTouchableHistory", static_cast<G4TouchableHistory * (G4Navigator::*)()  const>(&G4Navigator::CreateTouchableHistory));

    DEBUG_MSG("Adding wrapper for G4TouchableHistory * G4Navigator::CreateTouchableHistory(const G4NavigationHistory *) (" __HERE__ ")");
    // signature to use in the veto list: G4TouchableHistory * G4Navigator::CreateTouchableHistory(const G4NavigationHistory *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:198:32
    t.method("CreateTouchableHistory", static_cast<G4TouchableHistory * (G4Navigator::*)(const G4NavigationHistory *)  const>(&G4Navigator::CreateTouchableHistory));

    DEBUG_MSG("Adding wrapper for G4TouchableHandle G4Navigator::CreateTouchableHistoryHandle() (" __HERE__ ")");
    // signature to use in the veto list: G4TouchableHandle G4Navigator::CreateTouchableHistoryHandle()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:201:31
    t.method("CreateTouchableHistoryHandle", static_cast<G4TouchableHandle (G4Navigator::*)()  const>(&G4Navigator::CreateTouchableHistoryHandle));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Navigator::GetLocalExitNormal(G4bool *) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Navigator::GetLocalExitNormal(G4bool *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:204:27
    t.method("GetLocalExitNormal", static_cast<G4ThreeVector (G4Navigator::*)(G4bool *) >(&G4Navigator::GetLocalExitNormal));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Navigator::GetLocalExitNormalAndCheck(const G4ThreeVector &, G4bool *) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Navigator::GetLocalExitNormalAndCheck(const G4ThreeVector &, G4bool *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:205:27
    t.method("GetLocalExitNormalAndCheck", static_cast<G4ThreeVector (G4Navigator::*)(const G4ThreeVector &, G4bool *) >(&G4Navigator::GetLocalExitNormalAndCheck));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Navigator::GetGlobalExitNormal(const G4ThreeVector &, G4bool *) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Navigator::GetGlobalExitNormal(const G4ThreeVector &, G4bool *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:207:27
    t.method("GetGlobalExitNormal", static_cast<G4ThreeVector (G4Navigator::*)(const G4ThreeVector &, G4bool *) >(&G4Navigator::GetGlobalExitNormal));

    DEBUG_MSG("Adding wrapper for G4int G4Navigator::GetVerboseLevel() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Navigator::GetVerboseLevel()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:221:18
    t.method("GetVerboseLevel", static_cast<G4int (G4Navigator::*)()  const>(&G4Navigator::GetVerboseLevel));

    DEBUG_MSG("Adding wrapper for void G4Navigator::SetVerboseLevel(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::SetVerboseLevel(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:222:18
    t.method("SetVerboseLevel", static_cast<void (G4Navigator::*)(G4int) >(&G4Navigator::SetVerboseLevel));

    DEBUG_MSG("Adding wrapper for G4bool G4Navigator::IsActive() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Navigator::IsActive()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:226:19
    t.method("IsActive", static_cast<G4bool (G4Navigator::*)()  const>(&G4Navigator::IsActive));

    DEBUG_MSG("Adding wrapper for void G4Navigator::Activate(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::Activate(G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:228:18
    t.method("Activate", static_cast<void (G4Navigator::*)(G4bool) >(&G4Navigator::Activate));

    DEBUG_MSG("Adding wrapper for G4bool G4Navigator::EnteredDaughterVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Navigator::EnteredDaughterVolume()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:231:19
    t.method("EnteredDaughterVolume", static_cast<G4bool (G4Navigator::*)()  const>(&G4Navigator::EnteredDaughterVolume));

    DEBUG_MSG("Adding wrapper for G4bool G4Navigator::ExitedMotherVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Navigator::ExitedMotherVolume()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:240:19
    t.method("ExitedMotherVolume", static_cast<G4bool (G4Navigator::*)()  const>(&G4Navigator::ExitedMotherVolume));

    DEBUG_MSG("Adding wrapper for void G4Navigator::CheckMode(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::CheckMode(G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:243:17
    t.method("CheckMode", static_cast<void (G4Navigator::*)(G4bool) >(&G4Navigator::CheckMode));

    DEBUG_MSG("Adding wrapper for G4bool G4Navigator::IsCheckModeActive() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Navigator::IsCheckModeActive()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:247:19
    t.method("IsCheckModeActive", static_cast<G4bool (G4Navigator::*)()  const>(&G4Navigator::IsCheckModeActive));

    DEBUG_MSG("Adding wrapper for void G4Navigator::SetPushVerbosity(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::SetPushVerbosity(G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:248:17
    t.method("SetPushVerbosity", static_cast<void (G4Navigator::*)(G4bool) >(&G4Navigator::SetPushVerbosity));

    DEBUG_MSG("Adding wrapper for void G4Navigator::PrintState() (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::PrintState()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:251:10
    t.method("PrintState", static_cast<void (G4Navigator::*)()  const>(&G4Navigator::PrintState));

    DEBUG_MSG("Adding wrapper for const G4AffineTransform & G4Navigator::GetGlobalToLocalTransform() (" __HERE__ ")");
    // signature to use in the veto list: const G4AffineTransform & G4Navigator::GetGlobalToLocalTransform()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:255:37
    t.method("GetGlobalToLocalTransform", static_cast<const G4AffineTransform & (G4Navigator::*)()  const>(&G4Navigator::GetGlobalToLocalTransform));

    DEBUG_MSG("Adding wrapper for const G4AffineTransform G4Navigator::GetLocalToGlobalTransform() (" __HERE__ ")");
    // signature to use in the veto list: const G4AffineTransform G4Navigator::GetLocalToGlobalTransform()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:256:37
    t.method("GetLocalToGlobalTransform", static_cast<const G4AffineTransform (G4Navigator::*)()  const>(&G4Navigator::GetLocalToGlobalTransform));

    DEBUG_MSG("Adding wrapper for G4AffineTransform G4Navigator::GetMotherToDaughterTransform(G4VPhysicalVolume *, G4int, EVolume) (" __HERE__ ")");
    // signature to use in the veto list: G4AffineTransform G4Navigator::GetMotherToDaughterTransform(G4VPhysicalVolume *, G4int, EVolume)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:260:23
    t.method("GetMotherToDaughterTransform", static_cast<G4AffineTransform (G4Navigator::*)(G4VPhysicalVolume *, G4int, EVolume) >(&G4Navigator::GetMotherToDaughterTransform));

    DEBUG_MSG("Adding wrapper for void G4Navigator::ResetStackAndState() (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::ResetStackAndState()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:265:17
    t.method("ResetStackAndState", static_cast<void (G4Navigator::*)() >(&G4Navigator::ResetStackAndState));

    DEBUG_MSG("Adding wrapper for G4int G4Navigator::SeverityOfZeroStepping(G4int *) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Navigator::SeverityOfZeroStepping(G4int *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:270:18
    t.method("SeverityOfZeroStepping", static_cast<G4int (G4Navigator::*)(G4int *)  const>(&G4Navigator::SeverityOfZeroStepping));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Navigator::GetCurrentLocalCoordinate() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Navigator::GetCurrentLocalCoordinate()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:276:26
    t.method("GetCurrentLocalCoordinate", static_cast<G4ThreeVector (G4Navigator::*)()  const>(&G4Navigator::GetCurrentLocalCoordinate));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Navigator::NetTranslation() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Navigator::NetTranslation()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:281:26
    t.method("NetTranslation", static_cast<G4ThreeVector (G4Navigator::*)()  const>(&G4Navigator::NetTranslation));

    DEBUG_MSG("Adding wrapper for G4RotationMatrix G4Navigator::NetRotation() (" __HERE__ ")");
    // signature to use in the veto list: G4RotationMatrix G4Navigator::NetRotation()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:282:29
    t.method("NetRotation", static_cast<G4RotationMatrix (G4Navigator::*)()  const>(&G4Navigator::NetRotation));

    DEBUG_MSG("Adding wrapper for void G4Navigator::EnableBestSafety(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::EnableBestSafety(G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:285:17
    t.method("EnableBestSafety", static_cast<void (G4Navigator::*)(G4bool) >(&G4Navigator::EnableBestSafety));
    t.method("EnableBestSafety", [](G4Navigator& a)->void { a.EnableBestSafety(); });
    t.method("EnableBestSafety", [](G4Navigator* a)->void { a->EnableBestSafety(); });

    DEBUG_MSG("Adding wrapper for G4VExternalNavigation * G4Navigator::GetExternalNavigation() (" __HERE__ ")");
    // signature to use in the veto list: G4VExternalNavigation * G4Navigator::GetExternalNavigation()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:288:35
    t.method("GetExternalNavigation", static_cast<G4VExternalNavigation * (G4Navigator::*)()  const>(&G4Navigator::GetExternalNavigation));

    DEBUG_MSG("Adding wrapper for void G4Navigator::SetExternalNavigation(G4VExternalNavigation *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::SetExternalNavigation(G4VExternalNavigation *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:289:10
    t.method("SetExternalNavigation", static_cast<void (G4Navigator::*)(G4VExternalNavigation *) >(&G4Navigator::SetExternalNavigation));

    DEBUG_MSG("Adding wrapper for G4VoxelNavigation & G4Navigator::GetVoxelNavigator() (" __HERE__ ")");
    // signature to use in the veto list: G4VoxelNavigation & G4Navigator::GetVoxelNavigator()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:292:31
    t.method("GetVoxelNavigator", static_cast<G4VoxelNavigation & (G4Navigator::*)() >(&G4Navigator::GetVoxelNavigator));

    DEBUG_MSG("Adding wrapper for void G4Navigator::SetVoxelNavigation(G4VoxelNavigation *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::SetVoxelNavigation(G4VoxelNavigation *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:293:10
    t.method("SetVoxelNavigation", static_cast<void (G4Navigator::*)(G4VoxelNavigation *) >(&G4Navigator::SetVoxelNavigation));

    DEBUG_MSG("Adding wrapper for G4Navigator * G4Navigator::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4Navigator * G4Navigator::Clone()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:296:25
    t.method("Clone", static_cast<G4Navigator * (G4Navigator::*)()  const>(&G4Navigator::Clone));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Navigator::GetLastStepEndPoint() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Navigator::GetLastStepEndPoint()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:301:26
    t.method("GetLastStepEndPoint", static_cast<G4ThreeVector (G4Navigator::*)()  const>(&G4Navigator::GetLastStepEndPoint));

    DEBUG_MSG("Adding wrapper for void G4Navigator::InformLastStep(G4double, G4bool, G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4Navigator::InformLastStep(G4double, G4bool, G4bool)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Navigator.hh:304:10
    t.method("InformLastStep", static_cast<void (G4Navigator::*)(G4double, G4bool, G4bool) >(&G4Navigator::InformLastStep));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4Navigator>> type_;
};
std::shared_ptr<Wrapper> newJlG4Navigator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4Navigator(module));
}