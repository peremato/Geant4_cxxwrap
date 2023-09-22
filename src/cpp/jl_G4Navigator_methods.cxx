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
void add_methods_for_G4Navigator(jlcxx::Module& types, jlcxx::TypeWrapper<G4Navigator>& t140) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4Navigator
   */

  DEBUG_MSG("Adding wrapper for G4double G4Navigator::ComputeStep(const G4ThreeVector &, const G4ThreeVector &, const G4double, G4double &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Navigator::ComputeStep(const G4ThreeVector &, const G4ThreeVector &, const G4double, G4double &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:87:20
  t140.method("ComputeStep", static_cast<G4double (G4Navigator::*)(const G4ThreeVector &, const G4ThreeVector &, const G4double, G4double &) >(&G4Navigator::ComputeStep));

  DEBUG_MSG("Adding wrapper for G4double G4Navigator::CheckNextStep(const G4ThreeVector &, const G4ThreeVector &, const G4double, G4double &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Navigator::CheckNextStep(const G4ThreeVector &, const G4ThreeVector &, const G4double, G4double &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:104:12
  t140.method("CheckNextStep", static_cast<G4double (G4Navigator::*)(const G4ThreeVector &, const G4ThreeVector &, const G4double, G4double &) >(&G4Navigator::CheckNextStep));

  DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4Navigator::ResetHierarchyAndLocate(const G4ThreeVector &, const G4ThreeVector &, const G4TouchableHistory &) (" __HERE__ ")");
  // signature to use in the veto list: G4VPhysicalVolume * G4Navigator::ResetHierarchyAndLocate(const G4ThreeVector &, const G4ThreeVector &, const G4TouchableHistory &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:111:22
  t140.method("ResetHierarchyAndLocate", static_cast<G4VPhysicalVolume * (G4Navigator::*)(const G4ThreeVector &, const G4ThreeVector &, const G4TouchableHistory &) >(&G4Navigator::ResetHierarchyAndLocate));

  DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4Navigator::LocateGlobalPointAndSetup(const G4ThreeVector &, const G4ThreeVector *, const G4bool, const G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4VPhysicalVolume * G4Navigator::LocateGlobalPointAndSetup(const G4ThreeVector &, const G4ThreeVector *, const G4bool, const G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:125:22
  t140.method("LocateGlobalPointAndSetup", static_cast<G4VPhysicalVolume * (G4Navigator::*)(const G4ThreeVector &, const G4ThreeVector *, const G4bool, const G4bool) >(&G4Navigator::LocateGlobalPointAndSetup));
  t140.method("LocateGlobalPointAndSetup", [](G4Navigator& a, const G4ThreeVector & arg0)->G4VPhysicalVolume *{ return a.LocateGlobalPointAndSetup(arg0); });
  t140.method("LocateGlobalPointAndSetup", [](G4Navigator& a, const G4ThreeVector & arg0, const G4ThreeVector * arg1)->G4VPhysicalVolume *{ return a.LocateGlobalPointAndSetup(arg0, arg1); });
  t140.method("LocateGlobalPointAndSetup", [](G4Navigator& a, const G4ThreeVector & arg0, const G4ThreeVector * arg1, const G4bool arg2)->G4VPhysicalVolume *{ return a.LocateGlobalPointAndSetup(arg0, arg1, arg2); });
  t140.method("LocateGlobalPointAndSetup", [](G4Navigator* a, const G4ThreeVector & arg0)->G4VPhysicalVolume *{ return a->LocateGlobalPointAndSetup(arg0); });
  t140.method("LocateGlobalPointAndSetup", [](G4Navigator* a, const G4ThreeVector & arg0, const G4ThreeVector * arg1)->G4VPhysicalVolume *{ return a->LocateGlobalPointAndSetup(arg0, arg1); });
  t140.method("LocateGlobalPointAndSetup", [](G4Navigator* a, const G4ThreeVector & arg0, const G4ThreeVector * arg1, const G4bool arg2)->G4VPhysicalVolume *{ return a->LocateGlobalPointAndSetup(arg0, arg1, arg2); });

  DEBUG_MSG("Adding wrapper for void G4Navigator::LocateGlobalPointWithinVolume(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::LocateGlobalPointWithinVolume(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:146:8
  t140.method("LocateGlobalPointWithinVolume", static_cast<void (G4Navigator::*)(const G4ThreeVector &) >(&G4Navigator::LocateGlobalPointWithinVolume));

  DEBUG_MSG("Adding wrapper for void G4Navigator::LocateGlobalPointAndUpdateTouchable(const G4ThreeVector &, const G4ThreeVector &, G4VTouchable *, const G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::LocateGlobalPointAndUpdateTouchable(const G4ThreeVector &, const G4ThreeVector &, G4VTouchable *, const G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:164:15
  t140.method("LocateGlobalPointAndUpdateTouchable", static_cast<void (G4Navigator::*)(const G4ThreeVector &, const G4ThreeVector &, G4VTouchable *, const G4bool) >(&G4Navigator::LocateGlobalPointAndUpdateTouchable));
  t140.method("LocateGlobalPointAndUpdateTouchable", [](G4Navigator& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, G4VTouchable * arg2)->void{ a.LocateGlobalPointAndUpdateTouchable(arg0, arg1, arg2); });
  t140.method("LocateGlobalPointAndUpdateTouchable", [](G4Navigator* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, G4VTouchable * arg2)->void{ a->LocateGlobalPointAndUpdateTouchable(arg0, arg1, arg2); });

  DEBUG_MSG("Adding wrapper for void G4Navigator::LocateGlobalPointAndUpdateTouchable(const G4ThreeVector &, G4VTouchable *, const G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::LocateGlobalPointAndUpdateTouchable(const G4ThreeVector &, G4VTouchable *, const G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:173:15
  t140.method("LocateGlobalPointAndUpdateTouchable", static_cast<void (G4Navigator::*)(const G4ThreeVector &, G4VTouchable *, const G4bool) >(&G4Navigator::LocateGlobalPointAndUpdateTouchable));
  t140.method("LocateGlobalPointAndUpdateTouchable", [](G4Navigator& a, const G4ThreeVector & arg0, G4VTouchable * arg1)->void{ a.LocateGlobalPointAndUpdateTouchable(arg0, arg1); });
  t140.method("LocateGlobalPointAndUpdateTouchable", [](G4Navigator* a, const G4ThreeVector & arg0, G4VTouchable * arg1)->void{ a->LocateGlobalPointAndUpdateTouchable(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for void G4Navigator::SetGeometricallyLimitedStep() (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::SetGeometricallyLimitedStep()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:179:15
  t140.method("SetGeometricallyLimitedStep", static_cast<void (G4Navigator::*)() >(&G4Navigator::SetGeometricallyLimitedStep));

  DEBUG_MSG("Adding wrapper for G4double G4Navigator::ComputeSafety(const G4ThreeVector &, const G4double, const G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Navigator::ComputeSafety(const G4ThreeVector &, const G4double, const G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:183:20
  t140.method("ComputeSafety", static_cast<G4double (G4Navigator::*)(const G4ThreeVector &, const G4double, const G4bool) >(&G4Navigator::ComputeSafety));
  t140.method("ComputeSafety", [](G4Navigator& a, const G4ThreeVector & arg0)->G4double{ return a.ComputeSafety(arg0); });
  t140.method("ComputeSafety", [](G4Navigator& a, const G4ThreeVector & arg0, const G4double arg1)->G4double{ return a.ComputeSafety(arg0, arg1); });
  t140.method("ComputeSafety", [](G4Navigator* a, const G4ThreeVector & arg0)->G4double{ return a->ComputeSafety(arg0); });
  t140.method("ComputeSafety", [](G4Navigator* a, const G4ThreeVector & arg0, const G4double arg1)->G4double{ return a->ComputeSafety(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4Navigator::GetWorldVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4VPhysicalVolume * G4Navigator::GetWorldVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:195:29
  t140.method("GetWorldVolume", static_cast<G4VPhysicalVolume * (G4Navigator::*)()  const>(&G4Navigator::GetWorldVolume));

  DEBUG_MSG("Adding wrapper for void G4Navigator::SetWorldVolume(G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::SetWorldVolume(G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:198:15
  t140.method("SetWorldVolume", static_cast<void (G4Navigator::*)(G4VPhysicalVolume *) >(&G4Navigator::SetWorldVolume));

  DEBUG_MSG("Adding wrapper for G4GRSVolume * G4Navigator::CreateGRSVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4GRSVolume * G4Navigator::CreateGRSVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:202:23
  t140.method("CreateGRSVolume", static_cast<G4GRSVolume * (G4Navigator::*)()  const>(&G4Navigator::CreateGRSVolume));

  DEBUG_MSG("Adding wrapper for G4GRSSolid * G4Navigator::CreateGRSSolid() (" __HERE__ ")");
  // signature to use in the veto list: G4GRSSolid * G4Navigator::CreateGRSSolid()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:203:22
  t140.method("CreateGRSSolid", static_cast<G4GRSSolid * (G4Navigator::*)()  const>(&G4Navigator::CreateGRSSolid));

  DEBUG_MSG("Adding wrapper for G4TouchableHistory * G4Navigator::CreateTouchableHistory() (" __HERE__ ")");
  // signature to use in the veto list: G4TouchableHistory * G4Navigator::CreateTouchableHistory()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:204:30
  t140.method("CreateTouchableHistory", static_cast<G4TouchableHistory * (G4Navigator::*)()  const>(&G4Navigator::CreateTouchableHistory));

  DEBUG_MSG("Adding wrapper for G4TouchableHistory * G4Navigator::CreateTouchableHistory(const G4NavigationHistory *) (" __HERE__ ")");
  // signature to use in the veto list: G4TouchableHistory * G4Navigator::CreateTouchableHistory(const G4NavigationHistory *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:205:30
  t140.method("CreateTouchableHistory", static_cast<G4TouchableHistory * (G4Navigator::*)(const G4NavigationHistory *)  const>(&G4Navigator::CreateTouchableHistory));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Navigator::GetLocalExitNormal(G4bool *) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Navigator::GetLocalExitNormal(G4bool *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:211:25
  t140.method("GetLocalExitNormal", static_cast<G4ThreeVector (G4Navigator::*)(G4bool *) >(&G4Navigator::GetLocalExitNormal));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Navigator::GetLocalExitNormalAndCheck(const G4ThreeVector &, G4bool *) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Navigator::GetLocalExitNormalAndCheck(const G4ThreeVector &, G4bool *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:212:25
  t140.method("GetLocalExitNormalAndCheck", static_cast<G4ThreeVector (G4Navigator::*)(const G4ThreeVector &, G4bool *) >(&G4Navigator::GetLocalExitNormalAndCheck));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Navigator::GetGlobalExitNormal(const G4ThreeVector &, G4bool *) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Navigator::GetGlobalExitNormal(const G4ThreeVector &, G4bool *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:214:25
  t140.method("GetGlobalExitNormal", static_cast<G4ThreeVector (G4Navigator::*)(const G4ThreeVector &, G4bool *) >(&G4Navigator::GetGlobalExitNormal));

  DEBUG_MSG("Adding wrapper for G4int G4Navigator::GetVerboseLevel() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4Navigator::GetVerboseLevel()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:228:16
  t140.method("GetVerboseLevel", static_cast<G4int (G4Navigator::*)()  const>(&G4Navigator::GetVerboseLevel));

  DEBUG_MSG("Adding wrapper for void G4Navigator::SetVerboseLevel(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::SetVerboseLevel(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:229:16
  t140.method("SetVerboseLevel", static_cast<void (G4Navigator::*)(G4int) >(&G4Navigator::SetVerboseLevel));

  DEBUG_MSG("Adding wrapper for G4bool G4Navigator::IsActive() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4Navigator::IsActive()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:233:17
  t140.method("IsActive", static_cast<G4bool (G4Navigator::*)()  const>(&G4Navigator::IsActive));

  DEBUG_MSG("Adding wrapper for void G4Navigator::Activate(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::Activate(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:235:16
  t140.method("Activate", static_cast<void (G4Navigator::*)(G4bool) >(&G4Navigator::Activate));

  DEBUG_MSG("Adding wrapper for G4bool G4Navigator::EnteredDaughterVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4Navigator::EnteredDaughterVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:238:17
  t140.method("EnteredDaughterVolume", static_cast<G4bool (G4Navigator::*)()  const>(&G4Navigator::EnteredDaughterVolume));

  DEBUG_MSG("Adding wrapper for G4bool G4Navigator::ExitedMotherVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4Navigator::ExitedMotherVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:247:17
  t140.method("ExitedMotherVolume", static_cast<G4bool (G4Navigator::*)()  const>(&G4Navigator::ExitedMotherVolume));

  DEBUG_MSG("Adding wrapper for void G4Navigator::CheckMode(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::CheckMode(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:250:15
  t140.method("CheckMode", static_cast<void (G4Navigator::*)(G4bool) >(&G4Navigator::CheckMode));

  DEBUG_MSG("Adding wrapper for G4bool G4Navigator::IsCheckModeActive() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4Navigator::IsCheckModeActive()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:254:17
  t140.method("IsCheckModeActive", static_cast<G4bool (G4Navigator::*)()  const>(&G4Navigator::IsCheckModeActive));

  DEBUG_MSG("Adding wrapper for void G4Navigator::SetPushVerbosity(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::SetPushVerbosity(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:255:15
  t140.method("SetPushVerbosity", static_cast<void (G4Navigator::*)(G4bool) >(&G4Navigator::SetPushVerbosity));

  DEBUG_MSG("Adding wrapper for void G4Navigator::PrintState() (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::PrintState()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:258:8
  t140.method("PrintState", static_cast<void (G4Navigator::*)()  const>(&G4Navigator::PrintState));

  DEBUG_MSG("Adding wrapper for const G4AffineTransform & G4Navigator::GetGlobalToLocalTransform() (" __HERE__ ")");
  // signature to use in the veto list: const G4AffineTransform & G4Navigator::GetGlobalToLocalTransform()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:262:35
  t140.method("GetGlobalToLocalTransform", static_cast<const G4AffineTransform & (G4Navigator::*)()  const>(&G4Navigator::GetGlobalToLocalTransform));

  DEBUG_MSG("Adding wrapper for const G4AffineTransform G4Navigator::GetLocalToGlobalTransform() (" __HERE__ ")");
  // signature to use in the veto list: const G4AffineTransform G4Navigator::GetLocalToGlobalTransform()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:263:35
  t140.method("GetLocalToGlobalTransform", static_cast<const G4AffineTransform (G4Navigator::*)()  const>(&G4Navigator::GetLocalToGlobalTransform));

  DEBUG_MSG("Adding wrapper for G4AffineTransform G4Navigator::GetMotherToDaughterTransform(G4VPhysicalVolume *, G4int, EVolume) (" __HERE__ ")");
  // signature to use in the veto list: G4AffineTransform G4Navigator::GetMotherToDaughterTransform(G4VPhysicalVolume *, G4int, EVolume)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:267:21
  t140.method("GetMotherToDaughterTransform", static_cast<G4AffineTransform (G4Navigator::*)(G4VPhysicalVolume *, G4int, EVolume) >(&G4Navigator::GetMotherToDaughterTransform));

  DEBUG_MSG("Adding wrapper for void G4Navigator::ResetStackAndState() (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::ResetStackAndState()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:272:15
  t140.method("ResetStackAndState", static_cast<void (G4Navigator::*)() >(&G4Navigator::ResetStackAndState));

  DEBUG_MSG("Adding wrapper for G4int G4Navigator::SeverityOfZeroStepping(G4int *) (" __HERE__ ")");
  // signature to use in the veto list: G4int G4Navigator::SeverityOfZeroStepping(G4int *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:277:16
  t140.method("SeverityOfZeroStepping", static_cast<G4int (G4Navigator::*)(G4int *)  const>(&G4Navigator::SeverityOfZeroStepping));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Navigator::GetCurrentLocalCoordinate() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Navigator::GetCurrentLocalCoordinate()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:283:24
  t140.method("GetCurrentLocalCoordinate", static_cast<G4ThreeVector (G4Navigator::*)()  const>(&G4Navigator::GetCurrentLocalCoordinate));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Navigator::NetTranslation() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Navigator::NetTranslation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:288:24
  t140.method("NetTranslation", static_cast<G4ThreeVector (G4Navigator::*)()  const>(&G4Navigator::NetTranslation));

  DEBUG_MSG("Adding wrapper for G4RotationMatrix G4Navigator::NetRotation() (" __HERE__ ")");
  // signature to use in the veto list: G4RotationMatrix G4Navigator::NetRotation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:289:27
  t140.method("NetRotation", static_cast<G4RotationMatrix (G4Navigator::*)()  const>(&G4Navigator::NetRotation));

  DEBUG_MSG("Adding wrapper for void G4Navigator::EnableBestSafety(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::EnableBestSafety(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:292:15
  t140.method("EnableBestSafety", static_cast<void (G4Navigator::*)(G4bool) >(&G4Navigator::EnableBestSafety));
  t140.method("EnableBestSafety", [](G4Navigator& a)->void{ a.EnableBestSafety(); });
  t140.method("EnableBestSafety", [](G4Navigator* a)->void{ a->EnableBestSafety(); });

  DEBUG_MSG("Adding wrapper for G4VExternalNavigation * G4Navigator::GetExternalNavigation() (" __HERE__ ")");
  // signature to use in the veto list: G4VExternalNavigation * G4Navigator::GetExternalNavigation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:295:33
  t140.method("GetExternalNavigation", static_cast<G4VExternalNavigation * (G4Navigator::*)()  const>(&G4Navigator::GetExternalNavigation));

  DEBUG_MSG("Adding wrapper for void G4Navigator::SetExternalNavigation(G4VExternalNavigation *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::SetExternalNavigation(G4VExternalNavigation *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:296:15
  t140.method("SetExternalNavigation", static_cast<void (G4Navigator::*)(G4VExternalNavigation *) >(&G4Navigator::SetExternalNavigation));

  DEBUG_MSG("Adding wrapper for G4Navigator * G4Navigator::Clone() (" __HERE__ ")");
  // signature to use in the veto list: G4Navigator * G4Navigator::Clone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:299:23
  t140.method("Clone", static_cast<G4Navigator * (G4Navigator::*)()  const>(&G4Navigator::Clone));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Navigator::GetLastStepEndPoint() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Navigator::GetLastStepEndPoint()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:304:24
  t140.method("GetLastStepEndPoint", static_cast<G4ThreeVector (G4Navigator::*)()  const>(&G4Navigator::GetLastStepEndPoint));

  DEBUG_MSG("Adding wrapper for void G4Navigator::InformLastStep(G4double, G4bool, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::InformLastStep(G4double, G4bool, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:307:10
  t140.method("InformLastStep", static_cast<void (G4Navigator::*)(G4double, G4bool, G4bool) >(&G4Navigator::InformLastStep));

  DEBUG_MSG("Adding wrapper for void G4Navigator::SetVoxelNavigation(G4VoxelNavigation *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Navigator::SetVoxelNavigation(G4VoxelNavigation *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Navigator.hh:362:10
  t140.method("SetVoxelNavigation", static_cast<void (G4Navigator::*)(G4VoxelNavigation *) >(&G4Navigator::SetVoxelNavigation));

  /* End of G4Navigator class method wrappers
   **********************************************************************/

}
