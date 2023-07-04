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
void add_methods_for_G4LogicalVolume(jlcxx::Module& types, jlcxx::TypeWrapper<G4LogicalVolume>& t29) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4LogicalVolume
   */


  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::G4LogicalVolume(G4VSolid *, G4Material *, const G4String &, G4FieldManager *, G4VSensitiveDetector *, G4UserLimits *, G4bool) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:186:5
  t29.constructor<G4VSolid *, G4Material *, const G4String &>(/*finalize=*/false);
  t29.constructor<G4VSolid *, G4Material *, const G4String &, G4FieldManager *>(/*finalize=*/false);
  t29.constructor<G4VSolid *, G4Material *, const G4String &, G4FieldManager *, G4VSensitiveDetector *>(/*finalize=*/false);
  t29.constructor<G4VSolid *, G4Material *, const G4String &, G4FieldManager *, G4VSensitiveDetector *, G4UserLimits *>(/*finalize=*/false);
  t29.constructor<G4VSolid *, G4Material *, const G4String &, G4FieldManager *, G4VSensitiveDetector *, G4UserLimits *, G4bool>(/*finalize=*/false);

  DEBUG_MSG("Adding wrapper for const G4String & G4LogicalVolume::GetName() (" __HERE__ ")");
  // signature to use in the veto list: const G4String & G4LogicalVolume::GetName()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:209:28
  t29.method("GetName", static_cast<const G4String & (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetName));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetName(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetName(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:210:10
  t29.method("SetName", static_cast<void (G4LogicalVolume::*)(const G4String &) >(&G4LogicalVolume::SetName));

  DEBUG_MSG("Adding wrapper for size_t G4LogicalVolume::GetNoDaughters() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4LogicalVolume::GetNoDaughters()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:213:24
  t29.method("GetNoDaughters", static_cast<size_t (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetNoDaughters));

  DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4LogicalVolume::GetDaughter(const size_t) (" __HERE__ ")");
  // signature to use in the veto list: G4VPhysicalVolume * G4LogicalVolume::GetDaughter(const size_t)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:215:31
  t29.method("GetDaughter", static_cast<G4VPhysicalVolume * (G4LogicalVolume::*)(const size_t)  const>(&G4LogicalVolume::GetDaughter));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::AddDaughter(G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::AddDaughter(G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:218:10
  t29.method("AddDaughter", static_cast<void (G4LogicalVolume::*)(G4VPhysicalVolume *) >(&G4LogicalVolume::AddDaughter));

  DEBUG_MSG("Adding wrapper for G4bool G4LogicalVolume::IsDaughter(const G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4LogicalVolume::IsDaughter(const G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:220:19
  t29.method("IsDaughter", static_cast<G4bool (G4LogicalVolume::*)(const G4VPhysicalVolume *)  const>(&G4LogicalVolume::IsDaughter));

  DEBUG_MSG("Adding wrapper for G4bool G4LogicalVolume::IsAncestor(const G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4LogicalVolume::IsAncestor(const G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:223:12
  t29.method("IsAncestor", static_cast<G4bool (G4LogicalVolume::*)(const G4VPhysicalVolume *)  const>(&G4LogicalVolume::IsAncestor));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::RemoveDaughter(const G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::RemoveDaughter(const G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:227:10
  t29.method("RemoveDaughter", static_cast<void (G4LogicalVolume::*)(const G4VPhysicalVolume *) >(&G4LogicalVolume::RemoveDaughter));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::ClearDaughters() (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::ClearDaughters()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:230:10
  t29.method("ClearDaughters", static_cast<void (G4LogicalVolume::*)() >(&G4LogicalVolume::ClearDaughters));

  DEBUG_MSG("Adding wrapper for G4int G4LogicalVolume::TotalVolumeEntities() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4LogicalVolume::TotalVolumeEntities()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:233:11
  t29.method("TotalVolumeEntities", static_cast<G4int (G4LogicalVolume::*)()  const>(&G4LogicalVolume::TotalVolumeEntities));

  DEBUG_MSG("Adding wrapper for EVolume G4LogicalVolume::CharacteriseDaughters() (" __HERE__ ")");
  // signature to use in the veto list: EVolume G4LogicalVolume::CharacteriseDaughters()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:236:20
  t29.method("CharacteriseDaughters", static_cast<EVolume (G4LogicalVolume::*)()  const>(&G4LogicalVolume::CharacteriseDaughters));

  DEBUG_MSG("Adding wrapper for EVolume G4LogicalVolume::DeduceDaughtersType() (" __HERE__ ")");
  // signature to use in the veto list: EVolume G4LogicalVolume::DeduceDaughtersType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:238:20
  t29.method("DeduceDaughtersType", static_cast<EVolume (G4LogicalVolume::*)()  const>(&G4LogicalVolume::DeduceDaughtersType));

  DEBUG_MSG("Adding wrapper for G4VSolid * G4LogicalVolume::GetSolid() (" __HERE__ ")");
  // signature to use in the veto list: G4VSolid * G4LogicalVolume::GetSolid()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:241:15
  t29.method("GetSolid", static_cast<G4VSolid * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetSolid));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetSolid(G4VSolid *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetSolid(G4VSolid *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:242:10
  t29.method("SetSolid", static_cast<void (G4LogicalVolume::*)(G4VSolid *) >(&G4LogicalVolume::SetSolid));

  DEBUG_MSG("Adding wrapper for G4Material * G4LogicalVolume::GetMaterial() (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4LogicalVolume::GetMaterial()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:245:17
  t29.method("GetMaterial", static_cast<G4Material * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetMaterial));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetMaterial(G4Material *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetMaterial(G4Material *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:246:10
  t29.method("SetMaterial", static_cast<void (G4LogicalVolume::*)(G4Material *) >(&G4LogicalVolume::SetMaterial));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::UpdateMaterial(G4Material *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::UpdateMaterial(G4Material *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:248:10
  t29.method("UpdateMaterial", static_cast<void (G4LogicalVolume::*)(G4Material *) >(&G4LogicalVolume::UpdateMaterial));

  DEBUG_MSG("Adding wrapper for G4double G4LogicalVolume::GetMass(G4bool, G4bool, G4Material *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4LogicalVolume::GetMass(G4bool, G4bool, G4Material *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:252:14
  t29.method("GetMass", static_cast<G4double (G4LogicalVolume::*)(G4bool, G4bool, G4Material *) >(&G4LogicalVolume::GetMass));
  t29.method("GetMass", [](G4LogicalVolume& a)->G4double{ return a.GetMass(); });
  t29.method("GetMass", [](G4LogicalVolume& a, G4bool arg0)->G4double{ return a.GetMass(arg0); });
  t29.method("GetMass", [](G4LogicalVolume& a, G4bool arg0, G4bool arg1)->G4double{ return a.GetMass(arg0, arg1); });
  t29.method("GetMass", [](G4LogicalVolume* a)->G4double{ return a->GetMass(); });
  t29.method("GetMass", [](G4LogicalVolume* a, G4bool arg0)->G4double{ return a->GetMass(arg0); });
  t29.method("GetMass", [](G4LogicalVolume* a, G4bool arg0, G4bool arg1)->G4double{ return a->GetMass(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::ResetMass() (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::ResetMass()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:267:10
  t29.method("ResetMass", static_cast<void (G4LogicalVolume::*)() >(&G4LogicalVolume::ResetMass));

  DEBUG_MSG("Adding wrapper for G4FieldManager * G4LogicalVolume::GetFieldManager() (" __HERE__ ")");
  // signature to use in the veto list: G4FieldManager * G4LogicalVolume::GetFieldManager()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:272:21
  t29.method("GetFieldManager", static_cast<G4FieldManager * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetFieldManager));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetFieldManager(G4FieldManager *, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetFieldManager(G4FieldManager *, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:274:10
  t29.method("SetFieldManager", static_cast<void (G4LogicalVolume::*)(G4FieldManager *, G4bool) >(&G4LogicalVolume::SetFieldManager));

  DEBUG_MSG("Adding wrapper for G4VSensitiveDetector * G4LogicalVolume::GetSensitiveDetector() (" __HERE__ ")");
  // signature to use in the veto list: G4VSensitiveDetector * G4LogicalVolume::GetSensitiveDetector()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:281:27
  t29.method("GetSensitiveDetector", static_cast<G4VSensitiveDetector * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetSensitiveDetector));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetSensitiveDetector(G4VSensitiveDetector *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetSensitiveDetector(G4VSensitiveDetector *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:283:10
  t29.method("SetSensitiveDetector", static_cast<void (G4LogicalVolume::*)(G4VSensitiveDetector *) >(&G4LogicalVolume::SetSensitiveDetector));

  DEBUG_MSG("Adding wrapper for G4UserLimits * G4LogicalVolume::GetUserLimits() (" __HERE__ ")");
  // signature to use in the veto list: G4UserLimits * G4LogicalVolume::GetUserLimits()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:286:26
  t29.method("GetUserLimits", static_cast<G4UserLimits * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetUserLimits));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetUserLimits(G4UserLimits *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetUserLimits(G4UserLimits *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:287:17
  t29.method("SetUserLimits", static_cast<void (G4LogicalVolume::*)(G4UserLimits *) >(&G4LogicalVolume::SetUserLimits));

  DEBUG_MSG("Adding wrapper for G4SmartVoxelHeader * G4LogicalVolume::GetVoxelHeader() (" __HERE__ ")");
  // signature to use in the veto list: G4SmartVoxelHeader * G4LogicalVolume::GetVoxelHeader()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:290:32
  t29.method("GetVoxelHeader", static_cast<G4SmartVoxelHeader * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetVoxelHeader));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetVoxelHeader(G4SmartVoxelHeader *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetVoxelHeader(G4SmartVoxelHeader *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:291:17
  t29.method("SetVoxelHeader", static_cast<void (G4LogicalVolume::*)(G4SmartVoxelHeader *) >(&G4LogicalVolume::SetVoxelHeader));

  DEBUG_MSG("Adding wrapper for G4double G4LogicalVolume::GetSmartless() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4LogicalVolume::GetSmartless()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:294:21
  t29.method("GetSmartless", static_cast<G4double (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetSmartless));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetSmartless(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetSmartless(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:295:17
  t29.method("SetSmartless", static_cast<void (G4LogicalVolume::*)(G4double) >(&G4LogicalVolume::SetSmartless));

  DEBUG_MSG("Adding wrapper for G4bool G4LogicalVolume::IsToOptimise() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4LogicalVolume::IsToOptimise()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:298:19
  t29.method("IsToOptimise", static_cast<G4bool (G4LogicalVolume::*)()  const>(&G4LogicalVolume::IsToOptimise));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetOptimisation(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetOptimisation(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:301:17
  t29.method("SetOptimisation", static_cast<void (G4LogicalVolume::*)(G4bool) >(&G4LogicalVolume::SetOptimisation));

  DEBUG_MSG("Adding wrapper for G4bool G4LogicalVolume::IsRootRegion() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4LogicalVolume::IsRootRegion()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:306:19
  t29.method("IsRootRegion", static_cast<G4bool (G4LogicalVolume::*)()  const>(&G4LogicalVolume::IsRootRegion));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetRegionRootFlag(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetRegionRootFlag(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:308:17
  t29.method("SetRegionRootFlag", static_cast<void (G4LogicalVolume::*)(G4bool) >(&G4LogicalVolume::SetRegionRootFlag));

  DEBUG_MSG("Adding wrapper for G4bool G4LogicalVolume::IsRegion() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4LogicalVolume::IsRegion()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:310:19
  t29.method("IsRegion", static_cast<G4bool (G4LogicalVolume::*)()  const>(&G4LogicalVolume::IsRegion));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetRegion(G4Region *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetRegion(G4Region *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:312:17
  t29.method("SetRegion", static_cast<void (G4LogicalVolume::*)(G4Region *) >(&G4LogicalVolume::SetRegion));

  DEBUG_MSG("Adding wrapper for G4Region * G4LogicalVolume::GetRegion() (" __HERE__ ")");
  // signature to use in the veto list: G4Region * G4LogicalVolume::GetRegion()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:314:22
  t29.method("GetRegion", static_cast<G4Region * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetRegion));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::PropagateRegion() (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::PropagateRegion()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:316:17
  t29.method("PropagateRegion", static_cast<void (G4LogicalVolume::*)() >(&G4LogicalVolume::PropagateRegion));

  DEBUG_MSG("Adding wrapper for const G4MaterialCutsCouple * G4LogicalVolume::GetMaterialCutsCouple() (" __HERE__ ")");
  // signature to use in the veto list: const G4MaterialCutsCouple * G4LogicalVolume::GetMaterialCutsCouple()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:319:33
  t29.method("GetMaterialCutsCouple", static_cast<const G4MaterialCutsCouple * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetMaterialCutsCouple));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetMaterialCutsCouple(G4MaterialCutsCouple *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetMaterialCutsCouple(G4MaterialCutsCouple *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:320:10
  t29.method("SetMaterialCutsCouple", static_cast<void (G4LogicalVolume::*)(G4MaterialCutsCouple *) >(&G4LogicalVolume::SetMaterialCutsCouple));
  types.set_override_module(jl_base_module);

  DEBUG_MSG("Adding wrapper for G4bool G4LogicalVolume::operator==(const G4LogicalVolume &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4LogicalVolume::operator==(const G4LogicalVolume &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:323:12
  t29.method("==", static_cast<G4bool (G4LogicalVolume::*)(const G4LogicalVolume &)  const>(&G4LogicalVolume::operator==));

  types.unset_override_module();

  DEBUG_MSG("Adding wrapper for const G4VisAttributes * G4LogicalVolume::GetVisAttributes() (" __HERE__ ")");
  // signature to use in the veto list: const G4VisAttributes * G4LogicalVolume::GetVisAttributes()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:327:28
  t29.method("GetVisAttributes", static_cast<const G4VisAttributes * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetVisAttributes));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetVisAttributes(const G4VisAttributes *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetVisAttributes(const G4VisAttributes *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:328:10
  t29.method("SetVisAttributes", static_cast<void (G4LogicalVolume::*)(const G4VisAttributes *) >(&G4LogicalVolume::SetVisAttributes));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetVisAttributes(const G4VisAttributes &) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetVisAttributes(const G4VisAttributes &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:329:10
  t29.method("SetVisAttributes", static_cast<void (G4LogicalVolume::*)(const G4VisAttributes &) >(&G4LogicalVolume::SetVisAttributes));

  DEBUG_MSG("Adding wrapper for G4FastSimulationManager * G4LogicalVolume::GetFastSimulationManager() (" __HERE__ ")");
  // signature to use in the veto list: G4FastSimulationManager * G4LogicalVolume::GetFastSimulationManager()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:333:37
  t29.method("GetFastSimulationManager", static_cast<G4FastSimulationManager * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetFastSimulationManager));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetBiasWeight(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetBiasWeight(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:336:17
  t29.method("SetBiasWeight", static_cast<void (G4LogicalVolume::*)(G4double) >(&G4LogicalVolume::SetBiasWeight));

  DEBUG_MSG("Adding wrapper for G4double G4LogicalVolume::GetBiasWeight() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4LogicalVolume::GetBiasWeight()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:337:21
  t29.method("GetBiasWeight", static_cast<G4double (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetBiasWeight));

  DEBUG_MSG("Adding wrapper for G4bool G4LogicalVolume::IsExtended() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4LogicalVolume::IsExtended()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:347:20
  t29.method("IsExtended", static_cast<G4bool (G4LogicalVolume::*)()  const>(&G4LogicalVolume::IsExtended));

  DEBUG_MSG("Adding wrapper for G4FieldManager * G4LogicalVolume::GetMasterFieldManager() (" __HERE__ ")");
  // signature to use in the veto list: G4FieldManager * G4LogicalVolume::GetMasterFieldManager()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:350:28
  t29.method("GetMasterFieldManager", static_cast<G4FieldManager * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetMasterFieldManager));

  DEBUG_MSG("Adding wrapper for G4VSensitiveDetector * G4LogicalVolume::GetMasterSensitiveDetector() (" __HERE__ ")");
  // signature to use in the veto list: G4VSensitiveDetector * G4LogicalVolume::GetMasterSensitiveDetector()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:352:34
  t29.method("GetMasterSensitiveDetector", static_cast<G4VSensitiveDetector * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetMasterSensitiveDetector));

  DEBUG_MSG("Adding wrapper for G4VSolid * G4LogicalVolume::GetMasterSolid() (" __HERE__ ")");
  // signature to use in the veto list: G4VSolid * G4LogicalVolume::GetMasterSolid()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:354:22
  t29.method("GetMasterSolid", static_cast<G4VSolid * (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetMasterSolid));

  DEBUG_MSG("Adding wrapper for G4int G4LogicalVolume::GetInstanceID() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4LogicalVolume::GetInstanceID()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:357:18
  t29.method("GetInstanceID", static_cast<G4int (G4LogicalVolume::*)()  const>(&G4LogicalVolume::GetInstanceID));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::Clean() (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::Clean()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:363:17
  t29.method("G4LogicalVolume!Clean", static_cast<void (*)() >(&G4LogicalVolume::Clean));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::Lock() (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::Lock()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:366:17
  t29.method("Lock", static_cast<void (G4LogicalVolume::*)() >(&G4LogicalVolume::Lock));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::InitialiseWorker(G4LogicalVolume *, G4VSolid *, G4VSensitiveDetector *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::InitialiseWorker(G4LogicalVolume *, G4VSolid *, G4VSensitiveDetector *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:369:10
  t29.method("InitialiseWorker", static_cast<void (G4LogicalVolume::*)(G4LogicalVolume *, G4VSolid *, G4VSensitiveDetector *) >(&G4LogicalVolume::InitialiseWorker));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::TerminateWorker(G4LogicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::TerminateWorker(G4LogicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:374:10
  t29.method("TerminateWorker", static_cast<void (G4LogicalVolume::*)(G4LogicalVolume *) >(&G4LogicalVolume::TerminateWorker));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::AssignFieldManager(G4FieldManager *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::AssignFieldManager(G4FieldManager *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:378:10
  t29.method("AssignFieldManager", static_cast<void (G4LogicalVolume::*)(G4FieldManager *) >(&G4LogicalVolume::AssignFieldManager));

  DEBUG_MSG("Adding wrapper for G4VSolid * G4LogicalVolume::GetSolid(G4LVData &) (" __HERE__ ")");
  // signature to use in the veto list: G4VSolid * G4LogicalVolume::GetSolid(G4LVData &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:381:22
  t29.method("G4LogicalVolume!GetSolid", static_cast<G4VSolid * (*)(G4LVData &) >(&G4LogicalVolume::GetSolid));

  DEBUG_MSG("Adding wrapper for void G4LogicalVolume::SetSolid(G4LVData &, G4VSolid *) (" __HERE__ ")");
  // signature to use in the veto list: void G4LogicalVolume::SetSolid(G4LVData &, G4VSolid *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:382:17
  t29.method("G4LogicalVolume!SetSolid", static_cast<void (*)(G4LVData &, G4VSolid *) >(&G4LogicalVolume::SetSolid));

  DEBUG_MSG("Adding wrapper for G4bool G4LogicalVolume::ChangeDaughtersType(EVolume) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4LogicalVolume::ChangeDaughtersType(EVolume)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:385:12
  t29.method("ChangeDaughtersType", static_cast<G4bool (G4LogicalVolume::*)(EVolume) >(&G4LogicalVolume::ChangeDaughtersType));

  /* End of G4LogicalVolume class method wrappers
   **********************************************************************/

}