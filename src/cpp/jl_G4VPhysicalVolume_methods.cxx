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
void add_methods_for_G4VPhysicalVolume(jlcxx::Module& types, jlcxx::TypeWrapper<G4VPhysicalVolume>& t28) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4VPhysicalVolume
   */

  types.set_override_module(jl_base_module);

  DEBUG_MSG("Adding wrapper for G4bool G4VPhysicalVolume::operator==(const G4VPhysicalVolume &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VPhysicalVolume::operator==(const G4VPhysicalVolume &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:107:19
  t28.method("==", static_cast<G4bool (G4VPhysicalVolume::*)(const G4VPhysicalVolume &)  const>(&G4VPhysicalVolume::operator==));

  types.unset_override_module();

  DEBUG_MSG("Adding wrapper for G4RotationMatrix * G4VPhysicalVolume::GetObjectRotation() (" __HERE__ ")");
  // signature to use in the veto list: G4RotationMatrix * G4VPhysicalVolume::GetObjectRotation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:117:23
  t28.method("GetObjectRotation", static_cast<G4RotationMatrix * (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetObjectRotation));

  DEBUG_MSG("Adding wrapper for G4RotationMatrix G4VPhysicalVolume::GetObjectRotationValue() (" __HERE__ ")");
  // signature to use in the veto list: G4RotationMatrix G4VPhysicalVolume::GetObjectRotationValue()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:118:23
  t28.method("GetObjectRotationValue", static_cast<G4RotationMatrix (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetObjectRotationValue));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4VPhysicalVolume::GetObjectTranslation() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4VPhysicalVolume::GetObjectTranslation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:119:20
  t28.method("GetObjectTranslation", static_cast<G4ThreeVector (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetObjectTranslation));

  DEBUG_MSG("Adding wrapper for const G4RotationMatrix * G4VPhysicalVolume::GetFrameRotation() (" __HERE__ ")");
  // signature to use in the veto list: const G4RotationMatrix * G4VPhysicalVolume::GetFrameRotation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:121:29
  t28.method("GetFrameRotation", static_cast<const G4RotationMatrix * (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetFrameRotation));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4VPhysicalVolume::GetFrameTranslation() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4VPhysicalVolume::GetFrameTranslation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:122:19
  t28.method("GetFrameTranslation", static_cast<G4ThreeVector (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetFrameTranslation));

  DEBUG_MSG("Adding wrapper for const G4ThreeVector G4VPhysicalVolume::GetTranslation() (" __HERE__ ")");
  // signature to use in the veto list: const G4ThreeVector G4VPhysicalVolume::GetTranslation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:128:25
  t28.method("GetTranslation", static_cast<const G4ThreeVector (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetTranslation));

  DEBUG_MSG("Adding wrapper for const G4RotationMatrix * G4VPhysicalVolume::GetRotation() (" __HERE__ ")");
  // signature to use in the veto list: const G4RotationMatrix * G4VPhysicalVolume::GetRotation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:129:29
  t28.method("GetRotation", static_cast<const G4RotationMatrix * (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetRotation));

  DEBUG_MSG("Adding wrapper for void G4VPhysicalVolume::SetTranslation(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4VPhysicalVolume::SetTranslation(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:135:10
  t28.method("SetTranslation", static_cast<void (G4VPhysicalVolume::*)(const G4ThreeVector &) >(&G4VPhysicalVolume::SetTranslation));

  DEBUG_MSG("Adding wrapper for G4RotationMatrix * G4VPhysicalVolume::GetRotation() (" __HERE__ ")");
  // signature to use in the veto list: G4RotationMatrix * G4VPhysicalVolume::GetRotation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:136:23
  t28.method("GetRotation", static_cast<G4RotationMatrix * (G4VPhysicalVolume::*)() >(&G4VPhysicalVolume::GetRotation));

  DEBUG_MSG("Adding wrapper for void G4VPhysicalVolume::SetRotation(G4RotationMatrix *) (" __HERE__ ")");
  // signature to use in the veto list: void G4VPhysicalVolume::SetRotation(G4RotationMatrix *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:137:10
  t28.method("SetRotation", static_cast<void (G4VPhysicalVolume::*)(G4RotationMatrix *) >(&G4VPhysicalVolume::SetRotation));

  DEBUG_MSG("Adding wrapper for G4LogicalVolume * G4VPhysicalVolume::GetLogicalVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4LogicalVolume * G4VPhysicalVolume::GetLogicalVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:142:29
  t28.method("GetLogicalVolume", static_cast<G4LogicalVolume * (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetLogicalVolume));

  DEBUG_MSG("Adding wrapper for void G4VPhysicalVolume::SetLogicalVolume(G4LogicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4VPhysicalVolume::SetLogicalVolume(G4LogicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:144:17
  t28.method("SetLogicalVolume", static_cast<void (G4VPhysicalVolume::*)(G4LogicalVolume *) >(&G4VPhysicalVolume::SetLogicalVolume));

  DEBUG_MSG("Adding wrapper for G4LogicalVolume * G4VPhysicalVolume::GetMotherLogical() (" __HERE__ ")");
  // signature to use in the veto list: G4LogicalVolume * G4VPhysicalVolume::GetMotherLogical()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:147:29
  t28.method("GetMotherLogical", static_cast<G4LogicalVolume * (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetMotherLogical));

  DEBUG_MSG("Adding wrapper for void G4VPhysicalVolume::SetMotherLogical(G4LogicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4VPhysicalVolume::SetMotherLogical(G4LogicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:149:17
  t28.method("SetMotherLogical", static_cast<void (G4VPhysicalVolume::*)(G4LogicalVolume *) >(&G4VPhysicalVolume::SetMotherLogical));

  DEBUG_MSG("Adding wrapper for const G4String & G4VPhysicalVolume::GetName() (" __HERE__ ")");
  // signature to use in the veto list: const G4String & G4VPhysicalVolume::GetName()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:152:28
  t28.method("GetName", static_cast<const G4String & (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetName));

  DEBUG_MSG("Adding wrapper for void G4VPhysicalVolume::SetName(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4VPhysicalVolume::SetName(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:154:10
  t28.method("SetName", static_cast<void (G4VPhysicalVolume::*)(const G4String &) >(&G4VPhysicalVolume::SetName));

  DEBUG_MSG("Adding wrapper for G4int G4VPhysicalVolume::GetMultiplicity() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VPhysicalVolume::GetMultiplicity()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:157:19
  t28.method("GetMultiplicity", static_cast<G4int (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetMultiplicity));

  DEBUG_MSG("Adding wrapper for EVolume G4VPhysicalVolume::VolumeType() (" __HERE__ ")");
  // signature to use in the veto list: EVolume G4VPhysicalVolume::VolumeType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:163:21
  t28.method("VolumeType", static_cast<EVolume (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::VolumeType));

  DEBUG_MSG("Adding wrapper for G4bool G4VPhysicalVolume::IsMany() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VPhysicalVolume::IsMany()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:165:20
  t28.method("IsMany", static_cast<G4bool (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::IsMany));

  DEBUG_MSG("Adding wrapper for G4int G4VPhysicalVolume::GetCopyNo() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VPhysicalVolume::GetCopyNo()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:167:19
  t28.method("GetCopyNo", static_cast<G4int (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetCopyNo));

  DEBUG_MSG("Adding wrapper for void G4VPhysicalVolume::SetCopyNo(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4VPhysicalVolume::SetCopyNo(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:169:19
  t28.method("SetCopyNo", static_cast<void (G4VPhysicalVolume::*)(G4int) >(&G4VPhysicalVolume::SetCopyNo));

  DEBUG_MSG("Adding wrapper for G4bool G4VPhysicalVolume::IsReplicated() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VPhysicalVolume::IsReplicated()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:171:20
  t28.method("IsReplicated", static_cast<G4bool (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::IsReplicated));

  DEBUG_MSG("Adding wrapper for G4bool G4VPhysicalVolume::IsParameterised() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VPhysicalVolume::IsParameterised()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:174:20
  t28.method("IsParameterised", static_cast<G4bool (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::IsParameterised));

  DEBUG_MSG("Adding wrapper for G4VPVParameterisation * G4VPhysicalVolume::GetParameterisation() (" __HERE__ ")");
  // signature to use in the veto list: G4VPVParameterisation * G4VPhysicalVolume::GetParameterisation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:177:36
  t28.method("GetParameterisation", static_cast<G4VPVParameterisation * (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetParameterisation));

  DEBUG_MSG("Adding wrapper for void G4VPhysicalVolume::GetReplicationData(EAxis &, G4int &, G4double &, G4double &, G4bool &) (" __HERE__ ")");
  // signature to use in the veto list: void G4VPhysicalVolume::GetReplicationData(EAxis &, G4int &, G4double &, G4double &, G4bool &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:180:18
  t28.method("GetReplicationData", static_cast<void (G4VPhysicalVolume::*)(EAxis &, G4int &, G4double &, G4double &, G4bool &)  const>(&G4VPhysicalVolume::GetReplicationData));

  DEBUG_MSG("Adding wrapper for G4bool G4VPhysicalVolume::IsRegularStructure() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VPhysicalVolume::IsRegularStructure()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:186:20
  t28.method("IsRegularStructure", static_cast<G4bool (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::IsRegularStructure));

  DEBUG_MSG("Adding wrapper for G4int G4VPhysicalVolume::GetRegularStructureId() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VPhysicalVolume::GetRegularStructureId()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:188:19
  t28.method("GetRegularStructureId", static_cast<G4int (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetRegularStructureId));

  DEBUG_MSG("Adding wrapper for G4bool G4VPhysicalVolume::CheckOverlaps(G4int, G4double, G4bool, G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VPhysicalVolume::CheckOverlaps(G4int, G4double, G4bool, G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:193:20
  t28.method("CheckOverlaps", static_cast<G4bool (G4VPhysicalVolume::*)(G4int, G4double, G4bool, G4int) >(&G4VPhysicalVolume::CheckOverlaps));
  t28.method("CheckOverlaps", [](G4VPhysicalVolume& a)->G4bool{ return a.CheckOverlaps(); });
  t28.method("CheckOverlaps", [](G4VPhysicalVolume& a, G4int arg0)->G4bool{ return a.CheckOverlaps(arg0); });
  t28.method("CheckOverlaps", [](G4VPhysicalVolume& a, G4int arg0, G4double arg1)->G4bool{ return a.CheckOverlaps(arg0, arg1); });
  t28.method("CheckOverlaps", [](G4VPhysicalVolume& a, G4int arg0, G4double arg1, G4bool arg2)->G4bool{ return a.CheckOverlaps(arg0, arg1, arg2); });
  t28.method("CheckOverlaps", [](G4VPhysicalVolume* a)->G4bool{ return a->CheckOverlaps(); });
  t28.method("CheckOverlaps", [](G4VPhysicalVolume* a, G4int arg0)->G4bool{ return a->CheckOverlaps(arg0); });
  t28.method("CheckOverlaps", [](G4VPhysicalVolume* a, G4int arg0, G4double arg1)->G4bool{ return a->CheckOverlaps(arg0, arg1); });
  t28.method("CheckOverlaps", [](G4VPhysicalVolume* a, G4int arg0, G4double arg1, G4bool arg2)->G4bool{ return a->CheckOverlaps(arg0, arg1, arg2); });

  DEBUG_MSG("Adding wrapper for G4int G4VPhysicalVolume::GetInstanceID() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VPhysicalVolume::GetInstanceID()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:208:18
  t28.method("GetInstanceID", static_cast<G4int (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::GetInstanceID));

  DEBUG_MSG("Adding wrapper for void G4VPhysicalVolume::Clean() (" __HERE__ ")");
  // signature to use in the veto list: void G4VPhysicalVolume::Clean()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:214:17
  t28.method("G4VPhysicalVolume!Clean", static_cast<void (*)() >(&G4VPhysicalVolume::Clean));

  DEBUG_MSG("Adding wrapper for EVolume G4VPhysicalVolume::DeduceVolumeType() (" __HERE__ ")");
  // signature to use in the veto list: EVolume G4VPhysicalVolume::DeduceVolumeType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VPhysicalVolume.hh:217:20
  t28.method("DeduceVolumeType", static_cast<EVolume (G4VPhysicalVolume::*)()  const>(&G4VPhysicalVolume::DeduceVolumeType));

  /* End of G4VPhysicalVolume class method wrappers
   **********************************************************************/

}