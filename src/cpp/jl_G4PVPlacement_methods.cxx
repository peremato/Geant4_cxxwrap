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
void add_methods_for_G4PVPlacement(jlcxx::Module& types, jlcxx::TypeWrapper<G4PVPlacement>& t154) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4PVPlacement
   */


  DEBUG_MSG("Adding wrapper for void G4PVPlacement::G4PVPlacement(G4RotationMatrix *, const G4ThreeVector &, G4LogicalVolume *, const G4String &, G4LogicalVolume *, G4bool, G4int, G4bool) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:47:5
  t154.constructor<G4RotationMatrix *, const G4ThreeVector &, G4LogicalVolume *, const G4String &, G4LogicalVolume *, G4bool, G4int>(/*finalize=*/false);
  t154.constructor<G4RotationMatrix *, const G4ThreeVector &, G4LogicalVolume *, const G4String &, G4LogicalVolume *, G4bool, G4int, G4bool>(/*finalize=*/false);


  DEBUG_MSG("Adding wrapper for void G4PVPlacement::G4PVPlacement(const G4Transform3D &, G4LogicalVolume *, const G4String &, G4LogicalVolume *, G4bool, G4int, G4bool) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:69:5
  t154.constructor<const G4Transform3D &, G4LogicalVolume *, const G4String &, G4LogicalVolume *, G4bool, G4int>(/*finalize=*/false);
  t154.constructor<const G4Transform3D &, G4LogicalVolume *, const G4String &, G4LogicalVolume *, G4bool, G4int, G4bool>(/*finalize=*/false);


  DEBUG_MSG("Adding wrapper for void G4PVPlacement::G4PVPlacement(G4RotationMatrix *, const G4ThreeVector &, const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, G4bool, G4int, G4bool) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:90:5
  t154.constructor<G4RotationMatrix *, const G4ThreeVector &, const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, G4bool, G4int>(/*finalize=*/false);
  t154.constructor<G4RotationMatrix *, const G4ThreeVector &, const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, G4bool, G4int, G4bool>(/*finalize=*/false);


  DEBUG_MSG("Adding wrapper for void G4PVPlacement::G4PVPlacement(const G4Transform3D &, const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, G4bool, G4int, G4bool) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:102:5
  t154.constructor<const G4Transform3D &, const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, G4bool, G4int>(/*finalize=*/false);
  t154.constructor<const G4Transform3D &, const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, G4bool, G4int, G4bool>(/*finalize=*/false);

  DEBUG_MSG("Adding wrapper for G4int G4PVPlacement::GetCopyNo() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4PVPlacement::GetCopyNo()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:117:18
  t154.method("GetCopyNo", static_cast<G4int (G4PVPlacement::*)()  const>(&G4PVPlacement::GetCopyNo));

  DEBUG_MSG("Adding wrapper for void G4PVPlacement::SetCopyNo(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4PVPlacement::SetCopyNo(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:119:10
  t154.method("SetCopyNo", static_cast<void (G4PVPlacement::*)(G4int) >(&G4PVPlacement::SetCopyNo));

  DEBUG_MSG("Adding wrapper for G4bool G4PVPlacement::CheckOverlaps(G4int, G4double, G4bool, G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4PVPlacement::CheckOverlaps(G4int, G4double, G4bool, G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:122:12
  t154.method("CheckOverlaps", static_cast<G4bool (G4PVPlacement::*)(G4int, G4double, G4bool, G4int) >(&G4PVPlacement::CheckOverlaps));
  t154.method("CheckOverlaps", [](G4PVPlacement& a)->G4bool{ return a.CheckOverlaps(); });
  t154.method("CheckOverlaps", [](G4PVPlacement& a, G4int arg0)->G4bool{ return a.CheckOverlaps(arg0); });
  t154.method("CheckOverlaps", [](G4PVPlacement& a, G4int arg0, G4double arg1)->G4bool{ return a.CheckOverlaps(arg0, arg1); });
  t154.method("CheckOverlaps", [](G4PVPlacement& a, G4int arg0, G4double arg1, G4bool arg2)->G4bool{ return a.CheckOverlaps(arg0, arg1, arg2); });
  t154.method("CheckOverlaps", [](G4PVPlacement* a)->G4bool{ return a->CheckOverlaps(); });
  t154.method("CheckOverlaps", [](G4PVPlacement* a, G4int arg0)->G4bool{ return a->CheckOverlaps(arg0); });
  t154.method("CheckOverlaps", [](G4PVPlacement* a, G4int arg0, G4double arg1)->G4bool{ return a->CheckOverlaps(arg0, arg1); });
  t154.method("CheckOverlaps", [](G4PVPlacement* a, G4int arg0, G4double arg1, G4bool arg2)->G4bool{ return a->CheckOverlaps(arg0, arg1, arg2); });

  DEBUG_MSG("Adding wrapper for G4bool G4PVPlacement::IsMany() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4PVPlacement::IsMany()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:143:12
  t154.method("IsMany", static_cast<G4bool (G4PVPlacement::*)()  const>(&G4PVPlacement::IsMany));

  DEBUG_MSG("Adding wrapper for G4bool G4PVPlacement::IsReplicated() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4PVPlacement::IsReplicated()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:144:12
  t154.method("IsReplicated", static_cast<G4bool (G4PVPlacement::*)()  const>(&G4PVPlacement::IsReplicated));

  DEBUG_MSG("Adding wrapper for G4bool G4PVPlacement::IsParameterised() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4PVPlacement::IsParameterised()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:145:12
  t154.method("IsParameterised", static_cast<G4bool (G4PVPlacement::*)()  const>(&G4PVPlacement::IsParameterised));

  DEBUG_MSG("Adding wrapper for G4VPVParameterisation * G4PVPlacement::GetParameterisation() (" __HERE__ ")");
  // signature to use in the veto list: G4VPVParameterisation * G4PVPlacement::GetParameterisation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:146:28
  t154.method("GetParameterisation", static_cast<G4VPVParameterisation * (G4PVPlacement::*)()  const>(&G4PVPlacement::GetParameterisation));

  DEBUG_MSG("Adding wrapper for void G4PVPlacement::GetReplicationData(EAxis &, G4int &, G4double &, G4double &, G4bool &) (" __HERE__ ")");
  // signature to use in the veto list: void G4PVPlacement::GetReplicationData(EAxis &, G4int &, G4double &, G4double &, G4bool &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:147:10
  t154.method("GetReplicationData", static_cast<void (G4PVPlacement::*)(EAxis &, G4int &, G4double &, G4double &, G4bool &)  const>(&G4PVPlacement::GetReplicationData));

  DEBUG_MSG("Adding wrapper for G4bool G4PVPlacement::IsRegularStructure() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4PVPlacement::IsRegularStructure()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:152:12
  t154.method("IsRegularStructure", static_cast<G4bool (G4PVPlacement::*)()  const>(&G4PVPlacement::IsRegularStructure));

  DEBUG_MSG("Adding wrapper for G4int G4PVPlacement::GetRegularStructureId() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4PVPlacement::GetRegularStructureId()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:153:11
  t154.method("GetRegularStructureId", static_cast<G4int (G4PVPlacement::*)()  const>(&G4PVPlacement::GetRegularStructureId));

  DEBUG_MSG("Adding wrapper for EVolume G4PVPlacement::VolumeType() (" __HERE__ ")");
  // signature to use in the veto list: EVolume G4PVPlacement::VolumeType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVPlacement.hh:155:13
  t154.method("VolumeType", static_cast<EVolume (G4PVPlacement::*)()  const>(&G4PVPlacement::VolumeType));

  /* End of G4PVPlacement class method wrappers
   **********************************************************************/

}
