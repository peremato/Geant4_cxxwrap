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
void add_methods_for_G4UnionSolid(jlcxx::Module& types, jlcxx::TypeWrapper<G4UnionSolid>& t210) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4UnionSolid
   */


  DEBUG_MSG("Adding wrapper for void G4UnionSolid::G4UnionSolid(const G4String &, G4VSolid *, G4VSolid *) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:49:5
  t210.constructor<const G4String &, G4VSolid *, G4VSolid *>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void G4UnionSolid::G4UnionSolid(const G4String &, G4VSolid *, G4VSolid *, G4RotationMatrix *, const G4ThreeVector &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:53:5
  t210.constructor<const G4String &, G4VSolid *, G4VSolid *, G4RotationMatrix *, const G4ThreeVector &>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void G4UnionSolid::G4UnionSolid(const G4String &, G4VSolid *, G4VSolid *, const G4Transform3D &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:59:5
  t210.constructor<const G4String &, G4VSolid *, G4VSolid *, const G4Transform3D &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4GeometryType G4UnionSolid::GetEntityType() (" __HERE__ ")");
  // signature to use in the veto list: G4GeometryType G4UnionSolid::GetEntityType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:66:21
  t210.method("GetEntityType", static_cast<G4GeometryType (G4UnionSolid::*)()  const>(&G4UnionSolid::GetEntityType));

  DEBUG_MSG("Adding wrapper for G4VSolid * G4UnionSolid::Clone() (" __HERE__ ")");
  // signature to use in the veto list: G4VSolid * G4UnionSolid::Clone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:68:15
  t210.method("Clone", static_cast<G4VSolid * (G4UnionSolid::*)()  const>(&G4UnionSolid::Clone));


  DEBUG_MSG("Adding wrapper for void G4UnionSolid::G4UnionSolid(const G4UnionSolid &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:75:5
  t210.constructor<const G4UnionSolid &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4UnionSolid & G4UnionSolid::operator=(const G4UnionSolid &) (" __HERE__ ")");
  // signature to use in the veto list: G4UnionSolid & G4UnionSolid::operator=(const G4UnionSolid &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:76:19
  t210.method("assign", static_cast<G4UnionSolid & (G4UnionSolid::*)(const G4UnionSolid &) >(&G4UnionSolid::operator=));

  DEBUG_MSG("Adding wrapper for void G4UnionSolid::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4UnionSolid::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:79:10
  t210.method("BoundingLimits", static_cast<void (G4UnionSolid::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4UnionSolid::BoundingLimits));

  DEBUG_MSG("Adding wrapper for EInside G4UnionSolid::Inside(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: EInside G4UnionSolid::Inside(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:86:13
  t210.method("Inside", static_cast<EInside (G4UnionSolid::*)(const G4ThreeVector &)  const>(&G4UnionSolid::Inside));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4UnionSolid::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4UnionSolid::SurfaceNormal(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:88:19
  t210.method("SurfaceNormal", static_cast<G4ThreeVector (G4UnionSolid::*)(const G4ThreeVector &)  const>(&G4UnionSolid::SurfaceNormal));

  DEBUG_MSG("Adding wrapper for G4double G4UnionSolid::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4UnionSolid::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:90:14
  t210.method("DistanceToIn", static_cast<G4double (G4UnionSolid::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4UnionSolid::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4UnionSolid::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4UnionSolid::DistanceToIn(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:93:14
  t210.method("DistanceToIn", static_cast<G4double (G4UnionSolid::*)(const G4ThreeVector &)  const>(&G4UnionSolid::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4UnionSolid::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4UnionSolid::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:95:14
  t210.method("DistanceToOut", static_cast<G4double (G4UnionSolid::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4UnionSolid::DistanceToOut));
  t210.method("DistanceToOut", [](G4UnionSolid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a.DistanceToOut(arg0, arg1); });
  t210.method("DistanceToOut", [](G4UnionSolid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a.DistanceToOut(arg0, arg1, arg2); });
  t210.method("DistanceToOut", [](G4UnionSolid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a.DistanceToOut(arg0, arg1, arg2, arg3); });
  t210.method("DistanceToOut", [](G4UnionSolid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a->DistanceToOut(arg0, arg1); });
  t210.method("DistanceToOut", [](G4UnionSolid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a->DistanceToOut(arg0, arg1, arg2); });
  t210.method("DistanceToOut", [](G4UnionSolid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a->DistanceToOut(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for G4double G4UnionSolid::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4UnionSolid::DistanceToOut(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:101:14
  t210.method("DistanceToOut", static_cast<G4double (G4UnionSolid::*)(const G4ThreeVector &)  const>(&G4UnionSolid::DistanceToOut));

  DEBUG_MSG("Adding wrapper for void G4UnionSolid::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UnionSolid::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:104:10
  t210.method("ComputeDimensions", static_cast<void (G4UnionSolid::*)(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) >(&G4UnionSolid::ComputeDimensions));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4UnionSolid::CreatePolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4UnionSolid::CreatePolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:109:19
  t210.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4UnionSolid::*)()  const>(&G4UnionSolid::CreatePolyhedron));

  DEBUG_MSG("Adding wrapper for G4double G4UnionSolid::GetCubicVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4UnionSolid::GetCubicVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UnionSolid.hh:111:22
  t210.method("GetCubicVolume", static_cast<G4double (G4UnionSolid::*)() >(&G4UnionSolid::GetCubicVolume));

  /* End of G4UnionSolid class method wrappers
   **********************************************************************/

}
