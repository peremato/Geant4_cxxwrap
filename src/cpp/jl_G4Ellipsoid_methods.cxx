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
void add_methods_for_G4Ellipsoid(jlcxx::Module& types, jlcxx::TypeWrapper<G4Ellipsoid>& t168) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4Ellipsoid
   */


  DEBUG_MSG("Adding wrapper for void G4Ellipsoid::G4Ellipsoid(const G4String &, G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:68:5
  t168.constructor<const G4String &, G4double, G4double, G4double>(/*finalize=*/true);
  t168.constructor<const G4String &, G4double, G4double, G4double, G4double>(/*finalize=*/true);
  t168.constructor<const G4String &, G4double, G4double, G4double, G4double, G4double>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4double G4Ellipsoid::GetDx() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Ellipsoid::GetDx()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:79:21
  t168.method("GetDx", static_cast<G4double (G4Ellipsoid::*)()  const>(&G4Ellipsoid::GetDx));

  DEBUG_MSG("Adding wrapper for G4double G4Ellipsoid::GetDy() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Ellipsoid::GetDy()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:80:21
  t168.method("GetDy", static_cast<G4double (G4Ellipsoid::*)()  const>(&G4Ellipsoid::GetDy));

  DEBUG_MSG("Adding wrapper for G4double G4Ellipsoid::GetDz() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Ellipsoid::GetDz()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:81:21
  t168.method("GetDz", static_cast<G4double (G4Ellipsoid::*)()  const>(&G4Ellipsoid::GetDz));

  DEBUG_MSG("Adding wrapper for G4double G4Ellipsoid::GetSemiAxisMax(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Ellipsoid::GetSemiAxisMax(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:82:21
  t168.method("GetSemiAxisMax", static_cast<G4double (G4Ellipsoid::*)(G4int)  const>(&G4Ellipsoid::GetSemiAxisMax));

  DEBUG_MSG("Adding wrapper for G4double G4Ellipsoid::GetZBottomCut() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Ellipsoid::GetZBottomCut()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:83:21
  t168.method("GetZBottomCut", static_cast<G4double (G4Ellipsoid::*)()  const>(&G4Ellipsoid::GetZBottomCut));

  DEBUG_MSG("Adding wrapper for G4double G4Ellipsoid::GetZTopCut() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Ellipsoid::GetZTopCut()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:84:21
  t168.method("GetZTopCut", static_cast<G4double (G4Ellipsoid::*)()  const>(&G4Ellipsoid::GetZTopCut));

  DEBUG_MSG("Adding wrapper for void G4Ellipsoid::SetSemiAxis(G4double, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Ellipsoid::SetSemiAxis(G4double, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:87:17
  t168.method("SetSemiAxis", static_cast<void (G4Ellipsoid::*)(G4double, G4double, G4double) >(&G4Ellipsoid::SetSemiAxis));

  DEBUG_MSG("Adding wrapper for void G4Ellipsoid::SetZCuts(G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Ellipsoid::SetZCuts(G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:88:17
  t168.method("SetZCuts", static_cast<void (G4Ellipsoid::*)(G4double, G4double) >(&G4Ellipsoid::SetZCuts));

  DEBUG_MSG("Adding wrapper for void G4Ellipsoid::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Ellipsoid::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:91:10
  t168.method("ComputeDimensions", static_cast<void (G4Ellipsoid::*)(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) >(&G4Ellipsoid::ComputeDimensions));

  DEBUG_MSG("Adding wrapper for void G4Ellipsoid::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4Ellipsoid::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:95:10
  t168.method("BoundingLimits", static_cast<void (G4Ellipsoid::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4Ellipsoid::BoundingLimits));

  DEBUG_MSG("Adding wrapper for EInside G4Ellipsoid::Inside(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: EInside G4Ellipsoid::Inside(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:101:13
  t168.method("Inside", static_cast<EInside (G4Ellipsoid::*)(const G4ThreeVector &)  const>(&G4Ellipsoid::Inside));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Ellipsoid::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Ellipsoid::SurfaceNormal(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:102:19
  t168.method("SurfaceNormal", static_cast<G4ThreeVector (G4Ellipsoid::*)(const G4ThreeVector &)  const>(&G4Ellipsoid::SurfaceNormal));

  DEBUG_MSG("Adding wrapper for G4double G4Ellipsoid::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Ellipsoid::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:103:14
  t168.method("DistanceToIn", static_cast<G4double (G4Ellipsoid::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Ellipsoid::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4Ellipsoid::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Ellipsoid::DistanceToIn(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:105:14
  t168.method("DistanceToIn", static_cast<G4double (G4Ellipsoid::*)(const G4ThreeVector &)  const>(&G4Ellipsoid::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4Ellipsoid::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Ellipsoid::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:106:14
  t168.method("DistanceToOut", static_cast<G4double (G4Ellipsoid::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4Ellipsoid::DistanceToOut));
  t168.method("DistanceToOut", [](G4Ellipsoid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a.DistanceToOut(arg0, arg1); });
  t168.method("DistanceToOut", [](G4Ellipsoid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a.DistanceToOut(arg0, arg1, arg2); });
  t168.method("DistanceToOut", [](G4Ellipsoid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a.DistanceToOut(arg0, arg1, arg2, arg3); });
  t168.method("DistanceToOut", [](G4Ellipsoid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a->DistanceToOut(arg0, arg1); });
  t168.method("DistanceToOut", [](G4Ellipsoid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a->DistanceToOut(arg0, arg1, arg2); });
  t168.method("DistanceToOut", [](G4Ellipsoid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a->DistanceToOut(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for G4double G4Ellipsoid::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Ellipsoid::DistanceToOut(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:111:14
  t168.method("DistanceToOut", static_cast<G4double (G4Ellipsoid::*)(const G4ThreeVector &)  const>(&G4Ellipsoid::DistanceToOut));

  DEBUG_MSG("Adding wrapper for G4GeometryType G4Ellipsoid::GetEntityType() (" __HERE__ ")");
  // signature to use in the veto list: G4GeometryType G4Ellipsoid::GetEntityType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:113:20
  t168.method("GetEntityType", static_cast<G4GeometryType (G4Ellipsoid::*)()  const>(&G4Ellipsoid::GetEntityType));

  DEBUG_MSG("Adding wrapper for G4VSolid * G4Ellipsoid::Clone() (" __HERE__ ")");
  // signature to use in the veto list: G4VSolid * G4Ellipsoid::Clone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:115:15
  t168.method("Clone", static_cast<G4VSolid * (G4Ellipsoid::*)()  const>(&G4Ellipsoid::Clone));

  DEBUG_MSG("Adding wrapper for G4double G4Ellipsoid::GetCubicVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Ellipsoid::GetCubicVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:119:14
  t168.method("GetCubicVolume", static_cast<G4double (G4Ellipsoid::*)() >(&G4Ellipsoid::GetCubicVolume));

  DEBUG_MSG("Adding wrapper for G4double G4Ellipsoid::GetSurfaceArea() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Ellipsoid::GetSurfaceArea()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:120:14
  t168.method("GetSurfaceArea", static_cast<G4double (G4Ellipsoid::*)() >(&G4Ellipsoid::GetSurfaceArea));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Ellipsoid::GetPointOnSurface() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Ellipsoid::GetPointOnSurface()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:122:19
  t168.method("GetPointOnSurface", static_cast<G4ThreeVector (G4Ellipsoid::*)()  const>(&G4Ellipsoid::GetPointOnSurface));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Ellipsoid::CreatePolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4Ellipsoid::CreatePolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:127:19
  t168.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4Ellipsoid::*)()  const>(&G4Ellipsoid::CreatePolyhedron));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Ellipsoid::GetPolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4Ellipsoid::GetPolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:128:19
  t168.method("GetPolyhedron", static_cast<G4Polyhedron * (G4Ellipsoid::*)()  const>(&G4Ellipsoid::GetPolyhedron));


  DEBUG_MSG("Adding wrapper for void G4Ellipsoid::G4Ellipsoid(const G4Ellipsoid &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:138:5
  t168.constructor<const G4Ellipsoid &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4Ellipsoid & G4Ellipsoid::operator=(const G4Ellipsoid &) (" __HERE__ ")");
  // signature to use in the veto list: G4Ellipsoid & G4Ellipsoid::operator=(const G4Ellipsoid &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Ellipsoid.hh:141:18
  t168.method("assign", static_cast<G4Ellipsoid & (G4Ellipsoid::*)(const G4Ellipsoid &) >(&G4Ellipsoid::operator=));

  /* End of G4Ellipsoid class method wrappers
   **********************************************************************/

}
