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
void add_methods_for_G4Cons(jlcxx::Module& types, jlcxx::TypeWrapper<G4Cons>& t165) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4Cons
   */


  DEBUG_MSG("Adding wrapper for void G4Cons::G4Cons(const G4String &, G4double, G4double, G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:81:5
  t165.constructor<const G4String &, G4double, G4double, G4double, G4double, G4double, G4double, G4double>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetInnerRadiusMinusZ() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetInnerRadiusMinusZ()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:95:21
  t165.method("GetInnerRadiusMinusZ", static_cast<G4double (G4Cons::*)()  const>(&G4Cons::GetInnerRadiusMinusZ));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetOuterRadiusMinusZ() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetOuterRadiusMinusZ()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:96:21
  t165.method("GetOuterRadiusMinusZ", static_cast<G4double (G4Cons::*)()  const>(&G4Cons::GetOuterRadiusMinusZ));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetInnerRadiusPlusZ() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetInnerRadiusPlusZ()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:97:21
  t165.method("GetInnerRadiusPlusZ", static_cast<G4double (G4Cons::*)()  const>(&G4Cons::GetInnerRadiusPlusZ));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetOuterRadiusPlusZ() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetOuterRadiusPlusZ()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:98:21
  t165.method("GetOuterRadiusPlusZ", static_cast<G4double (G4Cons::*)()  const>(&G4Cons::GetOuterRadiusPlusZ));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetZHalfLength() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetZHalfLength()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:99:21
  t165.method("GetZHalfLength", static_cast<G4double (G4Cons::*)()  const>(&G4Cons::GetZHalfLength));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetStartPhiAngle() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetStartPhiAngle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:100:21
  t165.method("GetStartPhiAngle", static_cast<G4double (G4Cons::*)()  const>(&G4Cons::GetStartPhiAngle));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetDeltaPhiAngle() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetDeltaPhiAngle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:101:21
  t165.method("GetDeltaPhiAngle", static_cast<G4double (G4Cons::*)()  const>(&G4Cons::GetDeltaPhiAngle));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetSinStartPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetSinStartPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:102:21
  t165.method("GetSinStartPhi", static_cast<G4double (G4Cons::*)()  const>(&G4Cons::GetSinStartPhi));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetCosStartPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetCosStartPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:103:21
  t165.method("GetCosStartPhi", static_cast<G4double (G4Cons::*)()  const>(&G4Cons::GetCosStartPhi));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetSinEndPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetSinEndPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:104:21
  t165.method("GetSinEndPhi", static_cast<G4double (G4Cons::*)()  const>(&G4Cons::GetSinEndPhi));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetCosEndPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetCosEndPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:105:21
  t165.method("GetCosEndPhi", static_cast<G4double (G4Cons::*)()  const>(&G4Cons::GetCosEndPhi));

  DEBUG_MSG("Adding wrapper for void G4Cons::SetInnerRadiusMinusZ(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Cons::SetInnerRadiusMinusZ(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:109:17
  t165.method("SetInnerRadiusMinusZ", static_cast<void (G4Cons::*)(G4double) >(&G4Cons::SetInnerRadiusMinusZ));

  DEBUG_MSG("Adding wrapper for void G4Cons::SetOuterRadiusMinusZ(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Cons::SetOuterRadiusMinusZ(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:110:17
  t165.method("SetOuterRadiusMinusZ", static_cast<void (G4Cons::*)(G4double) >(&G4Cons::SetOuterRadiusMinusZ));

  DEBUG_MSG("Adding wrapper for void G4Cons::SetInnerRadiusPlusZ(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Cons::SetInnerRadiusPlusZ(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:111:17
  t165.method("SetInnerRadiusPlusZ", static_cast<void (G4Cons::*)(G4double) >(&G4Cons::SetInnerRadiusPlusZ));

  DEBUG_MSG("Adding wrapper for void G4Cons::SetOuterRadiusPlusZ(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Cons::SetOuterRadiusPlusZ(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:112:17
  t165.method("SetOuterRadiusPlusZ", static_cast<void (G4Cons::*)(G4double) >(&G4Cons::SetOuterRadiusPlusZ));

  DEBUG_MSG("Adding wrapper for void G4Cons::SetZHalfLength(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Cons::SetZHalfLength(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:113:17
  t165.method("SetZHalfLength", static_cast<void (G4Cons::*)(G4double) >(&G4Cons::SetZHalfLength));

  DEBUG_MSG("Adding wrapper for void G4Cons::SetStartPhiAngle(G4double, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4Cons::SetStartPhiAngle(G4double, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:114:17
  t165.method("SetStartPhiAngle", static_cast<void (G4Cons::*)(G4double, G4bool) >(&G4Cons::SetStartPhiAngle));
  t165.method("SetStartPhiAngle", [](G4Cons& a, G4double arg0)->void{ a.SetStartPhiAngle(arg0); });
  t165.method("SetStartPhiAngle", [](G4Cons* a, G4double arg0)->void{ a->SetStartPhiAngle(arg0); });

  DEBUG_MSG("Adding wrapper for void G4Cons::SetDeltaPhiAngle(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Cons::SetDeltaPhiAngle(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:115:17
  t165.method("SetDeltaPhiAngle", static_cast<void (G4Cons::*)(G4double) >(&G4Cons::SetDeltaPhiAngle));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetCubicVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetCubicVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:119:21
  t165.method("GetCubicVolume", static_cast<G4double (G4Cons::*)() >(&G4Cons::GetCubicVolume));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::GetSurfaceArea() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::GetSurfaceArea()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:120:21
  t165.method("GetSurfaceArea", static_cast<G4double (G4Cons::*)() >(&G4Cons::GetSurfaceArea));

  DEBUG_MSG("Adding wrapper for void G4Cons::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Cons::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:122:10
  t165.method("ComputeDimensions", static_cast<void (G4Cons::*)(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) >(&G4Cons::ComputeDimensions));

  DEBUG_MSG("Adding wrapper for void G4Cons::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4Cons::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:126:10
  t165.method("BoundingLimits", static_cast<void (G4Cons::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4Cons::BoundingLimits));

  DEBUG_MSG("Adding wrapper for EInside G4Cons::Inside(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: EInside G4Cons::Inside(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:133:13
  t165.method("Inside", static_cast<EInside (G4Cons::*)(const G4ThreeVector &)  const>(&G4Cons::Inside));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Cons::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Cons::SurfaceNormal(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:135:19
  t165.method("SurfaceNormal", static_cast<G4ThreeVector (G4Cons::*)(const G4ThreeVector &)  const>(&G4Cons::SurfaceNormal));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:137:14
  t165.method("DistanceToIn", static_cast<G4double (G4Cons::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Cons::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::DistanceToIn(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:139:14
  t165.method("DistanceToIn", static_cast<G4double (G4Cons::*)(const G4ThreeVector &)  const>(&G4Cons::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4Cons::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:140:14
  t165.method("DistanceToOut", static_cast<G4double (G4Cons::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4Cons::DistanceToOut));
  t165.method("DistanceToOut", [](G4Cons const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a.DistanceToOut(arg0, arg1); });
  t165.method("DistanceToOut", [](G4Cons const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a.DistanceToOut(arg0, arg1, arg2); });
  t165.method("DistanceToOut", [](G4Cons const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a.DistanceToOut(arg0, arg1, arg2, arg3); });
  t165.method("DistanceToOut", [](G4Cons const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a->DistanceToOut(arg0, arg1); });
  t165.method("DistanceToOut", [](G4Cons const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a->DistanceToOut(arg0, arg1, arg2); });
  t165.method("DistanceToOut", [](G4Cons const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a->DistanceToOut(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for G4double G4Cons::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Cons::DistanceToOut(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:145:14
  t165.method("DistanceToOut", static_cast<G4double (G4Cons::*)(const G4ThreeVector &)  const>(&G4Cons::DistanceToOut));

  DEBUG_MSG("Adding wrapper for G4GeometryType G4Cons::GetEntityType() (" __HERE__ ")");
  // signature to use in the veto list: G4GeometryType G4Cons::GetEntityType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:147:20
  t165.method("GetEntityType", static_cast<G4GeometryType (G4Cons::*)()  const>(&G4Cons::GetEntityType));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Cons::GetPointOnSurface() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Cons::GetPointOnSurface()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:149:19
  t165.method("GetPointOnSurface", static_cast<G4ThreeVector (G4Cons::*)()  const>(&G4Cons::GetPointOnSurface));

  DEBUG_MSG("Adding wrapper for G4VSolid * G4Cons::Clone() (" __HERE__ ")");
  // signature to use in the veto list: G4VSolid * G4Cons::Clone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:151:15
  t165.method("Clone", static_cast<G4VSolid * (G4Cons::*)()  const>(&G4Cons::Clone));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Cons::CreatePolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4Cons::CreatePolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:158:19
  t165.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4Cons::*)()  const>(&G4Cons::CreatePolyhedron));


  DEBUG_MSG("Adding wrapper for void G4Cons::G4Cons(const G4Cons &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:168:5
  t165.constructor<const G4Cons &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4Cons & G4Cons::operator=(const G4Cons &) (" __HERE__ ")");
  // signature to use in the veto list: G4Cons & G4Cons::operator=(const G4Cons &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Cons.hh:169:13
  t165.method("assign", static_cast<G4Cons & (G4Cons::*)(const G4Cons &) >(&G4Cons::operator=));

  /* End of G4Cons class method wrappers
   **********************************************************************/

}
