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
void add_methods_for_G4Torus(jlcxx::Module& types, jlcxx::TypeWrapper<G4Torus>& t194) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4Torus
   */


  DEBUG_MSG("Adding wrapper for void G4Torus::G4Torus(const G4String &, G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:96:5
  t194.constructor<const G4String &, G4double, G4double, G4double, G4double, G4double>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4double G4Torus::GetRmin() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::GetRmin()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:107:21
  t194.method("GetRmin", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetRmin));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::GetRmax() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::GetRmax()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:108:21
  t194.method("GetRmax", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetRmax));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::GetRtor() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::GetRtor()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:109:21
  t194.method("GetRtor", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetRtor));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::GetSPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::GetSPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:110:21
  t194.method("GetSPhi", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetSPhi));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::GetDPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::GetDPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:111:21
  t194.method("GetDPhi", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetDPhi));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::GetSinStartPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::GetSinStartPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:112:21
  t194.method("GetSinStartPhi", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetSinStartPhi));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::GetCosStartPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::GetCosStartPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:113:21
  t194.method("GetCosStartPhi", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetCosStartPhi));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::GetSinEndPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::GetSinEndPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:114:21
  t194.method("GetSinEndPhi", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetSinEndPhi));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::GetCosEndPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::GetCosEndPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:115:21
  t194.method("GetCosEndPhi", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetCosEndPhi));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::GetCubicVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::GetCubicVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:119:21
  t194.method("GetCubicVolume", static_cast<G4double (G4Torus::*)() >(&G4Torus::GetCubicVolume));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::GetSurfaceArea() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::GetSurfaceArea()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:120:21
  t194.method("GetSurfaceArea", static_cast<G4double (G4Torus::*)() >(&G4Torus::GetSurfaceArea));

  DEBUG_MSG("Adding wrapper for EInside G4Torus::Inside(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: EInside G4Torus::Inside(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:122:13
  t194.method("Inside", static_cast<EInside (G4Torus::*)(const G4ThreeVector &)  const>(&G4Torus::Inside));

  DEBUG_MSG("Adding wrapper for void G4Torus::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4Torus::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:123:10
  t194.method("BoundingLimits", static_cast<void (G4Torus::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4Torus::BoundingLimits));

  DEBUG_MSG("Adding wrapper for void G4Torus::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Torus::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:128:10
  t194.method("ComputeDimensions", static_cast<void (G4Torus::*)(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) >(&G4Torus::ComputeDimensions));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Torus::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Torus::SurfaceNormal(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:131:19
  t194.method("SurfaceNormal", static_cast<G4ThreeVector (G4Torus::*)(const G4ThreeVector &)  const>(&G4Torus::SurfaceNormal));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:132:14
  t194.method("DistanceToIn", static_cast<G4double (G4Torus::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Torus::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::DistanceToIn(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:133:14
  t194.method("DistanceToIn", static_cast<G4double (G4Torus::*)(const G4ThreeVector &)  const>(&G4Torus::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4Torus::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:134:14
  t194.method("DistanceToOut", static_cast<G4double (G4Torus::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4Torus::DistanceToOut));
  t194.method("DistanceToOut", [](G4Torus const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a.DistanceToOut(arg0, arg1); });
  t194.method("DistanceToOut", [](G4Torus const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a.DistanceToOut(arg0, arg1, arg2); });
  t194.method("DistanceToOut", [](G4Torus const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a.DistanceToOut(arg0, arg1, arg2, arg3); });
  t194.method("DistanceToOut", [](G4Torus const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a->DistanceToOut(arg0, arg1); });
  t194.method("DistanceToOut", [](G4Torus const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a->DistanceToOut(arg0, arg1, arg2); });
  t194.method("DistanceToOut", [](G4Torus const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a->DistanceToOut(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for G4double G4Torus::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Torus::DistanceToOut(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:138:14
  t194.method("DistanceToOut", static_cast<G4double (G4Torus::*)(const G4ThreeVector &)  const>(&G4Torus::DistanceToOut));

  DEBUG_MSG("Adding wrapper for G4GeometryType G4Torus::GetEntityType() (" __HERE__ ")");
  // signature to use in the veto list: G4GeometryType G4Torus::GetEntityType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:140:20
  t194.method("GetEntityType", static_cast<G4GeometryType (G4Torus::*)()  const>(&G4Torus::GetEntityType));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Torus::GetPointOnSurface() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Torus::GetPointOnSurface()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:142:19
  t194.method("GetPointOnSurface", static_cast<G4ThreeVector (G4Torus::*)()  const>(&G4Torus::GetPointOnSurface));

  DEBUG_MSG("Adding wrapper for G4VSolid * G4Torus::Clone() (" __HERE__ ")");
  // signature to use in the veto list: G4VSolid * G4Torus::Clone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:144:15
  t194.method("Clone", static_cast<G4VSolid * (G4Torus::*)()  const>(&G4Torus::Clone));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Torus::CreatePolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4Torus::CreatePolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:151:25
  t194.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4Torus::*)()  const>(&G4Torus::CreatePolyhedron));

  DEBUG_MSG("Adding wrapper for void G4Torus::SetAllParameters(G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Torus::SetAllParameters(G4double, G4double, G4double, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:155:10
  t194.method("SetAllParameters", static_cast<void (G4Torus::*)(G4double, G4double, G4double, G4double, G4double) >(&G4Torus::SetAllParameters));


  DEBUG_MSG("Adding wrapper for void G4Torus::G4Torus(const G4Torus &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:163:5
  t194.constructor<const G4Torus &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4Torus & G4Torus::operator=(const G4Torus &) (" __HERE__ ")");
  // signature to use in the veto list: G4Torus & G4Torus::operator=(const G4Torus &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Torus.hh:164:14
  t194.method("assign", static_cast<G4Torus & (G4Torus::*)(const G4Torus &) >(&G4Torus::operator=));

  /* End of G4Torus class method wrappers
   **********************************************************************/

}
