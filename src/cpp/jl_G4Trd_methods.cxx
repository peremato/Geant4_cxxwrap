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
void add_methods_for_G4Trd(jlcxx::Module& types, jlcxx::TypeWrapper<G4Trd>& t152) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4Trd
   */


  DEBUG_MSG("Adding wrapper for void G4Trd::G4Trd(const G4String &, G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:66:5
  t152.constructor<const G4String &, G4double, G4double, G4double, G4double, G4double>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4double G4Trd::GetXHalfLength1() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Trd::GetXHalfLength1()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:79:21
  t152.method("GetXHalfLength1", static_cast<G4double (G4Trd::*)()  const>(&G4Trd::GetXHalfLength1));

  DEBUG_MSG("Adding wrapper for G4double G4Trd::GetXHalfLength2() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Trd::GetXHalfLength2()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:80:21
  t152.method("GetXHalfLength2", static_cast<G4double (G4Trd::*)()  const>(&G4Trd::GetXHalfLength2));

  DEBUG_MSG("Adding wrapper for G4double G4Trd::GetYHalfLength1() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Trd::GetYHalfLength1()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:81:21
  t152.method("GetYHalfLength1", static_cast<G4double (G4Trd::*)()  const>(&G4Trd::GetYHalfLength1));

  DEBUG_MSG("Adding wrapper for G4double G4Trd::GetYHalfLength2() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Trd::GetYHalfLength2()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:82:21
  t152.method("GetYHalfLength2", static_cast<G4double (G4Trd::*)()  const>(&G4Trd::GetYHalfLength2));

  DEBUG_MSG("Adding wrapper for G4double G4Trd::GetZHalfLength() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Trd::GetZHalfLength()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:83:21
  t152.method("GetZHalfLength", static_cast<G4double (G4Trd::*)()  const>(&G4Trd::GetZHalfLength));

  DEBUG_MSG("Adding wrapper for void G4Trd::SetXHalfLength1(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Trd::SetXHalfLength1(G4double)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:87:17
  t152.method("SetXHalfLength1", static_cast<void (G4Trd::*)(G4double) >(&G4Trd::SetXHalfLength1));

  DEBUG_MSG("Adding wrapper for void G4Trd::SetXHalfLength2(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Trd::SetXHalfLength2(G4double)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:88:17
  t152.method("SetXHalfLength2", static_cast<void (G4Trd::*)(G4double) >(&G4Trd::SetXHalfLength2));

  DEBUG_MSG("Adding wrapper for void G4Trd::SetYHalfLength1(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Trd::SetYHalfLength1(G4double)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:89:17
  t152.method("SetYHalfLength1", static_cast<void (G4Trd::*)(G4double) >(&G4Trd::SetYHalfLength1));

  DEBUG_MSG("Adding wrapper for void G4Trd::SetYHalfLength2(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Trd::SetYHalfLength2(G4double)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:90:17
  t152.method("SetYHalfLength2", static_cast<void (G4Trd::*)(G4double) >(&G4Trd::SetYHalfLength2));

  DEBUG_MSG("Adding wrapper for void G4Trd::SetZHalfLength(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Trd::SetZHalfLength(G4double)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:91:17
  t152.method("SetZHalfLength", static_cast<void (G4Trd::*)(G4double) >(&G4Trd::SetZHalfLength));

  DEBUG_MSG("Adding wrapper for void G4Trd::SetAllParameters(G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Trd::SetAllParameters(G4double, G4double, G4double, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:93:10
  t152.method("SetAllParameters", static_cast<void (G4Trd::*)(G4double, G4double, G4double, G4double, G4double) >(&G4Trd::SetAllParameters));

  DEBUG_MSG("Adding wrapper for G4double G4Trd::GetCubicVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Trd::GetCubicVolume()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:99:14
  t152.method("GetCubicVolume", static_cast<G4double (G4Trd::*)() >(&G4Trd::GetCubicVolume));

  DEBUG_MSG("Adding wrapper for G4double G4Trd::GetSurfaceArea() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Trd::GetSurfaceArea()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:100:14
  t152.method("GetSurfaceArea", static_cast<G4double (G4Trd::*)() >(&G4Trd::GetSurfaceArea));

  DEBUG_MSG("Adding wrapper for void G4Trd::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Trd::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:102:10
  t152.method("ComputeDimensions", static_cast<void (G4Trd::*)(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) >(&G4Trd::ComputeDimensions));

  DEBUG_MSG("Adding wrapper for void G4Trd::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4Trd::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:106:10
  t152.method("BoundingLimits", static_cast<void (G4Trd::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4Trd::BoundingLimits));

  DEBUG_MSG("Adding wrapper for EInside G4Trd::Inside(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: EInside G4Trd::Inside(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:113:13
  t152.method("Inside", static_cast<EInside (G4Trd::*)(const G4ThreeVector &)  const>(&G4Trd::Inside));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Trd::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Trd::SurfaceNormal(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:115:19
  t152.method("SurfaceNormal", static_cast<G4ThreeVector (G4Trd::*)(const G4ThreeVector &)  const>(&G4Trd::SurfaceNormal));

  DEBUG_MSG("Adding wrapper for G4double G4Trd::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Trd::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:117:14
  t152.method("DistanceToIn", static_cast<G4double (G4Trd::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Trd::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4Trd::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Trd::DistanceToIn(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:120:14
  t152.method("DistanceToIn", static_cast<G4double (G4Trd::*)(const G4ThreeVector &)  const>(&G4Trd::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4Trd::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Trd::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:122:14
  t152.method("DistanceToOut", static_cast<G4double (G4Trd::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4Trd::DistanceToOut));
  t152.method("DistanceToOut", [](G4Trd const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a.DistanceToOut(arg0, arg1); });
  t152.method("DistanceToOut", [](G4Trd const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a.DistanceToOut(arg0, arg1, arg2); });
  t152.method("DistanceToOut", [](G4Trd const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a.DistanceToOut(arg0, arg1, arg2, arg3); });
  t152.method("DistanceToOut", [](G4Trd const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a->DistanceToOut(arg0, arg1); });
  t152.method("DistanceToOut", [](G4Trd const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a->DistanceToOut(arg0, arg1, arg2); });
  t152.method("DistanceToOut", [](G4Trd const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a->DistanceToOut(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for G4double G4Trd::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Trd::DistanceToOut(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:128:14
  t152.method("DistanceToOut", static_cast<G4double (G4Trd::*)(const G4ThreeVector &)  const>(&G4Trd::DistanceToOut));

  DEBUG_MSG("Adding wrapper for G4GeometryType G4Trd::GetEntityType() (" __HERE__ ")");
  // signature to use in the veto list: G4GeometryType G4Trd::GetEntityType()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:130:20
  t152.method("GetEntityType", static_cast<G4GeometryType (G4Trd::*)()  const>(&G4Trd::GetEntityType));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Trd::GetPointOnSurface() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Trd::GetPointOnSurface()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:132:19
  t152.method("GetPointOnSurface", static_cast<G4ThreeVector (G4Trd::*)()  const>(&G4Trd::GetPointOnSurface));

  DEBUG_MSG("Adding wrapper for G4VSolid * G4Trd::Clone() (" __HERE__ ")");
  // signature to use in the veto list: G4VSolid * G4Trd::Clone()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:134:15
  t152.method("Clone", static_cast<G4VSolid * (G4Trd::*)()  const>(&G4Trd::Clone));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Trd::CreatePolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4Trd::CreatePolyhedron()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:141:19
  t152.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4Trd::*)()  const>(&G4Trd::CreatePolyhedron));


  DEBUG_MSG("Adding wrapper for void G4Trd::G4Trd(const G4Trd &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:150:5
  t152.constructor<const G4Trd &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4Trd & G4Trd::operator=(const G4Trd &) (" __HERE__ ")");
  // signature to use in the veto list: G4Trd & G4Trd::operator=(const G4Trd &)
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4Trd.hh:151:12
  t152.method("assign", static_cast<G4Trd & (G4Trd::*)(const G4Trd &) >(&G4Trd::operator=));

  /* End of G4Trd class method wrappers
   **********************************************************************/

}
