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
void add_methods_for_G4EllipticalTube(jlcxx::Module& types, jlcxx::TypeWrapper<G4EllipticalTube>& t221) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4EllipticalTube
   */


  DEBUG_MSG("Adding wrapper for void G4EllipticalTube::G4EllipticalTube(const G4String &, G4double, G4double, G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:64:5
  t221.constructor<const G4String &, G4double, G4double, G4double>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for void G4EllipticalTube::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4EllipticalTube::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:73:10
  t221.method("BoundingLimits", static_cast<void (G4EllipticalTube::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4EllipticalTube::BoundingLimits));

  DEBUG_MSG("Adding wrapper for EInside G4EllipticalTube::Inside(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: EInside G4EllipticalTube::Inside(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:80:13
  t221.method("Inside", static_cast<EInside (G4EllipticalTube::*)(const G4ThreeVector &)  const>(&G4EllipticalTube::Inside));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4EllipticalTube::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4EllipticalTube::SurfaceNormal(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:82:19
  t221.method("SurfaceNormal", static_cast<G4ThreeVector (G4EllipticalTube::*)(const G4ThreeVector &)  const>(&G4EllipticalTube::SurfaceNormal));

  DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4EllipticalTube::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:84:14
  t221.method("DistanceToIn", static_cast<G4double (G4EllipticalTube::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4EllipticalTube::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4EllipticalTube::DistanceToIn(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:87:14
  t221.method("DistanceToIn", static_cast<G4double (G4EllipticalTube::*)(const G4ThreeVector &)  const>(&G4EllipticalTube::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4EllipticalTube::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:89:14
  t221.method("DistanceToOut", static_cast<G4double (G4EllipticalTube::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4EllipticalTube::DistanceToOut));
  t221.method("DistanceToOut", [](G4EllipticalTube const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a.DistanceToOut(arg0, arg1); });
  t221.method("DistanceToOut", [](G4EllipticalTube const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a.DistanceToOut(arg0, arg1, arg2); });
  t221.method("DistanceToOut", [](G4EllipticalTube const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a.DistanceToOut(arg0, arg1, arg2, arg3); });
  t221.method("DistanceToOut", [](G4EllipticalTube const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a->DistanceToOut(arg0, arg1); });
  t221.method("DistanceToOut", [](G4EllipticalTube const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a->DistanceToOut(arg0, arg1, arg2); });
  t221.method("DistanceToOut", [](G4EllipticalTube const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a->DistanceToOut(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4EllipticalTube::DistanceToOut(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:95:14
  t221.method("DistanceToOut", static_cast<G4double (G4EllipticalTube::*)(const G4ThreeVector &)  const>(&G4EllipticalTube::DistanceToOut));

  DEBUG_MSG("Adding wrapper for G4GeometryType G4EllipticalTube::GetEntityType() (" __HERE__ ")");
  // signature to use in the veto list: G4GeometryType G4EllipticalTube::GetEntityType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:97:20
  t221.method("GetEntityType", static_cast<G4GeometryType (G4EllipticalTube::*)()  const>(&G4EllipticalTube::GetEntityType));

  DEBUG_MSG("Adding wrapper for G4VSolid * G4EllipticalTube::Clone() (" __HERE__ ")");
  // signature to use in the veto list: G4VSolid * G4EllipticalTube::Clone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:99:15
  t221.method("Clone", static_cast<G4VSolid * (G4EllipticalTube::*)()  const>(&G4EllipticalTube::Clone));

  DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::GetCubicVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4EllipticalTube::GetCubicVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:103:14
  t221.method("GetCubicVolume", static_cast<G4double (G4EllipticalTube::*)() >(&G4EllipticalTube::GetCubicVolume));

  DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::GetSurfaceArea() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4EllipticalTube::GetSurfaceArea()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:104:14
  t221.method("GetSurfaceArea", static_cast<G4double (G4EllipticalTube::*)() >(&G4EllipticalTube::GetSurfaceArea));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4EllipticalTube::GetPointOnSurface() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4EllipticalTube::GetPointOnSurface()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:106:19
  t221.method("GetPointOnSurface", static_cast<G4ThreeVector (G4EllipticalTube::*)()  const>(&G4EllipticalTube::GetPointOnSurface));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4EllipticalTube::CreatePolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4EllipticalTube::CreatePolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:110:19
  t221.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4EllipticalTube::*)()  const>(&G4EllipticalTube::CreatePolyhedron));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4EllipticalTube::GetPolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4EllipticalTube::GetPolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:111:19
  t221.method("GetPolyhedron", static_cast<G4Polyhedron * (G4EllipticalTube::*)()  const>(&G4EllipticalTube::GetPolyhedron));

  DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::GetDx() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4EllipticalTube::GetDx()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:117:21
  t221.method("GetDx", static_cast<G4double (G4EllipticalTube::*)()  const>(&G4EllipticalTube::GetDx));

  DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::GetDy() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4EllipticalTube::GetDy()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:118:21
  t221.method("GetDy", static_cast<G4double (G4EllipticalTube::*)()  const>(&G4EllipticalTube::GetDy));

  DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::GetDz() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4EllipticalTube::GetDz()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:119:21
  t221.method("GetDz", static_cast<G4double (G4EllipticalTube::*)()  const>(&G4EllipticalTube::GetDz));

  DEBUG_MSG("Adding wrapper for void G4EllipticalTube::SetDx(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4EllipticalTube::SetDx(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:121:17
  t221.method("SetDx", static_cast<void (G4EllipticalTube::*)(G4double) >(&G4EllipticalTube::SetDx));

  DEBUG_MSG("Adding wrapper for void G4EllipticalTube::SetDy(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4EllipticalTube::SetDy(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:122:17
  t221.method("SetDy", static_cast<void (G4EllipticalTube::*)(G4double) >(&G4EllipticalTube::SetDy));

  DEBUG_MSG("Adding wrapper for void G4EllipticalTube::SetDz(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4EllipticalTube::SetDz(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:123:17
  t221.method("SetDz", static_cast<void (G4EllipticalTube::*)(G4double) >(&G4EllipticalTube::SetDz));


  DEBUG_MSG("Adding wrapper for void G4EllipticalTube::G4EllipticalTube(const G4EllipticalTube &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:132:5
  t221.constructor<const G4EllipticalTube &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4EllipticalTube & G4EllipticalTube::operator=(const G4EllipticalTube &) (" __HERE__ ")");
  // signature to use in the veto list: G4EllipticalTube & G4EllipticalTube::operator=(const G4EllipticalTube &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EllipticalTube.hh:133:23
  t221.method("assign", static_cast<G4EllipticalTube & (G4EllipticalTube::*)(const G4EllipticalTube &) >(&G4EllipticalTube::operator=));

  /* End of G4EllipticalTube class method wrappers
   **********************************************************************/

}
