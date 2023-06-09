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
void add_methods_for_G4Tet(jlcxx::Module& types, jlcxx::TypeWrapper<G4Tet>& t197) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4Tet
   */


  DEBUG_MSG("Adding wrapper for void G4Tet::G4Tet(const G4String &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, G4bool *) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:61:5
  t197.constructor<const G4String &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &>(/*finalize=*/true);
  t197.constructor<const G4String &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, G4bool *>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for void G4Tet::SetVertices(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, G4bool *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Tet::SetVertices(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, G4bool *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:72:10
  t197.method("SetVertices", static_cast<void (G4Tet::*)(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, G4bool *) >(&G4Tet::SetVertices));
  t197.method("SetVertices", [](G4Tet& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4ThreeVector & arg2, const G4ThreeVector & arg3)->void{ a.SetVertices(arg0, arg1, arg2, arg3); });
  t197.method("SetVertices", [](G4Tet* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4ThreeVector & arg2, const G4ThreeVector & arg3)->void{ a->SetVertices(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for void G4Tet::GetVertices(G4ThreeVector &, G4ThreeVector &, G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4Tet::GetVertices(G4ThreeVector &, G4ThreeVector &, G4ThreeVector &, G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:79:10
  t197.method("GetVertices", static_cast<void (G4Tet::*)(G4ThreeVector &, G4ThreeVector &, G4ThreeVector &, G4ThreeVector &)  const>(&G4Tet::GetVertices));

  DEBUG_MSG("Adding wrapper for std::vector<G4ThreeVector> G4Tet::GetVertices() (" __HERE__ ")");
  // signature to use in the veto list: std::vector<G4ThreeVector> G4Tet::GetVertices()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:83:32
  t197.method("GetVertices", static_cast<std::vector<G4ThreeVector> (G4Tet::*)()  const>(&G4Tet::GetVertices));

  DEBUG_MSG("Adding wrapper for void G4Tet::PrintWarnings(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4Tet::PrintWarnings(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:86:10
  t197.method("PrintWarnings", static_cast<void (G4Tet::*)(G4bool) >(&G4Tet::PrintWarnings));

  DEBUG_MSG("Adding wrapper for G4bool G4Tet::CheckDegeneracy(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4Tet::CheckDegeneracy(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:89:12
  t197.method("CheckDegeneracy", static_cast<G4bool (G4Tet::*)(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Tet::CheckDegeneracy));

  DEBUG_MSG("Adding wrapper for void G4Tet::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Tet::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:95:10
  t197.method("ComputeDimensions", static_cast<void (G4Tet::*)(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) >(&G4Tet::ComputeDimensions));

  DEBUG_MSG("Adding wrapper for void G4Tet::SetBoundingLimits(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4Tet::SetBoundingLimits(const G4ThreeVector &, const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:99:10
  t197.method("SetBoundingLimits", static_cast<void (G4Tet::*)(const G4ThreeVector &, const G4ThreeVector &) >(&G4Tet::SetBoundingLimits));

  DEBUG_MSG("Adding wrapper for void G4Tet::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4Tet::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:100:10
  t197.method("BoundingLimits", static_cast<void (G4Tet::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4Tet::BoundingLimits));

  DEBUG_MSG("Adding wrapper for EInside G4Tet::Inside(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: EInside G4Tet::Inside(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:106:13
  t197.method("Inside", static_cast<EInside (G4Tet::*)(const G4ThreeVector &)  const>(&G4Tet::Inside));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Tet::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Tet::SurfaceNormal(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:107:19
  t197.method("SurfaceNormal", static_cast<G4ThreeVector (G4Tet::*)(const G4ThreeVector &)  const>(&G4Tet::SurfaceNormal));

  DEBUG_MSG("Adding wrapper for G4double G4Tet::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Tet::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:108:14
  t197.method("DistanceToIn", static_cast<G4double (G4Tet::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Tet::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4Tet::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Tet::DistanceToIn(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:110:14
  t197.method("DistanceToIn", static_cast<G4double (G4Tet::*)(const G4ThreeVector &)  const>(&G4Tet::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4Tet::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Tet::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:111:14
  t197.method("DistanceToOut", static_cast<G4double (G4Tet::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4Tet::DistanceToOut));
  t197.method("DistanceToOut", [](G4Tet const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a.DistanceToOut(arg0, arg1); });
  t197.method("DistanceToOut", [](G4Tet const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a.DistanceToOut(arg0, arg1, arg2); });
  t197.method("DistanceToOut", [](G4Tet const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a.DistanceToOut(arg0, arg1, arg2, arg3); });
  t197.method("DistanceToOut", [](G4Tet const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a->DistanceToOut(arg0, arg1); });
  t197.method("DistanceToOut", [](G4Tet const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a->DistanceToOut(arg0, arg1, arg2); });
  t197.method("DistanceToOut", [](G4Tet const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a->DistanceToOut(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for G4double G4Tet::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Tet::DistanceToOut(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:116:14
  t197.method("DistanceToOut", static_cast<G4double (G4Tet::*)(const G4ThreeVector &)  const>(&G4Tet::DistanceToOut));

  DEBUG_MSG("Adding wrapper for G4GeometryType G4Tet::GetEntityType() (" __HERE__ ")");
  // signature to use in the veto list: G4GeometryType G4Tet::GetEntityType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:118:20
  t197.method("GetEntityType", static_cast<G4GeometryType (G4Tet::*)()  const>(&G4Tet::GetEntityType));

  DEBUG_MSG("Adding wrapper for G4VSolid * G4Tet::Clone() (" __HERE__ ")");
  // signature to use in the veto list: G4VSolid * G4Tet::Clone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:120:15
  t197.method("Clone", static_cast<G4VSolid * (G4Tet::*)()  const>(&G4Tet::Clone));

  DEBUG_MSG("Adding wrapper for G4double G4Tet::GetCubicVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Tet::GetCubicVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:124:14
  t197.method("GetCubicVolume", static_cast<G4double (G4Tet::*)() >(&G4Tet::GetCubicVolume));

  DEBUG_MSG("Adding wrapper for G4double G4Tet::GetSurfaceArea() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Tet::GetSurfaceArea()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:125:14
  t197.method("GetSurfaceArea", static_cast<G4double (G4Tet::*)() >(&G4Tet::GetSurfaceArea));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4Tet::GetPointOnSurface() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4Tet::GetPointOnSurface()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:127:19
  t197.method("GetPointOnSurface", static_cast<G4ThreeVector (G4Tet::*)()  const>(&G4Tet::GetPointOnSurface));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Tet::CreatePolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4Tet::CreatePolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:132:19
  t197.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4Tet::*)()  const>(&G4Tet::CreatePolyhedron));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Tet::GetPolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4Tet::GetPolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:133:19
  t197.method("GetPolyhedron", static_cast<G4Polyhedron * (G4Tet::*)()  const>(&G4Tet::GetPolyhedron));


  DEBUG_MSG("Adding wrapper for void G4Tet::G4Tet(const G4Tet &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:143:5
  t197.constructor<const G4Tet &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4Tet & G4Tet::operator=(const G4Tet &) (" __HERE__ ")");
  // signature to use in the veto list: G4Tet & G4Tet::operator=(const G4Tet &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Tet.hh:146:12
  t197.method("assign", static_cast<G4Tet & (G4Tet::*)(const G4Tet &) >(&G4Tet::operator=));

  /* End of G4Tet class method wrappers
   **********************************************************************/

}
