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
void add_methods_for_G4SubtractionSolid(jlcxx::Module& types, jlcxx::TypeWrapper<G4SubtractionSolid>& t210) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4SubtractionSolid
   */


  DEBUG_MSG("Adding wrapper for void G4SubtractionSolid::G4SubtractionSolid(const G4String &, G4VSolid *, G4VSolid *) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:50:5
  t210.constructor<const G4String &, G4VSolid *, G4VSolid *>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void G4SubtractionSolid::G4SubtractionSolid(const G4String &, G4VSolid *, G4VSolid *, G4RotationMatrix *, const G4ThreeVector &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:54:5
  t210.constructor<const G4String &, G4VSolid *, G4VSolid *, G4RotationMatrix *, const G4ThreeVector &>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void G4SubtractionSolid::G4SubtractionSolid(const G4String &, G4VSolid *, G4VSolid *, const G4Transform3D &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:60:5
  t210.constructor<const G4String &, G4VSolid *, G4VSolid *, const G4Transform3D &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4GeometryType G4SubtractionSolid::GetEntityType() (" __HERE__ ")");
  // signature to use in the veto list: G4GeometryType G4SubtractionSolid::GetEntityType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:67:21
  t210.method("GetEntityType", static_cast<G4GeometryType (G4SubtractionSolid::*)()  const>(&G4SubtractionSolid::GetEntityType));

  DEBUG_MSG("Adding wrapper for G4VSolid * G4SubtractionSolid::Clone() (" __HERE__ ")");
  // signature to use in the veto list: G4VSolid * G4SubtractionSolid::Clone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:69:15
  t210.method("Clone", static_cast<G4VSolid * (G4SubtractionSolid::*)()  const>(&G4SubtractionSolid::Clone));


  DEBUG_MSG("Adding wrapper for void G4SubtractionSolid::G4SubtractionSolid(const G4SubtractionSolid &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:76:5
  t210.constructor<const G4SubtractionSolid &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4SubtractionSolid & G4SubtractionSolid::operator=(const G4SubtractionSolid &) (" __HERE__ ")");
  // signature to use in the veto list: G4SubtractionSolid & G4SubtractionSolid::operator=(const G4SubtractionSolid &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:77:25
  t210.method("assign", static_cast<G4SubtractionSolid & (G4SubtractionSolid::*)(const G4SubtractionSolid &) >(&G4SubtractionSolid::operator=));

  DEBUG_MSG("Adding wrapper for void G4SubtractionSolid::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4SubtractionSolid::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:80:10
  t210.method("BoundingLimits", static_cast<void (G4SubtractionSolid::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4SubtractionSolid::BoundingLimits));

  DEBUG_MSG("Adding wrapper for EInside G4SubtractionSolid::Inside(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: EInside G4SubtractionSolid::Inside(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:87:13
  t210.method("Inside", static_cast<EInside (G4SubtractionSolid::*)(const G4ThreeVector &)  const>(&G4SubtractionSolid::Inside));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4SubtractionSolid::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4SubtractionSolid::SurfaceNormal(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:89:19
  t210.method("SurfaceNormal", static_cast<G4ThreeVector (G4SubtractionSolid::*)(const G4ThreeVector &)  const>(&G4SubtractionSolid::SurfaceNormal));

  DEBUG_MSG("Adding wrapper for G4double G4SubtractionSolid::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4SubtractionSolid::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:91:14
  t210.method("DistanceToIn", static_cast<G4double (G4SubtractionSolid::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4SubtractionSolid::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4SubtractionSolid::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4SubtractionSolid::DistanceToIn(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:94:14
  t210.method("DistanceToIn", static_cast<G4double (G4SubtractionSolid::*)(const G4ThreeVector &)  const>(&G4SubtractionSolid::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4SubtractionSolid::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4SubtractionSolid::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:96:14
  t210.method("DistanceToOut", static_cast<G4double (G4SubtractionSolid::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4SubtractionSolid::DistanceToOut));
  t210.method("DistanceToOut", [](G4SubtractionSolid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a.DistanceToOut(arg0, arg1); });
  t210.method("DistanceToOut", [](G4SubtractionSolid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a.DistanceToOut(arg0, arg1, arg2); });
  t210.method("DistanceToOut", [](G4SubtractionSolid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a.DistanceToOut(arg0, arg1, arg2, arg3); });
  t210.method("DistanceToOut", [](G4SubtractionSolid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a->DistanceToOut(arg0, arg1); });
  t210.method("DistanceToOut", [](G4SubtractionSolid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a->DistanceToOut(arg0, arg1, arg2); });
  t210.method("DistanceToOut", [](G4SubtractionSolid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a->DistanceToOut(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for G4double G4SubtractionSolid::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4SubtractionSolid::DistanceToOut(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:102:14
  t210.method("DistanceToOut", static_cast<G4double (G4SubtractionSolid::*)(const G4ThreeVector &)  const>(&G4SubtractionSolid::DistanceToOut));

  DEBUG_MSG("Adding wrapper for void G4SubtractionSolid::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4SubtractionSolid::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:105:10
  t210.method("ComputeDimensions", static_cast<void (G4SubtractionSolid::*)(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) >(&G4SubtractionSolid::ComputeDimensions));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4SubtractionSolid::CreatePolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4SubtractionSolid::CreatePolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:110:19
  t210.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4SubtractionSolid::*)()  const>(&G4SubtractionSolid::CreatePolyhedron));

  DEBUG_MSG("Adding wrapper for G4double G4SubtractionSolid::GetCubicVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4SubtractionSolid::GetCubicVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SubtractionSolid.hh:112:22
  t210.method("GetCubicVolume", static_cast<G4double (G4SubtractionSolid::*)() >(&G4SubtractionSolid::GetCubicVolume));

  /* End of G4SubtractionSolid class method wrappers
   **********************************************************************/

}
