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
void add_methods_for_G4CutTubs(jlcxx::Module& types, jlcxx::TypeWrapper<G4CutTubs>& t227) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4CutTubs
   */


  DEBUG_MSG("Adding wrapper for void G4CutTubs::G4CutTubs(const G4String &, G4double, G4double, G4double, G4double, G4double, G4ThreeVector, G4ThreeVector) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:63:5
  t227.constructor<const G4String &, G4double, G4double, G4double, G4double, G4double, G4ThreeVector, G4ThreeVector>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::GetInnerRadius() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::GetInnerRadius()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:80:21
  t227.method("GetInnerRadius", static_cast<G4double (G4CutTubs::*)()  const>(&G4CutTubs::GetInnerRadius));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::GetOuterRadius() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::GetOuterRadius()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:81:21
  t227.method("GetOuterRadius", static_cast<G4double (G4CutTubs::*)()  const>(&G4CutTubs::GetOuterRadius));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::GetZHalfLength() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::GetZHalfLength()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:82:21
  t227.method("GetZHalfLength", static_cast<G4double (G4CutTubs::*)()  const>(&G4CutTubs::GetZHalfLength));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::GetStartPhiAngle() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::GetStartPhiAngle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:83:21
  t227.method("GetStartPhiAngle", static_cast<G4double (G4CutTubs::*)()  const>(&G4CutTubs::GetStartPhiAngle));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::GetDeltaPhiAngle() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::GetDeltaPhiAngle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:84:21
  t227.method("GetDeltaPhiAngle", static_cast<G4double (G4CutTubs::*)()  const>(&G4CutTubs::GetDeltaPhiAngle));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::GetSinStartPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::GetSinStartPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:85:21
  t227.method("GetSinStartPhi", static_cast<G4double (G4CutTubs::*)()  const>(&G4CutTubs::GetSinStartPhi));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::GetCosStartPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::GetCosStartPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:86:21
  t227.method("GetCosStartPhi", static_cast<G4double (G4CutTubs::*)()  const>(&G4CutTubs::GetCosStartPhi));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::GetSinEndPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::GetSinEndPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:87:21
  t227.method("GetSinEndPhi", static_cast<G4double (G4CutTubs::*)()  const>(&G4CutTubs::GetSinEndPhi));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::GetCosEndPhi() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::GetCosEndPhi()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:88:21
  t227.method("GetCosEndPhi", static_cast<G4double (G4CutTubs::*)()  const>(&G4CutTubs::GetCosEndPhi));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4CutTubs::GetLowNorm() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4CutTubs::GetLowNorm()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:89:26
  t227.method("GetLowNorm", static_cast<G4ThreeVector (G4CutTubs::*)()  const>(&G4CutTubs::GetLowNorm));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4CutTubs::GetHighNorm() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4CutTubs::GetHighNorm()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:90:26
  t227.method("GetHighNorm", static_cast<G4ThreeVector (G4CutTubs::*)()  const>(&G4CutTubs::GetHighNorm));

  DEBUG_MSG("Adding wrapper for void G4CutTubs::SetInnerRadius(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4CutTubs::SetInnerRadius(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:94:17
  t227.method("SetInnerRadius", static_cast<void (G4CutTubs::*)(G4double) >(&G4CutTubs::SetInnerRadius));

  DEBUG_MSG("Adding wrapper for void G4CutTubs::SetOuterRadius(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4CutTubs::SetOuterRadius(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:95:17
  t227.method("SetOuterRadius", static_cast<void (G4CutTubs::*)(G4double) >(&G4CutTubs::SetOuterRadius));

  DEBUG_MSG("Adding wrapper for void G4CutTubs::SetZHalfLength(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4CutTubs::SetZHalfLength(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:96:17
  t227.method("SetZHalfLength", static_cast<void (G4CutTubs::*)(G4double) >(&G4CutTubs::SetZHalfLength));

  DEBUG_MSG("Adding wrapper for void G4CutTubs::SetStartPhiAngle(G4double, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4CutTubs::SetStartPhiAngle(G4double, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:97:17
  t227.method("SetStartPhiAngle", static_cast<void (G4CutTubs::*)(G4double, G4bool) >(&G4CutTubs::SetStartPhiAngle));
  t227.method("SetStartPhiAngle", [](G4CutTubs& a, G4double arg0)->void{ a.SetStartPhiAngle(arg0); });
  t227.method("SetStartPhiAngle", [](G4CutTubs* a, G4double arg0)->void{ a->SetStartPhiAngle(arg0); });

  DEBUG_MSG("Adding wrapper for void G4CutTubs::SetDeltaPhiAngle(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4CutTubs::SetDeltaPhiAngle(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:98:17
  t227.method("SetDeltaPhiAngle", static_cast<void (G4CutTubs::*)(G4double) >(&G4CutTubs::SetDeltaPhiAngle));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::GetCubicVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::GetCubicVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:102:14
  t227.method("GetCubicVolume", static_cast<G4double (G4CutTubs::*)() >(&G4CutTubs::GetCubicVolume));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::GetSurfaceArea() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::GetSurfaceArea()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:103:14
  t227.method("GetSurfaceArea", static_cast<G4double (G4CutTubs::*)() >(&G4CutTubs::GetSurfaceArea));

  DEBUG_MSG("Adding wrapper for void G4CutTubs::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4CutTubs::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:105:10
  t227.method("BoundingLimits", static_cast<void (G4CutTubs::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4CutTubs::BoundingLimits));

  DEBUG_MSG("Adding wrapper for EInside G4CutTubs::Inside(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: EInside G4CutTubs::Inside(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:112:13
  t227.method("Inside", static_cast<EInside (G4CutTubs::*)(const G4ThreeVector &)  const>(&G4CutTubs::Inside));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4CutTubs::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4CutTubs::SurfaceNormal(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:114:19
  t227.method("SurfaceNormal", static_cast<G4ThreeVector (G4CutTubs::*)(const G4ThreeVector &)  const>(&G4CutTubs::SurfaceNormal));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:116:14
  t227.method("DistanceToIn", static_cast<G4double (G4CutTubs::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4CutTubs::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::DistanceToIn(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:117:14
  t227.method("DistanceToIn", static_cast<G4double (G4CutTubs::*)(const G4ThreeVector &)  const>(&G4CutTubs::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:118:14
  t227.method("DistanceToOut", static_cast<G4double (G4CutTubs::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4CutTubs::DistanceToOut));
  t227.method("DistanceToOut", [](G4CutTubs const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a.DistanceToOut(arg0, arg1); });
  t227.method("DistanceToOut", [](G4CutTubs const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a.DistanceToOut(arg0, arg1, arg2); });
  t227.method("DistanceToOut", [](G4CutTubs const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a.DistanceToOut(arg0, arg1, arg2, arg3); });
  t227.method("DistanceToOut", [](G4CutTubs const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a->DistanceToOut(arg0, arg1); });
  t227.method("DistanceToOut", [](G4CutTubs const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a->DistanceToOut(arg0, arg1, arg2); });
  t227.method("DistanceToOut", [](G4CutTubs const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a->DistanceToOut(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for G4double G4CutTubs::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4CutTubs::DistanceToOut(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:122:14
  t227.method("DistanceToOut", static_cast<G4double (G4CutTubs::*)(const G4ThreeVector &)  const>(&G4CutTubs::DistanceToOut));

  DEBUG_MSG("Adding wrapper for G4GeometryType G4CutTubs::GetEntityType() (" __HERE__ ")");
  // signature to use in the veto list: G4GeometryType G4CutTubs::GetEntityType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:124:20
  t227.method("GetEntityType", static_cast<G4GeometryType (G4CutTubs::*)()  const>(&G4CutTubs::GetEntityType));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4CutTubs::GetPointOnSurface() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4CutTubs::GetPointOnSurface()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:126:19
  t227.method("GetPointOnSurface", static_cast<G4ThreeVector (G4CutTubs::*)()  const>(&G4CutTubs::GetPointOnSurface));

  DEBUG_MSG("Adding wrapper for G4VSolid * G4CutTubs::Clone() (" __HERE__ ")");
  // signature to use in the veto list: G4VSolid * G4CutTubs::Clone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:128:15
  t227.method("Clone", static_cast<G4VSolid * (G4CutTubs::*)()  const>(&G4CutTubs::Clone));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4CutTubs::CreatePolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4CutTubs::CreatePolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:135:25
  t227.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4CutTubs::*)()  const>(&G4CutTubs::CreatePolyhedron));


  DEBUG_MSG("Adding wrapper for void G4CutTubs::G4CutTubs(const G4CutTubs &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:145:5
  t227.constructor<const G4CutTubs &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4CutTubs & G4CutTubs::operator=(const G4CutTubs &) (" __HERE__ ")");
  // signature to use in the veto list: G4CutTubs & G4CutTubs::operator=(const G4CutTubs &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4CutTubs.hh:146:16
  t227.method("assign", static_cast<G4CutTubs & (G4CutTubs::*)(const G4CutTubs &) >(&G4CutTubs::operator=));

  /* End of G4CutTubs class method wrappers
   **********************************************************************/

}
