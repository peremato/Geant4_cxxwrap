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
void add_methods_for_G4VCSGfaceted(jlcxx::Module& types, jlcxx::TypeWrapper<G4VCSGfaceted>& t189) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4VCSGfaceted
   */



  DEBUG_MSG("Adding wrapper for G4VCSGfaceted & G4VCSGfaceted::operator=(const G4VCSGfaceted &) (" __HERE__ ")");
  // signature to use in the veto list: G4VCSGfaceted & G4VCSGfaceted::operator=(const G4VCSGfaceted &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:51:20
  t189.method("assign", static_cast<G4VCSGfaceted & (G4VCSGfaceted::*)(const G4VCSGfaceted &) >(&G4VCSGfaceted::operator=));

  DEBUG_MSG("Adding wrapper for EInside G4VCSGfaceted::Inside(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: EInside G4VCSGfaceted::Inside(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:58:21
  t189.method("Inside", static_cast<EInside (G4VCSGfaceted::*)(const G4ThreeVector &)  const>(&G4VCSGfaceted::Inside));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4VCSGfaceted::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4VCSGfaceted::SurfaceNormal(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:60:27
  t189.method("SurfaceNormal", static_cast<G4ThreeVector (G4VCSGfaceted::*)(const G4ThreeVector &)  const>(&G4VCSGfaceted::SurfaceNormal));

  DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VCSGfaceted::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:62:22
  t189.method("DistanceToIn", static_cast<G4double (G4VCSGfaceted::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4VCSGfaceted::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VCSGfaceted::DistanceToIn(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:64:22
  t189.method("DistanceToIn", static_cast<G4double (G4VCSGfaceted::*)(const G4ThreeVector &)  const>(&G4VCSGfaceted::DistanceToIn));

  DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VCSGfaceted::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:65:22
  t189.method("DistanceToOut", static_cast<G4double (G4VCSGfaceted::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4VCSGfaceted::DistanceToOut));
  t189.method("DistanceToOut", [](G4VCSGfaceted const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a.DistanceToOut(arg0, arg1); });
  t189.method("DistanceToOut", [](G4VCSGfaceted const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a.DistanceToOut(arg0, arg1, arg2); });
  t189.method("DistanceToOut", [](G4VCSGfaceted const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a.DistanceToOut(arg0, arg1, arg2, arg3); });
  t189.method("DistanceToOut", [](G4VCSGfaceted const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double{ return a->DistanceToOut(arg0, arg1); });
  t189.method("DistanceToOut", [](G4VCSGfaceted const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double{ return a->DistanceToOut(arg0, arg1, arg2); });
  t189.method("DistanceToOut", [](G4VCSGfaceted const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double{ return a->DistanceToOut(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VCSGfaceted::DistanceToOut(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:70:22
  t189.method("DistanceToOut", static_cast<G4double (G4VCSGfaceted::*)(const G4ThreeVector &)  const>(&G4VCSGfaceted::DistanceToOut));

  DEBUG_MSG("Adding wrapper for G4GeometryType G4VCSGfaceted::GetEntityType() (" __HERE__ ")");
  // signature to use in the veto list: G4GeometryType G4VCSGfaceted::GetEntityType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:72:28
  t189.method("GetEntityType", static_cast<G4GeometryType (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::GetEntityType));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4VCSGfaceted::CreatePolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4VCSGfaceted::CreatePolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:76:27
  t189.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::CreatePolyhedron));

  DEBUG_MSG("Adding wrapper for G4Polyhedron * G4VCSGfaceted::GetPolyhedron() (" __HERE__ ")");
  // signature to use in the veto list: G4Polyhedron * G4VCSGfaceted::GetPolyhedron()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:82:27
  t189.method("GetPolyhedron", static_cast<G4Polyhedron * (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::GetPolyhedron));

  DEBUG_MSG("Adding wrapper for G4int G4VCSGfaceted::GetCubVolStatistics() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VCSGfaceted::GetCubVolStatistics()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:84:11
  t189.method("GetCubVolStatistics", static_cast<G4int (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::GetCubVolStatistics));

  DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::GetCubVolEpsilon() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VCSGfaceted::GetCubVolEpsilon()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:85:14
  t189.method("GetCubVolEpsilon", static_cast<G4double (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::GetCubVolEpsilon));

  DEBUG_MSG("Adding wrapper for void G4VCSGfaceted::SetCubVolStatistics(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4VCSGfaceted::SetCubVolStatistics(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:86:10
  t189.method("SetCubVolStatistics", static_cast<void (G4VCSGfaceted::*)(G4int) >(&G4VCSGfaceted::SetCubVolStatistics));

  DEBUG_MSG("Adding wrapper for void G4VCSGfaceted::SetCubVolEpsilon(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4VCSGfaceted::SetCubVolEpsilon(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:87:10
  t189.method("SetCubVolEpsilon", static_cast<void (G4VCSGfaceted::*)(G4double) >(&G4VCSGfaceted::SetCubVolEpsilon));

  DEBUG_MSG("Adding wrapper for G4int G4VCSGfaceted::GetAreaStatistics() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VCSGfaceted::GetAreaStatistics()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:88:11
  t189.method("GetAreaStatistics", static_cast<G4int (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::GetAreaStatistics));

  DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::GetAreaAccuracy() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VCSGfaceted::GetAreaAccuracy()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:89:14
  t189.method("GetAreaAccuracy", static_cast<G4double (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::GetAreaAccuracy));

  DEBUG_MSG("Adding wrapper for void G4VCSGfaceted::SetAreaStatistics(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4VCSGfaceted::SetAreaStatistics(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:90:10
  t189.method("SetAreaStatistics", static_cast<void (G4VCSGfaceted::*)(G4int) >(&G4VCSGfaceted::SetAreaStatistics));

  DEBUG_MSG("Adding wrapper for void G4VCSGfaceted::SetAreaAccuracy(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4VCSGfaceted::SetAreaAccuracy(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:91:10
  t189.method("SetAreaAccuracy", static_cast<void (G4VCSGfaceted::*)(G4double) >(&G4VCSGfaceted::SetAreaAccuracy));

  DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::GetCubicVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VCSGfaceted::GetCubicVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:93:22
  t189.method("GetCubicVolume", static_cast<G4double (G4VCSGfaceted::*)() >(&G4VCSGfaceted::GetCubicVolume));

  DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::GetSurfaceArea() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VCSGfaceted::GetSurfaceArea()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VCSGfaceted.hh:96:22
  t189.method("GetSurfaceArea", static_cast<G4double (G4VCSGfaceted::*)() >(&G4VCSGfaceted::GetSurfaceArea));

  /* End of G4VCSGfaceted class method wrappers
   **********************************************************************/

}
