// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4VCSGfaceted> : std::false_type { };
  template<> struct DefaultConstructible<G4VCSGfaceted> : std::false_type { };
template<> struct SuperType<G4VCSGfaceted> { typedef G4VSolid type; };
}

// Class generating the wrapper for type G4VCSGfaceted
// signature to use in the veto file: G4VCSGfaceted
struct JlG4VCSGfaceted: public Wrapper {

  JlG4VCSGfaceted(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4VCSGfaceted (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:43:7
    jlcxx::TypeWrapper<G4VCSGfaceted>  t = jlModule.add_type<G4VCSGfaceted>("G4VCSGfaceted",
      jlcxx::julia_base_type<G4VSolid>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4VCSGfaceted>>(new jlcxx::TypeWrapper<G4VCSGfaceted>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;



    DEBUG_MSG("Adding wrapper for G4VCSGfaceted & G4VCSGfaceted::operator=(const G4VCSGfaceted &) (" __HERE__ ")");
    // signature to use in the veto list: G4VCSGfaceted & G4VCSGfaceted::operator=(const G4VCSGfaceted &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:51:20
    t.method("assign", static_cast<G4VCSGfaceted & (G4VCSGfaceted::*)(const G4VCSGfaceted &) >(&G4VCSGfaceted::operator=));

    DEBUG_MSG("Adding wrapper for EInside G4VCSGfaceted::Inside(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: EInside G4VCSGfaceted::Inside(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:58:13
    t.method("Inside", static_cast<EInside (G4VCSGfaceted::*)(const G4ThreeVector &)  const>(&G4VCSGfaceted::Inside));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4VCSGfaceted::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4VCSGfaceted::SurfaceNormal(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:60:19
    t.method("SurfaceNormal", static_cast<G4ThreeVector (G4VCSGfaceted::*)(const G4ThreeVector &)  const>(&G4VCSGfaceted::SurfaceNormal));

    DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VCSGfaceted::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:62:14
    t.method("DistanceToIn", static_cast<G4double (G4VCSGfaceted::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4VCSGfaceted::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VCSGfaceted::DistanceToIn(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:64:14
    t.method("DistanceToIn", static_cast<G4double (G4VCSGfaceted::*)(const G4ThreeVector &)  const>(&G4VCSGfaceted::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VCSGfaceted::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:65:14
    t.method("DistanceToOut", static_cast<G4double (G4VCSGfaceted::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4VCSGfaceted::DistanceToOut));
    t.method("DistanceToOut", [](G4VCSGfaceted const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a.DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4VCSGfaceted const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a.DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4VCSGfaceted const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a.DistanceToOut(arg0, arg1, arg2, arg3); });
    t.method("DistanceToOut", [](G4VCSGfaceted const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a->DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4VCSGfaceted const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a->DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4VCSGfaceted const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a->DistanceToOut(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VCSGfaceted::DistanceToOut(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:70:14
    t.method("DistanceToOut", static_cast<G4double (G4VCSGfaceted::*)(const G4ThreeVector &)  const>(&G4VCSGfaceted::DistanceToOut));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4VCSGfaceted::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4VCSGfaceted::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:72:20
    t.method("GetEntityType", static_cast<G4GeometryType (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::GetEntityType));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4VCSGfaceted::CreatePolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4VCSGfaceted::CreatePolyhedron()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:76:19
    t.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::CreatePolyhedron));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4VCSGfaceted::GetPolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4VCSGfaceted::GetPolyhedron()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:82:19
    t.method("GetPolyhedron", static_cast<G4Polyhedron * (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::GetPolyhedron));

    DEBUG_MSG("Adding wrapper for G4int G4VCSGfaceted::GetCubVolStatistics() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4VCSGfaceted::GetCubVolStatistics()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:84:11
    t.method("GetCubVolStatistics", static_cast<G4int (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::GetCubVolStatistics));

    DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::GetCubVolEpsilon() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VCSGfaceted::GetCubVolEpsilon()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:85:14
    t.method("GetCubVolEpsilon", static_cast<G4double (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::GetCubVolEpsilon));

    DEBUG_MSG("Adding wrapper for void G4VCSGfaceted::SetCubVolStatistics(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4VCSGfaceted::SetCubVolStatistics(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:86:10
    t.method("SetCubVolStatistics", static_cast<void (G4VCSGfaceted::*)(G4int) >(&G4VCSGfaceted::SetCubVolStatistics));

    DEBUG_MSG("Adding wrapper for void G4VCSGfaceted::SetCubVolEpsilon(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4VCSGfaceted::SetCubVolEpsilon(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:87:10
    t.method("SetCubVolEpsilon", static_cast<void (G4VCSGfaceted::*)(G4double) >(&G4VCSGfaceted::SetCubVolEpsilon));

    DEBUG_MSG("Adding wrapper for G4int G4VCSGfaceted::GetAreaStatistics() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4VCSGfaceted::GetAreaStatistics()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:88:11
    t.method("GetAreaStatistics", static_cast<G4int (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::GetAreaStatistics));

    DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::GetAreaAccuracy() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VCSGfaceted::GetAreaAccuracy()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:89:14
    t.method("GetAreaAccuracy", static_cast<G4double (G4VCSGfaceted::*)()  const>(&G4VCSGfaceted::GetAreaAccuracy));

    DEBUG_MSG("Adding wrapper for void G4VCSGfaceted::SetAreaStatistics(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4VCSGfaceted::SetAreaStatistics(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:90:10
    t.method("SetAreaStatistics", static_cast<void (G4VCSGfaceted::*)(G4int) >(&G4VCSGfaceted::SetAreaStatistics));

    DEBUG_MSG("Adding wrapper for void G4VCSGfaceted::SetAreaAccuracy(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4VCSGfaceted::SetAreaAccuracy(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:91:10
    t.method("SetAreaAccuracy", static_cast<void (G4VCSGfaceted::*)(G4double) >(&G4VCSGfaceted::SetAreaAccuracy));

    DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::GetCubicVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VCSGfaceted::GetCubicVolume()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:93:14
    t.method("GetCubicVolume", static_cast<G4double (G4VCSGfaceted::*)() >(&G4VCSGfaceted::GetCubicVolume));

    DEBUG_MSG("Adding wrapper for G4double G4VCSGfaceted::GetSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4VCSGfaceted::GetSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4VCSGfaceted.hh:96:14
    t.method("GetSurfaceArea", static_cast<G4double (G4VCSGfaceted::*)() >(&G4VCSGfaceted::GetSurfaceArea));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4VCSGfaceted>> type_;
};
std::shared_ptr<Wrapper> newJlG4VCSGfaceted(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4VCSGfaceted(module));
}
