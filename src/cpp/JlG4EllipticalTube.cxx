// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4EllipticalTube> : std::false_type { };
  template<> struct DefaultConstructible<G4EllipticalTube> : std::false_type { };
template<> struct SuperType<G4EllipticalTube> { typedef G4VSolid type; };
}

// Class generating the wrapper for type G4EllipticalTube
// signature to use in the veto file: G4EllipticalTube
struct JlG4EllipticalTube: public Wrapper {

  JlG4EllipticalTube(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4EllipticalTube (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:60:7
    jlcxx::TypeWrapper<G4EllipticalTube>  t = jlModule.add_type<G4EllipticalTube>("G4EllipticalTube",
      jlcxx::julia_base_type<G4VSolid>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4EllipticalTube>>(new jlcxx::TypeWrapper<G4EllipticalTube>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4EllipticalTube::G4EllipticalTube(const G4String &, G4double, G4double, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:64:5
    t.constructor<const G4String &, G4double, G4double, G4double>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void G4EllipticalTube::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4EllipticalTube::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:73:10
    t.method("BoundingLimits", static_cast<void (G4EllipticalTube::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4EllipticalTube::BoundingLimits));

    DEBUG_MSG("Adding wrapper for EInside G4EllipticalTube::Inside(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: EInside G4EllipticalTube::Inside(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:80:13
    t.method("Inside", static_cast<EInside (G4EllipticalTube::*)(const G4ThreeVector &)  const>(&G4EllipticalTube::Inside));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4EllipticalTube::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4EllipticalTube::SurfaceNormal(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:82:19
    t.method("SurfaceNormal", static_cast<G4ThreeVector (G4EllipticalTube::*)(const G4ThreeVector &)  const>(&G4EllipticalTube::SurfaceNormal));

    DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4EllipticalTube::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:84:14
    t.method("DistanceToIn", static_cast<G4double (G4EllipticalTube::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4EllipticalTube::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4EllipticalTube::DistanceToIn(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:87:14
    t.method("DistanceToIn", static_cast<G4double (G4EllipticalTube::*)(const G4ThreeVector &)  const>(&G4EllipticalTube::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4EllipticalTube::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:89:14
    t.method("DistanceToOut", static_cast<G4double (G4EllipticalTube::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4EllipticalTube::DistanceToOut));
    t.method("DistanceToOut", [](G4EllipticalTube const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a.DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4EllipticalTube const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a.DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4EllipticalTube const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a.DistanceToOut(arg0, arg1, arg2, arg3); });
    t.method("DistanceToOut", [](G4EllipticalTube const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a->DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4EllipticalTube const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a->DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4EllipticalTube const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a->DistanceToOut(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4EllipticalTube::DistanceToOut(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:95:14
    t.method("DistanceToOut", static_cast<G4double (G4EllipticalTube::*)(const G4ThreeVector &)  const>(&G4EllipticalTube::DistanceToOut));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4EllipticalTube::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4EllipticalTube::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:97:20
    t.method("GetEntityType", static_cast<G4GeometryType (G4EllipticalTube::*)()  const>(&G4EllipticalTube::GetEntityType));

    DEBUG_MSG("Adding wrapper for G4VSolid * G4EllipticalTube::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4VSolid * G4EllipticalTube::Clone()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:99:15
    t.method("Clone", static_cast<G4VSolid * (G4EllipticalTube::*)()  const>(&G4EllipticalTube::Clone));

    DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::GetCubicVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4EllipticalTube::GetCubicVolume()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:103:14
    t.method("GetCubicVolume", static_cast<G4double (G4EllipticalTube::*)() >(&G4EllipticalTube::GetCubicVolume));

    DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::GetSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4EllipticalTube::GetSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:104:14
    t.method("GetSurfaceArea", static_cast<G4double (G4EllipticalTube::*)() >(&G4EllipticalTube::GetSurfaceArea));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4EllipticalTube::GetPointOnSurface() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4EllipticalTube::GetPointOnSurface()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:106:19
    t.method("GetPointOnSurface", static_cast<G4ThreeVector (G4EllipticalTube::*)()  const>(&G4EllipticalTube::GetPointOnSurface));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4EllipticalTube::CreatePolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4EllipticalTube::CreatePolyhedron()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:110:19
    t.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4EllipticalTube::*)()  const>(&G4EllipticalTube::CreatePolyhedron));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4EllipticalTube::GetPolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4EllipticalTube::GetPolyhedron()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:111:19
    t.method("GetPolyhedron", static_cast<G4Polyhedron * (G4EllipticalTube::*)()  const>(&G4EllipticalTube::GetPolyhedron));

    DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::GetDx() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4EllipticalTube::GetDx()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:117:21
    t.method("GetDx", static_cast<G4double (G4EllipticalTube::*)()  const>(&G4EllipticalTube::GetDx));

    DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::GetDy() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4EllipticalTube::GetDy()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:118:21
    t.method("GetDy", static_cast<G4double (G4EllipticalTube::*)()  const>(&G4EllipticalTube::GetDy));

    DEBUG_MSG("Adding wrapper for G4double G4EllipticalTube::GetDz() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4EllipticalTube::GetDz()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:119:21
    t.method("GetDz", static_cast<G4double (G4EllipticalTube::*)()  const>(&G4EllipticalTube::GetDz));

    DEBUG_MSG("Adding wrapper for void G4EllipticalTube::SetDx(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4EllipticalTube::SetDx(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:121:17
    t.method("SetDx", static_cast<void (G4EllipticalTube::*)(G4double) >(&G4EllipticalTube::SetDx));

    DEBUG_MSG("Adding wrapper for void G4EllipticalTube::SetDy(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4EllipticalTube::SetDy(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:122:17
    t.method("SetDy", static_cast<void (G4EllipticalTube::*)(G4double) >(&G4EllipticalTube::SetDy));

    DEBUG_MSG("Adding wrapper for void G4EllipticalTube::SetDz(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4EllipticalTube::SetDz(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:123:17
    t.method("SetDz", static_cast<void (G4EllipticalTube::*)(G4double) >(&G4EllipticalTube::SetDz));


    DEBUG_MSG("Adding wrapper for void G4EllipticalTube::G4EllipticalTube(const G4EllipticalTube &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:130:5
    t.constructor<const G4EllipticalTube &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for G4EllipticalTube & G4EllipticalTube::operator=(const G4EllipticalTube &) (" __HERE__ ")");
    // signature to use in the veto list: G4EllipticalTube & G4EllipticalTube::operator=(const G4EllipticalTube &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4EllipticalTube.hh:131:23
    t.method("assign", static_cast<G4EllipticalTube & (G4EllipticalTube::*)(const G4EllipticalTube &) >(&G4EllipticalTube::operator=));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4EllipticalTube>> type_;
};
std::shared_ptr<Wrapper> newJlG4EllipticalTube(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4EllipticalTube(module));
}
