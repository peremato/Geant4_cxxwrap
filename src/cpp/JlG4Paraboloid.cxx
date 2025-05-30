// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4Paraboloid> : std::false_type { };
  template<> struct DefaultConstructible<G4Paraboloid> : std::false_type { };
template<> struct SuperType<G4Paraboloid> { typedef G4VSolid type; };
}

// Class generating the wrapper for type G4Paraboloid
// signature to use in the veto file: G4Paraboloid
struct JlG4Paraboloid: public Wrapper {

  JlG4Paraboloid(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4Paraboloid (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:67:7
    jlcxx::TypeWrapper<G4Paraboloid>  t = jlModule.add_type<G4Paraboloid>("G4Paraboloid",
      jlcxx::julia_base_type<G4VSolid>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4Paraboloid>>(new jlcxx::TypeWrapper<G4Paraboloid>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4Paraboloid::G4Paraboloid(const G4String &, G4double, G4double, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:71:5
    t.constructor<const G4String &, G4double, G4double, G4double>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for G4double G4Paraboloid::GetZHalfLength() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Paraboloid::GetZHalfLength()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:80:21
    t.method("GetZHalfLength", static_cast<G4double (G4Paraboloid::*)()  const>(&G4Paraboloid::GetZHalfLength));

    DEBUG_MSG("Adding wrapper for G4double G4Paraboloid::GetRadiusMinusZ() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Paraboloid::GetRadiusMinusZ()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:81:21
    t.method("GetRadiusMinusZ", static_cast<G4double (G4Paraboloid::*)()  const>(&G4Paraboloid::GetRadiusMinusZ));

    DEBUG_MSG("Adding wrapper for G4double G4Paraboloid::GetRadiusPlusZ() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Paraboloid::GetRadiusPlusZ()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:82:21
    t.method("GetRadiusPlusZ", static_cast<G4double (G4Paraboloid::*)()  const>(&G4Paraboloid::GetRadiusPlusZ));

    DEBUG_MSG("Adding wrapper for G4double G4Paraboloid::GetCubicVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Paraboloid::GetCubicVolume()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:84:21
    t.method("GetCubicVolume", static_cast<G4double (G4Paraboloid::*)() >(&G4Paraboloid::GetCubicVolume));

    DEBUG_MSG("Adding wrapper for G4double G4Paraboloid::GetSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Paraboloid::GetSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:85:21
    t.method("GetSurfaceArea", static_cast<G4double (G4Paraboloid::*)() >(&G4Paraboloid::GetSurfaceArea));

    DEBUG_MSG("Adding wrapper for G4double G4Paraboloid::CalculateSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Paraboloid::CalculateSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:86:21
    t.method("CalculateSurfaceArea", static_cast<G4double (G4Paraboloid::*)()  const>(&G4Paraboloid::CalculateSurfaceArea));

    DEBUG_MSG("Adding wrapper for void G4Paraboloid::SetZHalfLength(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Paraboloid::SetZHalfLength(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:90:17
    t.method("SetZHalfLength", static_cast<void (G4Paraboloid::*)(G4double) >(&G4Paraboloid::SetZHalfLength));

    DEBUG_MSG("Adding wrapper for void G4Paraboloid::SetRadiusMinusZ(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Paraboloid::SetRadiusMinusZ(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:91:17
    t.method("SetRadiusMinusZ", static_cast<void (G4Paraboloid::*)(G4double) >(&G4Paraboloid::SetRadiusMinusZ));

    DEBUG_MSG("Adding wrapper for void G4Paraboloid::SetRadiusPlusZ(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Paraboloid::SetRadiusPlusZ(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:92:17
    t.method("SetRadiusPlusZ", static_cast<void (G4Paraboloid::*)(G4double) >(&G4Paraboloid::SetRadiusPlusZ));

    DEBUG_MSG("Adding wrapper for void G4Paraboloid::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Paraboloid::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:96:10
    t.method("BoundingLimits", static_cast<void (G4Paraboloid::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4Paraboloid::BoundingLimits));

    DEBUG_MSG("Adding wrapper for EInside G4Paraboloid::Inside(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: EInside G4Paraboloid::Inside(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:101:13
    t.method("Inside", static_cast<EInside (G4Paraboloid::*)(const G4ThreeVector &)  const>(&G4Paraboloid::Inside));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Paraboloid::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Paraboloid::SurfaceNormal(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:102:19
    t.method("SurfaceNormal", static_cast<G4ThreeVector (G4Paraboloid::*)(const G4ThreeVector &)  const>(&G4Paraboloid::SurfaceNormal));

    DEBUG_MSG("Adding wrapper for G4double G4Paraboloid::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Paraboloid::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:103:14
    t.method("DistanceToIn", static_cast<G4double (G4Paraboloid::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Paraboloid::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4Paraboloid::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Paraboloid::DistanceToIn(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:105:14
    t.method("DistanceToIn", static_cast<G4double (G4Paraboloid::*)(const G4ThreeVector &)  const>(&G4Paraboloid::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4Paraboloid::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Paraboloid::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:106:14
    t.method("DistanceToOut", static_cast<G4double (G4Paraboloid::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4Paraboloid::DistanceToOut));
    t.method("DistanceToOut", [](G4Paraboloid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a.DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4Paraboloid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a.DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4Paraboloid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a.DistanceToOut(arg0, arg1, arg2, arg3); });
    t.method("DistanceToOut", [](G4Paraboloid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a->DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4Paraboloid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a->DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4Paraboloid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a->DistanceToOut(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for G4double G4Paraboloid::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Paraboloid::DistanceToOut(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:111:14
    t.method("DistanceToOut", static_cast<G4double (G4Paraboloid::*)(const G4ThreeVector &)  const>(&G4Paraboloid::DistanceToOut));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4Paraboloid::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4Paraboloid::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:113:20
    t.method("GetEntityType", static_cast<G4GeometryType (G4Paraboloid::*)()  const>(&G4Paraboloid::GetEntityType));

    DEBUG_MSG("Adding wrapper for G4VSolid * G4Paraboloid::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4VSolid * G4Paraboloid::Clone()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:115:15
    t.method("Clone", static_cast<G4VSolid * (G4Paraboloid::*)()  const>(&G4Paraboloid::Clone));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Paraboloid::GetPointOnSurface() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Paraboloid::GetPointOnSurface()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:119:19
    t.method("GetPointOnSurface", static_cast<G4ThreeVector (G4Paraboloid::*)()  const>(&G4Paraboloid::GetPointOnSurface));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Paraboloid::CreatePolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4Paraboloid::CreatePolyhedron()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:124:19
    t.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4Paraboloid::*)()  const>(&G4Paraboloid::CreatePolyhedron));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Paraboloid::GetPolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4Paraboloid::GetPolyhedron()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:125:19
    t.method("GetPolyhedron", static_cast<G4Polyhedron * (G4Paraboloid::*)()  const>(&G4Paraboloid::GetPolyhedron));


    DEBUG_MSG("Adding wrapper for void G4Paraboloid::G4Paraboloid(const G4Paraboloid &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:132:5
    t.constructor<const G4Paraboloid &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for G4Paraboloid & G4Paraboloid::operator=(const G4Paraboloid &) (" __HERE__ ")");
    // signature to use in the veto list: G4Paraboloid & G4Paraboloid::operator=(const G4Paraboloid &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Paraboloid.hh:133:19
    t.method("assign", static_cast<G4Paraboloid & (G4Paraboloid::*)(const G4Paraboloid &) >(&G4Paraboloid::operator=));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4Paraboloid>> type_;
};
std::shared_ptr<Wrapper> newJlG4Paraboloid(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4Paraboloid(module));
}
