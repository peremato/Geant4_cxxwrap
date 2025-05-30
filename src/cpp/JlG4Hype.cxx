// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4Hype> : std::false_type { };
  template<> struct DefaultConstructible<G4Hype> : std::false_type { };
template<> struct SuperType<G4Hype> { typedef G4VSolid type; };
}

// Class generating the wrapper for type G4Hype
// signature to use in the veto file: G4Hype
struct JlG4Hype: public Wrapper {

  JlG4Hype(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4Hype (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:68:7
    jlcxx::TypeWrapper<G4Hype>  t = jlModule.add_type<G4Hype>("G4Hype",
      jlcxx::julia_base_type<G4VSolid>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4Hype>>(new jlcxx::TypeWrapper<G4Hype>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4Hype::G4Hype(const G4String &, G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:72:5
    t.constructor<const G4String &, G4double, G4double, G4double, G4double, G4double>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void G4Hype::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Hype::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:81:10
    t.method("ComputeDimensions", static_cast<void (G4Hype::*)(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) >(&G4Hype::ComputeDimensions));

    DEBUG_MSG("Adding wrapper for void G4Hype::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Hype::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:85:10
    t.method("BoundingLimits", static_cast<void (G4Hype::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4Hype::BoundingLimits));

    DEBUG_MSG("Adding wrapper for G4double G4Hype::GetInnerRadius() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Hype::GetInnerRadius()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:92:21
    t.method("GetInnerRadius", static_cast<G4double (G4Hype::*)()  const>(&G4Hype::GetInnerRadius));

    DEBUG_MSG("Adding wrapper for G4double G4Hype::GetOuterRadius() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Hype::GetOuterRadius()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:93:21
    t.method("GetOuterRadius", static_cast<G4double (G4Hype::*)()  const>(&G4Hype::GetOuterRadius));

    DEBUG_MSG("Adding wrapper for G4double G4Hype::GetZHalfLength() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Hype::GetZHalfLength()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:94:21
    t.method("GetZHalfLength", static_cast<G4double (G4Hype::*)()  const>(&G4Hype::GetZHalfLength));

    DEBUG_MSG("Adding wrapper for G4double G4Hype::GetInnerStereo() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Hype::GetInnerStereo()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:95:21
    t.method("GetInnerStereo", static_cast<G4double (G4Hype::*)()  const>(&G4Hype::GetInnerStereo));

    DEBUG_MSG("Adding wrapper for G4double G4Hype::GetOuterStereo() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Hype::GetOuterStereo()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:96:21
    t.method("GetOuterStereo", static_cast<G4double (G4Hype::*)()  const>(&G4Hype::GetOuterStereo));

    DEBUG_MSG("Adding wrapper for void G4Hype::SetInnerRadius(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Hype::SetInnerRadius(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:98:17
    t.method("SetInnerRadius", static_cast<void (G4Hype::*)(G4double) >(&G4Hype::SetInnerRadius));

    DEBUG_MSG("Adding wrapper for void G4Hype::SetOuterRadius(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Hype::SetOuterRadius(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:99:17
    t.method("SetOuterRadius", static_cast<void (G4Hype::*)(G4double) >(&G4Hype::SetOuterRadius));

    DEBUG_MSG("Adding wrapper for void G4Hype::SetZHalfLength(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Hype::SetZHalfLength(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:100:17
    t.method("SetZHalfLength", static_cast<void (G4Hype::*)(G4double) >(&G4Hype::SetZHalfLength));

    DEBUG_MSG("Adding wrapper for void G4Hype::SetInnerStereo(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Hype::SetInnerStereo(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:101:17
    t.method("SetInnerStereo", static_cast<void (G4Hype::*)(G4double) >(&G4Hype::SetInnerStereo));

    DEBUG_MSG("Adding wrapper for void G4Hype::SetOuterStereo(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Hype::SetOuterStereo(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:102:17
    t.method("SetOuterStereo", static_cast<void (G4Hype::*)(G4double) >(&G4Hype::SetOuterStereo));

    DEBUG_MSG("Adding wrapper for EInside G4Hype::Inside(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: EInside G4Hype::Inside(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:104:13
    t.method("Inside", static_cast<EInside (G4Hype::*)(const G4ThreeVector &)  const>(&G4Hype::Inside));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Hype::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Hype::SurfaceNormal(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:106:19
    t.method("SurfaceNormal", static_cast<G4ThreeVector (G4Hype::*)(const G4ThreeVector &)  const>(&G4Hype::SurfaceNormal));

    DEBUG_MSG("Adding wrapper for G4double G4Hype::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Hype::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:108:14
    t.method("DistanceToIn", static_cast<G4double (G4Hype::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Hype::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4Hype::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Hype::DistanceToIn(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:110:14
    t.method("DistanceToIn", static_cast<G4double (G4Hype::*)(const G4ThreeVector &)  const>(&G4Hype::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4Hype::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Hype::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:111:14
    t.method("DistanceToOut", static_cast<G4double (G4Hype::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4Hype::DistanceToOut));
    t.method("DistanceToOut", [](G4Hype const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a.DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4Hype const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a.DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4Hype const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a.DistanceToOut(arg0, arg1, arg2, arg3); });
    t.method("DistanceToOut", [](G4Hype const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a->DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4Hype const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a->DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4Hype const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a->DistanceToOut(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for G4double G4Hype::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Hype::DistanceToOut(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:115:14
    t.method("DistanceToOut", static_cast<G4double (G4Hype::*)(const G4ThreeVector &)  const>(&G4Hype::DistanceToOut));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4Hype::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4Hype::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:117:21
    t.method("GetEntityType", static_cast<G4GeometryType (G4Hype::*)()  const>(&G4Hype::GetEntityType));

    DEBUG_MSG("Adding wrapper for G4VSolid * G4Hype::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4VSolid * G4Hype::Clone()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:119:15
    t.method("Clone", static_cast<G4VSolid * (G4Hype::*)()  const>(&G4Hype::Clone));

    DEBUG_MSG("Adding wrapper for G4double G4Hype::GetCubicVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Hype::GetCubicVolume()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:123:14
    t.method("GetCubicVolume", static_cast<G4double (G4Hype::*)() >(&G4Hype::GetCubicVolume));

    DEBUG_MSG("Adding wrapper for G4double G4Hype::GetSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Hype::GetSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:124:14
    t.method("GetSurfaceArea", static_cast<G4double (G4Hype::*)() >(&G4Hype::GetSurfaceArea));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Hype::GetPointOnSurface() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Hype::GetPointOnSurface()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:126:19
    t.method("GetPointOnSurface", static_cast<G4ThreeVector (G4Hype::*)()  const>(&G4Hype::GetPointOnSurface));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Hype::CreatePolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4Hype::CreatePolyhedron()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:130:19
    t.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4Hype::*)()  const>(&G4Hype::CreatePolyhedron));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Hype::GetPolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4Hype::GetPolyhedron()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:131:19
    t.method("GetPolyhedron", static_cast<G4Polyhedron * (G4Hype::*)()  const>(&G4Hype::GetPolyhedron));


    DEBUG_MSG("Adding wrapper for void G4Hype::G4Hype(const G4Hype &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:138:5
    t.constructor<const G4Hype &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for G4Hype & G4Hype::operator=(const G4Hype &) (" __HERE__ ")");
    // signature to use in the veto list: G4Hype & G4Hype::operator=(const G4Hype &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Hype.hh:139:13
    t.method("assign", static_cast<G4Hype & (G4Hype::*)(const G4Hype &) >(&G4Hype::operator=));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4Hype>> type_;
};
std::shared_ptr<Wrapper> newJlG4Hype(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4Hype(module));
}
