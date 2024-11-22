// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4Tet> : std::false_type { };
  template<> struct DefaultConstructible<G4Tet> : std::false_type { };
template<> struct SuperType<G4Tet> { typedef G4VSolid type; };
}

// Class generating the wrapper for type G4Tet
// signature to use in the veto file: G4Tet
struct JlG4Tet: public Wrapper {

  JlG4Tet(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4Tet (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:55:7
    jlcxx::TypeWrapper<G4Tet>  t = jlModule.add_type<G4Tet>("G4Tet",
      jlcxx::julia_base_type<G4VSolid>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4Tet>>(new jlcxx::TypeWrapper<G4Tet>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4Tet::G4Tet(const G4String &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, G4bool *) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:59:5
    t.constructor<const G4String &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const G4String &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, G4bool *>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void G4Tet::SetVertices(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, G4bool *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Tet::SetVertices(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, G4bool *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:68:10
    t.method("SetVertices", static_cast<void (G4Tet::*)(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, G4bool *) >(&G4Tet::SetVertices));
    t.method("SetVertices", [](G4Tet& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4ThreeVector & arg2, const G4ThreeVector & arg3)->void { a.SetVertices(arg0, arg1, arg2, arg3); });
    t.method("SetVertices", [](G4Tet* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4ThreeVector & arg2, const G4ThreeVector & arg3)->void { a->SetVertices(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for void G4Tet::GetVertices(G4ThreeVector &, G4ThreeVector &, G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Tet::GetVertices(G4ThreeVector &, G4ThreeVector &, G4ThreeVector &, G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:75:10
    t.method("GetVertices", static_cast<void (G4Tet::*)(G4ThreeVector &, G4ThreeVector &, G4ThreeVector &, G4ThreeVector &)  const>(&G4Tet::GetVertices));

    DEBUG_MSG("Adding wrapper for std::vector<G4ThreeVector> G4Tet::GetVertices() (" __HERE__ ")");
    // signature to use in the veto list: std::vector<G4ThreeVector> G4Tet::GetVertices()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:79:32
    t.method("GetVertices", static_cast<std::vector<G4ThreeVector> (G4Tet::*)()  const>(&G4Tet::GetVertices));

    DEBUG_MSG("Adding wrapper for void G4Tet::PrintWarnings(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4Tet::PrintWarnings(G4bool)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:82:17
    t.method("PrintWarnings", static_cast<void (G4Tet::*)(G4bool) >(&G4Tet::PrintWarnings));

    DEBUG_MSG("Adding wrapper for G4bool G4Tet::CheckDegeneracy(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Tet::CheckDegeneracy(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:85:12
    t.method("CheckDegeneracy", static_cast<G4bool (G4Tet::*)(const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Tet::CheckDegeneracy));

    DEBUG_MSG("Adding wrapper for void G4Tet::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Tet::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:91:10
    t.method("ComputeDimensions", static_cast<void (G4Tet::*)(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) >(&G4Tet::ComputeDimensions));

    DEBUG_MSG("Adding wrapper for void G4Tet::SetBoundingLimits(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Tet::SetBoundingLimits(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:95:10
    t.method("SetBoundingLimits", static_cast<void (G4Tet::*)(const G4ThreeVector &, const G4ThreeVector &) >(&G4Tet::SetBoundingLimits));

    DEBUG_MSG("Adding wrapper for void G4Tet::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Tet::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:96:10
    t.method("BoundingLimits", static_cast<void (G4Tet::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4Tet::BoundingLimits));

    DEBUG_MSG("Adding wrapper for EInside G4Tet::Inside(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: EInside G4Tet::Inside(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:102:13
    t.method("Inside", static_cast<EInside (G4Tet::*)(const G4ThreeVector &)  const>(&G4Tet::Inside));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Tet::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Tet::SurfaceNormal(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:103:19
    t.method("SurfaceNormal", static_cast<G4ThreeVector (G4Tet::*)(const G4ThreeVector &)  const>(&G4Tet::SurfaceNormal));

    DEBUG_MSG("Adding wrapper for G4double G4Tet::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Tet::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:104:14
    t.method("DistanceToIn", static_cast<G4double (G4Tet::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Tet::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4Tet::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Tet::DistanceToIn(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:106:14
    t.method("DistanceToIn", static_cast<G4double (G4Tet::*)(const G4ThreeVector &)  const>(&G4Tet::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4Tet::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Tet::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:107:14
    t.method("DistanceToOut", static_cast<G4double (G4Tet::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4Tet::DistanceToOut));
    t.method("DistanceToOut", [](G4Tet const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a.DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4Tet const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a.DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4Tet const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a.DistanceToOut(arg0, arg1, arg2, arg3); });
    t.method("DistanceToOut", [](G4Tet const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a->DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4Tet const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a->DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4Tet const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a->DistanceToOut(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for G4double G4Tet::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Tet::DistanceToOut(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:112:14
    t.method("DistanceToOut", static_cast<G4double (G4Tet::*)(const G4ThreeVector &)  const>(&G4Tet::DistanceToOut));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4Tet::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4Tet::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:114:20
    t.method("GetEntityType", static_cast<G4GeometryType (G4Tet::*)()  const>(&G4Tet::GetEntityType));

    DEBUG_MSG("Adding wrapper for G4VSolid * G4Tet::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4VSolid * G4Tet::Clone()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:116:15
    t.method("Clone", static_cast<G4VSolid * (G4Tet::*)()  const>(&G4Tet::Clone));

    DEBUG_MSG("Adding wrapper for G4double G4Tet::GetCubicVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Tet::GetCubicVolume()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:120:14
    t.method("GetCubicVolume", static_cast<G4double (G4Tet::*)() >(&G4Tet::GetCubicVolume));

    DEBUG_MSG("Adding wrapper for G4double G4Tet::GetSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Tet::GetSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:121:14
    t.method("GetSurfaceArea", static_cast<G4double (G4Tet::*)() >(&G4Tet::GetSurfaceArea));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Tet::GetPointOnSurface() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Tet::GetPointOnSurface()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:123:19
    t.method("GetPointOnSurface", static_cast<G4ThreeVector (G4Tet::*)()  const>(&G4Tet::GetPointOnSurface));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Tet::CreatePolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4Tet::CreatePolyhedron()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:128:19
    t.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4Tet::*)()  const>(&G4Tet::CreatePolyhedron));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Tet::GetPolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4Tet::GetPolyhedron()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:129:19
    t.method("GetPolyhedron", static_cast<G4Polyhedron * (G4Tet::*)()  const>(&G4Tet::GetPolyhedron));


    DEBUG_MSG("Adding wrapper for void G4Tet::G4Tet(const G4Tet &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:137:5
    t.constructor<const G4Tet &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for G4Tet & G4Tet::operator=(const G4Tet &) (" __HERE__ ")");
    // signature to use in the veto list: G4Tet & G4Tet::operator=(const G4Tet &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Tet.hh:140:12
    t.method("assign", static_cast<G4Tet & (G4Tet::*)(const G4Tet &) >(&G4Tet::operator=));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4Tet>> type_;
};
std::shared_ptr<Wrapper> newJlG4Tet(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4Tet(module));
}
