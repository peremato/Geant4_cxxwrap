// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4Polyhedra> : std::false_type { };
  template<> struct DefaultConstructible<G4Polyhedra> : std::false_type { };
template<> struct SuperType<G4Polyhedra> { typedef G4VCSGfaceted type; };
}

// Class generating the wrapper for type G4Polyhedra
// signature to use in the veto file: G4Polyhedra
struct JlG4Polyhedra: public Wrapper {

  JlG4Polyhedra(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4Polyhedra (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:74:7
    jlcxx::TypeWrapper<G4Polyhedra>  t = jlModule.add_type<G4Polyhedra>("G4Polyhedra",
      jlcxx::julia_base_type<G4VCSGfaceted>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4Polyhedra>>(new jlcxx::TypeWrapper<G4Polyhedra>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4Polyhedra::G4Polyhedra(const G4String &, G4double, G4double, G4int, G4int, const G4double [], const G4double [], const G4double []) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:78:5
    t.constructor<const G4String &, G4double, G4double, G4int, G4int, const G4double [], const G4double [], const G4double []>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void G4Polyhedra::G4Polyhedra(const G4String &, G4double, G4double, G4int, G4int, const G4double [], const G4double []) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:87:5
    t.constructor<const G4String &, G4double, G4double, G4int, G4int, const G4double [], const G4double []>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for EInside G4Polyhedra::Inside(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: EInside G4Polyhedra::Inside(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:97:13
    t.method("Inside", static_cast<EInside (G4Polyhedra::*)(const G4ThreeVector &)  const>(&G4Polyhedra::Inside));

    DEBUG_MSG("Adding wrapper for G4double G4Polyhedra::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Polyhedra::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:98:14
    t.method("DistanceToIn", static_cast<G4double (G4Polyhedra::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Polyhedra::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4Polyhedra::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Polyhedra::DistanceToIn(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:100:14
    t.method("DistanceToIn", static_cast<G4double (G4Polyhedra::*)(const G4ThreeVector &)  const>(&G4Polyhedra::DistanceToIn));

    DEBUG_MSG("Adding wrapper for void G4Polyhedra::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Polyhedra::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:102:10
    t.method("BoundingLimits", static_cast<void (G4Polyhedra::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4Polyhedra::BoundingLimits));

    DEBUG_MSG("Adding wrapper for void G4Polyhedra::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Polyhedra::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:108:10
    t.method("ComputeDimensions", static_cast<void (G4Polyhedra::*)(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) >(&G4Polyhedra::ComputeDimensions));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4Polyhedra::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4Polyhedra::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:112:21
    t.method("GetEntityType", static_cast<G4GeometryType (G4Polyhedra::*)()  const>(&G4Polyhedra::GetEntityType));

    DEBUG_MSG("Adding wrapper for G4VSolid * G4Polyhedra::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4VSolid * G4Polyhedra::Clone()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:114:15
    t.method("Clone", static_cast<G4VSolid * (G4Polyhedra::*)()  const>(&G4Polyhedra::Clone));

    DEBUG_MSG("Adding wrapper for G4double G4Polyhedra::GetCubicVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Polyhedra::GetCubicVolume()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:116:14
    t.method("GetCubicVolume", static_cast<G4double (G4Polyhedra::*)() >(&G4Polyhedra::GetCubicVolume));

    DEBUG_MSG("Adding wrapper for G4double G4Polyhedra::GetSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Polyhedra::GetSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:117:14
    t.method("GetSurfaceArea", static_cast<G4double (G4Polyhedra::*)() >(&G4Polyhedra::GetSurfaceArea));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Polyhedra::GetPointOnSurface() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Polyhedra::GetPointOnSurface()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:119:19
    t.method("GetPointOnSurface", static_cast<G4ThreeVector (G4Polyhedra::*)()  const>(&G4Polyhedra::GetPointOnSurface));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Polyhedra::CreatePolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4Polyhedra::CreatePolyhedron()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:123:19
    t.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4Polyhedra::*)()  const>(&G4Polyhedra::CreatePolyhedron));

    DEBUG_MSG("Adding wrapper for G4bool G4Polyhedra::Reset() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Polyhedra::Reset()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:125:12
    t.method("Reset", static_cast<G4bool (G4Polyhedra::*)() >(&G4Polyhedra::Reset));

    DEBUG_MSG("Adding wrapper for G4int G4Polyhedra::GetNumSide() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Polyhedra::GetNumSide()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:129:18
    t.method("GetNumSide", static_cast<G4int (G4Polyhedra::*)()  const>(&G4Polyhedra::GetNumSide));

    DEBUG_MSG("Adding wrapper for G4double G4Polyhedra::GetStartPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Polyhedra::GetStartPhi()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:130:21
    t.method("GetStartPhi", static_cast<G4double (G4Polyhedra::*)()  const>(&G4Polyhedra::GetStartPhi));

    DEBUG_MSG("Adding wrapper for G4double G4Polyhedra::GetEndPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Polyhedra::GetEndPhi()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:131:21
    t.method("GetEndPhi", static_cast<G4double (G4Polyhedra::*)()  const>(&G4Polyhedra::GetEndPhi));

    DEBUG_MSG("Adding wrapper for G4double G4Polyhedra::GetSinStartPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Polyhedra::GetSinStartPhi()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:132:21
    t.method("GetSinStartPhi", static_cast<G4double (G4Polyhedra::*)()  const>(&G4Polyhedra::GetSinStartPhi));

    DEBUG_MSG("Adding wrapper for G4double G4Polyhedra::GetCosStartPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Polyhedra::GetCosStartPhi()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:133:21
    t.method("GetCosStartPhi", static_cast<G4double (G4Polyhedra::*)()  const>(&G4Polyhedra::GetCosStartPhi));

    DEBUG_MSG("Adding wrapper for G4double G4Polyhedra::GetSinEndPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Polyhedra::GetSinEndPhi()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:134:21
    t.method("GetSinEndPhi", static_cast<G4double (G4Polyhedra::*)()  const>(&G4Polyhedra::GetSinEndPhi));

    DEBUG_MSG("Adding wrapper for G4double G4Polyhedra::GetCosEndPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Polyhedra::GetCosEndPhi()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:135:21
    t.method("GetCosEndPhi", static_cast<G4double (G4Polyhedra::*)()  const>(&G4Polyhedra::GetCosEndPhi));

    DEBUG_MSG("Adding wrapper for G4bool G4Polyhedra::IsOpen() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Polyhedra::IsOpen()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:136:19
    t.method("IsOpen", static_cast<G4bool (G4Polyhedra::*)()  const>(&G4Polyhedra::IsOpen));

    DEBUG_MSG("Adding wrapper for G4bool G4Polyhedra::IsGeneric() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Polyhedra::IsGeneric()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:137:19
    t.method("IsGeneric", static_cast<G4bool (G4Polyhedra::*)()  const>(&G4Polyhedra::IsGeneric));

    DEBUG_MSG("Adding wrapper for G4int G4Polyhedra::GetNumRZCorner() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Polyhedra::GetNumRZCorner()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:138:18
    t.method("GetNumRZCorner", static_cast<G4int (G4Polyhedra::*)()  const>(&G4Polyhedra::GetNumRZCorner));

    DEBUG_MSG("Adding wrapper for G4PolyhedraSideRZ G4Polyhedra::GetCorner(const G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4PolyhedraSideRZ G4Polyhedra::GetCorner(const G4int)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:139:30
    t.method("GetCorner", static_cast<G4PolyhedraSideRZ (G4Polyhedra::*)(const G4int)  const>(&G4Polyhedra::GetCorner));

    DEBUG_MSG("Adding wrapper for G4PolyhedraHistorical * G4Polyhedra::GetOriginalParameters() (" __HERE__ ")");
    // signature to use in the veto list: G4PolyhedraHistorical * G4Polyhedra::GetOriginalParameters()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:141:35
    t.method("GetOriginalParameters", static_cast<G4PolyhedraHistorical * (G4Polyhedra::*)()  const>(&G4Polyhedra::GetOriginalParameters));

    DEBUG_MSG("Adding wrapper for void G4Polyhedra::SetOriginalParameters(G4PolyhedraHistorical *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Polyhedra::SetOriginalParameters(G4PolyhedraHistorical *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:143:17
    t.method("SetOriginalParameters", static_cast<void (G4Polyhedra::*)(G4PolyhedraHistorical *) >(&G4Polyhedra::SetOriginalParameters));


    DEBUG_MSG("Adding wrapper for void G4Polyhedra::G4Polyhedra(const G4Polyhedra &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:153:5
    t.constructor<const G4Polyhedra &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for G4Polyhedra & G4Polyhedra::operator=(const G4Polyhedra &) (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedra & G4Polyhedra::operator=(const G4Polyhedra &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Polyhedra.hh:154:18
    t.method("assign", static_cast<G4Polyhedra & (G4Polyhedra::*)(const G4Polyhedra &) >(&G4Polyhedra::operator=));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4Polyhedra>> type_;
};
std::shared_ptr<Wrapper> newJlG4Polyhedra(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4Polyhedra(module));
}
