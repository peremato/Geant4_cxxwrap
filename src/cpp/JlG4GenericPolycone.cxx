// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4GenericPolycone> : std::false_type { };
  template<> struct DefaultConstructible<G4GenericPolycone> : std::false_type { };
template<> struct SuperType<G4GenericPolycone> { typedef G4VCSGfaceted type; };
}

// Class generating the wrapper for type G4GenericPolycone
// signature to use in the veto file: G4GenericPolycone
struct JlG4GenericPolycone: public Wrapper {

  JlG4GenericPolycone(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4GenericPolycone (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:65:7
    jlcxx::TypeWrapper<G4GenericPolycone>  t = jlModule.add_type<G4GenericPolycone>("G4GenericPolycone",
      jlcxx::julia_base_type<G4VCSGfaceted>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4GenericPolycone>>(new jlcxx::TypeWrapper<G4GenericPolycone>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4GenericPolycone::G4GenericPolycone(const G4String &, G4double, G4double, G4int, const G4double [], const G4double []) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:69:5
    t.constructor<const G4String &, G4double, G4double, G4int, const G4double [], const G4double []>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for EInside G4GenericPolycone::Inside(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: EInside G4GenericPolycone::Inside(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:78:13
    t.method("Inside", static_cast<EInside (G4GenericPolycone::*)(const G4ThreeVector &)  const>(&G4GenericPolycone::Inside));

    DEBUG_MSG("Adding wrapper for G4double G4GenericPolycone::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GenericPolycone::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:79:14
    t.method("DistanceToIn", static_cast<G4double (G4GenericPolycone::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4GenericPolycone::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4GenericPolycone::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GenericPolycone::DistanceToIn(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:81:14
    t.method("DistanceToIn", static_cast<G4double (G4GenericPolycone::*)(const G4ThreeVector &)  const>(&G4GenericPolycone::DistanceToIn));

    DEBUG_MSG("Adding wrapper for void G4GenericPolycone::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4GenericPolycone::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:83:10
    t.method("BoundingLimits", static_cast<void (G4GenericPolycone::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4GenericPolycone::BoundingLimits));

    DEBUG_MSG("Adding wrapper for G4double G4GenericPolycone::GetCubicVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GenericPolycone::GetCubicVolume()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:89:14
    t.method("GetCubicVolume", static_cast<G4double (G4GenericPolycone::*)() >(&G4GenericPolycone::GetCubicVolume));

    DEBUG_MSG("Adding wrapper for G4double G4GenericPolycone::GetSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GenericPolycone::GetSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:90:14
    t.method("GetSurfaceArea", static_cast<G4double (G4GenericPolycone::*)() >(&G4GenericPolycone::GetSurfaceArea));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4GenericPolycone::GetPointOnSurface() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4GenericPolycone::GetPointOnSurface()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:92:19
    t.method("GetPointOnSurface", static_cast<G4ThreeVector (G4GenericPolycone::*)()  const>(&G4GenericPolycone::GetPointOnSurface));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4GenericPolycone::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4GenericPolycone::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:94:20
    t.method("GetEntityType", static_cast<G4GeometryType (G4GenericPolycone::*)()  const>(&G4GenericPolycone::GetEntityType));

    DEBUG_MSG("Adding wrapper for G4VSolid * G4GenericPolycone::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4VSolid * G4GenericPolycone::Clone()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:96:15
    t.method("Clone", static_cast<G4VSolid * (G4GenericPolycone::*)()  const>(&G4GenericPolycone::Clone));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4GenericPolycone::CreatePolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4GenericPolycone::CreatePolyhedron()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:100:19
    t.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4GenericPolycone::*)()  const>(&G4GenericPolycone::CreatePolyhedron));

    DEBUG_MSG("Adding wrapper for G4bool G4GenericPolycone::Reset() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4GenericPolycone::Reset()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:102:12
    t.method("Reset", static_cast<G4bool (G4GenericPolycone::*)() >(&G4GenericPolycone::Reset));

    DEBUG_MSG("Adding wrapper for G4double G4GenericPolycone::GetStartPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GenericPolycone::GetStartPhi()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:106:21
    t.method("GetStartPhi", static_cast<G4double (G4GenericPolycone::*)()  const>(&G4GenericPolycone::GetStartPhi));

    DEBUG_MSG("Adding wrapper for G4double G4GenericPolycone::GetEndPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GenericPolycone::GetEndPhi()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:107:21
    t.method("GetEndPhi", static_cast<G4double (G4GenericPolycone::*)()  const>(&G4GenericPolycone::GetEndPhi));

    DEBUG_MSG("Adding wrapper for G4double G4GenericPolycone::GetSinStartPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GenericPolycone::GetSinStartPhi()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:108:21
    t.method("GetSinStartPhi", static_cast<G4double (G4GenericPolycone::*)()  const>(&G4GenericPolycone::GetSinStartPhi));

    DEBUG_MSG("Adding wrapper for G4double G4GenericPolycone::GetCosStartPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GenericPolycone::GetCosStartPhi()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:109:21
    t.method("GetCosStartPhi", static_cast<G4double (G4GenericPolycone::*)()  const>(&G4GenericPolycone::GetCosStartPhi));

    DEBUG_MSG("Adding wrapper for G4double G4GenericPolycone::GetSinEndPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GenericPolycone::GetSinEndPhi()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:110:21
    t.method("GetSinEndPhi", static_cast<G4double (G4GenericPolycone::*)()  const>(&G4GenericPolycone::GetSinEndPhi));

    DEBUG_MSG("Adding wrapper for G4double G4GenericPolycone::GetCosEndPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GenericPolycone::GetCosEndPhi()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:111:21
    t.method("GetCosEndPhi", static_cast<G4double (G4GenericPolycone::*)()  const>(&G4GenericPolycone::GetCosEndPhi));

    DEBUG_MSG("Adding wrapper for G4bool G4GenericPolycone::IsOpen() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4GenericPolycone::IsOpen()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:112:19
    t.method("IsOpen", static_cast<G4bool (G4GenericPolycone::*)()  const>(&G4GenericPolycone::IsOpen));

    DEBUG_MSG("Adding wrapper for G4int G4GenericPolycone::GetNumRZCorner() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4GenericPolycone::GetNumRZCorner()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:113:19
    t.method("GetNumRZCorner", static_cast<G4int (G4GenericPolycone::*)()  const>(&G4GenericPolycone::GetNumRZCorner));

    DEBUG_MSG("Adding wrapper for G4PolyconeSideRZ G4GenericPolycone::GetCorner(G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4PolyconeSideRZ G4GenericPolycone::GetCorner(G4int)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:114:29
    t.method("GetCorner", static_cast<G4PolyconeSideRZ (G4GenericPolycone::*)(G4int)  const>(&G4GenericPolycone::GetCorner));


    DEBUG_MSG("Adding wrapper for void G4GenericPolycone::G4GenericPolycone(const G4GenericPolycone &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:121:5
    t.constructor<const G4GenericPolycone &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for G4GenericPolycone & G4GenericPolycone::operator=(const G4GenericPolycone &) (" __HERE__ ")");
    // signature to use in the veto list: G4GenericPolycone & G4GenericPolycone::operator=(const G4GenericPolycone &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4GenericPolycone.hh:122:24
    t.method("assign", static_cast<G4GenericPolycone & (G4GenericPolycone::*)(const G4GenericPolycone &) >(&G4GenericPolycone::operator=));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4GenericPolycone>> type_;
};
std::shared_ptr<Wrapper> newJlG4GenericPolycone(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4GenericPolycone(module));
}